// -----------------------------------------------------------------------------
//  digitize —— hits(Geant4) -> GasDigital -> digi.root
//
//  对应 docs/digitize_flow.dot 的四段：
//    ① 读入文件    TChain 把所有 hits 树拼起来，读最大事件号
//    ② 循环列表    按 event 切组（树里同一事件的步是连续的）
//    ③ 处理一个事件 每个沉积步 -> 按 --substep 切成子段
//                   每子段 edep -> 泊松/Fano 抽电子数（不是 edep/W 的确定值）
//                   -> 电子分成电离簇团，每个电子有自己的产生坐标
//                   -> 逐个电子漂移到芯片（y = --chip-y），横向扩散 sigma 随
//                      漂移高度增长：sigma_d = D*sqrt(L)，D = --diffusion
//                   -> 再叠加读出分辨率，得到芯片落点
//    ④ 输出文件    默认每个电子一条写进 digi.root 的 digi 树；
//                   --per-cluster 则每个簇团一条（ne = 簇内电子数，位置 = 簇心）
//
//  为什么要在步内插值：Geant4 的一步不是"在起点瞬间沉积 edep"，而是粒子
//  从起点沿直线走到终点、一路上连续电离。只用起点把电荷全堆在步边界上，
//  位置就被量化到 10 um 步长的格点上，分布会出现台阶和假结构。
//
//  为什么每个电子单独摆：edep 落下来是一团团电离（簇团），不是一堆电子叠在
//  同一点。电子数本身有泊松/Fano 涨落，簇团有大有小，簇心沿径迹分开，
//  簇内电子又有自己的产生坐标——这些都由 GasDigital::Ionize 给出。
//
//  编译：cmake --build build -j8        （工程里叫 digitize）
//  运行：./build/digitize pixel.root -o digi.root
//
//  只处理部分事件（不扫全文件）：
//    --event 123            只做事件 123
//    --event-start 1000 --event-count 500   做事件 1000..1499
//    --part 3/16            事件按 16 等分，只做第 3 份（可多进程并行）
//  定位靠 hits 文件里附带的 events 树（event/I, n/L）：行数前缀和就是
//  该事件在 hits 里的起始行，直接 GetEntry 跳过去，不用从头遍历；
//  旧文件没有这棵树时退回二分查找，复杂度仍是 O(log N)。
//
//  单位：输入坐标 mm、edep MeV、t ns；输出坐标 mm。
//  分支类型：event / ne 是编号与计数 -> Int_t；x / z / edep / t 是物理量 -> Double_t。
// -----------------------------------------------------------------------------
#include "GasDigital.hh"

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
  enum class Mode { Latest, Points };

  std::vector<std::string> inputs;
  std::string output   = "digi.root";
  double w_eV          = 33.97;    // gauss 后端的 W（空气）
  double fano          = 0.2;      // 电离涨落 Var(N) = F·N；<0 关掉
  double sigma_mm      = 0.16;     // 160 um 读出分辨率
  double diffusion     = 0.085;    // 85 um/sqrt(cm) 漂移横向扩散
  double chip_y_mm     = -40.0;    // 芯片平面 y = 柱阵底面（80 mm 高，中心在原点）
  double substep_mm    = 0.002;    // 一步内部插值的子段长度（2 um；<=0 关掉）
  double cluster_mean  = 2.5;      // 一个电离簇团平均几个电子
  double cluster_sigma = 5e-5;     // 簇内产生点的空间涨落 [mm] = 50 nm
  // Latest：每个事件只用"最晚到达"的那一个点（旧行为）
  // Points：**所有能量沉积点**都处理（含次级粒子），每个点各自折算电子、各自抽样
  Mode   mode          = Mode::Points;
  bool   per_electron  = true;     // true: 每个电子一条（默认）；false: 每个簇团一条
  GasDigital::Backend backend = GasDigital::Backend::Garfield;
  unsigned long seed   = 12345;
  int    max_events    = -1;       // 只处理前若干个事件（调试用）

  // 只处理一段事件号（并行/调试用）。event_start < 0 表示不限制。
  int    event_start   = -1;
  int    event_count   = -1;       // <= 0 表示一直到最后一个事件
  int    part          = -1;       // --part i/N：按事件号 N 等分，只做第 i 份
  int    nparts        = 1;
};

void Usage(const char* prog)
{
  std::cout <<
    "用法: " << prog << " [选项] 输入文件...\n"
    "\n"
    "  -i, --input  <f>     hits 文件（可给多个；第一个也可以直接写在最后）\n"
    "  -o, --output <f>     输出文件（默认 digi.root）\n"
    "      --backend <b>    gauss | garfield（默认 garfield；没编进来会自动退回 gauss）\n"
    "      --w <eV>         gauss 后端的 W 值（默认 33.97，干燥空气）\n"
    "      --fano <f>       Fano 因子（默认 0.2；<= 0 表示不加涨落）\n"
    "      --sigma <mm>     芯片位置分辨率（默认 0.16 = 160 um）\n"
    "      --diffusion <d> 漂移横向扩散系数 [um/sqrt(cm)]（默认 85）；\n"
    "                       sigma_d = d * sqrt(漂移距离[cm])，产生得越高扩散越大；\n"
    "                       <= 0 表示不考虑扩散\n"
    "      --chip-y <mm>    芯片所在平面 y 坐标（默认 -40 = 柱阵底面）\n"
    "      --substep <mm>   一步内部插值的子段长度（默认 0.002 = 2 um；\n"
    "                       <= 0 表示不插值，每步只有一个点）\n"
    "      --cluster-mean <n>   一个电离簇团平均几个电子（默认 2.5；<= 1 表示一电子一簇）\n"
    "      --cluster-sigma <mm> 簇内电子产生点的空间涨落（默认 5e-5 = 50 nm）\n"
    "      --mode <m>       points | latest（默认 points：所有沉积点都处理；\n"
    "                       latest：每个事件只用最晚到达的一个点）\n"
    "      --per-electron   每个电子单独一条（默认）：位置就是它自己的产生点\n"
    "      --per-cluster    每个电离簇团一条（ne = 簇内电子数，位置 = 簇心）\n"
    "      --seed <n>       随机数种子（默认 12345）\n"
    "      --max-events <n> 只处理前 n 个事件（调试用）\n"
    "      --event <id>     只处理这一个事件号\n"
    "      --event-start <id> --event-count <k>\n"
    "                       只处理事件号 [id, id+k) 这一段\n"
    "      --part <i/N>     把事件号按 N 等分，只做第 i 份（i 从 0 起），\n"
    "                       用来多进程并行数字化：每份一个进程，最后 hadd 合并\n"
    "                       注：上面三个选项靠 hits 里附带的 events 树直接跳到\n"
    "                       目标事件，不扫描前面的数据；旧文件则用二分查找\n";
}

bool ParseArgs(int argc, char** argv, Options& opt)
{
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](const char* what) -> const char* {
      if (i + 1 >= argc) {
        std::cerr << "缺参数: " << what << "\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (a == "-h" || a == "--help") { Usage(argv[0]); std::exit(0); }
    else if (a == "-i" || a == "--input")  opt.inputs.emplace_back(next("--input"));
    else if (a == "-o" || a == "--output") opt.output = next("--output");
    else if (a == "--w")                   opt.w_eV = std::atof(next("--w"));
    else if (a == "--fano")                opt.fano = std::atof(next("--fano"));
    else if (a == "--sigma")               opt.sigma_mm = std::atof(next("--sigma"));
    else if (a == "--diffusion")           opt.diffusion = std::atof(next("--diffusion"));
    else if (a == "--chip-y")              opt.chip_y_mm = std::atof(next("--chip-y"));
    else if (a == "--substep")             opt.substep_mm = std::atof(next("--substep"));
    else if (a == "--cluster-mean")        opt.cluster_mean = std::atof(next("--cluster-mean"));
    else if (a == "--cluster-sigma")       opt.cluster_sigma = std::atof(next("--cluster-sigma"));
    else if (a == "--seed")                opt.seed = std::strtoul(next("--seed"), nullptr, 10);
    else if (a == "--max-events")          opt.max_events = std::atoi(next("--max-events"));
    else if (a == "--event") {
      opt.event_start = std::atoi(next("--event"));
      opt.event_count = 1;
    }
    else if (a == "--event-start")         opt.event_start = std::atoi(next("--event-start"));
    else if (a == "--event-count")         opt.event_count = std::atoi(next("--event-count"));
    else if (a == "--part") {
      const char* s = next("--part");
      if (std::sscanf(s, "%d/%d", &opt.part, &opt.nparts) != 2 ||
          opt.nparts <= 0 || opt.part < 0 || opt.part >= opt.nparts) {
        std::cerr << "--part 需要 i/N 形式，且 0 <= i < N：" << s << "\n";
        return false;
      }
    }
    else if (a == "--per-electron")        opt.per_electron = true;
    else if (a == "--per-cluster")         opt.per_electron = false;
    else if (a == "--mode") {
      const std::string m = next("--mode");
      if (m == "points")      opt.mode = Options::Mode::Points;
      else if (m == "latest") opt.mode = Options::Mode::Latest;
      else { std::cerr << "未知 mode: " << m << "\n"; return false; }
    }
    else if (a == "--backend") {
      const std::string b = next("--backend");
      if (b == "gauss")         opt.backend = GasDigital::Backend::Gauss;
      else if (b == "garfield") opt.backend = GasDigital::Backend::Garfield;
      else { std::cerr << "未知 backend: " << b << "\n"; return false; }
    }
    else if (!a.empty() && a[0] == '-') { std::cerr << "未知选项: " << a << "\n"; return false; }
    else opt.inputs.push_back(a);
  }
  return !opt.inputs.empty();
}

}  // namespace

int main(int argc, char** argv)
{
  Options opt;
  if (!ParseArgs(argc, argv, opt)) { Usage(argv[0]); return 1; }

  // ---------------- ① 读入文件 ----------------
  TChain chain("hits");
  for (const auto& f : opt.inputs) chain.Add(f.c_str());
  if (chain.GetNtrees() == 0) {
    std::cerr << "读不到 hits 树，检查输入文件\n";
    return 1;
  }

  int    ev = 0;
  double hx = 0, hy = 0, hz = 0, edep = 0, ht = 0;
  chain.SetBranchAddress("event", &ev);
  chain.SetBranchAddress("x", &hx);
  chain.SetBranchAddress("y", &hy);
  chain.SetBranchAddress("z", &hz);
  chain.SetBranchAddress("edep", &edep);

  // 老文件没有 t 分支：退化成"同一事件组的最后一条"，并明确警告
  const bool hasT = (chain.GetBranch("t") != nullptr);
  if (hasT) chain.SetBranchAddress("t", &ht);
  else std::cerr << "[警告] hits 里没有 t 分支（旧文件），"
                    "改用每组最后一条作为\"最晚到达\"的点\n";

  // 步的终点（x2,y2,z2）与落在气体里的比例（fin）：有了才能做步内插值。
  // 老文件没有这几个分支，就退化成"每步只用一个点（起点）"。
  double hx2 = 0, hy2 = 0, hz2 = 0, ht2 = 0, hfin = 1.0;
  const bool hasSeg = (chain.GetBranch("x2") != nullptr) &&
                      (chain.GetBranch("y2") != nullptr) &&
                      (chain.GetBranch("z2") != nullptr);
  if (hasSeg) {
    chain.SetBranchAddress("x2", &hx2);
    chain.SetBranchAddress("y2", &hy2);
    chain.SetBranchAddress("z2", &hz2);
    if (chain.GetBranch("t2"))  chain.SetBranchAddress("t2", &ht2);
    if (chain.GetBranch("fin")) chain.SetBranchAddress("fin", &hfin);
  } else {
    std::cerr << "[警告] hits 里没有 x2/y2/z2（步终点），无法在步内插值，"
                 "退化成每步只取起点；重新跑一遍 pixel 就有\n";
  }

  // 事件总数 = 最大事件号 + 1（事件号从 0 起）
  chain.GetEntry(chain.GetEntries() - 1);
  const int maxEvent = ev;
  const Long64_t nEntries = chain.GetEntries();
  std::cout << "[输入] " << opt.inputs.size() << " 个文件, " << nEntries
            << " 条命中, 事件号 0.." << maxEvent << "（共 " << maxEvent + 1
            << " 个事件）\n";

  // ---------------- ② 事件号 -> 起始行 ----------------
  // hits 文件里附带一棵 `events` 树（event/I, n/L）：n 的前缀和就是该事件在
  // hits 里的起始行。hadd 不改数据，合并后重新 cumsum 依然对。
  TChain evchain("events");
  for (const auto& f : opt.inputs) evchain.Add(f.c_str());

  std::vector<int>      idxEvent;   // 事件号（升序）
  std::vector<Long64_t> idxStart;   // 该事件第一行在 hits 里的行号
  bool haveIndex = false;
  if (evchain.GetNtrees() > 0 && evchain.GetBranch("n") != nullptr) {
    int      ie = 0;
    Long64_t in = 0;
    evchain.SetBranchAddress("event", &ie);
    evchain.SetBranchAddress("n", &in);
    const Long64_t nEv = evchain.GetEntries();
    idxEvent.reserve(nEv + 1);
    idxStart.reserve(nEv + 1);
    Long64_t acc = 0;
    for (Long64_t i = 0; i < nEv; ++i) {
      evchain.GetEntry(i);
      idxEvent.push_back(ie);
      idxStart.push_back(acc);
      acc += in;
    }
    idxEvent.push_back(maxEvent + 1);   // 哨兵：查"最后一个事件之后"落在哪
    idxStart.push_back(acc);
    haveIndex = (acc == nEntries);
    if (!haveIndex)
      std::cerr << "[警告] events 树的行数合计 " << acc << " 与 hits 的 "
                << nEntries << " 行不一致，索引作废，改用二分查找\n";
    else
      std::cout << "[索引] 用 events 树直接定位（" << nEv
                << " 个事件，不必扫描全文件）\n";
  } else {
    std::cout << "[索引] 文件里没有 events 树（旧文件），改用二分查找定位\n";
  }

  // 第一个 event >= target 的 hits 行号
  auto LowerBound = [&](int target) -> Long64_t {
    if (haveIndex) {
      const auto it = std::lower_bound(idxEvent.begin(), idxEvent.end(), target);
      return idxStart[static_cast<std::size_t>(it - idxEvent.begin())];
    }
    Long64_t lo = 0, hi = nEntries;
    while (lo < hi) {
      const Long64_t mid = lo + (hi - lo) / 2;
      chain.GetEntry(mid);              // 只为了读 event 分支
      if (ev < target) lo = mid + 1;
      else             hi = mid;
    }
    return lo;
  };

  // ---- 决定要处理的事件号区间 [evStart, evStop) ----
  if (opt.part >= 0) {
    const int total = maxEvent + 1;
    const int base  = total / opt.nparts;
    const int rem   = total % opt.nparts;
    int start = 0;
    for (int j = 0; j < opt.part; ++j) start += base + (j < rem ? 1 : 0);
    opt.event_start = start;
    opt.event_count = base + (opt.part < rem ? 1 : 0);
    std::cout << "[分片] " << opt.part << "/" << opt.nparts << " : 事件 "
              << opt.event_start << ".."
              << (opt.event_start + opt.event_count - 1) << "\n";
  }
  const int evStart = (opt.event_start > 0) ? opt.event_start : 0;
  const int evStop  = (opt.event_count > 0)
                        ? evStart + opt.event_count : maxEvent + 1;

  const Long64_t i0 = LowerBound(evStart);
  const Long64_t i1 = LowerBound(evStop);
  std::cout << "[输入] 处理事件 " << evStart << ".." << (evStop - 1)
            << "  ->  hits 行 " << i0 << ".." << (i1 - 1)
            << "（共 " << (i1 - i0) << " 行，直接定位，不扫描前面的数据）\n";

  // ---------------- 数字化引擎 ----------------
  GasDigital gas(opt.w_eV, opt.fano, opt.sigma_mm, opt.seed);
  gas.SetSubStep(opt.substep_mm);
  gas.SetClusterMean(opt.cluster_mean);
  gas.SetClusterSigma(opt.cluster_sigma);
  gas.SetDiffusion(opt.diffusion * 1e-3);   // 选项单位 um/sqrt(cm)，内部用 mm/sqrt(cm)
  gas.SetChipY(opt.chip_y_mm);
  gas.SetBackend(opt.backend);
  std::cout << "[数字化] 后端 = " << gas.BackendName()
            << ", W = " << gas.GetW() << " eV, Fano = " << gas.GetFano()
            << ", 读出 sigma = " << gas.GetSigma() * 1000 << " um\n"
            << "[数字化] 漂移扩散 D = " << gas.GetDiffusion() * 1000
            << " um/sqrt(cm), 芯片 y = " << gas.GetChipY() << " mm"
            << "（扩散 sigma = D*sqrt(漂移距离)）\n"
            << "[数字化] 步内插值子段 = "
            << (gas.GetSubStep() > 0 ? gas.GetSubStep() * 1000 : 0.0)
            << " um" << (gas.GetSubStep() > 0 ? "" : "（关闭：每步只取起点）")
            << "\n"
            << "[数字化] 电离簇团：平均 " << gas.GetClusterMean()
            << " 电子/簇，簇内产生点 sigma = " << gas.GetClusterSigma() * 1e6
            << " nm" << (gas.GetClusterSigma() > 0 ? "" : "（关闭：一簇一个点）")
            << "\n"
            << "[数字化] 输出粒度："
            << (opt.per_electron ? "每个电子一条" : "每个电离簇团一条")
            << "\n";

  // ---------------- ④ 输出（先建树） ----------------
  TFile fout(opt.output.c_str(), "RECREATE");
  TTree tree("digi", "芯片上的一个落点");

  int    bEvent = 0, bNe = 0, bNstep = 0;
  double bX = 0, bZ = 0, bEdep = 0, bT = 0, bYsrc = 0;

  tree.Branch("event", &bEvent, "event/I");   // 事件号（编号 -> 整型）
  tree.Branch("x",     &bX,     "x/D");       // 芯片 x [mm]
  tree.Branch("z",     &bZ,     "z/D");       // 芯片 z [mm]
  tree.Branch("ne",    &bNe,    "ne/I");      // 电子数（计数 -> 整型；每电子一条时为 1）
  tree.Branch("edep",  &bEdep,  "edep/D");    // 该事件总沉积 [MeV]
  tree.Branch("t",     &bT,     "t/D");       // 这个电子/簇团的产生时间 [ns]
  tree.Branch("ysrc",  &bYsrc,  "ysrc/D");    // 产生点的 y [mm]（算漂移距离用）
  tree.Branch("nstep", &bNstep, "nstep/I");   // 该事件在柱子里有多少步

  // ---------------- ② 循环列表：按 event 切组 ----------------
  // 一个沉积步：起点 (x,y,z,t)、终点 (x2,y2,z2,t2)、这一步的 edep。
  // 终点已经按 fin 裁到"还留在气体里"的位置，所以插值只在这段上做。
  struct Deposit { double x, y, z, t, x2, y2, z2, t2, edep; };

  std::vector<Deposit> buf;      // 当前事件的所有沉积步（含次级粒子）
  int  curEvent = -1;
  long nEventOut = 0, nElectron = 0;
  long nPoint = 0;               // 写出的行数
  long nDeposit = 0;             // 处理过的沉积步数
  long nSubSeg = 0;              // 切出来的子段总数
  long nCluster = 0;             // 电离出来的簇团总数
  double sumX = 0, sumX2 = 0, sumZ = 0, sumZ2 = 0;
  bool warnedOrder = false;

  // 把当前读到的这一条变成 Deposit：终点按 fin 裁到"还在气体里"的那一段。
  // 没有 x2/y2/z2 的老文件、或 fin == 1 的普通步，终点就是终点本身。
  auto MakeDeposit = [&]() {
    const double t0 = hasT ? ht : 0.0;
    double x2 = hx, y2 = hy, z2 = hz, t2 = t0;
    if (hasSeg) {
      double f = (hfin > 0.0) ? hfin : 0.0;
      if (f > 1.0) f = 1.0;
      x2 = hx + f * (hx2 - hx);
      y2 = hy + f * (hy2 - hy);
      z2 = hz + f * (hz2 - hz);
      t2 = t0 + f * (ht2 - t0);
    }
    return Deposit{hx, hy, hz, t0, x2, y2, z2, t2, edep};
  };

  // 写一行：按"读出分辨率 + 该高度的漂移扩散"抽一次落点后落盘
  auto writeRow = [&](int event, double x, double y, double z, double t,
                      double edepEvent, long ne, int nstep) {
    const auto s = gas.SmearAt(x, y, z);
    bEvent = event; bX = s[0]; bZ = s[1]; bNe = static_cast<int>(ne);
    bEdep = edepEvent; bT = t; bYsrc = y; bNstep = nstep;
    tree.Fill();
    sumX += s[0]; sumX2 += s[0] * s[0];
    sumZ += s[1]; sumZ2 += s[1] * s[1];
    ++nPoint;
    nElectron += ne;
  };

  // 把一团电离出来的电子写出去。
  //   per-electron：每个电子一条，位置就是它自己的产生坐标，各自抽一次读出错
  //   per-cluster ：每个簇团一条，位置取簇内电子产生坐标的平均（= 簇心），
  //                 ne = 簇内电子数
  auto writeCloud = [&](int event, const std::vector<GasDigital::Electron>& cloud,
                        double edepEvent, int nstep) {
    if (cloud.empty()) return;

    // 簇团统计：同一片段内的 cluster 编号是递增的，数一下跳变就是簇团数
    ++nCluster;
    for (std::size_t i = 1; i < cloud.size(); ++i)
      if (cloud[i].cluster != cloud[i - 1].cluster) ++nCluster;

    if (opt.per_electron) {
      for (const auto& e : cloud)
        writeRow(event, e.x, e.y, e.z, e.t, edepEvent, 1, nstep);
      return;
    }

    std::size_t i = 0;
    while (i < cloud.size()) {
      const int cid = cloud[i].cluster;
      double sx = 0, sy = 0, sz = 0, st = 0;
      long k = 0;
      for (; i < cloud.size() && cloud[i].cluster == cid; ++i) {
        sx += cloud[i].x; sy += cloud[i].y; sz += cloud[i].z; st += cloud[i].t;
        ++k;
      }
      writeRow(event, sx / k, sy / k, sz / k, st / k, edepEvent, k, nstep);
    }
  };

  // 把一个沉积步数字化：
  //   ① 在步的起点和终点之间切成子段（粒子一路走一路电离，不是一个点）
  //   ② 每个子段的 edep 抽电子数（泊松/Fano 涨落），分成电离簇团，
  //      每个电子有自己的产生坐标
  //   ③ 逐个电子（或按簇团）在芯片上抽样落点
  auto emit = [&](int event, const Deposit& d, double edepEvent, int nstep) {
    const auto segs = gas.SplitSegments(d.x, d.y, d.z, d.t,
                                        d.x2, d.y2, d.z2, d.t2, d.edep);
    nSubSeg += static_cast<long>(segs.size());
    for (const auto& s : segs) {
      const auto cloud = gas.Ionize(s);
      writeCloud(event, cloud, edepEvent, nstep);
    }
  };

  auto flush = [&](int event) {
    if (buf.empty()) return;
    const int nstep = static_cast<int>(buf.size());

    double edepSum = 0;
    for (const auto& p : buf) edepSum += p.edep;

    if (opt.mode == Options::Mode::Latest) {
      // 旧行为：一个事件只取"最晚到达"的那个点，电子数用整个事件的总沉积。
      // 电子仍然按新模型摆：整段沉积退化成一个零长度的子段（簇团叠在那点上）
      const Deposit* best = &buf.front();
      for (const auto& p : buf)
        if (!hasT || p.t > best->t) best = &p;
      const GasDigital::Segment seg{best->x, best->y, best->z, best->t,
                                    best->x, best->y, best->z, best->t, edepSum};
      const auto cloud = gas.Ionize(seg);
      writeCloud(event, cloud, edepSum, nstep);
    } else {
      // 每个沉积步都处理：含次级粒子（δ 射线、散射电子），
      // 步内先切子段、再抽电子数——不挑"主粒子那一个点"
      for (const auto& p : buf) emit(event, p, edepSum, nstep);
    }

    nDeposit += static_cast<long>(buf.size());
    buf.clear();
    ++nEventOut;
  };

  bool stopped = false;
  for (Long64_t i = i0; i < i1; ++i) {
    chain.GetEntry(i);

    if (ev != curEvent) {
      if (curEvent >= 0) {
        flush(curEvent);
        if (opt.max_events > 0 && nEventOut >= opt.max_events) {
          stopped = true;
          break;
        }
      }
      if (ev < curEvent && !warnedOrder) {
        std::cerr << "[警告] 事件号回退（" << curEvent << " -> " << ev
                  << "），树不是按事件排序的，结果会错\n";
        warnedOrder = true;
      }
      curEvent = ev;
    }

    buf.push_back(MakeDeposit());
  }
  if (!stopped) flush(curEvent);

  // ---------------- 收尾 ----------------
  tree.Write();
  const Long64_t nRows = tree.GetEntries();
  fout.Close();

  const double meanX = nPoint ? sumX / nPoint : 0.0;
  const double meanZ = nPoint ? sumZ / nPoint : 0.0;
  const double rmsX = nPoint > 1 ? std::sqrt(sumX2 / nPoint - meanX * meanX) : 0.0;
  const double rmsZ = nPoint > 1 ? std::sqrt(sumZ2 / nPoint - meanZ * meanZ) : 0.0;

  std::cout << "[输出] " << opt.output << " : " << nRows << " 条, "
            << nEventOut << " 个事件, " << nDeposit << " 个沉积步"
            << (opt.mode == Options::Mode::Points ? "（逐点）" : "（最晚点）")
            << (opt.per_electron ? "，每电子一条" : "，每簇团一条") << "\n"
            << "[统计] 步内插值：沉积步 " << nDeposit << " -> 子段 " << nSubSeg
            << "（每步平均 "
            << (nDeposit ? double(nSubSeg) / nDeposit : 0.0) << " 个，子段 "
            << (gas.GetSubStep() > 0 ? gas.GetSubStep() * 1000 : 0.0) << " um）\n"
            << "[统计] 电离：" << nCluster << " 个簇团 -> " << nElectron
            << " 个电子（平均 "
            << (nCluster ? double(nElectron) / nCluster : 0.0) << " 电子/簇）\n"
            << "[统计] 电子总数 " << nElectron << "，平均每事件 "
            << (nEventOut ? double(nElectron) / nEventOut : 0.0)
            << "，每条记录 " << (nPoint ? double(nElectron) / nPoint : 0.0)
            << " 个电子\n"
            << "[统计] 芯片落点 x: mean " << meanX << " mm, RMS " << rmsX
            << " mm; z: mean " << meanZ << " mm, RMS " << rmsZ << " mm\n"
            << "[统计] RMS 里含束斑展宽与漂移扩散，不会被读出 sigma="
            << gas.GetSigma() * 1000 << " um 单独决定\n";
  return 0;
}
