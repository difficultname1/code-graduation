// -----------------------------------------------------------------------------
//  digitize —— hits(Geant4) -> GasDigital -> digi.root
//
//  对应 docs/digitize_flow.dot 的四段：
//    ① 读入文件    TChain 把所有 hits 树拼起来，读最大事件号
//    ② 循环列表    按 event 切组（树里同一事件的步是连续的）
//    ③ 处理一个事件 取到达时间最晚的点 -> 投影到芯片 (x, z)
//                   Ne = edep / W（带 Fano 涨落）
//                   按 N(0, sigma) 打到芯片上（sigma 默认 160 um，读出分辨率）
//    ④ 输出文件    每个事件一条写进 digi.root 的 digi 树
//
//  编译：cmake --build build -j8        （工程里叫 digitize）
//  运行：./build/digitize pixel.root -o digi.root
//
//  单位：输入坐标 mm、edep MeV、t ns；输出坐标 mm。
//  分支类型：event / ne 是编号与计数 -> Int_t；x / z / edep / t 是物理量 -> Double_t。
// -----------------------------------------------------------------------------
#include "GasDigital.hh"

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
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
  double chip_y_mm     = 0.0;      // 芯片平面 y = chip_y（只影响漂移距离的记录）
  // Latest：每个事件只用"最晚到达"的那一个点（旧行为）
  // Points：**所有能量沉积点**都处理（含次级粒子），每个点各自折算电子、各自抽样
  Mode   mode          = Mode::Points;
  bool   per_electron  = false;    // true: 每个电子一条；false: 每个沉积点一条
  GasDigital::Backend backend = GasDigital::Backend::Garfield;
  unsigned long seed   = 12345;
  int    max_events    = -1;       // 只处理前若干个事件（调试用）
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
    "      --chip-y <mm>    芯片所在平面 y 坐标（默认 0）\n"
    "      --mode <m>       points | latest（默认 points：所有沉积点都处理；\n"
    "                       latest：每个事件只用最晚到达的一个点）\n"
    "      --per-electron   每个电子单独一条（否则每个沉积点一条，ne 为电子数）\n"
    "      --seed <n>       随机数种子（默认 12345）\n"
    "      --max-events <n> 只处理前 n 个事件（调试用）\n";
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
    else if (a == "--chip-y")              opt.chip_y_mm = std::atof(next("--chip-y"));
    else if (a == "--seed")                opt.seed = std::strtoul(next("--seed"), nullptr, 10);
    else if (a == "--max-events")          opt.max_events = std::atoi(next("--max-events"));
    else if (a == "--per-electron")        opt.per_electron = true;
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

  // 事件总数 = 最大事件号 + 1（事件号从 0 起）
  chain.GetEntry(chain.GetEntries() - 1);
  const int maxEvent = ev;
  const Long64_t nEntries = chain.GetEntries();
  std::cout << "[输入] " << opt.inputs.size() << " 个文件, " << nEntries
            << " 条命中, 事件号 0.." << maxEvent << "（共 " << maxEvent + 1
            << " 个事件）\n";

  // ---------------- 数字化引擎 ----------------
  GasDigital gas(opt.w_eV, opt.fano, opt.sigma_mm, opt.seed);
  gas.SetBackend(opt.backend);
  std::cout << "[数字化] 后端 = " << gas.BackendName()
            << ", W = " << gas.GetW() << " eV, Fano = " << gas.GetFano()
            << ", sigma = " << gas.GetSigma() * 1000 << " um\n";

  // ---------------- ④ 输出（先建树） ----------------
  TFile fout(opt.output.c_str(), "RECREATE");
  TTree tree("digi", "芯片上的一个落点");

  int    bEvent = 0, bNe = 0, bNstep = 0;
  double bX = 0, bZ = 0, bEdep = 0, bT = 0, bYsrc = 0;

  tree.Branch("event", &bEvent, "event/I");   // 事件号（编号 -> 整型）
  tree.Branch("x",     &bX,     "x/D");       // 芯片 x [mm]
  tree.Branch("z",     &bZ,     "z/D");       // 芯片 z [mm]
  tree.Branch("ne",    &bNe,    "ne/I");      // 电子数（计数 -> 整型）
  tree.Branch("edep",  &bEdep,  "edep/D");    // 该事件总沉积 [MeV]
  tree.Branch("t",     &bT,     "t/D");       // 被选中那一步的时间 [ns]
  tree.Branch("ysrc",  &bYsrc,  "ysrc/D");    // 被选中那一步的 y [mm]
  tree.Branch("nstep", &bNstep, "nstep/I");   // 该事件在柱子里有多少步

  // ---------------- ② 循环列表：按 event 切组 ----------------
  struct Deposit { double x, y, z, edep, t; };   // 一个能量沉积点

  std::vector<Deposit> buf;      // 当前事件的所有沉积点（含次级粒子）
  int  curEvent = -1;
  long nEventOut = 0, nElectron = 0;
  long nPoint = 0;               // 写出的行数
  long nDeposit = 0;             // 处理过的沉积点数
  double sumX = 0, sumX2 = 0, sumZ = 0, sumZ2 = 0;
  bool warnedOrder = false;

  // 写一行：把这个点的位置按 sigma 抽样后落盘
  auto writeRow = [&](int event, double x, double y, double z, double t,
                      double edepEvent, long ne, int nstep) {
    const auto s = gas.Smear(x, z);
    bEvent = event; bX = s[0]; bZ = s[1]; bNe = static_cast<int>(ne);
    bEdep = edepEvent; bT = t; bYsrc = y; bNstep = nstep;
    tree.Fill();
    sumX += s[0]; sumX2 += s[0] * s[0];
    sumZ += s[1]; sumZ2 += s[1] * s[1];
    ++nPoint;
    nElectron += ne;
  };

  auto flush = [&](int event) {
    if (buf.empty()) return;
    const int nstep = static_cast<int>(buf.size());

    double edepSum = 0;
    for (const auto& p : buf) edepSum += p.edep;

    if (opt.mode == Options::Mode::Latest) {
      // 旧行为：一个事件只取"最晚到达"的那个点，电子数用整个事件的总沉积
      const Deposit* best = &buf.front();
      for (const auto& p : buf)
        if (!hasT || p.t > best->t) best = &p;
      const long ne = gas.SampleElectrons(edepSum);
      if (opt.per_electron) {
        for (long k = 0; k < ne; ++k)
          writeRow(event, best->x, best->y, best->z, best->t, edepSum, 1, nstep);
      } else if (ne > 0) {
        writeRow(event, best->x, best->y, best->z, best->t, edepSum, ne, nstep);
      }
    } else {
      // 每个沉积点都处理：含次级粒子（δ 射线、散射电子），
      // 每点自己折算电子、自己在芯片上抽样位置——不挑"主粒子那一个点"
      for (const auto& p : buf) {
        const long ne = gas.SampleElectrons(p.edep);
        if (ne <= 0) continue;      // 这一步没电离出电子，芯片上就没有信号
        if (opt.per_electron) {
          for (long k = 0; k < ne; ++k)
            writeRow(event, p.x, p.y, p.z, p.t, edepSum, 1, nstep);
        } else {
          writeRow(event, p.x, p.y, p.z, p.t, edepSum, ne, nstep);
        }
      }
    }

    nDeposit += static_cast<long>(buf.size());
    buf.clear();
    ++nEventOut;
  };

  bool stopped = false;
  for (Long64_t i = 0; i < nEntries; ++i) {
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

    buf.push_back({hx, hy, hz, edep, hasT ? ht : 0.0});
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
            << nEventOut << " 个事件, " << nDeposit << " 个沉积点"
            << (opt.mode == Options::Mode::Points ? "（逐点）" : "（最晚点）")
            << (opt.per_electron ? "，每电子一条" : "") << "\n"
            << "[统计] 电子总数 " << nElectron
            << ", 平均每事件 " << (nEventOut ? double(nElectron) / nEventOut : 0.0)
            << ", 平均每条 " << (nPoint ? double(nElectron) / nPoint : 0.0) << "\n"
            << "[统计] 芯片落点 x: mean " << meanX << " mm, RMS " << rmsX
            << " mm; z: mean " << meanZ << " mm, RMS " << rmsZ << " mm\n"
            << "[统计] RMS 里含束斑展宽，不会被 sigma=" << gas.GetSigma() * 1000
            << " um 单独决定（固定枪时 RMS ~ sigma）\n";
  return 0;
}
