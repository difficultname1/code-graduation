#include "GasDigital.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>

#ifdef HAVE_GARFIELD
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#endif

// Garfield 的对象只在 .cc 里出现，头文件保持干净（不装 Garfield 也能编）
struct GasDigital::GasImpl {
#ifdef HAVE_GARFIELD
  Garfield::MediumMagboltz gas;
  Garfield::ComponentAnalyticField cmp;   // 只为把介质挂到 Sensor 上
  Garfield::Sensor sensor;
  Garfield::TrackHeed track;              // 用来取 Heed 的 W / Fano
#endif
};

GasDigital::GasDigital(double w_eV, double fano, double sigma_mm,
                       unsigned long seed)
  : fW(w_eV), fFano(fano), fSigma(sigma_mm), fRng(seed)
{
}

GasDigital::~GasDigital() = default;

bool GasDigital::BackendBuiltIn(Backend b)
{
  switch (b) {
    case Backend::Gauss:
      return true;
    case Backend::Garfield:
#ifdef HAVE_GARFIELD
      return true;
#else
      return false;
#endif
  }
  return false;
}

const char* GasDigital::BackendName() const
{
  switch (fBackend) {
    case Backend::Gauss:    return "gauss";
    case Backend::Garfield: return "garfield";
  }
  return "?";
}

void GasDigital::SetBackend(Backend b)
{
  if (b == Backend::Garfield && !BackendBuiltIn(b)) {
    std::cerr << "[GasDigital] 这个 build 没把 Garfield++ 编进来"
                 "（cmake -DUSE_GARFIELD=ON），后端保持 gauss\n";
    fBackend = Backend::Gauss;
    return;
  }
  fBackend = b;
  if (fBackend == Backend::Garfield) BuildGas();
}

// 只在 garfield 后端用到：建一份干燥空气，并从 Heed 取 W 与 Fano 因子
void GasDigital::BuildGas()
{
#ifdef HAVE_GARFIELD
  if (fGas) return;

  // Heed 的原初电离表靠环境变量找（老版本 Garfield 没有 SetDatabasePath）。
  // 编译时把安装目录写进 GARFIELD_HEED_DB，运行时自动补上，省得用户手设。
#ifdef GARFIELD_HEED_DB
  if (std::getenv("HEED_DATABASE") == nullptr)
    setenv("HEED_DATABASE", GARFIELD_HEED_DB, 1);
  std::cout << "[GasDigital] Heed 数据库: " << std::getenv("HEED_DATABASE")
            << std::endl;
#endif

  fGas = std::make_unique<GasImpl>();

  // 干燥空气（体积比），293.15 K，760 Torr
  auto& gas = fGas->gas;
  gas.SetComposition("N2", 78.08, "O2", 20.95, "Ar", 0.93, "CO2", 0.04);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(false);

  // Heed 从 Sensor 的组件里拿介质；两块板只是让 AnalyticField 的 CellCheck 过关
  fGas->cmp.SetMedium(&gas);
  fGas->cmp.AddPlaneY(-0.5, 0., "bottom");
  fGas->cmp.AddPlaneY(0.5, -1000., "top");

  fGas->sensor.SetArea(-1., -1., -1., 1., 1., 1.);
  fGas->sensor.AddComponent(&fGas->cmp);

  // 模拟一条 100 MeV 电子的径迹：只是为了触发 Heed 的物性计算，
  // 不取它的簇团信息（这里的电子数 / 簇团结构由 GasDigital 自己抽）
  auto& track = fGas->track;
  track.SetSensor(&fGas->sensor);
  track.SetParticle("e-");
  track.SetEnergy(100.e6);

  if (!track.NewTrack(0., 0., 0., 0., 0., 0., 1.)) {
    std::cerr << "[GasDigital] Heed 建径迹失败，W / Fano 保持默认值\n";
    return;
  }

  const double w = track.GetW();
  const double fano = track.GetFanoFactor();
  if (w > 0.) { fW = w; fWFromGas = true; }
  if (fano > 0.) fFano = fano;

  std::cout << "[GasDigital] Garfield 后端: 空气 N2/O2/Ar/CO2 = 78.08/20.95/0.93/0.04,\n"
            << "             W = " << fW << " eV, Fano = " << fFano
            << (fWFromGas ? "（W/Fano 来自 Heed）" : "（Heed 没给值，用默认）")
            << std::endl;
#endif
}

double GasDigital::GetW() const
{
  return fW;
}

// 电离涨落：电子数多时 Var(N) = F·N（Fano）+ 高斯近似；
// 电子数很少时（单个沉积步只有几 eV，期望值远小于 1）高斯的含义就不清楚了，
// 那种情形用泊松更合适——它本来描述的就是"稀有独立事件"。
long GasDigital::SampleElectrons(double edep_MeV)
{
  const double mean = MeanElectrons(edep_MeV);
  if (!(mean > 0.0)) return 0;                       // 也挡住 NaN

  if (mean < 10.0) {
    std::poisson_distribution<long> p(mean);
    return p(fRng);
  }

  if (fFano <= 0.0) return std::llround(mean);       // 关掉涨落

  std::normal_distribution<double> g(mean, std::sqrt(fFano * mean));
  const long n = std::llround(g(fRng));
  return (n < 0) ? 0 : n;
}

// 横向扩散：电子从产生点沿 Y 漂到芯片（y = fChipY），漂移长度 L = |y - fChipY|，
// 扩散 sigma = D * sqrt(L)。D 的单位是 um/sqrt(cm)，这里换算成 mm/sqrt(cm)：
//   L[cm] = |y - chip|[mm] / 10
// 这样"产生得越高（离芯片越远）→ 扩散越大"就自动进来了。
double GasDigital::DiffusionSigma(double y) const
{
  if (!(fDiffusion > 0.0)) return 0.0;
  double L = y - fChipY;                 // mm
  if (L < 0.0) L = -L;                   // 电荷总是朝芯片漂，距离取绝对值
  return fDiffusion * std::sqrt(L / 10.0);
}

// 芯片上的落点：读出分辨率 + 漂移扩散，两个独立的展宽按平方和合成。
// x、z 各抽一次，互相独立（扩散在横向是各向同性的）。
std::array<double, 2> GasDigital::SmearAt(double x, double y, double z)
{
  const double sd = DiffusionSigma(y);
  const double sigma = std::sqrt(fSigma * fSigma + sd * sd);
  if (!(sigma > 0.0)) return { x, z };

  std::normal_distribution<double> g(0.0, sigma);
  return { x + g(fRng), z + g(fRng) };
}

// 没有产生点高度时（旧接口）：把产生点当成就在芯片平面上，只剩读出分辨率
std::array<double, 2> GasDigital::Smear(double x, double z)
{
  return SmearAt(x, fChipY, z);
}

void GasDigital::SmearInPlace(double& x, double& z)
{
  const auto p = Smear(x, z);
  x = p[0];
  z = p[1];
}

std::vector<std::array<double, 2>> GasDigital::SmearMany(double x, double y,
                                                         double z, long n)
{
  std::vector<std::array<double, 2>> out;
  if (n <= 0) return out;
  out.reserve(static_cast<std::size_t>(n));
  for (long i = 0; i < n; ++i) out.push_back(SmearAt(x, y, z));
  return out;
}

// 一步不是"一个点"：电子沿这段路径一路走一路碰，能量是连续沉积的。
// 把这一段按 fSubStep 切成等长小段，每段的能量平均分。
//
// 子段的两个端点都要保留：电离出来的电子要沿这段子段分开摆（簇心沿段随机），
// 只留一个"代表点"就做不到"每个电子有自己的产生坐标"。
std::vector<GasDigital::Segment> GasDigital::SplitSegments(
    double x0, double y0, double z0, double t0,
    double x1, double y1, double z1, double t1, double edep)
{
  std::vector<Segment> out;
  if (!(edep > 0.0)) return out;                 // 也挡住 NaN

  const double dx = x1 - x0;
  const double dy = y1 - y0;
  const double dz = z1 - z0;
  const double len = std::sqrt(dx * dx + dy * dy + dz * dz);

  // 不插值（或这一步长度为 0）：老行为，整步的能量都算在起点上
  if (!(fSubStep > 0.0) || !(len > 0.0)) {
    out.push_back({x0, y0, z0, t0, x0, y0, z0, t0, edep});
    return out;
  }

  long n = static_cast<long>(std::ceil(len / fSubStep));
  if (n < 1) n = 1;
  const long kMaxSub = 64;                       // 防极端长步把内存撑爆
  if (n > kMaxSub) n = kMaxSub;

  out.reserve(static_cast<std::size_t>(n));
  for (long k = 0; k < n; ++k) {
    const double f0 = static_cast<double>(k) / n;
    const double f1 = static_cast<double>(k + 1) / n;
    out.push_back({x0 + f0 * dx, y0 + f0 * dy, z0 + f0 * dz,
                   t0 + f0 * (t1 - t0),
                   x0 + f1 * dx, y0 + f1 * dy, z0 + f1 * dz,
                   t0 + f1 * (t1 - t0), edep / n});
  }
  return out;
}

// 旧接口：每子段只给一个点。位置在子段内均匀随机——取中点的话子点会落在
// 2 um 间隔的固定格点上（固定枪、步长恒定时尤其明显），细直方图上能看见
// 一排"梳齿"。真实碰撞位置本来就是随机的。
std::vector<GasDigital::SubHit> GasDigital::SplitStep(
    double x0, double y0, double z0, double t0,
    double x1, double y1, double z1, double t1, double edep)
{
  std::vector<SubHit> out;
  const auto segs = SplitSegments(x0, y0, z0, t0, x1, y1, z1, t1, edep);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  out.reserve(segs.size());
  for (const auto& s : segs) {
    const double f = uni(fRng);                  // 子段内均匀一点
    out.push_back({s.x0 + f * (s.x1 - s.x0), s.y0 + f * (s.y1 - s.y0),
                   s.z0 + f * (s.z1 - s.z0), s.t0 + f * (s.t1 - s.t0), s.edep});
  }
  return out;
}

namespace {
// 给单位方向 u 造一对与之垂直的单位基矢 (a, b)，用来在"垂直于径迹"的平面里
// 抽簇内电子的横向偏移。
void TransverseBasis(double ux, double uy, double uz,
                     double& ax, double& ay, double& az,
                     double& bx, double& by, double& bz)
{
  // 参考轴：取与 u 夹角较大的那个坐标轴，避免叉乘退化
  double rx = 0.0, ry = 0.0, rz = 1.0;
  if (std::fabs(uz) > 0.9) { rx = 1.0; rz = 0.0; }

  ax = ry * uz - rz * uy;
  ay = rz * ux - rx * uz;
  az = rx * uy - ry * ux;
  const double na = std::sqrt(ax * ax + ay * ay + az * az);
  if (na > 0.0) { ax /= na; ay /= na; az /= na; }
  else          { ax = 1.0; ay = 0.0; az = 0.0; }

  bx = uy * az - uz * ay;                        // b = u × a，已经是单位向量
  by = uz * ax - ux * az;
  bz = ux * ay - uy * ax;
}
}  // namespace

// 一个电离簇团里的电子数：几何分布，均值 fClusterMean，最小 1。
// （气体里原初电离的簇大小分布尾部比几何分布还重，Heed 用 1/n² 型的分布；
//   这里取几何分布是因为它简单、只看一个均值参数，足够描述"簇团不是一个点"。）
long GasDigital::SampleClusterSize()
{
  if (!(fClusterMean > 1.0)) return 1;           // 均值 <= 1：一个电子一簇
  const double p = 1.0 / fClusterMean;           // 1 + Geom(p) 的均值 = 1/p
  std::geometric_distribution<long> geo(p);
  return 1 + geo(fRng);
}

// 把一段沉积电离成逐个电子。
//
// 能量 → 电子数不是确定换算：edep/W 只是期望值，实际个数按泊松/Fano 涨落
// （见 SampleElectrons）。拿到 n 个电子之后，再把它们分成若干簇团：
//   · 簇心沿子段均匀分布（碰撞位置本来就是沿程随机的），时间随位置插值；
//   · 同一簇团里的电子绕簇心做三维各向同性涨落，sigma = fClusterSigma；
//   · 不同簇团沿径迹分开，所以"同一段 edep 的电子"不会全叠在一个点上。
std::vector<GasDigital::Electron> GasDigital::Ionize(const Segment& seg)
{
  std::vector<Electron> out;

  const long n = SampleElectrons(seg.edep);
  if (n <= 0) return out;

  const double dx = seg.x1 - seg.x0;
  const double dy = seg.y1 - seg.y0;
  const double dz = seg.z1 - seg.z0;
  const double len = std::sqrt(dx * dx + dy * dy + dz * dz);

  double ux = 0.0, uy = 0.0, uz = 1.0;           // 径迹方向（零长度时随便给一个）
  if (len > 0.0) { ux = dx / len; uy = dy / len; uz = dz / len; }

  double ax, ay, az, bx, by, bz;
  TransverseBasis(ux, uy, uz, ax, ay, az, bx, by, bz);

  std::uniform_real_distribution<double> uni(0.0, 1.0);
  std::normal_distribution<double> gauss(0.0, 1.0);

  out.reserve(static_cast<std::size_t>(n));
  long remaining = n;
  int  cluster = 0;
  while (remaining > 0) {
    long k = SampleClusterSize();
    if (k > remaining) k = remaining;            // 最后一簇吃掉余数
    remaining -= k;

    // 簇心：沿子段均匀抽一个位置，时间按同一比例插值
    const double f  = uni(fRng);
    const double cx = seg.x0 + f * dx;
    const double cy = seg.y0 + f * dy;
    const double cz = seg.z0 + f * dz;
    const double ct = seg.t0 + f * (seg.t1 - seg.t0);

    for (long j = 0; j < k; ++j) {
      // 簇内产生点：三个方向各自高斯涨落 → 三维各向同性
      const double u = fClusterSigma * gauss(fRng);   // 横向 a
      const double v = fClusterSigma * gauss(fRng);   // 横向 b
      const double w = fClusterSigma * gauss(fRng);   // 纵向（沿径迹）
      out.push_back({cx + u * ax + v * bx + w * ux,
                     cy + u * ay + v * by + w * uy,
                     cz + u * az + v * bz + w * uz,
                     ct, cluster});
    }
    ++cluster;
  }
  return out;
}

// 一步 → 子段 → 电离。簇团编号在整步内连续，方便调用者按簇团分组写记录。
std::vector<GasDigital::Electron> GasDigital::IonizeStep(
    double x0, double y0, double z0, double t0,
    double x1, double y1, double z1, double t1, double edep)
{
  std::vector<Electron> out;
  int offset = 0;
  for (const auto& s : SplitSegments(x0, y0, z0, t0, x1, y1, z1, t1, edep)) {
    const auto e = Ionize(s);
    if (e.empty()) continue;
    out.reserve(out.size() + e.size());
    for (const auto& h : e) out.push_back({h.x, h.y, h.z, h.t, h.cluster + offset});
    offset += e.back().cluster + 1;
  }
  return out;
}
