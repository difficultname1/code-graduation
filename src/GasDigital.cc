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
  // 不取它的簇团信息（这套流程用 edep/W 数电子）
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

// 芯片的落点分辨率：x、z 各抽一次，互相独立
std::array<double, 2> GasDigital::Smear(double x, double z)
{
  std::normal_distribution<double> g(0.0, fSigma);
  return { x + g(fRng), z + g(fRng) };
}

void GasDigital::SmearInPlace(double& x, double& z)
{
  const auto p = Smear(x, z);
  x = p[0];
  z = p[1];
}

std::vector<std::array<double, 2>> GasDigital::SmearMany(double x, double z,
                                                         long n)
{
  std::vector<std::array<double, 2>> out;
  if (n <= 0) return out;
  out.reserve(static_cast<std::size_t>(n));
  for (long i = 0; i < n; ++i) out.push_back(Smear(x, z));
  return out;
}
