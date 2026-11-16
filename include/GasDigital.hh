// GasDigital —— 数字化核心：能量沉积 → 电子数 → 芯片平面上的高斯抽样
//
// 这套流程**不做漂移、不做扩散**，所以 Garfield++ 在这里只负责一件事：
// 把气体定义清楚、给出平均电离能 W（Magboltz 的气体表里带）。
// 电子数和落点仍然是 "edep/W + N(0, sigma)"，sigma 是芯片读出分辨率。
//
// 两个后端：
//   Gauss     —— 自带常数 W（默认 33.97 eV，干燥空气），零依赖；
//   Garfield  —— 用 Garfield::MediumMagboltz 建空气、取 W
//                （需要编译时 -DUSE_GARFIELD，见 CMakeLists.txt）。
//
// 单位约定：edep 用 MeV，W 用 eV，sigma / 坐标用 mm（输出文件的单位）。
#ifndef GasDigital_h
#define GasDigital_h 1

#include <array>
#include <memory>
#include <random>
#include <vector>

class GasDigital
{
  public:
    enum class Backend { Gauss, Garfield };

    // w_eV     : gauss 后端的平均电离能 W [eV]；garfield 后端会用气体表里的值覆盖
    // fano     : Fano 因子，描述电离涨落 Var(N) = F·N；<= 0 表示不加涨落
    // sigma_mm : 芯片位置分辨率 [mm]，160 um → 0.16
    explicit GasDigital(double w_eV     = 33.97,
                        double fano     = 0.2,
                        double sigma_mm = 0.16,
                        unsigned long seed = 12345);
    ~GasDigital();

    GasDigital(const GasDigital&) = delete;
    GasDigital& operator=(const GasDigital&) = delete;

    // ---- 参数 ----
    void SetFano(double fano)       { fFano = fano; }
    void SetSigma(double sigma_mm)  { fSigma = sigma_mm; }
    void SetSeed(unsigned long s)   { fRng.seed(s); }

    // 切后端；garfield 后端没编进来时打印警告并保持 gauss
    void SetBackend(Backend b);
    static bool BackendBuiltIn(Backend b);   // 该后端是否编译进来了

    const char* BackendName() const;
    double GetW() const;      // garfield 后端返回气体表里的 W
    double GetFano() const  { return fFano; }
    double GetSigma() const { return fSigma; }

    // ---- (1) 平均电子数：一段能量沉积平均能电离出多少个电子 ----
    double MeanElectrons(double edep_MeV) const { return edep_MeV * 1e6 / fW; }

    // ---- (2) 带涨落的电子数（整数；这是"计数"，不是物理量）----
    long SampleElectrons(double edep_MeV);

    // ---- (3) 位置抽样：(x, z) 是被测到的落点，返回芯片上实际的落点 ----
    std::array<double, 2> Smear(double x, double z);
    void SmearInPlace(double& x, double& z);

    // 每个电子各自独立抽一次（要"逐个电子打到芯片上"时用）
    std::vector<std::array<double, 2>> SmearMany(double x, double z, long n);

  private:
    void BuildGas();               // 建 garfield 后端的气体（懒加载）

    double fW     = 33.97;   // eV
    double fFano  = 0.2;
    double fSigma = 0.16;    // mm
    bool   fWFromGas = false;

    Backend fBackend = Backend::Gauss;
    std::mt19937_64 fRng;

    struct GasImpl;           // 内含 Garfield::MediumMagboltz，只在 .cc 里可见
    std::unique_ptr<GasImpl> fGas;
};

#endif
