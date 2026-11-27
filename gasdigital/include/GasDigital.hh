// GasDigital —— 数字化核心：能量沉积 → 电离簇团 → 逐个电子的产生点 → 漂移扩散 → 芯片落点
//
// 漂移只算**横向扩散**：电子沿 Y 从产生点漂到芯片（y = fChipY），横向偏移服从
//   sigma_d(y) = D * sqrt(L)，  L = |y - fChipY| 是该电子要漂的距离 [cm]
// 其中 D 是约化扩散系数 [um/sqrt(cm)]。所以离芯片越远（高度越高），落点越模糊——
// 这是之前缺的那一块。纵向扩散沿漂移方向，不影响芯片上的 (x, z)，这里不模拟。
// Garfield++ 在这套流程里只负责一件事：把气体定义清楚、给出平均电离能 W
// 和 Fano 因子（Magboltz / Heed 的气体表里带）。
//
// 物理模型分三层：
//   (1) 计数涨落：一段 edep 平均电离 edep/W 个电子，**实际个数不是这个确定值**，
//       而是抽出来的（期望值小于 10 时用泊松，否则用 Fano 修正的高斯：
//       Var(N) = F·N）。
//   (2) 簇团结构：Geant4 的一步是沿路径连续沉积的。先把一步按 fSubStep 切成子段，
//       子段分到的 edep 再电离成若干**簇团**：簇心沿子段随机、簇大小取几何分布
//       （均值 fClusterMean），同一个簇团里的电子产生点彼此靠近（fClusterSigma）。
//       所以电子不是"全堆在子段的一个点上"，每个电子都有自己的产生坐标。
//   (3) 读出：每个电子各自按 N(0, sigma) 抽一次落点，sigma 是芯片读出分辨率。
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
    // sigma_mm : 电离倍增小孔的定位宽度 [mm]，80 um → 0.08（芯片本身不加 sigma）
    explicit GasDigital(double w_eV     = 33.97,
                        double fano     = 0.2,
                        double sigma_mm = 0.08,
                        unsigned long seed = 12345);
    ~GasDigital();

    GasDigital(const GasDigital&) = delete;
    GasDigital& operator=(const GasDigital&) = delete;

    // ---- 参数 ----
    void SetFano(double fano)       { fFano = fano; }
    void SetSigma(double sigma_mm)  { fSigma = sigma_mm; }
    void SetSeed(unsigned long s)   { fRng.seed(s); }

    // 漂移横向扩散：约化扩散系数 D [mm/sqrt(cm)]（85 um/sqrt(cm) → 0.085），
    // 芯片平面 y = chip_mm。sigma_d = D * sqrt(|y - chip|/10)
    void   SetDiffusion(double D_mm_sqrtcm) { fDiffusion = (D_mm_sqrtcm > 0.0) ? D_mm_sqrtcm : 0.0; }
    void   SetChipY(double chip_mm)         { fChipY = chip_mm; }
    double GetDiffusion() const             { return fDiffusion; }
    double GetChipY() const                 { return fChipY; }

    // 某个 y 处的横向扩散 sigma [mm]（不含读出分辨率）
    double DiffusionSigma(double y) const;

    // 电离簇团：一个簇团平均几个电子、簇内产生点的空间涨落 [mm]
    void   SetClusterMean(double n)   { fClusterMean = (n > 1.0) ? n : 1.0; }
    void   SetClusterSigma(double mm) { fClusterSigma = (mm > 0.0) ? mm : 0.0; }
    double GetClusterMean() const     { return fClusterMean; }
    double GetClusterSigma() const    { return fClusterSigma; }

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

    // 一个电离簇团里的电子数（几何分布，均值 fClusterMean，最小 1）
    long SampleClusterSize();

    // ---- (2b) 一个电子的产生点（气体里，数字化前的坐标）----
    // cluster 是它所属的簇团编号（同一次 Ionize/IonizeStep 调用内从 0 起，
    // 同一个编号的电子在空间上是挨在一起的）
    struct Electron { double x, y, z, t; int cluster; };

    // 一段直线：Geant4 的一步被切成若干这样的子段，每段有自己的 edep
    struct Segment {
      double x0, y0, z0, t0;   // 子段起点
      double x1, y1, z1, t1;   // 子段终点
      double edep;             // 这段分到的沉积能量 [MeV]
    };

    // 把一步按 fSubStep 切成子段，edep 平分（fSubStep <= 0 或长度 0 时退化成
    // 一个"零长度"子段，即老行为：整步的能量都算在起点上）
    std::vector<Segment> SplitSegments(double x0, double y0, double z0, double t0,
                                       double x1, double y1, double z1, double t1,
                                       double edep);

    // 把一段沉积电离成逐个电子：
    //   ① 电子数按泊松/Fano 涨落抽（不是 edep/W 的确定值）
    //   ② 电子按簇团分组，簇心沿子段随机，簇内产生点按 fClusterSigma 散开
    // 返回的每个电子带自己的产生坐标和时间，调用者逐个处理（各自抽样落点）。
    std::vector<Electron> Ionize(const Segment& seg);

    // SplitSegments + Ionize 的合并版：直接给一步的两端和 edep
    std::vector<Electron> IonizeStep(double x0, double y0, double z0, double t0,
                                     double x1, double y1, double z1, double t1,
                                     double edep);

    // ---- (3) 位置抽样：(x, z) 是被测到的落点，返回芯片上实际的落点 ----
    // SmearAt 用产生点的 y 算漂移距离，sigma = sqrt(sigma_read^2 + sigma_d(y)^2)；
    // 两维独立抽样。Smear(x, z) 等价于把产生点放在芯片平面上（无扩散）。
    std::array<double, 2> SmearAt(double x, double y, double z);
    std::array<double, 2> Smear(double x, double z);
    void SmearInPlace(double& x, double& z);

    // 每个电子各自独立抽一次（要"逐个电子打到芯片上"时用）
    std::vector<std::array<double, 2>> SmearMany(double x, double y, double z, long n);

    // ---- (4) 一步内部插值（旧接口：每子段只给一个点，没有簇团结构）----
    // 一段直线上的一个子点：位置、时间、分到的沉积能量
    struct SubHit { double x, y, z, t, edep; };

    // 把 edep 沿 (x0,y0,z0,t0) -> (x1,y1,z1,t1) 这段直线均匀摊开，
    // 子段长度不超过 fSubStep；每个子点的位置在子段内均匀随机（电离碰撞
    // 本来就是沿程随机发生的，固定取中点会让位置落在规则的格点上）。
    // fSubStep <= 0、或这一步长度为 0 时只返回一个点：起点（老行为）。
    // 新代码请用 Ionize/IonizeStep：那里每个电子有自己的产生坐标。
    std::vector<SubHit> SplitStep(double x0, double y0, double z0, double t0,
                                  double x1, double y1, double z1, double t1,
                                  double edep);

    // 插值的子段长度 [mm]；<= 0 表示不插值（每步只取一个点）
    void   SetSubStep(double mm) { fSubStep = mm; }
    double GetSubStep() const    { return fSubStep; }

  private:
    void BuildGas();               // 建 garfield 后端的气体（懒加载）

    double fW     = 33.97;   // eV
    double fFano  = 0.2;
    double fSigma = 0.08;    // mm：电离倍增小孔的定位宽度（芯片无额外 sigma）
    double fDiffusion = 0.085; // mm/sqrt(cm) = 85 um/sqrt(cm)：漂移横向扩散
    double fChipY = 0.0;       // mm：芯片所在平面
    double fSubStep = 0.002; // mm：一步内部插值的子段长度（2 um）
    // 电离簇团：空气里原初电离簇的平均大小约 2~3 个电子（Heed 默认分布均值
    // 同量级）；簇内产生点在气体里热化尺度上是很小的（几十 nm），相对 80 um
    // 小孔定位宽度可以忽略——想要看到它的效果就把 sigma 调小、或把这里改大。
    double fClusterMean  = 2.5;      // 电子/簇团
    double fClusterSigma = 5e-5;     // mm = 50 nm
    bool   fWFromGas = false;

    Backend fBackend = Backend::Gauss;
    std::mt19937_64 fRng;

    struct GasImpl;           // 内含 Garfield::MediumMagboltz，只在 .cc 里可见
    std::unique_ptr<GasImpl> fGas;
};

#endif
