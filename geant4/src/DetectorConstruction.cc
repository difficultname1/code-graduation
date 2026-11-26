#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"

#include <cmath>

// ==================== 几何参数：要改就改这里 ====================
namespace {
  const G4double kSide    = 8.3e-5 * m;   // 柱子底面边长（正方形，83 um）
  const G4double kHeight  = 80.0 * mm;    // 柱子高度 = 漂移长度，沿 Y 轴
  const G4int    kN       = 100;          // 每边柱子数，共 kN*kN = 10000 根
  const G4double kMaxStep = 10.0 * um;    // 柱内带电粒子的最大步长
  const G4double kWorld   = 80.0 * mm;    // 世界：边长 80 mm 的立方体
}
// ================================================================
//
// 坐标约定：
//   柱阵沿 Y 轴放置，高度 80 mm；柱阵中心 = 世界立方体的中心 = 原点。
//   电荷沿 Y 漂移到柱子的一端，芯片就贴在那个端面上：y = -kHeight/2。
//   所以电离点离芯片的漂移距离 = y - chipY ∈ [0, 80] mm，
//   扩散 sigma 随这个高度增长（见 digitize 的 --diffusion）。

// World（真空）
//   └─ Array       100x100 根柱子拼起来的整体
//        └─ Slice     沿 Z 复制 kN 片          copyNo = iz
//             └─ Column  沿 X 复制 kN 根柱子    copyNo = ix   ← 柱子本体
//
// 用 G4PVReplica 而不是写 10000 次 G4PVPlacement：
// 复制体自动紧挨着排，导航也快（100 个事件不到 1 秒）。
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  const G4double halfXZ = 0.5 * kN * kSide;   // 阵列在 X、Z 方向的半宽
  const G4double halfH  = 0.5 * kHeight;      // 柱子半高

  // 这几个量别的类要用，存下来
  fN = kN;
  fHalfXZ = halfXZ;
  fHalfHeight = halfH;
  fChipY = -halfH;              // 芯片所在平面（柱阵底面，紧贴柱阵）

  // 电子枪：放在 XOY 平面的中心、紧贴阵列入射面（z = -halfXZ），朝 +Z 打。
  // 注意：kN = 100 是偶数，x = 0 正好是两根柱子（ix=49 / ix=50）的交界面，
  // 电子会沿着柱壁飞，结果可能不稳定（详见 README）。
  fGunPos = G4ThreeVector(0.0, 0.0, -halfXZ);

  G4cout << "[几何] 世界: 边长 " << kWorld / mm << " mm 的真空立方体, 中心在原点\n"
         << "[几何] " << kN << " x " << kN << " 根柱子, 底面 " << kSide / um
         << " um, 高 " << kHeight / mm << " mm, 阵列 " << 2 * halfXZ / um
         << " um x " << 2 * halfXZ / um << " um, 居中于立方体\n"
         << "[几何] 芯片贴在下端面 y = " << fChipY / mm
         << " mm, 漂移长度 0.." << kHeight / mm << " mm\n"
         << "[几何] 电子枪固定在 (" << fGunPos.x() / um << ", "
         << fGunPos.y() / um << ", " << fGunPos.z() / um
         << ") um, 方向 +Z" << G4endl;

  auto* nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

  // 世界：边长 80 mm 的真空立方体（柱阵正好落在它的中心，Y 向占满）
  auto* worldS = new G4Box("World", 0.5 * kWorld, 0.5 * kWorld, 0.5 * kWorld);
  auto* worldL = new G4LogicalVolume(worldS, vacuum, "World");
  auto* worldP = new G4PVPlacement(nullptr, G4ThreeVector(), worldL, "World",
                                   nullptr, false, 0);

  // 阵列容器：大小正好是 100x100 根柱子拼起来
  auto* arrayS = new G4Box("Array", halfXZ, halfH, halfXZ);
  auto* arrayL = new G4LogicalVolume(arrayS, air, "Array");
  new G4PVPlacement(nullptr, G4ThreeVector(), arrayL, "Array", worldL, false, 0);

  // 先沿 Z 切成 kN 片
  auto* sliceS = new G4Box("Slice", halfXZ, halfH, kSide / 2);
  auto* sliceL = new G4LogicalVolume(sliceS, air, "Slice");
  new G4PVReplica("Slice", sliceL, arrayL, kZAxis, kN, kSide);

  // 每片里再沿 X 复制 kN 根柱子
  auto* columnS = new G4Box("Column", kSide / 2, halfH, kSide / 2);
  fColumn = new G4LogicalVolume(columnS, air, "Column");

  // 限制步长，否则电子一步能跨几百 mm，柱子里的"点"就没意义了。
  // 注意：必须配合 main.cc 里的 G4StepLimiterPhysics 才生效。
  fColumn->SetUserLimits(new G4UserLimits(kMaxStep));

  new G4PVReplica("Column", fColumn, sliceL, kXAxis, kN, kSide);

  return worldP;
}

// 柱阵在几何上就是一个盒子（|x|,|z| <= halfXZ，|y| <= halfH），盒子里全是空气柱，
// 盒子外是真空。盒子是凸的，所以一条线段最多只会离开它一次：只要算出离开的那个
// 参数 t，就是这一段留在气体里的比例。
G4double DetectorConstruction::FractionInsideArray(const G4ThreeVector& p0,
                                                   const G4ThreeVector& p1) const
{
  const G4double lo[3] = { -fHalfXZ, -fHalfHeight, -fHalfXZ };
  const G4double hi[3] = {  fHalfXZ,  fHalfHeight,  fHalfXZ };

  G4double tExit = 1.0;
  for (int i = 0; i < 3; ++i) {
    const G4double a = p0[i];
    const G4double d = p1[i] - a;
    if (std::abs(d) < 1e-12) continue;            // 这个方向不动，不会从这儿出去
    const G4double t = (d > 0.0 ? hi[i] - a : lo[i] - a) / d;
    if (t < tExit) tExit = t;                     // 最早离开的那个面
  }

  if (!(tExit > 0.0)) return 0.0;                 // 起点就在边界外/数值异常
  return (tExit > 1.0) ? 1.0 : tExit;
}
