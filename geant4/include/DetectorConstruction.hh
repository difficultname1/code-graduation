// 探测器几何：100 x 100 根气体像素柱（间距 830 um），整体 83 x 40 x 83 mm。
// 具体尺寸、材料都在 src/DetectorConstruction.cc 里，这里只有接口。
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;

    G4VPhysicalVolume* Construct() override;

    // 下面几个是"构造完之后才能取"的几何量，给电子源和逐步记录用
    G4LogicalVolume* GetColumn()     const { return fColumn; }      // 柱子本体
    G4int            GetNPerSide()   const { return fN; }           // 每边多少根
    G4ThreeVector    GetGunPos()     const { return fGunPos; }      // 电子枪位置
    G4ThreeVector    GetGunDir()     const { return fGunDir; }      // 电子枪方向（单位矢量）
    G4double         GetHalfXZ()     const { return fHalfXZ; }      // 阵列半宽（X、Z）
    G4double         GetHalfHeight() const { return fHalfHeight; }  // 柱子半高（Y）
    G4double         GetChipY()      const { return fChipY; }       // 芯片平面 y（柱阵底面）

    // 一条直线段 p0->p1 落在柱阵盒子内部的比例（0..1）。
    // 步的终点还在柱子里时就是 1；最后一步飞出阵列（进真空）时小于 1，
    // 用来把这一步的沉积能量只摊到真正在气体里的那一段上。
    // 前提：p0 在柱阵里面（调用方已按"步起点在柱子里"筛过）。
    G4double FractionInsideArray(const G4ThreeVector& p0,
                                 const G4ThreeVector& p1) const;

  private:
    G4LogicalVolume* fColumn = nullptr;
    G4int    fN = 0;
    G4ThreeVector fGunPos;
    G4ThreeVector fGunDir = G4ThreeVector(0, 0, 1);
    G4double fHalfXZ = 0.0;
    G4double fHalfHeight = 0.0;
    G4double fChipY = 0.0;
};

#endif
