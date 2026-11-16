// 探测器几何：100 x 100 根气体像素柱。
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
    G4double         GetHalfXZ()     const { return fHalfXZ; }      // 阵列半宽（X、Z）
    G4double         GetHalfHeight() const { return fHalfHeight; }  // 柱子半高（Y）

  private:
    G4LogicalVolume* fColumn = nullptr;
    G4int    fN = 0;
    G4ThreeVector fGunPos;
    G4double fHalfXZ = 0.0;
    G4double fHalfHeight = 0.0;
};

#endif
