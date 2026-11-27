// 电子源：每个事件一个电子，沿几何里定好的 45° 方向（相对 +Z）飞。
//   默认：从固定的枪位（y = 14 mm、左下角外 1 mm）发射
//   /gun/randomXY true：起点在阵列截面上均匀随机（x ∈ ±半宽，y ∈ ±半高）
//   /gun/randomY  true：x 固定在枪位，y 沿整个柱高均匀（扫漂移高度）
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Types.hh"

class G4ParticleGun;
class G4Event;
class G4GenericMessenger;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    explicit PrimaryGeneratorAction(const DetectorConstruction* detector);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

  private:
    G4ParticleGun* fGun = nullptr;
    const DetectorConstruction* fDetector = nullptr;
    G4GenericMessenger* fMessenger = nullptr;
    G4bool fRandomXY = false;      // true = 起点在阵列截面上随机
    G4bool fRandomY  = false;      // true = x 固定、y 随机（扫高度）
};

#endif
