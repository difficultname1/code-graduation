// 电子源：每个事件一个电子，沿 +Z 飞。
//   默认：从固定的枪位（阵列截面中心、紧贴入射面）发射
//   /gun/randomXY true：起点在阵列截面上均匀随机（x ∈ ±半宽，y ∈ ±半高）
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
};

#endif
