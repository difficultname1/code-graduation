#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RandomTools.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(const DetectorConstruction* detector)
  : fDetector(detector)
{
  fGun = new G4ParticleGun(1);
  fGun->SetParticleDefinition(
      G4ParticleTable::GetParticleTable()->FindParticle("e-"));
  fGun->SetParticleEnergy(100.0 * MeV);              // 宏里可用 /gun/energy 改
  fGun->SetParticleMomentumDirection(fDetector->GetGunDir());   // 45° 斜入射

  fMessenger = new G4GenericMessenger(this, "/gun/", "电子枪");
  fMessenger->DeclareProperty("randomXY", fRandomXY,
      "true = 起点在阵列截面上均匀随机；false = 固定在枪位（默认）");
  fMessenger->DeclareProperty("randomY", fRandomY,
      "true = x 固定在枪位、y 在柱子全长内均匀随机（扫漂移高度用）");
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fMessenger;
  delete fGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // 默认：固定的枪位（y = 14 mm、左下角外 1 mm），方向固定为 45° 斜入射
  G4ThreeVector pos = fDetector->GetGunPos();

  if (fRandomXY) {
    // 全截面随机：x 在整个阵列宽度内、y 在整根柱子高度内均匀
    pos.setX((2 * G4UniformRand() - 1) * fDetector->GetHalfXZ());
    pos.setY((2 * G4UniformRand() - 1) * fDetector->GetHalfHeight());
  } else if (fRandomY) {
    // 笔形束扫高度：x 留在枪位（阵列正中），y 在整个柱高内均匀
    pos.setY((2 * G4UniformRand() - 1) * fDetector->GetHalfHeight());
  }

  fGun->SetParticlePosition(pos);
  fGun->SetParticleMomentumDirection(fDetector->GetGunDir());
  fGun->GeneratePrimaryVertex(event);
}
