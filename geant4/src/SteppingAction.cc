#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return;                       // 没有能量沉积：不记

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();

  // 只在柱子内部记录（世界和容器里的沉积不算）；不再取柱子编号
  if (pre->GetTouchableHandle()->GetVolume()->GetLogicalVolume()
      != fDetector->GetColumn()) return;

  // 一步里有能量沉积，不代表能量都堆在"起点"上：粒子沿着这一整段连续走，
  // 碰撞沿程随机发生。所以起点和终点都要交给 RunAction，数字化时在两者之间插值。
  const G4ThreeVector p0 = pre->GetPosition();
  const G4ThreeVector p1 = post->GetPosition();

  // 终点还在柱子里 -> 整段都在气体中，比例 1；
  // 最后一步飞出阵列（进真空）时，只算真正留在气体里的那一段。
  G4double fracInside = 1.0;
  const G4VPhysicalVolume* postVol = post->GetTouchableHandle()->GetVolume();
  if (postVol == nullptr ||
      postVol->GetLogicalVolume() != fDetector->GetColumn()) {
    fracInside = fDetector->FractionInsideArray(p0, p1);
  }

  fEventAction->AddHit(p0, p1, fracInside, edep,
                       pre->GetGlobalTime(), post->GetGlobalTime(),
                       step->GetTrack()->GetTrackID());
}
