#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include "G4VTouchable.hh"

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return;                       // 没有能量沉积：不记

  const G4StepPoint* pre = step->GetPreStepPoint();
  G4TouchableHandle touch = pre->GetTouchableHandle();

  // 只在柱子内部记录（世界和容器里的沉积不算）；不再取柱子编号
  if (touch->GetVolume()->GetLogicalVolume() != fDetector->GetColumn()) return;

  // 这一步发生的位置：步起点就是"碰撞/能损开始"的地方
  const G4ThreeVector pos = pre->GetPosition();

  // 到达时间也用步起点，和上面的坐标配对；数字化时用它挑"最晚到达"的点
  fEventAction->AddHit(pos.x(), pos.y(), pos.z(), edep, pre->GetGlobalTime());
}
