#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::Build() const
{
  auto* runAction = new RunAction();
  SetUserAction(runAction);                        // pixel.root

  SetUserAction(new PrimaryGeneratorAction(fDetector));     // 电子

  auto* eventAction = new EventAction(runAction);  // 每根被打到的柱子、每事件一条
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(eventAction, fDetector));   // 逐步记录
}
