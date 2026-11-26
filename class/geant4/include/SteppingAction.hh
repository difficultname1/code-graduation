// 每个 step 看一眼：能量是不是沉积在柱子里？是就交给 EventAction 记下来。
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class G4Step;
class EventAction;
class DetectorConstruction;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(EventAction* eventAction, const DetectorConstruction* detector)
      : fEventAction(eventAction), fDetector(detector) {}

    void UserSteppingAction(const G4Step* step) override;

  private:
    EventAction* fEventAction = nullptr;
    const DetectorConstruction* fDetector = nullptr;
};

#endif
