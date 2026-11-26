// 只负责记住当前事件号：每有一步在柱子里沉积能量，就立刻转给 RunAction 写一条，
// 不做按柱汇总、不分配任何编号。
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4Event;
class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    explicit EventAction(RunAction* run) : fRun(run) {}

    void BeginOfEventAction(const G4Event* event) override;

    // 由 SteppingAction 调用：这一步从 p0 走到 p1（Geant4 单位），
    // 沉积了 edep；fracInside 是这一段里真正落在气体柱里的比例（0..1）。
    // t0 / t1 是两端的全局时间。
    void AddHit(const G4ThreeVector& p0, const G4ThreeVector& p1,
                G4double fracInside, G4double edep,
                G4double t0, G4double t1);

  private:
    RunAction* fRun = nullptr;
    G4int fEventID = 0;
};

#endif
