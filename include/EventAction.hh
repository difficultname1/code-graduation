// 只负责记住当前事件号：每有一步在柱子里沉积能量，就立刻转给 RunAction 写一条，
// 不做按柱汇总、不分配任何编号。
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Types.hh"

class G4Event;
class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    explicit EventAction(RunAction* run) : fRun(run) {}

    void BeginOfEventAction(const G4Event* event) override;

    // 由 SteppingAction 调用：这一步的起点在 (x, y, z)、沉积了 edep、
    // 到达时间 t
    void AddHit(G4double x, G4double y, G4double z, G4double edep, G4double t);

  private:
    RunAction* fRun = nullptr;
    G4int fEventID = 0;
};

#endif
