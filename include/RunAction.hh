// 输出：打开 pixel.root，写一棵 hits 树
//
//   hits  每个"在柱子里沉积了能量"的步一条，不再做任何按柱汇总
//              event  事件号（整数，从 0 起）
//              x      这个步起点的 X 坐标 [mm]
//              y      这个步起点的 Y 坐标 [mm]
//              z      这个步起点的 Z 坐标 [mm]
//              edep   这个步在该柱子内沉积的能量 [MeV]
//              t      这个步起点的全局时间 [ns]（数字化时用来挑"最晚到达"的点）
//
//   约定：只有 event 是编号，用整数（"event/I"）；
//         坐标、能量和时间都是物理量，用浮点
//         （"x/D"、"y/D"、"z/D"、"edep/D"、"t/D"）。
#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Types.hh"
#include "G4String.hh"

#include <type_traits>

class G4Run;
class G4GenericMessenger;
class TFile;
class TTree;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run* run) override;

    // 记一步：pos 是这一步起点的坐标（Geant4 单位），edep 是沉积能量，
    // t 是这一步起点的全局时间（Geant4 单位）
    void FillHit(G4int event, G4double x, G4double y, G4double z,
                 G4double edep, G4double t);

  private:
    G4String fFileName = "pixel.root";     // 宏里用 /output/fileName 改
    G4GenericMessenger* fMessenger = nullptr;
    TFile* fFile = nullptr;
    TTree* fTree = nullptr;

    // 写树用的缓冲区
    G4int    bEvent = 0;
    G4double bX = 0.0;
    G4double bY = 0.0;
    G4double bZ = 0.0;
    G4double bEdep = 0.0;
    G4double bT = 0.0;

    // 编号必须是整数：谁把它改成 double，这里就编译不过
    static_assert(std::is_integral_v<decltype(bEvent)>, "event 编号必须是整型");
};

#endif
