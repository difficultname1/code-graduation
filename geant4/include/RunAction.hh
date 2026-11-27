// 输出：打开 pixel.root，写一棵 hits 树
//
//   hits  每个"在柱子里沉积了能量"的步一条，不再做任何按柱汇总
//              event  事件号（整数，从 0 起）
//              mcTrackID  Geant4 的 MC 径迹编号（初级 = 1，次级 > 1）；
//                         这是真值，不是分类结果——分类后的 trackID 由 class 产生
//              x,y,z  这一步起点的坐标 [mm]
//              x2,y2,z2  这一步终点的坐标 [mm]（一步内部的真实走向）
//              edep   这个步在该柱子内沉积的能量 [MeV]
//              t      这一步起点的全局时间 [ns]
//              t2     这一步终点的全局时间 [ns]
//              fin    这一段的直线路径落在柱阵（= 气体）里的比例，0..1；
//                     正常步为 1，飞出阵列的那一步小于 1
//
//   为什么存两个端点：一步里能量是沿路径连续沉积的，不是全堆在起点。
//   数字化时在两端之间切子段、电离成簇团，每个电子有各自的产生坐标
//   （见 GasDigital::Ionize），才不会让所有电荷都挤在步边界上。
//
//   约定：只有 event 是编号，用整数（"event/I"）；
//         坐标、能量、时间、比例都是物理量，用浮点（"x/D" ...）。
//
//   另外还会写一棵很小的 **events** 树，每个有沉积的事件一行：
//              event  事件号（和 hits 里的一致，含 /output/eventOffset 偏移）
//              n      这个事件在 hits 里占多少行（Long64_t）
//   行数按事件号的先后顺序累加就等于 hits 里的行号，所以下游可以
//   "查表 -> cumsum -> 直接 GetEntry 跳到该事件第一行"，不用扫全文件。
//   存行数而不是存绝对行号，是因为 hadd 合并时不会修改数据：
//   合并后按文件顺序重新 cumsum，偏移依然正确。
#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"

#include "RtypesCore.h"   // Long64_t：events 树里的行数是 64 位计数

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

    // 记一步：p0 / p1 是这一步的起点与终点（Geant4 单位），
    // edep 是沉积能量，t0 / t1 是两端时间，fracInside 是落在气体里的比例
    void FillHit(G4int event, const G4ThreeVector& p0, const G4ThreeVector& p1,
                 G4double fracInside, G4double edep,
                 G4double t0, G4double t1, G4int mcTrackID);

  private:
    // 把上一个事件的行数落进 events 树（事件切换 / 收尾时调用）
    void FlushEvent();

    G4String fFileName = "pixel.root";     // 宏里用 /output/fileName 改
    // 事件号偏移：写进树里的 event = G4EventID + fEventOffset。
    // 并行跑时第 i 个进程填 i * 每进程事件数，hadd 合并后事件号才全局唯一。
    G4int fEventOffset = 0;
    G4GenericMessenger* fMessenger = nullptr;
    TFile* fFile = nullptr;
    TTree* fTree = nullptr;
    TTree* fEvents = nullptr;              // 每个事件一行：事件号 + 该事件的行数
    G4int fCurEvent = -1;                  // 当前正在计入的事件号（未偏移）
    Long64_t fCurCount = 0;                // 当前事件已经写了多少行

    // 写树用的缓冲区
    G4int    bEvent = 0;
    G4int    bMcTrackID = 0;
    G4double bX = 0.0;
    G4double bY = 0.0;
    G4double bZ = 0.0;
    G4double bX2 = 0.0;
    G4double bY2 = 0.0;
    G4double bZ2 = 0.0;
    G4double bEdep = 0.0;
    G4double bT = 0.0;
    G4double bT2 = 0.0;
    G4double bFin = 1.0;

    G4int    bEvtEvent = 0;                // events 树的缓冲区
    Long64_t bEvtN = 0;

    // 编号必须是整数：谁把它改成 double，这里就编译不过
    static_assert(std::is_integral_v<decltype(bEvent)>, "event 编号必须是整型");
    static_assert(std::is_integral_v<decltype(bMcTrackID)>, "mcTrackID 必须是整型");
    static_assert(std::is_integral_v<decltype(fEventOffset)>, "事件号偏移也必须是整型");
    static_assert(std::is_integral_v<decltype(bEvtEvent)>, "events 树的事件号必须是整型");
    static_assert(std::is_integral_v<decltype(bEvtN)>, "事件行数是计数，必须是整型");
};

#endif
