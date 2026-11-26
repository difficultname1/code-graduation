#include "RunAction.hh"

#include "G4GenericMessenger.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "TFile.h"
#include "TTree.h"

RunAction::RunAction()
{
  fMessenger = new G4GenericMessenger(this, "/output/", "ROOT 输出");
  fMessenger->DeclareProperty("fileName", fFileName,
      "输出文件名（默认 pixel.root）。注意每次 beamOn 都会重建该文件。");
  fMessenger->DeclareProperty("eventOffset", fEventOffset,
      "事件号偏移（默认 0）。并行跑时第 i 个进程填 i*每进程事件数，"
      "hadd 合并后事件号才全局唯一。");
}

RunAction::~RunAction()
{
  delete fMessenger;
}

// Geant4 默认带 ROOT 分析后端（G4AnalysisManager）时可以直接用，
// 本机这个 build 没有，所以用 ROOT 自己的 TFile / TTree。
void RunAction::BeginOfRunAction(const G4Run*)
{
  fFile = new TFile(fFileName, "RECREATE");

  // 一步一条：记这一步的两端坐标、沉积能量、两端时间和气体占比。
  // 不存柱子编号，也不在这边做"摊到起点"的简化——插值留给 digitize。
  fTree = new TTree("hits", "每个有沉积的步一条");
  fTree->Branch("event", &bEvent, "event/I");
  fTree->Branch("x",     &bX,     "x/D");      // 起点 [mm]
  fTree->Branch("y",     &bY,     "y/D");
  fTree->Branch("z",     &bZ,     "z/D");
  fTree->Branch("x2",    &bX2,    "x2/D");     // 终点 [mm]
  fTree->Branch("y2",    &bY2,    "y2/D");
  fTree->Branch("z2",    &bZ2,    "z2/D");
  fTree->Branch("edep",  &bEdep,  "edep/D");   // MeV
  fTree->Branch("t",     &bT,     "t/D");      // 起点时间 [ns]
  fTree->Branch("t2",    &bT2,    "t2/D");     // 终点时间 [ns]
  fTree->Branch("fin",   &bFin,   "fin/D");    // 落在气体里的比例 0..1

  // 附带一棵很小的索引树：每个事件一行（事件号 + 该事件的行数）。
  // 下游按 cumsum 就能直接跳到某个事件的第一行，不用扫整个文件。
  // 存"行数"不存"绝对行号"：hadd 不改数据，合并后重新 cumsum 就对了。
  fEvents = new TTree("events", "每个事件一行：事件号 + 该事件在 hits 里的行数");
  fEvents->Branch("event", &bEvtEvent, "event/I");
  fEvents->Branch("n",     &bEvtN,     "n/L");

  fCurEvent = -1;
  fCurCount = 0;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  // 先把最后一个事件的计数落进去，再取条数（Close() 之后树就没了）
  FlushEvent();
  const Long64_t nHits = fTree->GetEntries();
  const Long64_t nEventsWithHits = fEvents->GetEntries();

  fFile->cd();
  fTree->Write();
  fEvents->Write();
  fFile->Close();
  fFile = nullptr;
  fTree = nullptr;
  fEvents = nullptr;

  G4cout << "[输出] " << fFileName << " : " << nHits << " 条记录, "
         << run->GetNumberOfEvent() << " 个事件";
  if (fEventOffset != 0)
    G4cout << "（事件号 " << fEventOffset << ".."
           << fEventOffset + run->GetNumberOfEvent() - 1 << "）";
  G4cout << "，其中 " << nEventsWithHits
         << " 个有沉积（events 树里 " << nEventsWithHits << " 行）" << G4endl;
}

void RunAction::FlushEvent()
{
  if (!fEvents || fCurEvent < 0) return;
  bEvtEvent = fCurEvent + fEventOffset;   // 和 hits 里的 event 保持同一个数
  bEvtN = fCurCount;
  fEvents->Fill();
  fCurCount = 0;
}

void RunAction::FillHit(G4int event, const G4ThreeVector& p0,
                        const G4ThreeVector& p1, G4double fracInside,
                        G4double edep, G4double t0, G4double t1)
{
  if (!fTree) return;
  // 事件切换：先把上一个事件的行数落进 events 树
  if (event != fCurEvent) {
    FlushEvent();
    fCurEvent = event;
  }
  ++fCurCount;

  // 事件号 = 本次运行的编号 + 偏移；类型仍是 G4int，行内物理量不动
  bEvent = event + fEventOffset;
  bX  = p0.x() / mm;
  bY  = p0.y() / mm;
  bZ  = p0.z() / mm;
  bX2 = p1.x() / mm;
  bY2 = p1.y() / mm;
  bZ2 = p1.z() / mm;
  bEdep = edep / MeV;
  bT  = t0 / ns;
  bT2 = t1 / ns;
  bFin = fracInside;
  fTree->Fill();
}
