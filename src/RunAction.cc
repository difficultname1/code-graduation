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

  // 一步一条：只记坐标、沉积能量和事件号，不做按柱编号/汇总
  fTree = new TTree("hits", "每个有沉积的步一条");
  fTree->Branch("event", &bEvent, "event/I");
  fTree->Branch("x",     &bX,     "x/D");      // mm
  fTree->Branch("y",     &bY,     "y/D");      // mm
  fTree->Branch("z",     &bZ,     "z/D");      // mm
  fTree->Branch("edep",  &bEdep,  "edep/D");   // MeV
  fTree->Branch("t",     &bT,     "t/D");      // ns
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  // 先取条数，Close() 之后树就没了
  const Long64_t nHits = fTree->GetEntries();

  fFile->cd();
  fTree->Write();
  fFile->Close();
  fFile = nullptr;
  fTree = nullptr;

  G4cout << "[输出] " << fFileName << " : " << nHits << " 条记录, "
         << run->GetNumberOfEvent() << " 个事件" << G4endl;
}

void RunAction::FillHit(G4int event, G4double x, G4double y, G4double z,
                        G4double edep, G4double t)
{
  if (!fTree) return;
  bEvent = event;
  bX = x / mm;
  bY = y / mm;
  bZ = z / mm;
  bEdep = edep / MeV;
  bT = t / ns;
  fTree->Fill();
}
