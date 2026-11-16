// -----------------------------------------------------------------------------
//  100 x 100 气体像素探测器，一个电子沿 +Z 穿过，结果写到 pixel.root
//
//  编译: cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/hep/miniforge/envs/hep
//        cmake --build build -j8
//  运行: ./build/pixel macros/run.mac
// -----------------------------------------------------------------------------
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv)
{
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  // 1) 几何（参数在 include/DetectorConstruction.hh 顶部的 namespace geo 里）
  auto* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  // 2) 物理列表：标准电磁 + 强子。
  //    G4StepLimiterPhysics 让柱子的步长限制真生效，别去掉。
  G4PhysListFactory factory;
  G4VModularPhysicsList* physics = factory.GetReferencePhysList("FTFP_BERT");
  physics->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physics);

  // 3) 用户动作（必须在物理列表之后创建，否则 Geant4 报 Run0041）
  runManager->SetUserInitialization(new ActionInitialization(detector));

  // 4) 执行宏
  const char* macro = (argc > 1) ? argv[1] : "macros/run.mac";
  G4UImanager::GetUIpointer()->ApplyCommand(G4String("/control/execute ") + macro);

  delete runManager;
  return 0;
}
