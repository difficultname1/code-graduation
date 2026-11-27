# geant4

Geant4 仿真：一个电子沿 X–Z 平面的 45° 对角线穿过 1000 × 1000 根气体像素柱
（整体 83 × 40 × 83 mm，漂移高度 40 mm），把柱内每个沉积步连同 MC
真值 `mcTrackID` 写成 `pixel.root` 的 `hits` 树（外加一棵 `events` 索引树）。

```
include/          类的声明（DetectorConstruction、RunAction、SteppingAction ...）
src/              类的实现；main.cc 是入口（物理列表 + 执行宏）
macros/           run.mac / run_random.mac / run_height.mac
scripts/          run_parallel.sh：多进程并行跑 + hadd 合并
```

```bash
./build/pixel geant4/macros/run.mac
./geant4/scripts/run_parallel.sh geant4/macros/run_height.mac 20000 15 pixel_height.root
```

几何参数在 `src/DetectorConstruction.cc` 顶部的常量里（`kSide` / `kHeight` / `kN` /
`kMaxStep` / `kGunHeight` / `kGunAngle`）；宏里只放能量、事件数和输出名。
