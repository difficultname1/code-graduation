# geant4

Geant4 仿真：一个电子沿 +Z 穿过 100 × 100 根气体像素柱，把柱内每个沉积步写成
`pixel.root` 的 `hits` 树（外加一棵 `events` 索引树）。

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
`kMaxStep`）；宏里只放能量、事件数和输出名。
