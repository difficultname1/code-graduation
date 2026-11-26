# gasdigital

气体数字化：`hits`（Geant4 沉积步）-> 电离 -> 漂移扩散 -> 芯片落点，写 `digi.root`。

```
include/, src/   GasDigital 类：步内切子段、抽电子数（泊松/Fano）、电离簇团、逐电子落点
app/digitize.cc  命令行入口（读 hits，按事件切组，写 digi 树）
scripts/         digitize_parallel.sh：按事件号分段并行数字化 + hadd 合并
```

```bash
./build/digitize pixel.root -o digi.root
./gasdigital/scripts/digitize_parallel.sh pixel_height.root digi_height.root 8 --chip-y -40
```

可选的 Garfield++ 后端（`-DUSE_GARFIELD=ON`）只用来取 Heed 的 W 与 Fano 因子；
气体输运用的是解析模型 $\sigma_d = D\sqrt{L}$。
