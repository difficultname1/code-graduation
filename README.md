# pixel — 气体像素探测器仿真

一个电子沿 +Z 穿过 100 × 100 根气体像素柱。只记录**落在像素柱内的能量沉积**：
每个有沉积的步一条，内容是这一步起点的坐标 (x, y, z)、这一步沉积的能量，
以及整数事件号。不存柱子编号，多次事件的结果用 `plot.py` 把坐标折回柱号再求平均。

## 几何

| 项 | 值 |
| --- | --- |
| 底面 | 边长 83 µm 的正方形 |
| 柱高 | 1 m，沿 Y 轴（底面法向量 = +Y） |
| 数量 | 100 × 100 根，排在 X–Z 平面，彼此紧贴 |
| 阵列 | 8.3 mm (X) × 1000 mm (Y) × 8.3 mm (Z) |
| 入射面 | 电子沿 +Z 飞来，正对的是 **83 µm (X) × 1 m (Y)** 那个侧面 |
| 电子枪 | 固定于 **XOY 平面中心 (0, 0)**、紧贴入射面 (z = −415 µm)，方向 +Z |
| 柱子材料 | `G4_AIR`；外面是真空（否则电子在阵列外就散了） |

> `kN = 100` 是偶数，x = 0 恰好是 ix=49 / ix=50 两根柱子的交界面，
> 电子沿柱壁飞时“落在哪根”由数值判定决定（实测两边的记录数约各一半）。
> 想让信号稳定落在正中一根柱子上，把枪的 x 改成 `0.5 * kSide`（偏 41.5 µm）。

**要改参数，只改 `src/DetectorConstruction.cc` 文件开头的常量**
（`kSide` / `kHeight` / `kN` / `kMaxStep`），宏文件里不用写任何几何命令。
电子能量在宏里用 `/gun/energy` 改。

头文件（`include/`）里只有类的声明，没有参数。别的类需要几何量时向
`DetectorConstruction` 对象要，比如 `detector->GetGunPos()`。

## 编译 / 运行

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/hep/miniforge/envs/hep
cmake --build build -j8

./build/pixel macros/run.mac          # 固定枪 -> pixel.root
python3 plot.py pixel.root            # -> edep_map*.png / y_map.png / columns_mean.csv

./build/pixel macros/run_random.mac   # 随机枪 -> pixel_random.root
python3 plot.py pixel_random.root     # -> random_*.png / random_columns_mean.csv
```

两个宏互不干扰，各写各的文件；`plot.py` 的输出前缀从 root 文件名自动推出
（`pixel.root` -> 无前缀，`pixel_random.root` -> `random_`），也可以手动指定：
`python3 plot.py a.root myprefix_`。

宏里可用的命令（都在 `/run/initialize` **之后**）：

| 命令 | 作用 |
| --- | --- |
| `/output/fileName pixel_random.root` | 输出文件名（默认 `pixel.root`） |
| `/gun/randomXY true` | 起点在阵列截面上均匀随机；`false`（默认）用固定枪位 |
| `/gun/energy 100 MeV` | 电子能量 |

`/gun/randomXY true` 的随机范围是一切面：x ∈ ±415 µm（阵列全宽）、
y ∈ ±500 mm（柱子全长）。

## 输出 `pixel.root`

一棵树 **`hits`**：柱子里**每个沉积了能量的步**一条，不再按柱子汇总，也不存柱子编号。

| 分支 | 含义 |
| --- | --- |
| `event` | 事件号（整数，0 起） |
| `x`, `y`, `z` | 这一步**起点**的坐标 [mm] |
| `edep` | 这一步在该柱子里沉积的能量 [MeV] |
| `t` | 这一步起点的全局时间 [ns]，数字化时用来挑"最晚到达"的点 |

柱内步长限制是 10 µm，所以一个事件约 900 行（100 根柱子 × 约 8 步），
5000 个事件约 450 万行，`pixel.root` 约 110 MB、`pixel_random.root` 约 92 MB。

**类型约定**：

| 分支 | ROOT 叶子类型 | 分支写法 |
| --- | --- | --- |
| `event` | `Int_t` | `"event/I"` |
| `x`、`y`、`z`、`edep`、`t` | `Double_t` | `"x/D"`、`"y/D"`、`"z/D"`、`"edep/D"`、`"t/D"` |

只有事件号是编号（`static_assert` 守着它必须是整型），坐标和能量都是物理量，浮点。

柱子编号不进文件，统计也直接用物理坐标：`plot.py` 把每根柱子当一个 bin，
bin 宽 `PITCH`、范围 ±`HALF`（`HALF = PITCH * NCOL / 2`），图上的横纵轴就是
**x / z 的 mm 值**，全程不出现整数索引。`PITCH` / `NCOL` 与
`DetectorConstruction.cc` 里的 `kSide` / `kN` 对应，可用第 3、4 个参数覆盖：

```
python3 plot.py pixel.root "" 0.083 100
```

按柱子求平均后写出 `columns_mean.csv`（`x_mm` / `z_mm` 是 bin 中心坐标，
`y_mean_mm` 是柱内所有步的 Y 平均值）：

```
x_mm,z_mm,nStep,edep_sum_MeV,edep_mean_MeV,y_mean_mm
```

## 代码结构

```
main.cc                        物理列表 + 执行宏
digitize.cc                    数字化：hits -> GasDigital -> digi.root
plot_digi.py                   读 digi.root，画芯片热力图（按 83 µm 像素分 bin）
include/, src/                 一个类一个文件
  DetectorConstruction         几何（两级 G4PVReplica 搭 10000 根柱子）
  PrimaryGeneratorAction       电子源
  RunAction                    pixel.root 和 hits 树（一步一行）
  EventAction                  只记当前事件号，把每步原样转给 RunAction
  SteppingAction               逐步判断"在不在柱子里"，交出自己的坐标、edep 和时间
  GasDigital                   数字化核心：edep -> 电子数 -> 芯片上的高斯抽样
  ActionInitialization         把上面几个挂到 Geant4
plot.py                        读 root 文件、按柱(坐标分 bin)求平均、画图/写 CSV
macros/run.mac                 固定枪的运行宏
macros/run_random.mac          随机枪的运行宏
```

## 数字化：hits → digi.root

两种模式（`--mode`）：

| 模式 | 做什么 | 何时用 |
| --- | --- | --- |
| `points`（**默认**） | **每个能量沉积点**都折算电子、各自在芯片上抽样位置（含次级粒子、δ 射线） | 想要真实的电荷分布 |
| `latest` | 每个事件只取到达时间最晚的那个点，电子数用事件总沉积 | 按 `docs/digitize_flow.dot` 的旧流程 |

```bash
./build/digitize pixel.root -o digi.root

# 常用参数
./build/digitize pixel.root -o digi.root \
  --mode points     # 默认；--mode latest 回到"每事件一个点"
  --sigma 0.16      # 芯片位置分辨率 [mm]，默认 160 um
  --backend gauss   # 不用 Garfield，用常数 W = 33.97 eV
  --per-electron    # 每个电子单独一条（否则每个沉积点一条，ne 为电子数）
  --max-events 100  # 只跑前 100 个事件
```

逐点模式里每个沉积点的 edep 只有约 2 eV，`2/33.77 ≈ 0.06` 个电子，所以用泊松抽样：
大多数点抽到 0（不出电子，跳过），少数点抽到 1~2——总量仍然守恒，平均 53.7 个电子/事件。
450 万个沉积点 → 约 23.3 万条有信号的记录。

输出 `digi.root`，一棵树 **`digi`**：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号 |
| `x`, `z` | `Double_t` | 芯片上的落点 [mm]（已加 N(0, σ) 抽样） |
| `ne` | `Int_t` | 这一条对应的电子数（计数，整型） |
| `edep` | `Double_t` | 所在事件的总沉积 [MeV] |
| `t` | `Double_t` | 这个点的到达时间 [ns] |
| `ysrc` | `Double_t` | 这个点的 y [mm]（算漂移距离用） |
| `nstep` | `Int_t` | 所在事件在柱子里共有多少步 |

### Garfield++ 后端（可选）

默认后端是 `garfield`：用 Garfield++ 的 `MediumMagboltz` 建干燥空气，
再用 `TrackHeed` 从 Heed 的物性表取 **W 与 Fano 因子**（实测 W = 33.7744 eV、Fano = 0.19，
对比常数 33.97 eV / 0.2）。**不做漂移、不做扩散**——`sigma = 160 µm` 是读出分辨率，
跟气体输运无关，所以这部分仍然是纯高斯抽样。

编译时需要打开开关并给出 Garfield++ 安装前缀：

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/hep/miniforge/envs/hep \
      -DUSE_GARFIELD=ON -DGARFIELD_ROOT=/home/hep/garfieldpp-install
cmake --build build -j8
```

找不到 Garfield++ 时 configure 只报 warning，`digitize` 自动退回自带的 `gauss` 后端，
所以不装 Garfield 也能跑（`-DUSE_GARFIELD=OFF`，默认就是 OFF）。

### 看一眼芯片上的落点分布

```bash
python3 plot_digi.py                     # -> chip_heatmap.png（左：整块芯片，右：束斑缩放）
python3 plot_digi.py --weight -o ne.png  # 按电子数 ne 加权
python3 plot_digi.py --zoom-half 0.3     # 缩放窗口改成 ±0.3 mm
```

一个柱截面就是一个像素（83 µm），所以热力图不做任何插值，直接就是像素图；
前两幅是 log 色阶，第三幅是 $x$ 方向的一维投影 + 高斯拟合（直接给出分辨率）。
逐点模式下束斑是沿 $z$ 拉长的一条带（粒子整个穿过阵列都在沉积能量），
所以左图看"打在哪"、中图看 $x$ 向剖面、右图把 $\sigma_x$ 拟出来（实测 161.5 µm，输入 160 µm）。

## 几个容易踩的点

- **步长限制必须配 `G4StepLimiterPhysics`**（`main.cc` 已注册）。
  否则电子在 1 m 空气里一步能跨几百 mm，柱子内只能落一两个点，Y 值就没意义了。
  现在 10 µm，每根 83 µm 的柱子约 8 步，穿过整个阵列约 900 步。
- `/gun/...` 命令要写在 `/run/initialize` 之后，否则命令还不存在。
- 每次 `/run/beamOn` 会**重新生成** `pixel.root`（不是追加）。
- 量级自查：空气 MIP ≈ 0.22 keV/mm，8.3 mm 阵列 → 每事件约 1.8 keV；
  10 µm 步长 → 每步约 2 eV，每根柱子约 16 eV。
- 阵列在束流方向只有 8.3 mm，电子基本笔直穿过同一行 100 根柱子，
  所以每个事件约 900 步、落在同一行的 100 根柱子里；`y` 几乎不变（µm 级抖动）。
- `hits` 一步一行。用 `plot.py` 之类的方式做统计时别在
  Python 里逐条循环，交给 ROOT 的 `TTree::Draw` / `RDataFrame`。
