# pixel — 气体像素探测器仿真与径迹/顶点重建

## 物理目标

重建 **真缪子偶素（true muonium，μ⁺μ⁻ 束缚态）衰变产生的正负电子对**：

- 真缪子偶素沿 +Z 方向前进，随后衰变为 e⁺e⁻；
- 正负电子穿过 83 × 83 × 40 mm 的气体像素探测器，各留下一条电离径迹；
- 目标是从芯片上的电子沉积出发，**把两条径迹分开并分别拟合**，
  再用两条径迹的公共交点重建 **衰变顶点**（衰变发生在探测器外时，就外推求交点）。

最终分析链：

```text
Geant4 (μ⁺μ⁻ → e⁺e⁻)
  → hits（气体沉积步）
  → electrons（芯片电子落点，只有 MC 真值 mcTrackID）
  → class（分到两个 trackID：e⁺ / e⁻）
  → track fit（每条径迹一条直线）
  → vertex fit（两条直线的公共顶点）
```

## 当前阶段

目前先用 **单个 100 MeV 电子沿 X–Z 对角线 45° 入射** 验证整条流程：
几何、数字化、命名、单径迹拟合已经打通。下一步把源换成 μ⁺μ⁻ → e⁺e⁻，
让一个事件里同时出现两条径迹，`class` 分成两个 `trackID`，`vertex fit`
再求两条直线的交点。

> 单径迹测试的目的是把每一步的数据结构和误差模型定死。最终的顶点分辨率
> 由两条径迹的夹角、探测器有效长度（40 mm 漂移）和单点分辨率
> （80 µm 电离倍增小孔 + 漂移扩散）共同决定。

## 仿真约定

`hits` 每个"在像素柱里沉积了能量"的步一条，记录这一步的起点和终点
`(x,y,z)` / `(x2,y2,z2)`、沉积能量 `edep`、事件号 `event` 与 MC 真值 `mcTrackID`。
一步不是"一个点"：粒子沿着这一步连续走、一路电离。存两个端点，数字化时
在两端之间插值（见下文 `--substep`），电荷才不会全堆在 10 µm 的步边界上。

## 探测器几何（当前单径迹测试配置）

| 项 | 值 |
| --- | --- |
| 世界 | 边长 **120 mm** 的真空立方体，中心在原点 |
| 底面 | 边长 830 µm 的正方形 |
| 柱高 | **40 mm**，沿 Y 轴（底面法向量 = +Y）= 最大漂移长度 |
| 数量 | 100 × 100 根，排在 X–Z 平面，彼此紧贴 |
| 阵列 | **83 mm (X) × 40 mm (Y) × 83 mm (Z)**，居中于立方体 |
| 芯片 | 贴在柱阵底面 **y = −20 mm**；电离点到芯片的漂移距离 y + 20 mm ∈ [0, 40] mm |
| 入射面（当前测试） | 单个电子从左下角 (−41.5, −41.5) 斜 45° 射到右上角 (+41.5, +41.5) |
| 电子枪（当前测试） | 固定于 **y = 14 mm**、左下角外 1 mm；方向 (sin45°, 0, cos45°) |
| 柱子材料 | `G4_AIR`；世界是真空（否则电子在阵列外就散了） |

> 45° 斜穿时不会沿着柱壁长距离飞行；左下角 (−41.5, −41.5) 是四根柱子的
> 共点边界，只有进入阵列的第一小段贴着边界，对整条径迹的统计影响可以忽略。
> 换成 μ⁺μ⁻ → e⁺e⁻ 后，两条径迹会从同一个顶点出发，几何本身不变。

**要改参数，只改 `geant4/src/DetectorConstruction.cc` 文件开头的常量**
（`kSide` / `kHeight` / `kN` / `kMaxStep` / `kGunHeight` / `kGunAngle`），
宏文件里不用写任何几何命令。电子能量在宏里用 `/gun/energy` 改。

头文件（`geant4/include/`）里只有类的声明，没有参数。别的类需要几何量时向
`DetectorConstruction` 对象要，比如 `detector->GetGunPos()`。

## 编译 / 运行

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/hep/miniforge/envs/hep
cmake --build build -j8

./build/pixel geant4/macros/run.mac            # 固定枪 -> pixel.root
python3 "picture draw/app/plot.py" pixel.root      # -> edep_map*.png / y_map.png / columns_mean.csv

./build/pixel geant4/macros/run_random.mac     # 随机枪 -> pixel_random.root
python3 "picture draw/app/plot.py" pixel_random.root
```

### 完整分析链

```bash
# 当前单径迹测试链；最终把第 1 步换成 μ⁺μ⁻ → e⁺e⁻ 的事件源
./build/pixel geant4/macros/run.mac                          # 1. Geant4 -> pixel.root
./build/digitize pixel.root -o digi.root                     # 2. 数字化 -> electrons
python3 "class/app/run_class.py" digi.root -o points.root    # 3. 分类 -> trackID
python3 "track fit/app/run_track_fit.py" points.root -o tracks.root   # 4. 直线拟合
python3 "vertex fit/app/run_vertex_fit.py" digi.root -o vertex.root   # 5. 顶点（目前是骨架）
```

`points.root` 的 `points` 树是每条径迹的落点（`event / trackID / x / z / ne / ysrc`）；
`tracks.root` 的 `tracks` 树是每条径迹的直线参数（`A / B`，满足 `x = A*z + B`），
外加 `sigmaA / sigmaB / covAB / corrAB`（误差传递）、`chi2 / ndf / chi2ndf / pvalue`
（卡方评估）和实际扣掉的本底 `bkg`。拟合按电子数加权，默认扣掉该径迹的平均
电荷本底（`--background mean`；每电子一条时自动退回不扣）；单点误差取
80 µm 小孔宽度与漂移扩散的平方和。
当前 `class` 假设分类完全正确，直接把 digitize 里的 MC 真值 `mcTrackID` 当作分类结果。
`trackID` 是 class 之后才出现的量，上游的 `hits` / `electrons` 里没有它。

### 最终目标：两条径迹 + 顶点

单电子测试跑通后，物理事件变成 μ⁺μ⁻ → e⁺e⁻：

| 阶段 | 单电子测试（现在） | μ⁺μ⁻ → e⁺e⁻（目标） |
| --- | --- | --- |
| Geant4 | 一个 100 MeV e⁻，45° 入射 | 一个真缪子偶素衰变出 e⁺e⁻ |
| `electrons` | 一条径迹的电子 | 两条径迹的电子混在一起 |
| `class` | 用 `mcTrackID` 做完美分类 | 分成 `trackID = 1`（e⁻）/ `2`（e⁺） |
| `track fit` | 一条直线 `x = A*z + B` | 两条直线，各自的 A / B / 误差 |
| `vertex fit` | 不用 | 两条直线求交，得到衰变顶点 |

当前 `track fit` 先做芯片平面 X–Z 的投影直线；顶点重建需要两条三维直线，
下一步补上 y(z)（产生高度 `ysrc` 对 z 的依赖）或者直接拿 3D 点拟合。
两条三维直线在有限分辨率下不一定严格相交，`vertex fit` 实际求它们的
最近接近点（或做两条直线的联合最小二乘），这个点就是衰变顶点。
`vertex fit/` 目前还是骨架，等两条径迹能稳定分开后再填。

两个宏互不干扰，各写各的文件；`plot.py` 的输出前缀从 root 文件名自动推出
（`pixel.root` -> 无前缀，`pixel_random.root` -> `random_`），也可以手动指定：
`python3 "picture draw/app/plot.py" a.root myprefix_`。

宏里可用的命令（都在 `/run/initialize` **之后**）：

| 命令 | 作用 |
| --- | --- |
| `/output/fileName pixel_random.root` | 输出文件名（默认 `pixel.root`） |
| `/output/eventOffset 5000` | 写进树里的 `event` = G4 事件号 + 这个偏移（默认 0）；并行跑时每个进程填自己的区间起点 |
| `/gun/randomXY true` | 起点在阵列截面上均匀随机；`false`（默认）用固定枪位 |
| `/gun/randomY true` | x 固定在枪位、只把 y 沿柱高随机（扫漂移高度用） |
| `/gun/energy 100 MeV` | 电子能量 |

### 并行跑（多进程 + `hadd`）

本机的 Geant4 是 serial build（`geant4-config --has-feature multithreaded` → `no`），
用不了 `G4MTRunManager`；并行靠**多进程**：每个进程写自己的文件，最后用 ROOT 的
`hadd` 拼起来。分支名、分支类型、每行的物理量都不变，`hadd` 只是把各文件的
TTree 桶按顺序串接（实测合并后仍是 `event/I`、`x/D` …）。

```bash
./geant4/scripts/run_parallel.sh geant4/macros/run_height.mac 20000 15 pixel_height.root
./build/digitize pixel_height.root -o digi_height.root --chip-y -20
```

数字化也能并行（按事件号切，靠 `events` 树直接跳）：

```bash
./gasdigital/scripts/digitize_parallel.sh pixel_height.root digi_height.root 8 --chip-y -20
```

每个进程跑 `digitize --part i/N`，种子用 `SEED+i` 错开，最后 `hadd`。`electrons` 的
`event` 来自 `hits`，本来就是全局唯一的，所以合并后不需要再偏移。

脚本按 worker 编号分配三样东西：独立输出文件 `.parallel/part_00i.root`、
独立种子 `(s1+i, s2+i)`、独立事件号区间（宏里写 `/output/eventOffset i*M`）。
跑完自动 `hadd` 并按 worker 顺序拼接，事件号因此整体递增。

**为什么必须加偏移**：事件号是 `G4RunManager` 的**进程内**计数器，每个进程都从
0 开始数。不偏移的话，合并后的编号序列是 `0..99, 0..99`：`hits` 里 200 个物理上
不同的事件只占 100 个号码。后果是 `plot.py` 的 `nEvents = GetMaximum("event") + 1`
少算一半（每事件平均能量翻倍），`plot_digi.py` 用 `bincount` 按事件号分组算
$x - x_{\text{event}}$，会把两个不同事件的入射 $x$ 混在一起，束斑宽度被人为撑宽。
加了偏移后事件号全局唯一，且 `digitize` 要求的"事件号连续递增"仍然成立
（实测合并后回退次数为 0）。

量级参考：16 核机器上 `run.mac` 单进程约 190 事件/秒（其中约 1 s 是几何/物理初始化），
`run_height.mac` 的 20000 事件串行约 1.7 min，分给 15 个进程后每个约 8 s。
中间文件总量随事件数线性增长（5000 事件约 200 MB），`KEEP=1` 可保留中间文件。

`/gun/randomXY true` 的随机范围是入射面：x ∈ ±41.5 mm（阵列全宽）、
y ∈ ±20 mm（柱阵半高）。

固定枪的电荷几乎都产生在同一个高度上（见「容易踩的点」），要验证漂移扩散随高度的变化，用第三个宏：

```bash
./build/pixel geant4/macros/run_height.mac   # 笔形束扫高度 -> pixel_height.root（末义 20k 事件）
./build/digitize pixel_height.root -o digi_height.root --chip-y -20
python3 "picture draw/app/plot_digi.py" -i digi_height.root -o chip_heatmap_height.png --chip-y -20
```

## 输出 `pixel.root`

一棵树 **`hits`**：柱子里**每个沉积了能量的步**一条，不再按柱子汇总，也不存柱子编号。

| 分支 | 含义 |
| --- | --- |
| `event` | 事件号（整数，0 起） |
| `mcTrackID` | MC 真值径迹号（整数；初级 = 1，次级 > 1），**不是分类结果** |
| `x`, `y`, `z` | 这一步**起点**的坐标 [mm] |
| `x2`, `y2`, `z2` | 这一步**终点**的坐标 [mm] |
| `edep` | 这一步在该柱子里沉积的能量 [MeV] |
| `t` | 这一步起点的全局时间 [ns] |
| `t2` | 这一步终点的全局时间 [ns] |
| `fin` | 这一步的直线路径落在柱阵（= 气体）里的比例，0..1。正常步是 1；
最后飞出阵列那一步小于 1，用它把 edep 只摊到气体里的那一段上 |

柱内步长限制是 10 µm，45° 斜穿阵列的路径约 117 mm，所以一个事件约 1.2 万行，
5000 个事件约 6000 万行，`pixel.root` 约 2 GB（跑大样本要留意磁盘）。
（多了两端的坐标/时间分支，比只存起点大一圈。）

**类型约定**：

| 分支 | ROOT 叶子类型 | 分支写法 |
| --- | --- | --- |
| `event`、`mcTrackID` | `Int_t` | `"event/I"`、`"mcTrackID/I"` |
| `x`、`y`、`z`、`x2`、`y2`、`z2`、`edep`、`t`、`t2`、`fin` | `Double_t` | `"x/D"` … `"fin/D"` |

同一个文件里还有一棵很小的 **`events`** 树，每个有沉积的事件一行：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号（与 `hits` 里的一致，含 `/output/eventOffset` 偏移） |
| `n` | `Long64_t` | 这个事件在 `hits` 里占多少行 |

把 `n` 按顺序累加就是该事件在 `hits` 里的起始行号，所以下游要处理某个事件时可以
**直接跳过去，不用扫描整个文件**。存"行数"而不存"绝对行号"是为 `hadd` 准备的：
合并不会修改任何数据，按文件顺序重新 cumsum 依然正确。5000 事件时这棵树只有
5000 行、几十 KB。

> 为什么不用 `TTree::BuildIndex`：实测 `hadd` **不会**把 `TTreeIndex` 带到合并文件里
> （合并后 `GetTreeIndex()` 返回 0），而且索引要求键唯一，一个事件 900 多行会不停
> 报 duplicate key 警告。自己写一棵索引树两件事都解决了。

只有事件号是编号（`static_assert` 守着它必须是整型），坐标和能量都是物理量，浮点。

柱子编号不进文件，统计也直接用物理坐标：`plot.py` 把每根柱子当一个 bin，
bin 宽 `PITCH`、范围 ±`HALF`（`HALF = PITCH * NCOL / 2`），图上的横纵轴就是
**x / z 的 mm 值**，全程不出现整数索引。`PITCH` / `NCOL` 与
`DetectorConstruction.cc` 里的 `kSide` / `kN` 对应，可用第 3、4 个参数覆盖：

```
python3 "picture draw/app/plot.py" pixel.root "" 0.083 100
```

按柱子求平均后写出 `columns_mean.csv`（`x_mm` / `z_mm` 是 bin 中心坐标，
`y_mean_mm` 是柱内所有步的 Y 平均值）：

```
x_mm,z_mm,nStep,edep_sum_MeV,edep_mean_MeV,y_mean_mm
```

## 代码结构

```
CMakeLists.txt                 两个可执行文件的构建（pixel / digitize）
geant4/                        Geant4 仿真：电子沿 45° 对角线穿过 100 × 100 根气体像素柱
  include/                     类声明
  src/                         类实现；main.cc 是入口
  macros/                      run.mac / run_random.mac / run_height.mac
  scripts/run_parallel.sh      多进程并行跑 + hadd（给每个进程分配事件号偏移）
gasdigital/                    气体数字化：hits -> GasDigital -> digi.root
  include/、src/               GasDigital 类的声明与实现
  app/digitize.cc              命令行入口
  scripts/digitize_parallel.sh 按事件号分段并行数字化 + hadd
class/                         分类：electrons -> 带 trackID 的 points（当前用真值跳过分类器）
  app/run_class.py             命令行入口
  classifier/                  包：io.py 读写、classify.py 分类钩子
track fit/                     径迹拟合：class 的 points -> tracks
  app/run_track_fit.py         命令行入口
  track_fit/                   包：io.py 读写、fit.py 加权最小二乘直线
vertex fit/                    入射位置重建：digi -> vertices
  app/run_vertex_fit.py        命令行入口
  vertex_fit/                  包：io.py 读写、fit.py 拟合
picture draw/                  画图
  app/plot.py                  读 hits，按柱求平均、画图、写 CSV
  app/plot_digi.py             读 digi，画芯片热力图
docs/                          数字化流程图
```

约定：**代码只放在分类目录的子目录里**（`app/`、`include/`、`src/`、包名），
分类根目录只放 `README.md`。

## 数字化：hits → digi.root

两种模式（`--mode`）：

| 模式 | 做什么 | 何时用 |
| --- | --- | --- |
| `points`（**默认**） | **每个能量沉积步**先按 `--substep` 切成子段；每子段的 edep 抽电子数（泊松/Fano 涨落）、分成**电离簇团**、每个电子有自己的产生坐标，再逐个在芯片上抽样位置 | 想要真实的电荷分布 |
| `latest` | 每个事件只取到达时间最晚的那个点，电子数用事件总沉积 | 按 `docs/digitize_flow.dot` 的旧流程 |

只做一部分事件（并行、调试都用得上），不用从头扫描文件：

| 选项 | 作用 |
| --- | --- |
| `--event 123` | 只处理事件 123 |
| `--event-start 1000 --event-count 500` | 处理事件 1000..1499 |
| `--part 3/16` | 事件号按 16 等分，只做第 3 份（i 从 0 起） |

定位优先用 `hits` 里附带的 `events` 树（前缀和 -> 行号，一次 `GetEntry` 就到位）；
旧文件没有这棵树时退回二分查找，复杂度仍是 $O(\log N)$，同样不遍历全表。
实测同一事件在"全量模式"和"`--event` 模式"下逐行完全一致（`|Δx| = 0`）。

电子数不是 `edep/W` 的确定值：期望值小于 10 个时用泊松，大量电离时用
Fano 修正的高斯（`Var(N) = F·N`，空气实测 F = 0.19）。拿到电子数后再把它们
分成簇团（默认均 2.5 电子/簇），簇心沿径迹随机、簇内产生点按 `--cluster-sigma`
散开——所以同一段沉积出来的电子不会全叠在一个点上。

```bash
./build/digitize pixel.root -o digi.root

# 常用参数
./build/digitize pixel.root -o digi.root \
  --mode points     # 默认；--mode latest 回到"每事件一个点"
  --sigma 0.08      # 电离倍增小孔的定位宽度 [mm]，默认 80 um（芯片本身不加 sigma）
  --diffusion 85    # 漂移横向扩散系数 [um/sqrt(cm)]，默认 85；0 = 关掉
  --chip-y -20      # 芯片平面 y [mm]，默认 -20 = 40 mm 柱阵的底面
  --substep 0.002   # 步内插值的子段长度 [mm]，默认 2 um；0 = 不插值
  --cluster-mean 2.5    # 一个电离簇团平均几个电子（<= 1 表示一电子一簇）
  --cluster-sigma 5e-5  # 簇内电子产生点的空间涨落 [mm]，默认 50 nm
  --backend gauss   # 不用 Garfield，用常数 W = 33.97 eV
  --per-electron    # 每个电子一条（默认）：位置就是它自己的产生坐标
  --per-cluster     # 每个电离簇团一条，ne = 簇内电子数、位置 = 簇心
  --max-events 100  # 只跑前 100 个事件
```

**漂移扩散随高度变化**：电子产生在 $y$ 处，沿柱轴漂到芯片（$y$ = `--chip-y`），
漂移长度 $L = |y - y_{\text{chip}}|$。横向扩散

$$\sigma_d(L) = D\,\sqrt{L},\qquad D = 85\ \mu\text{m}/\sqrt{\text{cm}}$$

和读出分辨率按平方和合成：$\sigma = \sqrt{\sigma_{\text{read}}^2 + \sigma_d^2}$。
柱高 40 mm 时 $L$ 从 0 到 40 mm，所以 $\sigma$ 从 80 µm 长到约 190 µm ——
**同一事件里产生得越高的电荷，落到芯片上越模糊**。纵向扩散沿漂移方向，
不影响芯片上的 $(x,z)$，这里没模拟。

**步内插值（`--substep`）**：Geant4 的一步是"从起点走到终点、一路上一直在电离"，
不是一个点上的瞬时沉积。数字化时把这一步的两端 `(x,y,z,t)` → `(x2,y2,z2,t2)`
按 `--substep` 切成若干子段，`edep` 平分给各段。
最后一步飞出阵列时按 `fin` 裁掉跑到真空里的那一截。

不插值时（`--substep 0`）所有电荷都堆在步起点上，固定枪 + 恒定步长下能直接在
细 bin 直方图上看到一条 10 µm 间隔的"梳齿"——那是步长格点，不是物理。
每步 10 µm、子段 2 µm 时，一个步变成 5 个子段。

逐点模式里每个子段的 edep 只有约 0.4 eV（一步 2 eV / 5 段），`0.4/33.77 ≈ 0.01` 个电子，
所以用泊松抽样：绝大多数子段抽到 0（不出电子，跳过），少数抽到 1——总量仍然守恒，
平均约 760 个电子/事件。一个事件约 1.2 万个沉积步 → 约 6 万个子段 → 约 760 条
电子记录（20 事件实测：25.6 万步 → 120 万子段 → 1.5 万电子）。

输出 `digi.root`，一棵树 **`electrons`**（芯片上的电子落点，不叫 hits）：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号 |
| `mcTrackID` | `Int_t` | MC 真值径迹号（初级 = 1），class 阶段用它做完美分类的捷径 |
| `x`, `z` | `Double_t` | 芯片上的落点 [mm]（已叠加 80 µm 小孔宽度与漂移扩散） |
| `ne` | `Int_t` | 这一条对应的电子数（计数，整型；每电子一条时为 1） |
| `edep` | `Double_t` | 所在事件的总沉积 [MeV] |
| `t` | `Double_t` | 这个电子（或簇团）的产生时间 [ns] |
| `ysrc` | `Double_t` | 这个电子（或簇团）产生点的 y [mm]（算漂移距离用） |
| `nstep` | `Int_t` | 所在事件在柱子里共有多少步 |

### Garfield++ 后端（可选）

默认后端是 `garfield`：用 Garfield++ 的 `MediumMagboltz` 建干燥空气，
再用 `TrackHeed` 从 Heed 的物性表取 **W 与 Fano 因子**（实测 W = 33.7744 eV、Fano = 0.19，
对比常数 33.97 eV / 0.2）。气体输运（漂移扩散）**不用 Garfield 的输运模块**，
而是用上面那个解析模型 $\sigma_d = D\sqrt{L}$，只依赖一个系数 $D$。

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
python3 "picture draw/app/plot_digi.py"                    # -> chip_heatmap.png（左：整块芯片，右：束斑缩放）
python3 "picture draw/app/plot_digi.py" --weight -o ne.png # 按电子数 ne 加权（每电子一条时 ne 恒为 1，
                                                            # 这条只在 --per-cluster 模式下有意义）
python3 "picture draw/app/plot_digi.py" --zoom-half 0.3    # 缩放窗口改成 ±0.3 mm
```

一个柱截面就是一个像素（830 µm），所以热力图不做任何插值，直接就是像素图；
前两幅是 log 色阶，第三幅是"扣掉该事件自身入射 $x$ 之后"的 $x$ 剖面 + 高斯拟合
（随机枪也能得到有意义的宽度）。
另外会写一张 `<输出名>_vs_drift.png`（四联）：左上 dx 随漂移长度 $L$ 的楔形展开
（虚线是模型 $\pm\sqrt{\sigma_{\text{read}}^2+(D\sqrt{L})^2}$），右上每个 bin
$\langle x-x_{\text{event}}\rangle \pm 1\sigma$（点是平均落点，误差棒上下各 1σ，棒尖正好落在模型曲线上；
同时拟合出 $\sigma_0$ 与 $D$），下排是低/高漂移长度两块芯片分布（共用同一色阶，直接比宽窄）。
一维剖面的误差棒是每个 bin 的 Poisson $\pm\sqrt{N}$。
固定枪的电荷都集中在同一个高度上，覆盖不到不同 $L$，这时这张图会跳过。
同时还会写一张 `<输出名>_log.png`：左边 $\sigma_x(L)$（log 纵轴 + 输入/拟合两条模型），
右边低 $L$ 与高 $L$ 的 $x$ 剖面叠在一起（log 纵轴，尾部差一个量级，宽度差别在这里最明显）。
逐点模式下束斑是沿 $z$ 拉长的一条带（粒子整个穿过阵列都在沉积能量），
所以左图看"打在哪"、中图看 $x$ 向剖面、右图把 $\sigma_x$ 拟出来。

## 几个容易踩的点

- **固定枪只在一个高度上**：y = 14 mm，整条径迹的漂移长度都是 34 mm，
  测不到扩散随高度的变化。要看高度依赖用 `run_height.mac`（x 固定、y 均匀）。
- **步长限制必须配 `G4StepLimiterPhysics`**（`main.cc` 已注册）。
  否则电子在 83 mm 的空气里几步就穿出去了，柱子内只能落一两个点。
  现在 10 µm，45° 斜穿阵列的路径约 117 mm，一个事件约 1.2 万步。
- **一步不是"一个点"**：步的 edep 沿整段路径连续释放，不是堆在起点。
  `hits` 存了两个端点，`digitize --substep` 在中间插值；不插值就会看到
  10 µm 间隔的位置格点（`--substep 2 um` 时变成 2 µm，且段内随机后彻底消失）。
  能量总数、电子总数不受影响，变的是**电荷在路径上的分布**。
- `/gun/...` 命令要写在 `/run/initialize` 之后，否则命令还不存在。
- 每次 `/run/beamOn` 会**重新生成** `pixel.root`（不是追加）。
- **事件号是进程内计数器**：每个进程都从 0 开始。并行跑必须用
  `/output/eventOffset` 错开区间，否则 `hadd` 之后重号，按事件分组的统计会串。
- 按事件切片（`digitize --part`）时每个进程的随机数起点不同，结果与串行**不逐位相同**
  （电子数的泊松/Fano 涨落本就随机），但统计等价；要逐位对比就用 `--event 0`，
  它和全量模式跑的完全是同一串随机数。
- 量级自查：100 MeV 电子在空气里约 0.25 keV/mm，83 mm 阵列 → 每事件约 20 keV
  （45° 路径 117 mm 时约 29 keV）；10 µm 步长 → 每步约 2.5 eV，
  总电子数约 800/事件。
- 45° 斜穿时 x 和 z 都从 −41.5 mm 变到 +41.5 mm，覆盖 100 根柱子；
  y 保持不变（方向矢量的 y 分量为 0）。
- `hits` 一步一行。用 `plot.py` 之类的方式做统计时别在
  Python 里逐条循环，交给 ROOT 的 `TTree::Draw` / `RDataFrame`。
