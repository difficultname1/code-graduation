# class

分类：把 `digitize` 的芯片沉积（`electrons`）按径迹分组，输出带 `trackID`
的 `points` 树。**`trackID` 是 class 的产物**，上游没有这个量。

```
classifier/io.py        读 electrons（event / mcTrackID / x / z / ne / ysrc），写 points
classifier/classify.py  classify() —— 分类钩子；perfect() —— 真值捷径
classifier/hits.py      按像素把电子合并成命中
app/run_class.py        命令行入口
```

```bash
python3 "class/app/run_class.py" digi.root -o points.root
```

## 输入：electrons

`digitize` 输出的 `electrons` 树，每个电子落点一条：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号 |
| `mcTrackID` | `Int_t` | MC 真值（初级 = 1，δ 射线 > 1）；**不是分类结果** |
| `x`, `z` | `Double_t` | 芯片上的落点 [mm] |
| `ne` | `Int_t` | 这一条代表的电子数 |
| `ysrc` | `Double_t` | 产生点 y [mm]（算漂移扩散用） |

## 输出：points

默认一个电子一条，`x / z / ne` 就是电子自己的落点。加 `--pitch 0.083` 则先按
83 µm 的像素把同一条径迹的电子合并，每个像素命中一条（`x / z` 是像素内电子的
电荷加权质心，`ne` 是该像素的总电荷）：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号 |
| `trackID` | `Int_t` | **class 分出来的径迹号**（同一事件内唯一） |
| `x`, `z` | `Double_t` | 落点（`--pitch 0`）或像素内电子的电荷加权质心 [mm] |
| `ne` | `Int_t` | 电子数（逐电子时恒为 1；按像素合并时是像素电荷） |
| `ysrc` | `Double_t` | 产生高度（按像素合并时是平均高度）[mm] |

## 现状：跳过分类器

真实分类器应该把电子按径迹聚类，再给每条径迹一个 `trackID`。当前阶段
`classify.classify()` 返回 `None`，入口脚本改走 `classify.perfect()`：
直接把 MC 真值 `mcTrackID` 当作 `trackID`，表示“分类完全正确”。
换成真实分类器时，只要让 `classify()` 返回带 `trackID` 的点列表，
`io.py`、像素合并和入口脚本都不用动。

按仓库约定：代码只放在子目录里，这个文件的旁边只放说明。
