# vertex fit

从 `electrons`（芯片上的电子落点）重建每个事件的入射位置。

## 现状

只有骨架：**读取和输出是通的**，重建算法还没写。

```
vertex_fit/io.py    读 electrons、按事件切组、建并写 vertices 树
vertex_fit/fit.py   fit_vertex() —— 算法在这里，现在直接 return None
app/run_vertex_fit.py   命令行入口
```

跑一遍：

```bash
python3 "vertex fit/app/run_vertex_fit.py" digi.root -o vertex_fit.root   # digi.root 里的 electrons 树
```

会写出一棵**结构正确但为空**的 `vertices` 树：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号 |
| `x`, `z` | `Double_t` | 芯片平面上的位置 [mm] |
| `y` | `Double_t` | 产生点的平均高度 [mm] |
| `t` | `Double_t` | 顶点时间 [ns] |
| `n` | `Int_t` | 参与重建的落点数 |

## 要填的地方

`vertex_fit/fit.py::fit_vertex(event, hits)`：

- 入参：事件号，以及这个事件在芯片上的落点 `x / z / ne / ysrc / t`；
- 返回值：`dict`（键对应上表）或 `None`（跳过）。

`io.py` 与命令行入口不用动。
