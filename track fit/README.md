# track fit

从 `hits`（Geant4 每个沉积步一条）拟合每个事件的径迹。

## 现状

只有骨架：**读取和输出是通的**，拟合算法还没写。

```
track_fit/io.py    读 hits、按事件切组、建并写 tracks 树
track_fit/fit.py   fit_track() —— 算法在这里，现在直接 return None
run_track_fit.py   命令行入口，串起 "读 -> 拟合 -> 写"
```

跑一遍：

```bash
python3 "track fit/run_track_fit.py" pixel.root -o track_fit.root
```

会读完整棵树（打印读了几个事件），写出一个**结构正确但为空**的 `tracks` 树：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号 |
| `nx`, `ny` | `Double_t` | 方向余弦 |
| `x0`, `y0` | `Double_t` | 参考平面上的截距 [mm] |
| `chi2` | `Double_t` | 拟合优度 |
| `ndf` | `Int_t` | 自由度 |

## 要填的地方

`track_fit/fit.py::fit_track(event, points)`：

- 入参：事件号，以及这个事件所有沉积点的 `x / y / z / edep`（mm / MeV）；
- 返回值：`dict`（键对应上表，缺的字段写 0）或 `None`（跳过这个事件）。

改完直接重跑上面那条命令，`io.py` 和命令行入口都不用动。
