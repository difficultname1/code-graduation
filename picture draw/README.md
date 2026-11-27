# picture draw

画图。入口脚本都在 `app/` 下，输入是前面几步产生的 ROOT 文件。

```bash
python3 "picture draw/app/plot.py" pixel.root            # hits：柱平均、edep 分布、y 分布
python3 "picture draw/app/plot_digi.py" -i digi.root     # electrons：芯片热力图 + x 剖面
```

| 脚本 | 读什么 | 画什么 |
| --- | --- | --- |
| `app/plot.py` | `hits` 树 | 每根柱子当 bin 的 x–z 能量图、edep / y 分布、`columns_mean.csv` |
| `app/plot_digi.py` | `electrons` 树 | 芯片像素热力图（830 µm 分 bin）、x 剖面拟合、随漂移长度的展宽 |

两个脚本都只依赖 ROOT（`plot_digi.py` 有 numpy 时走 `RDataFrame`，没有就退回逐条读）。
输出前缀从输入文件名自动推，也可以用参数改。
