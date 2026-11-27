# track fit

从 `class` 输出的 `points`（每条落点带 class 产生的 `trackID`）拟合每条径迹的直线。
上游 `electrons` 树里没有 `trackID`，只有 MC 真值 `mcTrackID`。

## 做法

芯片平面是 X–Z 平面，径迹参数化成

$$x = A\,z + B$$

用**加权最小二乘**求 A、B。每个落点有三样东西进场：

| 量 | 用法 |
| --- | --- |
| 电子数 `ne` | 电荷加权：q 个电子的质心误差是单电子的 $1/\sqrt q$ |
| 产生点高度 `ysrc` | 漂移扩散误差：$\sigma_d = D\sqrt{L}$，$L=\lvert y_{\rm src}-y_{\rm chip}\rvert$ |
| 本底 `bkg` | 每个点扣掉的本底电荷，默认取该径迹的平均 `ne` |

单电子位置方差

$$\sigma_i^2 = \sigma_{\rm read}^2 + \left(D\sqrt{L_i}\right)^2,$$

扣完本底后的信号电荷 $q_i=\max(ne_i-b,0)$，权重

$$w_i = \frac{q_i}{\sigma_i^2}.$$

拟合公式还是普通加权最小二乘：

$$
A = \frac{\sum w \sum w z x - \sum w z \sum w x}
         {\sum w \sum w z^2 - (\sum w z)^2},\qquad
B = \frac{\sum w x - A \sum w z}{\sum w}.
$$

`q_i <= 0` 的点不参与拟合。每电子一条时 `ne` 恒为 1，平均本底就是 1，
扣完一个点都不剩，程序会自动退回 `bkg = 0`（等价于普通最小二乘），
并在命令行给出提示；`digitize --per-cluster` 时 `ne` 有涨落，本底扣除才
真正起作用。如果扣完本底剩不到 3 个点，也会退回不扣本底。

### 卡方

x 和 z 都有测量误差，残差方差是 $\sigma_i^2(1+A^2)$，所以

$$
\chi^2 = \sum_i \frac{q_i\,r_i^2}{\sigma_i^2\,(1+A^2)},
\qquad r_i = x_i - A z_i - B,
\qquad \mathrm{ndf}=N_{\rm used}-2.
$$

输出 `chi2`、`ndf`、`chi2ndf` 和 `pvalue`（`ROOT.TMath.Prob`）。

### 误差传递

拟合用的权重实际是 $w_i/(1+A^2)$，所以参数协方差矩阵

$$
\operatorname{Cov} = (1+A^2)\,\frac{\chi^2}{\rm ndf}\,
\left(X^T W X\right)^{-1},
\qquad X = \begin{bmatrix} z & 1 \end{bmatrix},
$$

即

$$
\sigma_A^2 = (1+A^2)\,\frac{\chi^2}{\rm ndf}\frac{S_0}{\det},\quad
\sigma_B^2 = (1+A^2)\,\frac{\chi^2}{\rm ndf}\frac{S_2}{\det},\quad
\operatorname{Cov}(A,B) = -(1+A^2)\,\frac{\chi^2}{\rm ndf}\frac{S_1}{\det},
$$

其中 $S_0=\sum w$、$S_1=\sum w z$、$S_2=\sum w z^2$、
$\det=S_0S_2-S_1^2$。输出 `sigmaA`、`sigmaB`、`covAB`、`corrAB`。

```
track_fit/io.py    读 points、按 (event, trackID) 切组、写 tracks 树
track_fit/fit.py   fit_track() —— 电荷加权 + 本底扣除 + 误差传递
app/run_track_fit.py   命令行入口
```

跑一遍：

```bash
python3 "track fit/app/run_track_fit.py" points.root -o track_fit.root

# 误差模型 / 本底选项
python3 "track fit/app/run_track_fit.py" points.root -o track_fit.root \
  --sigma-read 0.08   # 电离倍增小孔 [mm]
  --diffusion 85      # [um/sqrt(cm)]
  --chip-y -20        # 芯片平面 [mm]
  --background mean   # mean（默认）/ none / 固定数值
```

输出 `tracks` 树，每条径迹一行：

| 分支 | 类型 | 含义 |
| --- | --- | --- |
| `event` | `Int_t` | 事件号 |
| `trackID` | `Int_t` | 径迹号 |
| `A` | `Double_t` | 直线斜率（`x = A*z + B`） |
| `B` | `Double_t` | 截距 [mm] |
| `sigmaA`, `sigmaB` | `Double_t` | A、B 的误差（误差传递） |
| `covAB`, `corrAB` | `Double_t` | 协方差、相关系数 |
| `chi2`, `chi2ndf`, `pvalue` | `Double_t` | 卡方、卡方/ndf、p 值 |
| `ndf` | `Int_t` | 自由度（有效点数 − 2） |
| `rms` | `Double_t` | 加权残差 RMS [mm] |
| `bkg` | `Double_t` | 实际扣掉的本底电荷（退回时为 0） |
| `n` | `Int_t` | 参与拟合的点数 |

少于 3 个有效点或所有 `z` 相同的径迹会被跳过（`fit_track()` 返回 `None`）。
当前场景是 45° 入射，主径迹（`trackID = 1`）的 A 应接近 1、B 接近 0，
`chi2/ndf` 应在 1 附近（20 事件实测：A 中位数 0.9988，$\sigma_A\approx4\times10^{-4}$，
`chi2/ndf` 中位数 1.00）。
