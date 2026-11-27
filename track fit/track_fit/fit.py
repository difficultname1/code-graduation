"""track fit 的算法层：芯片平面 X–Z 上的加权最小二乘直线拟合。

径迹参数化成

    x = A*z + B

每个落点带电荷数 ne（电子数）。拟合时做三件事：

1. **电荷加权**：一个落点由 q 个电子组成，位置误差是单电子误差的
   1/sqrt(q)，所以权重里带 q。每电子一条时 q = 1。
2. **本底扣除**：每个点先减去本底 b（默认取该径迹的平均电子数，
   也可以给固定值或关掉）。q_sig = ne - b <= 0 的点不参与拟合。
3. **误差传递**：单电子位置方差取

       sigma_i^2 = sigma_read^2 + (D * sqrt(L_i))^2,

   L_i 是产生点离芯片的漂移距离。拟合 x 对 z 时 z 也有同样的测量误差，
   残差方差要乘 (1 + A^2)，所以

       chi2 = sum q_i * r_i^2 / [sigma_i^2 * (1 + A^2)]

   拟合权重实际是 w/(1+A^2)，参数协方差为

       Cov = (1 + A^2) * (chi2/ndf) * (X^T W X)^(-1),

   从中给出 sigmaA、sigmaB、covAB、corrAB。
"""

import math

try:
    import ROOT
except ImportError:                 # 没有 PyROOT 也能做拟合，只是没有 p 值
    ROOT = None


def _background_value(mode, charges):
    """把 --background 的取值解析成一个数。

    "mean" -> 这条径迹的平均电子数；"none"/"off" -> 0；数字字符串 -> 固定本底。
    """
    if mode is None:
        return 0.0
    if isinstance(mode, str):
        m = mode.strip().lower()
        if m in ("none", "off"):
            return 0.0
        if m == "mean":
            return sum(charges) / len(charges) if charges else 0.0
        return float(m)
    return float(mode)


def fit_track(event, track_id, points, sigma_read=0.08, diffusion=85.0,
              chip_y=-20.0, background="mean"):
    """对一条径迹做直线拟合 x = A*z + B。

    参数
    ----
    event : int
        事件号（只用于调试，不参与计算）。
    track_id : int
        class 分出来的径迹号。
    points : list[dict]
        这条径迹在芯片平面上的落点，每项有 x / z（mm）、ne（电子数），
        可选 ysrc（产生点 y，mm；没有就按芯片平面算，扩散项为 0）。
    sigma_read : float
        电离倍增小孔的定位宽度 [mm]，默认 0.08。
    diffusion : float
        漂移横向扩散系数 [um/sqrt(cm)]，默认 85（内部换算成 mm/sqrt(cm)）。
    chip_y : float
        芯片平面 y [mm]，默认 -20。
    background : str or float
        每个点扣除的本底电荷。默认 "mean" = 该径迹的平均 ne；
        "none" 不扣；数字按固定值扣。

    返回
    ----
    dict 或 None
        None 表示点数不够或退化（所有 z 相同），调用方跳过。
        否则返回 A / B / sigmaA / sigmaB / covAB / corrAB / chi2 / ndf /
        chi2ndf / pvalue / rms / bkg / n。
    """
    if len(points) < 3:            # 两点必然过直线，没有自由度
        return None

    charges = [float(p.get("ne", 1)) for p in points]
    bkg = _background_value(background, charges)

    # 扣本底后的有效电荷；q <= 0 的点没有信号，不参与拟合
    rows = []
    for p, q0 in zip(points, charges):
        q = q0 - bkg
        if not (q > 0.0):
            continue
        y = p.get("ysrc")
        var = sigma_read * sigma_read
        if y is not None and math.isfinite(float(y)) and diffusion > 0.0:
            L = abs(float(y) - chip_y)          # mm
            sd = diffusion * 1e-3 * math.sqrt(L / 10.0)   # mm
            var += sd * sd
        if not (var > 0.0):
            var = 1e-12
        rows.append((float(p["x"]), float(p["z"]), q, var))

    # 本底减过头（比如每电子一条、ne 恒为 1，平均本底就是 1）：
    # 退回不扣本底，保证每电子一条时仍然能做普通最小二乘。
    if len(rows) < 3:
        if bkg != 0.0:
            return fit_track(event, track_id, points, sigma_read, diffusion,
                             chip_y, background=0.0)
        return None

    # 加权最小二乘：
    #   A = (Σw Σwzx - Σwz Σwx) / (Σw Σwz² - (Σwz)²)
    #   B = (Σwx - A Σwz) / Σw
    # 权重 w = q / sigma²
    S0 = S1 = S2 = T0 = T1 = 0.0
    for x, z, q, var in rows:
        w = q / var
        S0 += w
        S1 += w * z
        S2 += w * z * z
        T0 += w * x
        T1 += w * z * x

    det = S0 * S2 - S1 * S1
    # det 为 0（所有 z 一样）或数值上退化时不做拟合
    if not (det > 0.0) or det < 1e-12 * S0 * S2:
        return None

    A = (S0 * T1 - S1 * T0) / det
    B = (T0 - A * S1) / S0

    # 卡方：x 和 z 都有测量误差，残差方差 = sigma² * (1 + A²)
    scale = 1.0 + A * A
    chi2 = 0.0
    resid2 = 0.0
    for x, z, q, var in rows:
        r = x - (A * z + B)
        w = q / var
        chi2 += w * r * r / scale
        resid2 += w * r * r

    ndf = len(rows) - 2
    if ndf <= 0:
        return None

    # 误差传递：拟合权重实际是 w/(1+A^2)，所以
    #   Cov = (1+A^2) * (chi2/ndf) * (X^T W X)^(-1)
    #   X = [z, 1]，X^T W X = [[S2, S1], [S1, S0]]
    s2 = chi2 / ndf
    cov_scale = s2 * scale          # = sum(w r^2) / ndf
    sigmaA = math.sqrt(max(cov_scale * S0 / det, 0.0))
    sigmaB = math.sqrt(max(cov_scale * S2 / det, 0.0))
    covAB = -cov_scale * S1 / det
    corrAB = covAB / (sigmaA * sigmaB) if sigmaA > 0.0 and sigmaB > 0.0 else 0.0

    pvalue = 0.0
    if ROOT is not None:
        pvalue = float(ROOT.TMath.Prob(chi2, ndf))

    rms = math.sqrt(resid2 / S0) if S0 > 0.0 else 0.0

    return dict(A=A, B=B, sigmaA=sigmaA, sigmaB=sigmaB,
                covAB=covAB, corrAB=corrAB,
                chi2=chi2, ndf=ndf, chi2ndf=s2, pvalue=pvalue,
                rms=rms, bkg=bkg, n=len(rows))
