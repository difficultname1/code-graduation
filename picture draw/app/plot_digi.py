#!/usr/bin/env python3
"""芯片上的热力图：读 digi.root 里 electrons 树的 (x, z) 落点，画成像素热力图。

芯片平面就是 X–Z 平面，一个柱子的截面就是"一个像素"（默认 83 µm 见方），
所以按 83 µm 分 bin 画出来的就是真实的像素图，不做任何插值。

  python3 plot_digi.py                    # 读 digi.root，写 chip_heatmap.png
  python3 plot_digi.py -i a.root -o b.png
  python3 plot_digi.py --pitch 0.083 --half 41.5

左图：整块芯片（±83/2 mm）；右图：自动缩放到 99% 落点所在区域。
两图都用 log 色阶（单个像素的计数跳好几个量级）。
另外还会写一张 <输出名>_height.png：左上 dx 随漂移长度 L 的展开（楔形），
右上 sigma_x(L) 与模型，下排是低/高漂移长度两块芯片分布。
"""
import argparse
import math
import sys

import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="digi.root")
    ap.add_argument("-o", "--output", default="chip_heatmap.png")
    ap.add_argument("--tree", default="electrons")
    ap.add_argument("--pitch", type=float, default=0.083, help="像素间距 [mm]")
    ap.add_argument("--half", type=float, default=41.5, help="芯片半宽 [mm]")
    ap.add_argument("--zoom-half", type=float, default=0.6,
                    help="右图缩放半宽 [mm]（默认 ±0.6 mm）")
    ap.add_argument("--pad", type=float, default=1.0,
                    help="左图在芯片外多画出的边距 [mm]，用来看溢出")
    ap.add_argument("--weight", action="store_true", help="按电子数 ne 加权")
    ap.add_argument("--chip-y", type=float, default=-20.0,
                    help="芯片平面 y [mm]（默认 -20 = 40 mm 柱阵的底面）")
    ap.add_argument("--diffusion", type=float, default=85.0,
                    help="漂移横向扩散系数 [um/sqrt(cm)]（默认 85，画模型曲线用）")
    ap.add_argument("--sigma-read", type=float, default=0.08,
                    help="电离倍增小孔定位宽度 sigma [mm]（默认 0.08，画模型曲线用）")
    args = ap.parse_args()

    fin = ROOT.TFile.Open(args.input)
    if not fin or fin.IsZombie():
        sys.exit(f"打不开 {args.input}")
    tree = fin.Get(args.tree)
    if not tree:
        sys.exit(f"{args.input} 里没有 {args.tree} 树")

    n_event = int(tree.GetEntries())
    if n_event == 0:
        sys.exit("树是空的")

    # ---- 读数据 ----
    try:
        import numpy as np                       # noqa: PLC0415
    except ImportError:
        np = None

    if np is not None:
        # RDataFrame 一次读成 numpy 数组，比逐条 GetEntry 快得多
        arr = ROOT.RDataFrame(args.tree, args.input).AsNumpy(
            ["event", "x", "z", "ne", "ysrc"])
        xs, zs, ws = arr["x"], arr["z"], arr["ne"].astype(float)
        ysrc = arr["ysrc"]
        evs = arr["event"]
        cnt = np.bincount(evs)
        ssum = np.bincount(evs, weights=xs)
        dx = xs - (ssum / np.maximum(cnt, 1))[evs]
    else:
        xs, zs, ws = [], [], []
        ysrc, evs = [], []
        for i in range(n_event):
            tree.GetEntry(i)
            xs.append(tree.x)
            zs.append(tree.z)
            ws.append(float(tree.ne))
            ysrc.append(tree.ysrc)
            evs.append(int(tree.event))
        cnt, ssum = {}, {}
        for x, e in zip(xs, evs):
            cnt[e] = cnt.get(e, 0) + 1
            ssum[e] = ssum.get(e, 0.0) + x
        dx = [x - ssum[e] / cnt[e] for x, e in zip(xs, evs)]

    n = len(xs)

    def quantile(v, q):
        s = sorted(v)
        return s[min(len(s) - 1, max(0, int(q * (len(s) - 1))))]

    margin = 0.5  # mm，别贴边
    zx = (quantile(xs, 0.005) - margin, quantile(xs, 0.995) + margin)
    zz = (quantile(zs, 0.005) - margin, quantile(zs, 0.995) + margin)

    # 缩放窗口：x 用"中位数 ± zoom_half"，z 用 0.5%–99.5% 分位（逐点模式下
    # 落点沿束流方向铺满整个芯片，不能再按峰值附近截一个方块）
    mx = quantile(xs, 0.5)
    zh = args.zoom_half
    zx = (mx - zh, mx + zh)
    zz = (quantile(zs, 0.005) - 0.1, quantile(zs, 0.995) + 0.1)

    pitch = args.pitch
    span = args.half + args.pad          # 左图范围（含芯片外的溢出区）
    nfull = max(1, int(round(2 * span / pitch)))
    nx = max(4, int(round((zx[1] - zx[0]) / pitch)))
    nz = max(4, int(round((zz[1] - zz[0]) / pitch)))

    h_full = ROOT.TH2D("h_full", ";x [mm];z [mm]",
                       nfull, -span, span, nfull, -span, span)
    h_zoom = ROOT.TH2D("h_zoom", ";x [mm];z [mm]",
                       nx, zx[0], zx[1], nz, zz[0], zz[1])

    # 落到芯片外面的比例（芯片真实尺寸只有 ±half）
    n_out = sum(1 for i in range(n)
                if abs(xs[i]) > args.half or abs(zs[i]) > args.half)

    # 手工 Fill：bin 就是 83 µm 的像素，不做插值
    for i in range(n):
        w = ws[i] if args.weight else 1.0
        h_full.Fill(xs[i], zs[i], w)
        h_zoom.Fill(xs[i], zs[i], w)

    for h in (h_full, h_zoom):
        h.SetMinimum(1.0)
        h.GetZaxis().SetTitle("hits / pixel")

    canvas = ROOT.TCanvas("c", "chip heatmap", 1900, 700)
    canvas.Divide(3, 1, 0.002, 0.002)  # 注意第 5 个参数是颜色，不是 bool

    for i, h in enumerate((h_full, h_zoom), 1):
        canvas.cd(i)
        ROOT.gPad.SetLeftMargin(0.10)
        ROOT.gPad.SetRightMargin(0.16)
        ROOT.gPad.SetBottomMargin(0.11)
        ROOT.gPad.SetLogz(True)
        h.Draw("COLZ")
        ROOT.gPad.Update()
    # ---- 第三幅：x 方向的一维投影 + 高斯拟合（直接读出分辨率）----
    canvas.cd(3)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    h_prof = ROOT.TH1D("h_prof", ";x - x_{event} [mm];hits", 160, -0.8, 0.8)
    for i in range(n):
        h_prof.Fill(dx[i], ws[i] if args.weight else 1.0)
    h_prof.SetLineColor(ROOT.kAzure + 1)
    h_prof.SetLineWidth(2)
    h_prof.SetMarkerStyle(20)
    h_prof.SetMarkerSize(0.5)
    h_prof.SetStats(0)
    h_prof.SetMaximum(h_prof.GetMaximum() * 1.25)   # 给拟合信息留出位置
    h_prof.Draw("E1")                               # 误差棒 = 每个 bin 的 Poisson ±1σ

    fit = ROOT.TF1("fit_gaus", "gaus", -0.45, 0.45)
    fit.SetLineColor(ROOT.kRed + 1)
    fit.SetLineWidth(2)
    fit.SetParameters(h_prof.GetMaximum(), 0.0, 0.16)
    h_prof.Fit(fit, "RQS")           # R=只在区间里拟合，Q=安静，S=存进直方图

    f_fit = h_prof.GetFunction("fit_gaus")   # Fit 存的是副本，取它来画
    if f_fit:
        f_fit.SetLineColor(ROOT.kRed + 1)
        f_fit.SetLineWidth(2)
        f_fit.Draw("SAME")
    canv_txt = ROOT.TLatex()
    canv_txt.SetNDC()
    canv_txt.SetTextSize(0.034)
    canv_txt.DrawLatex(0.15, 0.90, "x profile (per-event) + gauss fit")
    canv_txt.SetTextSize(0.038)
    canv_txt.DrawLatex(0.15, 0.82,
                       f"#sigma = {f_fit.GetParameter(2) * 1000:.1f} um"
                       f"   (#mu = {f_fit.GetParameter(1):.3f} mm)")

    # 左图里把芯片真实边界画出来：框外就是"本该丢掉"的溢出
    canvas.cd(1)
    box = ROOT.TBox(-args.half, -args.half, args.half, args.half)
    box.SetFillStyle(0)
    box.SetLineColor(ROOT.kBlack)
    box.SetLineWidth(2)
    box.Draw()

    label = f"{n} rows, pitch = {pitch * 1000:.0f} um"
    if args.weight:
        label += ", weighted by ne"
    t = ROOT.TLatex()
    t.SetNDC()
    t.SetTextSize(0.032)
    canvas.cd(1)
    t.DrawLatex(0.13, 0.92, "full chip (" + label + ")")
    t.SetTextSize(0.028)
    t.DrawLatex(0.13, 0.875,
                f"box = chip +-{args.half:.2f} mm; "
                f"{100.0 * n_out / max(1, n):.1f}% of hits fall outside")
    t.SetTextSize(0.032)
    canvas.cd(2)
    t.DrawLatex(0.13, 0.92, f"zoom: x +-{zh:.2f} mm, all z")

    canvas.SaveAs(args.output)
    print(f"[输出] {args.output}  ({n} 条, pitch {pitch * 1000:.0f} um, "
          f"weight={'ne' if args.weight else 'hits'})")
    print(f"[缩放] x {zx[0]:.2f}..{zx[1]:.2f} mm, z {zz[0]:.2f}..{zz[1]:.2f} mm")
    print(f"[溢出] 落在芯片(±{args.half} mm)外: {n_out} 条 = "
          f"{100.0 * n_out / max(1, n):.1f}%")
    print(f"[拟合] sigma_x = {f_fit.GetParameter(2) * 1000:.1f} um, "
          f"mu_x = {f_fit.GetParameter(1):.3f} mm")

    # ---- 第四张图：漂移高度的影响 ----
    # 用"扣除该事件自身入射 x"后的偏差 dx，否则随机枪本身的 x 分布会把扩散淹掉
    make_height_figure(xs, zs, dx, ysrc, args)


def make_height_figure(xs, zs, dx, ysrc, args):
    """四联图：把"产生得越高越模糊"直接画出来。

      左上  dx 随漂移长度 L 的分布（楔形张角 = D sqrt(L)）+ 模型 ±sigma 包络
      右上  sigma_x(L) 与模型曲线
      左下  L 小的那部分电子在芯片上的分布（窄）
      右下  L 大的那部分电子在芯片上的分布（宽）
    """
    drift = [y - args.chip_y for y in ysrc]
    if not drift:
        return
    dmax = max(drift)
    if dmax <= 0.0:
        return

    # 固定枪的电荷几乎都产生在同一个高度上，只有一两个 bin 有统计，而且远端的
    # 少量 δ 射线会带出假的"宽度变大"。看 5%~95% 分位覆盖多宽，覆盖不到就不出图。
    ds = sorted(drift)
    q05 = ds[int(0.05 * (len(ds) - 1))]
    q95 = ds[int(0.95 * (len(ds) - 1))]
    if q95 - q05 < 0.3 * dmax:
        print(f"[扩散] 漂移距离的 5%~95% 只覆盖 {q05:.0f}..{q95:.0f} mm"
              f"（枪在同一高度上），跳过高度图；用 run_height.mac 或随机枪")
        return

    def model_sigma(d):          # um
        sd = args.diffusion * 1e-3 * math.sqrt(d / 10.0)     # um/sqrt(cm) -> mm
        return 1000.0 * math.sqrt(args.sigma_read ** 2 + sd * sd)

    def robust_sigma(v):
        v = sorted(v)
        def q(p):
            return v[min(len(v) - 1, max(0, int(p * (len(v) - 1))))]
        return 0.5 * (q(0.8413) - q(0.1587))      # (q84 - q16)/2 ~ sigma

    # ---------- 左上：dx vs L 楔形 ----------
    nbL = max(8, int(round(dmax / 2.0)))
    h_wedge = ROOT.TH2D("h_wedge", ";drift length L [mm];x - x_{event} [um]",
                        nbL, 0.0, dmax, 120, -450.0, 450.0)
    for i, d in enumerate(drift):
        h_wedge.Fill(d, dx[i] * 1000.0)

    up, dn = ROOT.TGraph(), ROOT.TGraph()
    for i in range(200):
        d = dmax * i / 199.0
        s = model_sigma(d)
        up.SetPoint(i, d, +s)
        dn.SetPoint(i, d, -s)
    for g in (up, dn):
        g.SetLineColor(ROOT.kRed + 1)
        g.SetLineWidth(2)
        g.SetLineStyle(2)

    # ---------- 右上：每个 bin 的平均落点 ±1σ（K 线式的误差棒）----------
    nb = max(4, int(round(dmax / 10.0)))       # 约 10 mm 一组
    hw = ROOT.TH1D("h_sigma_vs_y",
                   ";drift length L [mm];#LTx - x_{event}#GT [um]",
                   nb, 0.0, dmax)
    samples = [[] for _ in range(nb)]
    for i, d in enumerate(drift):
        k = min(max(hw.GetXaxis().FindBin(d) - 1, 0), nb - 1)
        samples[k].append(dx[i])

    gr = ROOT.TGraphErrors()                   # 只拿 sigma(L) 去做拟合
    ip = 0
    smax = 1.0
    for k, v in enumerate(samples):
        if len(v) < max(100, len(drift) // (nb * 20)):
            continue
        s = robust_sigma(v)                    # 这个 bin 的 1 sigma
        mean = sum(v) / len(v)                 # 平均落点（应贴着 0）
        xc = hw.GetXaxis().GetBinCenter(k + 1)
        # 点 = 平均落点，误差棒 = 点 ± 1 sigma（棒的总长 = 2 sigma）
        hw.SetBinContent(k + 1, mean * 1000.0)
        hw.SetBinError(k + 1, s * 1000.0)
        gr.SetPoint(ip, xc, s * 1000.0)
        gr.SetPointError(ip, 0.0, s * 1000.0)
        ip += 1
        smax = max(smax, s * 1000.0)
        print(f"[高度] L = {xc:5.1f} mm: <x> = {mean * 1000:+6.1f} um, "
              f"±1sigma = {s * 1000:6.1f} um"
              f"（输入模型 {model_sigma(xc):6.1f}），n = {len(v)}")
    hw.SetMarkerStyle(20)
    hw.SetMarkerSize(1.2)
    hw.SetLineWidth(2)
    hw.SetLineColor(ROOT.kAzure + 1)
    hw.SetMarkerColor(ROOT.kAzure + 1)
    hw.SetStats(0)
    hw.SetMinimum(-1.15 * smax)
    hw.SetMaximum(+1.15 * smax)

    # 用数据自己拟合 sigma^2 = sigma0^2 + D^2 * L，看拟合出的 D 和输入的 85
    # 差多少（sigma_0 里会吃掉气体本身的物理展宽，比 80 um 略大）
    fitf = ROOT.TF1("f_diff", "sqrt([0]*[0]+[1]*x)", 0.0, dmax)
    fitf.SetParameters(args.sigma_read * 1000.0, 722.5)   # D=85 -> 85^2/10
    fitf.SetLineColor(ROOT.kGreen + 2)
    fitf.SetLineWidth(3)
    a0, bfit = 0.0, 0.0
    if ip >= 2:
        gr.Fit(fitf, "RQN")
        a0 = abs(fitf.GetParameter(0))       # 参数以 a0^2 出现，取正
        bfit = fitf.GetParameter(1)
        print(f"[拟合] sigma_x(L) = sqrt({a0:.1f}^2 + {max(bfit, 0.0):.1f}*L)，"
              f"即 D = {math.sqrt(max(bfit, 0.0) * 10):.1f} um/sqrt(cm)")

    # 两条 ±sigma(L) 包络：输入模型（红虚线）和拟合结果（绿实线）
    def envelope(fn, color, style, width):
        gu, gd = ROOT.TGraph(), ROOT.TGraph()
        for i in range(200):
            d = dmax * i / 199.0
            s = fn(d)
            gu.SetPoint(i, d, +s)
            gd.SetPoint(i, d, -s)
        for g in (gu, gd):
            g.SetLineColor(color)
            g.SetLineStyle(style)
            g.SetLineWidth(width)
        return gu, gd

    up_mod, dn_mod = envelope(model_sigma, ROOT.kRed + 1, 2, 2)
    if ip >= 2:
        up_fit, dn_fit = envelope(
            lambda d: 1000.0 * math.sqrt(max(a0 * a0 + bfit * d, 0.0)),
            ROOT.kGreen + 2, 1, 2)

    # ---------- 下排：低/高漂移长度的芯片分布 ----------
    xhalf = max(0.4, 3.0 * smax / 1000.0)      # 两块图用同一个 x 范围

    def chip_map(name, idx):
        h = ROOT.TH2D(name, ";x - x_{event} [mm];z [mm]",
                      2 * int(round(xhalf / 0.01)), -xhalf, xhalf,
                      84, -4.2, 4.2)
        for i in idx:
            h.Fill(dx[i], zs[i])
        return h

    i_lo = [i for i, d in enumerate(drift) if d < 0.25 * dmax]
    i_hi = [i for i, d in enumerate(drift) if d > 0.75 * dmax]
    h_lo = chip_map("h_low", i_lo)
    h_hi = chip_map("h_high", i_hi)
    # 两块图用同一个色阶，才能直接比"窄/宽"
    ztop = max(h_lo.GetMaximum(), h_hi.GetMaximum(), 1.0)
    for h in (h_lo, h_hi):
        h.GetZaxis().SetRangeUser(1.0, ztop)
    s_lo = robust_sigma([dx[i] for i in i_lo]) * 1000.0 if i_lo else 0.0
    s_hi = robust_sigma([dx[i] for i in i_hi]) * 1000.0 if i_hi else 0.0

    c = ROOT.TCanvas("c_height", "drift height effect", 1800, 1200)
    c.Divide(2, 2, 0.002, 0.002)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetLogz(True)
    h_wedge.SetMinimum(1.0)
    h_wedge.GetZaxis().SetTitle("electrons / bin")
    h_wedge.Draw("COLZ")
    up.Draw("L SAME")
    dn.Draw("L SAME")
    t1 = ROOT.TLatex()
    t1.SetNDC()
    t1.SetTextSize(0.032)
    t1.DrawLatex(0.15, 0.965, "charge spread vs drift length")
    t1.SetTextSize(0.026)
    t1.DrawLatex(0.15, 0.905, "dashed: #pm#sqrt{#sigma_{read}^{2}+(D#sqrt{L})^{2}}")

    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.16)
    hw.Draw("E1")
    up_mod.Draw("L SAME")
    dn_mod.Draw("L SAME")
    if ip >= 2:
        up_fit.Draw("L SAME")
        dn_fit.Draw("L SAME")
    t2 = ROOT.TLatex()
    t2.SetNDC()
    t2.SetTextSize(0.031)
    t2.DrawLatex(0.16, 0.965, "hit position #pm 1#sigma vs drift length")
    t2.SetTextSize(0.023)
    t2.DrawLatex(0.16, 0.905, "point = mean, bar = point #pm 1#sigma (bin)")
    if ip >= 2:
        t2.DrawLatex(0.16, 0.855, "red dashed: input #pm#sqrt{80^{2}+(85#sqrt{L})^{2}}  |  "
                                 "green: fit #pm#sqrt{%0.0f^{2}+(%0.0f#sqrt{L})^{2}}"
                     % (a0, math.sqrt(max(bfit, 0.0) * 10)))

    sig_lines = []                       # 必须跨 pad 持有引用，否则被 Python 回收
    for pad, h, title, sig, dlo, dhi in (
            (3, h_lo, f"L < {0.25 * dmax:.0f} mm", s_lo, 0.0, 0.25 * dmax),
            (4, h_hi, f"L > {0.75 * dmax:.0f} mm", s_hi, 0.75 * dmax, dmax)):
        c.cd(pad)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTopMargin(0.16)
        ROOT.gPad.SetLogz(True)
        h.SetMinimum(1.0)
        h.GetZaxis().SetTitle("electrons / pixel")
        h.Draw("COLZ")
        t = ROOT.TLatex()
        t.SetNDC()
        t.SetTextSize(0.032)
        t.DrawLatex(0.14, 0.965, f"chip plane, {title}")
        t.SetTextSize(0.026)
        t.DrawLatex(0.14, 0.905, f"#sigma_{{x}} = {sig:.0f} um (dashed = #pm1#sigma)")
        for sx in (-sig / 1000.0, +sig / 1000.0):
            line = ROOT.TLine(sx, -4.2, sx, 4.2)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("SAME")
            sig_lines.append(line)

    c.SaveAs(args.output.replace(".png", "_vs_drift.png"))

    # ---------- 对数版：左 σx(L)，右 低/高漂移的 x 剖面 ----------
    cl = ROOT.TCanvas("c_log", "log view", 1600, 700)
    cl.Divide(2, 1, 0.002, 0.002)

    cl.cd(1)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetLogy(True)
    gsig = ROOT.TGraph()
    for i in range(gr.GetN()):
        gsig.SetPoint(i, gr.GetX()[i], gr.GetY()[i])
    gsig.SetMarkerStyle(20)
    gsig.SetMarkerSize(1.2)
    gsig.SetMarkerColor(ROOT.kAzure + 1)
    gsig.SetLineWidth(2)
    gsig.SetMinimum(80.0)
    gsig.SetMaximum(600.0)
    gsig.SetTitle(";drift length L [mm];#sigma_{x} [um]")
    gsig.Draw("AP")
    gline_in = ROOT.TF1("f_model_in", "sqrt([0]*[0]+[1]*x)", 0.0, dmax)
    gline_in.SetParameters(args.sigma_read * 1000.0,
                           (args.diffusion ** 2) / 10.0)   # (um/sqrt(cm))^2/10 = um^2/mm
    gline_in.SetLineColor(ROOT.kRed + 1)
    gline_in.SetLineStyle(2)
    gline_in.SetLineWidth(2)
    gline_in.Draw("L SAME")
    if ip >= 2:
        fitf.SetLineColor(ROOT.kGreen + 2)
        fitf.SetLineWidth(3)
        fitf.Draw("L SAME")
    leg1 = ROOT.TLegend(0.42, 0.14, 0.92, 0.32)
    leg1.AddEntry(gsig, "data (bin #sigma_{x})", "p")
    leg1.AddEntry(gline_in, "input #sqrt{80^{2}+(85#sqrt{L})^{2}}", "l")
    if ip >= 2:
        leg1.AddEntry(fitf, "fit #sqrt{%0.0f^{2}+(%0.0f#sqrt{L})^{2}}"
                      % (a0, math.sqrt(max(bfit, 0.0) * 10)), "l")
    leg1.Draw()
    tl = ROOT.TLatex()
    tl.SetNDC()
    tl.SetTextSize(0.030)
    tl.DrawLatex(0.15, 0.965, "#sigma_{x}(L), log scale")

    cl.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetLogy(True)
    hp_lo = ROOT.TH1D("hp_low", ";x - x_{event} [um];electrons / bin", 80, -800, 800)
    hp_hi = ROOT.TH1D("hp_high", "", 80, -800, 800)
    for i in i_lo:
        hp_lo.Fill(dx[i] * 1000.0)
    for i in i_hi:
        hp_hi.Fill(dx[i] * 1000.0)
    hp_lo.SetLineColor(ROOT.kAzure + 1)
    hp_hi.SetLineColor(ROOT.kRed + 1)
    for hp in (hp_lo, hp_hi):
        hp.SetLineWidth(2)
        hp.SetStats(0)
        hp.SetMinimum(0.5)
    hp_lo.SetMaximum(max(hp_lo.GetMaximum(), hp_hi.GetMaximum()) * 3.0)
    hp_lo.Draw("HIST")
    hp_hi.Draw("HIST SAME")
    leg2 = ROOT.TLegend(0.55, 0.60, 0.92, 0.76)
    leg2.AddEntry(hp_lo, f"L < {0.25 * dmax:.0f} mm: {s_lo:.0f} um", "l")
    leg2.AddEntry(hp_hi, f"L > {0.75 * dmax:.0f} mm: {s_hi:.0f} um", "l")
    leg2.Draw()
    tl2 = ROOT.TLatex()
    tl2.SetNDC()
    tl2.SetTextSize(0.030)
    tl2.DrawLatex(0.15, 0.965, "x profile at low / high drift, log scale")

    cl.SaveAs(args.output.replace(".png", "_log.png"))

    dmax_s = model_sigma(dmax)
    print(f"[扩散] D = {args.diffusion:.0f} um/sqrt(cm): "
          f"L=0 -> {model_sigma(0.0):.0f} um, "
          f"L={dmax:.0f} mm -> {dmax_s:.0f} um；"
          f"实测 {s_lo:.0f} um（低 L）-> {s_hi:.0f} um（高 L）")


if __name__ == "__main__":
    main()
