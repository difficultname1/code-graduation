#!/usr/bin/env python3
"""芯片上的热力图：读 digi.root 的 (x, z) 落点，画成像素热力图。

芯片平面就是 X–Z 平面，一个柱子的截面就是"一个像素"（默认 83 µm 见方），
所以按 83 µm 分 bin 画出来的就是真实的像素图，不做任何插值。

  python3 plot_digi.py                    # 读 digi.root，写 chip_heatmap.png
  python3 plot_digi.py -i a.root -o b.png
  python3 plot_digi.py --pitch 0.083 --half 4.15

左图：整块芯片（±8.3/2 mm）；右图：自动缩放到 99% 落点所在区域。
两图都用 log 色阶（单个像素的计数跨越好几个量级）。
"""
import argparse
import sys

import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="digi.root")
    ap.add_argument("-o", "--output", default="chip_heatmap.png")
    ap.add_argument("--tree", default="digi")
    ap.add_argument("--pitch", type=float, default=0.083, help="像素间距 [mm]")
    ap.add_argument("--half", type=float, default=4.15, help="芯片半宽 [mm]")
    ap.add_argument("--zoom-half", type=float, default=0.6,
                    help="右图缩放半宽 [mm]（默认 ±0.6 mm）")
    ap.add_argument("--pad", type=float, default=1.0,
                    help="左图在芯片外多画出的边距 [mm]，用来看溢出")
    ap.add_argument("--weight", action="store_true", help="按电子数 ne 加权")
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
        arr = ROOT.RDataFrame(args.tree, args.input).AsNumpy(["x", "z", "ne"])
        xs, zs, ws = arr["x"], arr["z"], arr["ne"].astype(float)
    else:
        xs, zs, ws = [], [], []
        for i in range(n_event):
            tree.GetEntry(i)
            xs.append(tree.x)
            zs.append(tree.z)
            ws.append(float(tree.ne))

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
    h_prof = ROOT.TH1D("h_prof", ";x [mm];hits", 160, mx - 0.8, mx + 0.8)
    for i in range(n):
        h_prof.Fill(xs[i], ws[i] if args.weight else 1.0)
    h_prof.SetLineColor(ROOT.kAzure + 1)
    h_prof.SetLineWidth(2)
    h_prof.SetStats(0)
    h_prof.SetMaximum(h_prof.GetMaximum() * 1.25)   # 给拟合信息留出位置
    h_prof.Draw("HIST")

    fit = ROOT.TF1("fit_gaus", "gaus", mx - 0.45, mx + 0.45)
    fit.SetLineColor(ROOT.kRed + 1)
    fit.SetLineWidth(2)
    fit.SetParameters(h_prof.GetMaximum(), mx, 0.16)
    h_prof.Fit(fit, "RQS")           # R=只在区间里拟合，Q=安静，S=存进直方图

    f_fit = h_prof.GetFunction("fit_gaus")   # Fit 存的是副本，取它来画
    if f_fit:
        f_fit.SetLineColor(ROOT.kRed + 1)
        f_fit.SetLineWidth(2)
        f_fit.Draw("SAME")
    canv_txt = ROOT.TLatex()
    canv_txt.SetNDC()
    canv_txt.SetTextSize(0.034)
    canv_txt.DrawLatex(0.15, 0.90, "x profile + gauss fit")
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


if __name__ == "__main__":
    main()
