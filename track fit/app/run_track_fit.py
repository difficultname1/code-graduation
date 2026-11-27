#!/usr/bin/env python3
"""track fit 命令行入口：读 class 输出 -> 直线拟合 -> 写 tracks 树。

输入是 class 的 points 树（event / trackID / x / z / ne / ysrc），
输出每条径迹一行：x = A*z + B 的 A、B、误差传递结果，以及 chi2 / ndf。
注意 trackID 是 class 的产物；上游 electrons 树里只有 MC 真值 mcTrackID。

用法:
    python3 "track fit/app/run_track_fit.py" class_points.root -o track_fit.root
"""
import argparse
import os
import statistics
import sys

# 本文件在 <分类>/app/ 下，包（track_fit/）在上一级，所以把上级目录加进 sys.path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import ROOT  # noqa: E402  （PyROOT 要在 sys.path 之后导入）

from track_fit import fit, io  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="class points -> 直线径迹拟合")
    ap.add_argument("input", help="class 的输出文件（class_points.root 等）")
    ap.add_argument("-o", "--output", default="track_fit.root", help="输出文件")
    ap.add_argument("--tree", default="points", help="输入树名（默认 points）")
    ap.add_argument("--sigma-read", type=float, default=0.08,
                    help="电离倍增小孔定位宽度 [mm]（默认 0.08）")
    ap.add_argument("--diffusion", type=float, default=85.0,
                    help="漂移横向扩散系数 [um/sqrt(cm)]（默认 85）")
    ap.add_argument("--chip-y", type=float, default=-20.0,
                    help="芯片平面 y [mm]（默认 -20）")
    ap.add_argument("--background", default="mean",
                    help="每个点扣除的本底电荷：mean（默认）/ none / 固定数值")
    args = ap.parse_args()

    _, points_tree = io.load_points(args.input, args.tree)
    fout, tracks, buf = io.open_output(args.output)

    cur_event = -1
    n_event = 0
    n_track = 0
    n_fit = 0
    n_bkg_fallback = 0
    chi2ndf = []
    chi2ndf_primary = []
    a_primary = []
    sig_a_primary = []
    bkgs = []
    for event, track_id, points in io.tracks(points_tree):
        if event != cur_event:      # io.tracks 保证同一事件连续
            cur_event = event
            n_event += 1
        n_track += 1
        result = fit.fit_track(event, track_id, points,
                               sigma_read=args.sigma_read,
                               diffusion=args.diffusion,
                               chip_y=args.chip_y,
                               background=args.background)
        if result is None:          # 点数不够 / 退化
            continue
        io.fill(tracks, buf, event, track_id, result)
        n_fit += 1
        if (result["bkg"] == 0.0
                and str(args.background).strip().lower() not in ("none", "off", "0")):
            n_bkg_fallback += 1
        bkgs.append(result["bkg"])
        if result["n"] >= 10:       # δ 射线只有几个点，不参与摘要统计
            chi2ndf.append(result["chi2ndf"])
        if track_id == 1:           # 初级电子：本场景的主径迹
            chi2ndf_primary.append(result["chi2ndf"])
            a_primary.append(result["A"])
            sig_a_primary.append(result["sigmaA"])

    tracks.Write()
    fout.Close()

    print(f"[track fit] {args.input} -> {args.output}")
    print(f"[track fit] 误差模型：小孔 sigma = {args.sigma_read * 1000:.0f} um, "
          f"D = {args.diffusion:.0f} um/sqrt(cm), 芯片 y = {args.chip_y} mm, "
          f"本底 = {args.background}")
    print(f"[track fit] 读了 {n_event} 个事件、{n_track} 条径迹，"
          f"写出 {n_fit} 条直线拟合")
    if chi2ndf:
        print(f"[track fit] n>=10 的径迹：chi2/ndf 中位数 = "
              f"{statistics.median(chi2ndf):.3f}")
    if chi2ndf_primary:
        print(f"[track fit] 主径迹 trackID=1：A 中位数 = "
              f"{statistics.median(a_primary):.6f} +- "
              f"{statistics.median(sig_a_primary):.2e}，chi2/ndf 中位数 = "
              f"{statistics.median(chi2ndf_primary):.3f}")
    if bkgs:
        print(f"[track fit] 实际本底中位数 = {statistics.median(bkgs):.3f}")
    if n_bkg_fallback == n_fit and n_fit > 0:
        print("[track fit] 提示：扣完本底后剩不到 3 个点（每电子一条时 ne 恒为 1，"
              "平均本底就是 1），已自动退回不扣本底", file=sys.stderr)
    if n_fit == 0:
        print("[track fit] 一条也没拟合出来：检查 points 树里有没有足够的点",
              file=sys.stderr)


if __name__ == "__main__":
    main()
