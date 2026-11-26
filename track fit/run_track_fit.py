#!/usr/bin/env python3
"""track fit 命令行入口：读 hits -> 拟合径迹 -> 写 tracks 树。

现在只有读写骨架：`track_fit/fit.py::fit_track()` 还返回 None，
所以跑出来的 tracks 树是空的（结构、分支类型都是最终形态）。
要真正出数，把 fit_track() 的 TODO 填掉即可，这里不用改。

用法:
    python3 "track fit/run_track_fit.py" pixel.root -o track_fit.root
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ROOT  # noqa: E402  （PyROOT 要在 sys.path 之后导入）

from track_fit import fit, io  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="hits -> 径迹拟合（骨架）")
    ap.add_argument("input", help="hits 文件（pixel.root 等）")
    ap.add_argument("-o", "--output", default="track_fit.root", help="输出文件")
    ap.add_argument("--tree", default="hits", help="输入树名（默认 hits）")
    args = ap.parse_args()

    _, hits = io.load_hits(args.input, args.tree)
    fout, tracks, buf = io.open_output(args.output)

    n_event = 0
    n_fit = 0
    for event, points in io.events(hits):
        n_event += 1
        result = fit.fit_track(event, points)
        if result is None:          # 还没实现 / 这个事件不拟合
            continue
        io.fill(tracks, buf, event, result)
        n_fit += 1

    tracks.Write()
    fout.Close()

    print(f"[track fit] {args.input} -> {args.output}")
    print(f"[track fit] 读了 {n_event} 个事件，写出 {n_fit} 条径迹")
    if n_fit == 0:
        print("[track fit] 一条也没写出来：fit_track() 还是 TODO", file=sys.stderr)


if __name__ == "__main__":
    main()
