#!/usr/bin/env python3
"""class 命令行入口：读 digi -> 按径迹分类 -> 合并成像素命中 -> 写 points 树。

本阶段假设分类完全正确：classify.classify() 返回 None，
于是用 digitize 记录的 MC 真值 mcTrackID 生成分类后的 trackID。
上游的 electrons 树里没有 trackID，只有真值 mcTrackID。
默认一个电子一条（points 里就是带电子的 trackID）；
加 --pitch 0.083 则按 83 µm 的像素把电子合并成命中，
这样每个点的 ne 才是该像素收集到的电荷，下游可以按电荷加权。

用法:
    python3 "class/app/run_class.py" digi.root -o class_points.root
    python3 "class/app/run_class.py" digi.root -o points.root --pitch 0  # 不合并
"""
import argparse
import os
import sys

# 本文件在 class/app/ 下，包（classifier/）在上一级，所以把上级目录加进 sys.path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import ROOT  # noqa: E402

from classifier import classify, hits, io  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="digi -> 按径迹分类的落点")
    ap.add_argument("input", help="digitize 的输出文件（digi.root 等）")
    ap.add_argument("-o", "--output", default="class_points.root", help="输出文件")
    ap.add_argument("--tree", default="electrons", help="输入树名（默认 electrons）")
    ap.add_argument("--pitch", type=float, default=0.0,
                    help="像素间距 [mm]；默认 0 = 不合并（逐电子），0.083 = 按像素合并")
    args = ap.parse_args()

    _, digi = io.load_digi(args.input, args.tree)
    fout, points, buf = io.open_output(args.output)

    n_event = 0
    n_raw = 0
    n_point = 0
    for event, rows in io.events(digi):
        n_event += 1
        n_raw += len(rows)
        result = classify.classify(event, rows)
        if result is None:          # 跳过分类器：完美分类，用 MC 真值当 trackID
            result = classify.perfect(event, rows)
        result = hits.pixel_hits(result, args.pitch)
        for p in result:
            io.fill(points, buf, event, p)
            n_point += 1

    points.Write()
    fout.Close()

    print(f"[class] {args.input} -> {args.output}")
    if args.pitch > 0.0:
        print(f"[class] 读了 {n_event} 个事件、{n_raw} 个电子落点，"
              f"按 {args.pitch:g} mm 像素合并成 {n_point} 个命中")
    else:
        print(f"[class] 读了 {n_event} 个事件、{n_raw} 个电子落点，"
              f"逐电子输出 {n_point} 个点（未合并）")


if __name__ == "__main__":
    main()
