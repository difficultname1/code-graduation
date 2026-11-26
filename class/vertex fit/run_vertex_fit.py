#!/usr/bin/env python3
"""vertex fit 命令行入口：读 digi -> 重建入射位置 -> 写 vertices 树。

现在只有读写骨架：`vertex_fit/fit.py::fit_vertex()` 还返回 None，
跑出来的 vertices 树是空的（结构、分支类型已经是最终形态）。
要真正出数，把 fit_vertex() 的 TODO 填掉，这里不用改。

用法:
    python3 "vertex fit/run_vertex_fit.py" digi.root -o vertex_fit.root
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ROOT  # noqa: E402

from vertex_fit import fit, io  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="digi -> 入射位置重建（骨架）")
    ap.add_argument("input", help="digi 文件（digitize 的输出）")
    ap.add_argument("-o", "--output", default="vertex_fit.root", help="输出文件")
    ap.add_argument("--tree", default="digi", help="输入树名（默认 digi）")
    args = ap.parse_args()

    _, hits = io.load_digi(args.input, args.tree)
    fout, vertices, buf = io.open_output(args.output)

    n_event = 0
    n_ok = 0
    for event, points in io.events(hits):
        n_event += 1
        result = fit.fit_vertex(event, points)
        if result is None:          # 还没实现 / 这个事件不重建
            continue
        io.fill(vertices, buf, event, result)
        n_ok += 1

    vertices.Write()
    fout.Close()

    print(f"[vertex fit] {args.input} -> {args.output}")
    print(f"[vertex fit] 读了 {n_event} 个事件，写出 {n_ok} 个顶点")
    if n_ok == 0:
        print("[vertex fit] 一个也没写出来：fit_vertex() 还是 TODO", file=sys.stderr)


if __name__ == "__main__":
    main()
