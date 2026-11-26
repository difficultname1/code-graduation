"""vertex fit 的读写层。

只做 I/O：把 digi 读成一串事件（芯片上的落点），把重建结果写成 vertices 树。
算法在 fit.py，这里把接口先定死。
"""
import array

import ROOT

_KEEP = []          # PyROOT 对象留引用，避免被 GC 回收后指针悬空


def load_digi(path, tree="digi"):
    """打开 digi 文件，返回 (TFile, TTree)。"""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"打不开 {path}")
    t = f.Get(tree)
    if not t:
        raise SystemExit(f"{path} 里没有 {tree} 树")
    _KEEP.extend([f, t])
    return f, t


def events(tree):
    """按事件号切组，逐个产出 (event, [落点 dict, ...])。

    每个落点带 x / z（芯片上的 mm）、ne（电子数）、ysrc（产生点的 y，mm）、
    t（产生时间，ns）。行是连续的，和 digitize 写出来的一致。
    """
    n = int(tree.GetEntries())
    cur, buf = None, []
    for i in range(n):
        tree.GetEntry(i)
        ev = int(tree.event)
        if ev != cur:
            if cur is not None:
                yield cur, buf
            cur, buf = ev, []
        buf.append(dict(x=tree.x, z=tree.z, ne=int(tree.ne),
                        ysrc=tree.ysrc, t=tree.t))
    if cur is not None:
        yield cur, buf


def open_output(path):
    """建输出文件与 vertices 树，返回 (TFile, TTree, 缓冲区)。

    分支先按"一个事件一个顶点"定：
      event  事件号（编号 -> Int_t）
      x,z    顶点在芯片平面上的位置 [mm]（物理量 -> Double_t）
      y      产生点的平均高度 [mm]
      t      顶点时间 [ns]
      n      参与重建的落点数（计数 -> Int_t）
    """
    f = ROOT.TFile.Open(path, "RECREATE")
    t = ROOT.TTree("vertices", "每事件一条：重建出的入射位置")

    buf = {
        "event": array.array("i", [0]),
        "x": array.array("d", [0.0]),
        "y": array.array("d", [0.0]),
        "z": array.array("d", [0.0]),
        "t": array.array("d", [0.0]),
        "n": array.array("i", [0]),
    }
    t.Branch("event", buf["event"], "event/I")
    t.Branch("x", buf["x"], "x/D")
    t.Branch("y", buf["y"], "y/D")
    t.Branch("z", buf["z"], "z/D")
    t.Branch("t", buf["t"], "t/D")
    t.Branch("n", buf["n"], "n/I")

    _KEEP.extend([f, t, buf])
    return f, t, buf


def fill(tree, buf, event, result):
    """把 fit.fit_vertex() 返回的 dict 写进树里（缺的字段留 0）。"""
    buf["event"][0] = int(event)
    for key, var in buf.items():
        if key == "event" or key not in result:
            continue
        var[0] = result[key]
    tree.Fill()
