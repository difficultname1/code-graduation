"""track fit 的读写层。

只做 I/O：把 hits 读成一串事件，把拟合结果写成 tracks 树。
算法在 fit.py，这里把接口先定死，后面往里填实现就行。
"""
import array

import ROOT

# PyROOT 对象必须留引用，否则被 Python 回收后树指针就悬空了
_KEEP = []


def load_hits(path, tree="hits"):
    """打开 hits 文件，返回 (TFile, TTree)。"""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"打不开 {path}")
    t = f.Get(tree)
    if not t:
        raise SystemExit(f"{path} 里没有 {tree} 树")
    _KEEP.extend([f, t])
    return f, t


def events(tree):
    """按事件号切组，逐个产出 (event, [点 dict, ...])。

    要求同一个事件的行在树里是连续的（pixel 写出来的就是这样）。
    每个点带 x / y / z / edep，单位 mm 与 MeV。
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
        buf.append(dict(x=tree.x, y=tree.y, z=tree.z, edep=tree.edep))
    if cur is not None:
        yield cur, buf


def open_output(path):
    """建输出文件与 tracks 树，返回 (TFile, TTree, 缓冲区)。

    tracks 每个事件最多一条，字段先按"直线径迹"定好：
      event  事件号（编号 -> Int_t）
      nx,ny  方向余弦（物理量 -> Double_t）
      x0,y0  参考平面上的截距 [mm]
      chi2   拟合优度
      ndf    自由度（计数 -> Int_t）
    真要用别的参数化（比如 helix），改这里和 fill() 即可。
    """
    f = ROOT.TFile.Open(path, "RECREATE")
    t = ROOT.TTree("tracks", "每事件一条：径迹拟合结果")

    buf = {
        "event": array.array("i", [0]),
        "nx": array.array("d", [0.0]),
        "ny": array.array("d", [0.0]),
        "x0": array.array("d", [0.0]),
        "y0": array.array("d", [0.0]),
        "chi2": array.array("d", [0.0]),
        "ndf": array.array("i", [0]),
    }
    t.Branch("event", buf["event"], "event/I")
    t.Branch("nx", buf["nx"], "nx/D")
    t.Branch("ny", buf["ny"], "ny/D")
    t.Branch("x0", buf["x0"], "x0/D")
    t.Branch("y0", buf["y0"], "y0/D")
    t.Branch("chi2", buf["chi2"], "chi2/D")
    t.Branch("ndf", buf["ndf"], "ndf/I")

    _KEEP.extend([f, t, buf])
    return f, t, buf


def fill(tree, buf, event, result):
    """把 fit.fit_track() 返回的 dict 写进树里（缺的字段留 0）。"""
    buf["event"][0] = int(event)
    for key, var in buf.items():
        if key == "event" or key not in result:
            continue
        var[0] = result[key]
    tree.Fill()
