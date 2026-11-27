"""class 的读写层。

读 digitize 的输出（electrons 树），写出按 trackID 分好组的点（points 树）。
真正的分类算法在 classify.py；本文件只负责搬运。
"""
import array

import ROOT

_KEEP = []          # PyROOT 对象留引用，避免被 GC 回收后指针悬空


def load_digi(path, tree="electrons"):
    """打开 digitize 的输出文件，返回 (TFile, TTree)。"""
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

    每个点带 mcTrackID / x / z / ne / ysrc：
      mcTrackID  digitize 从 Geant4 带下来的 MC 真值（初级 = 1），
                 不是分类结果——分类后的 trackID 由 class 负责产生
      x, z       芯片平面上的落点 [mm]
      ne         这个落点代表的电子数（每电子一条时为 1，按簇输出时为簇大小）
      ysrc       产生点的 y [mm]（算漂移扩散误差用；旧文件没有就写 NaN）
    同一个事件的行在 digitize 输出里是连续的。
    """
    n = int(tree.GetEntries())
    has_truth = tree.GetBranch("mcTrackID") is not None
    has_track = tree.GetBranch("trackID") is not None     # 兼容旧文件
    has_y = tree.GetBranch("ysrc") is not None
    cur, buf = None, []
    for i in range(n):
        tree.GetEntry(i)
        ev = int(tree.event)
        if ev != cur:
            if cur is not None:
                yield cur, buf
            cur, buf = ev, []
        if has_truth:
            truth = int(tree.mcTrackID)
        elif has_track:
            truth = int(tree.trackID)
        else:
            truth = 0
        buf.append(dict(mcTrackID=truth,
                        x=tree.x, z=tree.z, ne=int(tree.ne),
                        ysrc=float(tree.ysrc) if has_y else float("nan")))
    if cur is not None:
        yield cur, buf


def open_output(path):
    """建输出文件与 points 树，返回 (TFile, TTree, 缓冲区)。

    points 每个落点一条，字段就是下游径迹拟合需要的：
      event    事件号（编号 -> Int_t）
      trackID  分类后的径迹号（编号 -> Int_t；同一事件内唯一）
      x, z     芯片平面上的坐标 [mm]（物理量 -> Double_t）
      ne       电子数（计数 -> Int_t，拟合时当权重）
      ysrc     产生点的 y [mm]（物理量 -> Double_t，算漂移扩散误差用）
    """
    f = ROOT.TFile.Open(path, "RECREATE")
    t = ROOT.TTree("points", "每个落点一条：事件号 + 径迹号 + 坐标 + 电子数")

    buf = {
        "event": array.array("i", [0]),
        "trackID": array.array("i", [0]),
        "x": array.array("d", [0.0]),
        "z": array.array("d", [0.0]),
        "ne": array.array("i", [1]),
        "ysrc": array.array("d", [0.0]),
    }
    t.Branch("event", buf["event"], "event/I")
    t.Branch("trackID", buf["trackID"], "trackID/I")
    t.Branch("x", buf["x"], "x/D")
    t.Branch("z", buf["z"], "z/D")
    t.Branch("ne", buf["ne"], "ne/I")
    t.Branch("ysrc", buf["ysrc"], "ysrc/D")

    _KEEP.extend([f, t, buf])
    return f, t, buf


def fill(tree, buf, event, point):
    """把一个点写进 points 树。"""
    buf["event"][0] = int(event)
    buf["trackID"][0] = int(point["trackID"])
    buf["x"][0] = float(point["x"])
    buf["z"][0] = float(point["z"])
    buf["ne"][0] = int(point.get("ne", 1))
    buf["ysrc"][0] = float(point.get("ysrc", 0.0))
    tree.Fill()
