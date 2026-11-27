"""track fit 的读写层。

读 class 的输出（points 树：event / trackID / x / z / ne / ysrc），按
(event, trackID) 切组；把直线拟合结果写成 tracks 树。
算法在 fit.py。
"""
import array

import ROOT

# PyROOT 对象必须留引用，否则被 Python 回收后树指针就悬空了
_KEEP = []


def load_points(path, tree="points"):
    """打开 class 的输出文件，返回 (TFile, TTree)。"""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"打不开 {path}")
    t = f.Get(tree)
    if not t:
        raise SystemExit(f"{path} 里没有 {tree} 树")
    _KEEP.extend([f, t])
    return f, t


def tracks(tree):
    """按 (event, trackID) 切组，逐个产出 (event, trackID, [点 dict, ...])。

    同一个事件里不同 trackID 的行可能交错（digitize 按时间顺序写电子落点），
    所以这里用 dict 暂存，读完一个事件再统一吐出来。
    每个点带 x / z（mm）、ne（电子数，拟合权重）与 ysrc（产生点高度，
    用来算漂移扩散误差；旧文件没有这个分支时记为 None）。
    """
    n = int(tree.GetEntries())
    has_y = tree.GetBranch("ysrc") is not None
    cur_event = None
    groups = {}
    for i in range(n):
        tree.GetEntry(i)
        ev = int(tree.event)
        if ev != cur_event:
            if cur_event is not None:
                for tid in sorted(groups):
                    yield cur_event, tid, groups[tid]
            cur_event = ev
            groups = {}
        tid = int(tree.trackID)
        groups.setdefault(tid, []).append(
            dict(x=tree.x, z=tree.z, ne=int(tree.ne),
                 ysrc=float(tree.ysrc) if has_y else None))
    if cur_event is not None:
        for tid in sorted(groups):
            yield cur_event, tid, groups[tid]


def open_output(path):
    """建输出文件与 tracks 树，返回 (TFile, TTree, 缓冲区)。

    每条径迹一行，直线参数化成 x = A*z + B（芯片平面 X–Z）：
      event    事件号（编号 -> Int_t）
      trackID  径迹号（编号 -> Int_t）
      A        斜率（物理量 -> Double_t）
      B        截距 [mm]（物理量 -> Double_t）
      sigmaA   A 的误差（误差传递 -> Double_t）
      sigmaB   B 的误差 [mm]
      covAB    A、B 的协方差
      corrAB   A、B 的相关系数
      chi2     加权残差平方和
      ndf      自由度（点数 - 2，计数 -> Int_t）
      chi2ndf  chi2 / ndf
      pvalue   chi2 的 p 值（ROOT::TMath::Prob）
      rms      加权残差 RMS [mm]
      bkg      本底（每个点扣除的电荷数）
      n        参与拟合的点数（计数 -> Int_t）
    以后要换参数化（比如 3D 直线 / helix），改这里和 fill()。
    """
    f = ROOT.TFile.Open(path, "RECREATE")
    t = ROOT.TTree("tracks", "每条径迹一条：x = A*z + B")

    buf = {
        "event": array.array("i", [0]),
        "trackID": array.array("i", [0]),
        "A": array.array("d", [0.0]),
        "B": array.array("d", [0.0]),
        "sigmaA": array.array("d", [0.0]),
        "sigmaB": array.array("d", [0.0]),
        "covAB": array.array("d", [0.0]),
        "corrAB": array.array("d", [0.0]),
        "chi2": array.array("d", [0.0]),
        "ndf": array.array("i", [0]),
        "chi2ndf": array.array("d", [0.0]),
        "pvalue": array.array("d", [0.0]),
        "rms": array.array("d", [0.0]),
        "bkg": array.array("d", [0.0]),
        "n": array.array("i", [0]),
    }
    t.Branch("event", buf["event"], "event/I")
    t.Branch("trackID", buf["trackID"], "trackID/I")
    t.Branch("A", buf["A"], "A/D")
    t.Branch("B", buf["B"], "B/D")
    t.Branch("sigmaA", buf["sigmaA"], "sigmaA/D")
    t.Branch("sigmaB", buf["sigmaB"], "sigmaB/D")
    t.Branch("covAB", buf["covAB"], "covAB/D")
    t.Branch("corrAB", buf["corrAB"], "corrAB/D")
    t.Branch("chi2", buf["chi2"], "chi2/D")
    t.Branch("ndf", buf["ndf"], "ndf/I")
    t.Branch("chi2ndf", buf["chi2ndf"], "chi2ndf/D")
    t.Branch("pvalue", buf["pvalue"], "pvalue/D")
    t.Branch("rms", buf["rms"], "rms/D")
    t.Branch("bkg", buf["bkg"], "bkg/D")
    t.Branch("n", buf["n"], "n/I")

    _KEEP.extend([f, t, buf])
    return f, t, buf


def fill(tree, buf, event, track_id, result):
    """把 fit.fit_track() 返回的 dict 写进树里（缺的字段留 0）。"""
    buf["event"][0] = int(event)
    buf["trackID"][0] = int(track_id)
    for key in ("A", "B", "sigmaA", "sigmaB", "covAB", "corrAB",
                "chi2", "ndf", "chi2ndf", "pvalue", "rms", "bkg", "n"):
        if key in result:
            buf[key][0] = result[key]
    tree.Fill()
