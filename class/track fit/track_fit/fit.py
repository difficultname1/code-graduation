"""track fit 的算法层。

接口已经定好，具体过程留空。要填的就是 fit_track()。
"""


def fit_track(event, points):
    """对一个事件里的沉积点做径迹拟合。

    参数
    ----
    event : int
        事件号。
    points : list[dict]
        这个事件的所有沉积点，每项有 x / y / z（mm）与 edep（MeV）。
        行数、顺序就是 pixel 里 Geant4 记录的顺序（同一事件内按时间排列），
        但这里不保证按 z 排过序。

    返回
    ----
    dict 或 None
        返回 None 表示这个事件不做拟合 / 拟合失败，调用方会跳过。
        返回 dict 时，可以带的键（与 io.open_output() 的 tracks 分支对应）：
            nx, ny   方向余弦
            x0, y0   参考平面上的截距 [mm]
            chi2, ndf
    """
    # TODO: 具体的拟合过程
    #   1) 挑点：可选剔除 δ 射线的离群点
    #   2) 参数化：直线（或 helix）
    #   3) 最小化 chi2：可以用 ROOT.TGraph2D + Fit，也可以自己写迭代
    #   4) 填 chi2 / ndf
    return None
