"""vertex fit 的算法层。

接口已经定好，具体过程留空。要填的就是 fit_vertex()。
"""


def fit_vertex(event, hits):
    """从一个事件的芯片落点重建入射位置（顶点）。

    参数
    ----
    event : int
        事件号。
    hits : list[dict]
        这个事件在芯片上的落点，每项有 x / z（mm）、ne（电子数）、
        ysrc（产生点的 y，mm）、t（ns）。

    返回
    ----
    dict 或 None
        None 表示这个事件不重建 / 重建失败，调用方跳过。
        返回 dict 时可用键（与 io.open_output() 的 vertices 分支对应）：
            x, y, z, t, n
    """
    # TODO: 具体的重建过程
    #   1) 选权重：按 ne 加权，或者先砍掉远离主簇的孤立点
    #   2) 位置：加权平均 / 中位数 / 拟合簇心
    #   3) 高度 y 与时间 t 从 ysrc / t 汇总
    #   4) n 填参与重建的落点数
    return None
