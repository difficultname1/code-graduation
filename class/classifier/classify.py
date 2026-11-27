"""class 的算法层。

当前阶段假设分类**完全正确**：直接采用 digitize 从 Geant4 带下来的
MC 真值 mcTrackID 生成 trackID，所以函数体可以是空的。
"""


def classify(event, points):
    """给一个事件里的落点分配 trackID。

    参数
    ----
    event : int
        事件号。
    points : list[dict]
        这个事件的所有落点，每项有 mcTrackID / x / z / ne。
        mcTrackID 是 MC 真值（初级 = 1，次级 > 1），不是分类结果。

    返回
    ----
    list[dict] 或 None
        None 表示"还没有分类器"——调用方改走 perfect() 的真值捷径，
        用 mcTrackID 生成 trackID。这是当前阶段（完美分类）的行为。
        以后换成真实分类器时，返回分组后的点列表（每项同样带
        trackID / x / z / ne），或者直接原地修改 points 再返回。
    """
    # 真实分类器的挂载点。当前返回 None，由 perfect() 走真值捷径。
    return None


def perfect(event, points):
    """临时捷径：假设分类完全正确，直接把 MC 真值 mcTrackID 当作 trackID。

    正常流程里 trackID 是 class 的输出；这里只是为了跳过分类器，
    先用 Geant4 的真值顶上。换成真实分类器时删掉这个调用即可。
    """
    out = []
    for p in points:
        out.append(dict(trackID=int(p.get("mcTrackID", 0)),
                        x=p["x"], z=p["z"], ne=int(p.get("ne", 1)),
                        ysrc=p.get("ysrc", float("nan"))))
    return out
