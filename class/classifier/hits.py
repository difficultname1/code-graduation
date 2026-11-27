"""把逐电子落点合并成像素命中（hit）。

digitize 默认每个电子一条。真实像素读出是把落在同一个像素里的电荷加在一起，
所以这里按 (trackID, ix, iz) 分组，输出：

    x, z    该像素内电子的电荷加权质心 [mm]（不是像素中心，保留亚像素信息）
    ne      该像素收集到的电子总数
    ysrc    这些电子的平均产生高度 [mm]

分类仍然假设完全正确：分组键里的 trackID 是 class 已经产生的分类结果
（当前阶段由 MC 真值 mcTrackID 生成，见 classify.perfect()）。
"""
import math


def pixel_hits(points, pitch):
    """按 pitch [mm] 的像素把点合并；pitch <= 0 时原样返回。"""
    if not (pitch > 0.0) or not points:
        return points

    bins = {}
    for p in points:
        x = float(p["x"])
        z = float(p["z"])
        q = float(p.get("ne", 1))
        if not (q > 0.0):
            continue
        key = (int(p["trackID"]),
               int(math.floor(x / pitch)),
               int(math.floor(z / pitch)))
        acc = bins.get(key)
        if acc is None:
            bins[key] = [q, q * x, q * z, q * float(p.get("ysrc", float("nan")))]
        else:
            acc[0] += q
            acc[1] += q * x
            acc[2] += q * z
            acc[3] += q * float(p.get("ysrc", float("nan")))

    out = []
    for (tid, _ix, _iz), (q, sx, sz, sy) in bins.items():
        out.append(dict(trackID=tid, x=sx / q, z=sz / q,
                        ne=int(round(q)), ysrc=sy / q))
    return out
