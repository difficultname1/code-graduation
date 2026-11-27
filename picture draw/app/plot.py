#!/usr/bin/env python3
"""读 pixel.root / pixel_random.root 的 hits 树，按柱子归并求平均，画图并写 ROOT。

    python3 plot.py [root 文件] [输出前缀]

pixel.root        -> edep_map.png / edep_map_log.png / y_map.png / columns_mean.root
pixel_random.root -> random_ 开头的同一套文件（前缀自动从文件名推出来）

数据输出统一是 ROOT：columns_mean.root 里一棵 columns 树，每个被打到的柱子一行。
图片只是给人看的，不算数据产物。

树 hits：柱子里每个沉积了能量的步一条，
分支为 event(整数) / x,y,z(该步起点坐标, mm) / edep(该步沉积能量, MeV)。

文件里不存柱子编号，画图也直接用坐标：一根柱子一个 bin，bin 宽 PITCH、
范围 ±HALF，横纵轴就是 mm，全程不出现整数索引。
PITCH / NCOL 要和 DetectorConstruction.cc 里的 kSide / kN 一致。
"""
import array
import os
import sys
import ROOT

path = sys.argv[1] if len(sys.argv) > 1 else "pixel.root"

# 输出前缀：给了就用给的，否则从文件名推（pixel.root -> ""，pixel_random.root -> "random_"）
if len(sys.argv) > 2:
    prefix = sys.argv[2]
else:
    stem = os.path.basename(path).replace(".root", "")
    prefix = "" if stem == "pixel" else stem.replace("pixel_", "") + "_"

# 几何：默认对应 kSide = 83 um、kN = 1000，可用第 3、4 个参数覆盖
PITCH = float(sys.argv[3]) if len(sys.argv) > 3 else 0.083   # mm，柱子边长 = bin 宽
NCOL = int(sys.argv[4]) if len(sys.argv) > 4 else 1000       # 每边柱子数
HALF = PITCH * NCOL / 2.0                                     # 阵列半宽 [mm]

f = ROOT.TFile(path)
ROOT.gROOT.cd()          # 直方图建在内存里，不往只读文件里塞
tree = f.Get("hits")

n = NCOL
nEvents = int(tree.GetMaximum("event")) + 1     # 事件数（编号 0 起）
print("文件              :", path)
print("记录数（步）      :", tree.GetEntries())
print("事件数            :", nEvents)
print("每边柱子数 / bin宽: %d 根 / %.4f mm" % (n, PITCH))

# 一根柱子一个 bin，bin 的坐标范围就是物理坐标，轴标签直接是 mm
def newMap(name, title):
    return ROOT.TH2D(name, title, n, -HALF, HALF, n, -HALF, HALF)

hEdepSum = newMap("edep_sum", "sum edep;x [mm];z [mm]")     # [MeV]
hYSum    = newMap("y_sum",    "sum y;x [mm];z [mm]")        # [mm]
hCnt     = newMap("n_step",   "steps;x [mm];z [mm]")

# 四千多万条记录交给 ROOT 在 C++ 里扫，别在 Python 里循环
tree.Draw("z:x>>edep_sum", "edep", "goff")    # 权重 = edep
tree.Draw("z:x>>y_sum",    "y",    "goff")    # 权重 = y
tree.Draw("z:x>>n_step",   "1",    "goff")    # 只数步数

nStep = tree.GetEntries()
totalEdep = hEdepSum.Integral(0, n + 1, 0, n + 1)
nHitCol = sum(1 for bx in range(1, n + 1) for bz in range(1, n + 1)
              if hCnt.GetBinContent(bx, bz) > 0)

print("被打到过的柱子    :", nHitCol, "/ %d" % (n * n))
print("平均每事件步数    : %.0f" % (nStep / max(1, nEvents)))
print("平均每步沉积能量  : %.3f eV" % (totalEdep * 1e6 / max(1, nStep)))
print("总沉积能量        : %.4f MeV  (%.4f keV/事件)"
      % (totalEdep, totalEdep * 1e3 / max(1, nEvents)))

# ---- 每根柱子的事件平均值 ----
hEdep = hEdepSum.Clone("edep_mean")
hEdep.SetDirectory(0)
hEdep.Scale(1e6 / max(1, nEvents))              # MeV -> eV/事件
hEdep.SetTitle("mean edep per column (per event) [eV];x [mm];z [mm]")

hY = hYSum.Clone("y_mean")
hY.SetDirectory(0)
hY.Divide(hCnt)                                 # 每柱平均 y [mm]
hY.SetTitle("mean y per column;x [mm];z [mm];y [mm]")
hY.GetZaxis().SetTitleOffset(1.6)

ax, az = hEdepSum.GetXaxis(), hEdepSum.GetYaxis()
rows = []
for bx in range(1, n + 1):
    for bz in range(1, n + 1):
        cnt = hCnt.GetBinContent(bx, bz)
        if cnt <= 0:
            continue
        eSum = hEdepSum.GetBinContent(bx, bz)
        rows.append((ax.GetBinCenter(bx), az.GetBinCenter(bz), int(round(cnt)),
                     eSum, eSum / max(1, nEvents),
                     hYSum.GetBinContent(bx, bz) / cnt))

fout = ROOT.TFile(prefix + "columns_mean.root", "RECREATE")
tcols = ROOT.TTree("columns", "per-column means: one row per hit column")
bx = array.array("d", [0.0])
bz = array.array("d", [0.0])
bn = array.array("i", [0])
be = array.array("d", [0.0])
bm = array.array("d", [0.0])
by = array.array("d", [0.0])
tcols.Branch("x_mm", bx, "x_mm/D")
tcols.Branch("z_mm", bz, "z_mm/D")
tcols.Branch("nStep", bn, "nStep/I")
tcols.Branch("edep_sum_MeV", be, "edep_sum_MeV/D")
tcols.Branch("edep_mean_MeV", bm, "edep_mean_MeV/D")
tcols.Branch("y_mean_mm", by, "y_mean_mm/D")
for r in rows:
    bx[0], bz[0], bn[0], be[0], bm[0], by[0] = r
    tcols.Fill()
tcols.Write()
fout.Close()
ROOT.gROOT.cd()          # 回到内存目录，后面的画布/直方图不依赖已关闭的文件

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleSize(0.045, "xyz")
ROOT.gStyle.SetLabelSize(0.040, "xyz")
ROOT.gStyle.SetTitleOffset(1.1, "x")
ROOT.gStyle.SetTitleOffset(1.6, "y")
ROOT.gStyle.SetTitleOffset(2.6, "z")

W, H = 900, 760

def newCanvas(name):
    c = ROOT.TCanvas(name, "", W, H)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.20)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.09)
    return c

c1 = newCanvas("c1")
hEdep.Draw("COLZ")
c1.SaveAs(prefix + "edep_map.png")

# 同一张图改成对数色标：主电子经过的那两列（~9 eV）和 δ 射线偶然打到的列（~0.1 eV）
# 能同时看清
c1b = newCanvas("c1b")
c1b.SetLogz(1)
c1b.SetGrid()
hEdep.Draw("COLZ")
c1b.SaveAs(prefix + "edep_map_log.png")

c2 = newCanvas("c2")
hY.Draw("COLZ")
c2.SaveAs(prefix + "y_map.png")

print("已生成 %sedep_map.png / %sedep_map_log.png / %sy_map.png / %scolumns_mean.root"
      % (prefix, prefix, prefix, prefix))
