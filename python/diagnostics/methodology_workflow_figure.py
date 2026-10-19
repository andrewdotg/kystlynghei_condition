"""Draw Methodology Figure 1 (../kystlynghei_word_documents/methodology_figures/methodology_workflow.png)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path(__file__).resolve().parents[3] / "kystlynghei_word_documents" / "methodology_figures" / "methodology_workflow.png"
INPUT_FC, INPUT_EC = "#F0EBE0", "#8C7B5A"
CORE_FC, OUT_FC = "#315F3C", "#1F3D29"
MOD_FC, MOD_EC = "#E3EDE4", "#315F3C"
INK, INK2, LINE = "#1F2420", "#5B6157", "#4A5450"

fig, ax = plt.subplots(figsize=(11.5, 8.6))
ax.set_xlim(0, 115)
ax.set_ylim(0, 86)
ax.axis("off")


def box(x, y, w, h, title, sub, fc, ec=None, light=False, tsize=11.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                fc=fc, ec=ec or fc, lw=1.6))
    tc, sc = ("white", "#E3EDE4") if light else (INK, INK2)
    ax.text(x + w / 2, y + h * 0.68, title, ha="center", va="center", fontsize=tsize, weight="bold", color=tc,
            linespacing=1.15)
    ax.text(x + w / 2, y + h * 0.30, sub, ha="center", va="center", fontsize=8.4, color=sc, linespacing=1.3)


def arrow(x0, y0, x1, y1, dashed=False):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.4, shrinkA=0, shrinkB=0,
                                linestyle=(0, (3, 2)) if dashed else "-", mutation_scale=12))


def line(xs, ys):
    ax.plot(xs, ys, color=LINE, lw=1.4, solid_capstyle="butt")


for y, lab in ((76, "INPUT DATA"), (55, "CORE MODEL"), (31, "ANALYSIS MODULES"), (7, "OUTPUT")):
    ax.text(2.5, y, lab, rotation=90, ha="center", va="center", fontsize=10, weight="bold", color=INK2)

# inputs
box(8, 70, 30, 12, "Field-mapped heathland\npolygons (MI)", "7RA-SJ stage · grazing class (BT)\ncondition (tilstand) · survey year",
    INPUT_FC, INPUT_EC, tsize=10.5)
box(44, 70, 30, 12, "AlphaEarth Foundations\nembeddings", "annual · 64 dimensions · 10 m\n2017–2025", INPUT_FC, INPUT_EC, tsize=10.5)
box(80, 70, 32, 12, "Fire and burn-scar data", "NASA FIRMS active fire\n(VIIRS, MODIS) · Sentinel-2 L2A",
    INPUT_FC, INPUT_EC, tsize=10.5)
# core
box(8, 49, 30, 11, "Central model", "covariance-adjusted AE-SJ axis\ncalibrated on 201,215 field samples", CORE_FC, light=True, tsize=10.5)
box(44, 49, 30, 11, "Annual per-pixel scoring", "AE-SJ score and SJ-like class\n14,330 polygons × 8 years", CORE_FC, light=True, tsize=10.5)
arrow(23, 70, 23, 60)
arrow(52, 70, 34, 60)
arrow(59, 70, 59, 60)
arrow(38, 54.5, 44, 54.5)

# modules
mods = [("Loss to\ndevelopment", "Dynamic World\nbuilt-area change"),
        ("Trajectories\nand triage", "Theil–Sen trend,\nfive classes"),
        ("Grazing\neffect", "field grazing\nclass (BT)"),
        ("Agreement with\nfield data", "field 7RA-SJ,\nsame year"),
        ("Scrub, tree and\nspruce clearance", "pixel transitions,\nspruce mask"),
        ("Burn detection\nand fire type", "FIRMS trigger,\nSentinel-2 dNBR\npatches")]
w, gap, x0, y0, h = 16.2, 1.4, 7, 22, 18
xs = [x0 + i * (w + gap) for i in range(len(mods))]
for x, (t, s) in zip(xs, mods):
    box(x, y0, w, h, t, s, MOD_FC, MOD_EC, tsize=9.3)
cx = [x + w / 2 for x in xs]
# AE-SJ scoring feeds every module except burn detection
line([59, 59], [49, 44])
line([cx[0], cx[4]], [44, 44])
for i in (0, 1, 2, 3, 4):
    arrow(cx[i], 44, cx[i], y0 + h)
# fire data feed burn detection directly, not through AE-SJ
line([96, 96], [70, 47])
line([cx[5], 96], [47, 47])
arrow(cx[5], 47, cx[5], y0 + h)
# results passed between modules
arrow(xs[0] + w, y0 + h / 2, xs[1], y0 + h / 2, dashed=True)
arrow(xs[1] + w, y0 + h / 2, xs[2], y0 + h / 2, dashed=True)
arrow(xs[5], y0 + h / 2, xs[4] + w, y0 + h / 2, dashed=True)
# burn effect: burn layer and AE-SJ scores
ax.text(cx[5], y0 - 2.4, "+ AE-SJ response\nof burnt pixels", ha="center", va="center", fontsize=8,
        color=INK2, style="italic")

# output
box(30, 2, 55, 10, "Unified site summary", "one record per polygon · table, polygon and point layers", OUT_FC, light=True)
line([cx[1], cx[5]], [17, 17])
for i in (1, 4, 5):
    line([cx[i], cx[i]], [y0 - (4.6 if i == 5 else 0), 17])
arrow(57.5, 17, 57.5, 12)

ax.plot([8, 13], [-1.2, -1.2], color=LINE, lw=1.4, clip_on=False)
ax.text(14, -1.2, "data flow", va="center", fontsize=9, color=INK2)
ax.plot([30, 35], [-1.2, -1.2], color=LINE, lw=1.4, ls=(0, (3, 2)), clip_on=False)
ax.text(36, -1.2, "results passed between modules (development exclusion, trend, burn mask)",
        va="center", fontsize=9, color=INK2)
fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
print(OUT)
