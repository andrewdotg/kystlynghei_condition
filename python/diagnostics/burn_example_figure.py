"""Report figure of one mapped burn: Sentinel-2 true colour before and after
the fire with the burn outlined, and a third panel showing what followed.

usage: python diagnostics/burn_example_figure.py FID YEAR [recovery|aesj] [OUT.png]
  recovery  third panel: the same season in RECOVERY_YEAR (Sentinel-2)
  aesj      third panel: AE-SJ change across the site, first to last AlphaEarth year
Reads python/mi_dnbr_burn_results/ (run mi_dnbr_burn.py map first).
"""
import asyncio
import sys
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio.features
import rasterio.warp
from matplotlib.colors import TwoSlopeNorm
from odc.geo.geobox import GeoBox

import mi_dnbr_burn as burn
from aef_geobox import aoi_geobox

RECOVERY_YEAR = 2025
AESJ_FROM, AESJ_TO = 2018, 2025
BURN_EDGE = "#F05A1A"          # red-orange outline
SITE_EDGE = "white"
PAD_M = 200
INK, INK2 = "#22261f", "#5b6157"

fid, year = int(sys.argv[1]), int(sys.argv[2])
mode = sys.argv[3] if len(sys.argv) > 3 else "recovery"
out = Path(sys.argv[4]) if len(sys.argv) > 4 else burn.REPO / "report" / f"fig_burn_example_{fid}_{year}_{mode}.png"

R = burn.OUTPUT_DIR
se = pd.read_csv(R / "MI_burn_site_events.csv")
eid = se[(se.fid == fid) & (se.year == year)].sort_values("burn_area_ha").iloc[-1]["event_id"]
ev = pd.read_csv(R / "MI_firms_events.csv").set_index("event_id").loc[eid]
L, tf = burn.read_layers(R / "events" / f"{eid}.tif")
shape = L["burn"].shape
target = GeoBox(shape, tf, burn.MI_CRS)
site = burn.load_mi().set_index("fid").geometry.loc[fid]

# view: the site plus a margin, on the event grid
x0, y0, x1, y1 = site.bounds
c0, c1 = max(int((x0 - PAD_M - tf.c) // 10), 0), min(int((x1 + PAD_M - tf.c) // 10), shape[1])
r0, r1 = max(int((tf.f - (y1 + PAD_M)) // 10), 0), min(int((tf.f - (y0 - PAD_M)) // 10), shape[0])
crop = (slice(r0, r1), slice(c0, c1))
ext = [tf.c + c0 * 10, tf.c + c1 * 10, tf.f - r1 * 10, tf.f - r0 * 10]
poly = rasterio.features.geometry_mask([site], out_shape=shape, transform=tf, invert=True)[crop]
burnt = (L["burn"][crop] > 0) & poly

cat = burn.open_catalog()
bbox_ll = rasterio.warp.transform_bounds(burn.MI_CRS, "EPSG:4326", *target.boundingbox)


def read_rgb(items):
    return {b: burn.read_mosaic(items, b, target)[crop] for b in ("B02", "B03", "B04", "SCL")}


def scene_on(d):
    g = [g for g in burn.search_by_date(cat, bbox_ll, d, d + timedelta(days=1)) if g[0].datetime.date() == d]
    return read_rgb(g[0]), d


def clearest_near(d, days=40):
    """Least cloudy scene over the site within +-days of d."""
    best, best_d, best_clear = None, None, -1
    for g in burn.search_by_date(cat, bbox_ll, d - timedelta(days=days), d + timedelta(days=days)):
        scl = burn.read_mosaic(g, "SCL", target)[crop]
        clear = np.mean(np.isin(scl[poly], [4, 5]))
        if clear > best_clear:
            best, best_d, best_clear = g, g[0].datetime.date(), clear
    return read_rgb(best), best_d


pre, pre_d = scene_on(date.fromisoformat(str(ev["fig_pre_date"])))
post, post_d = scene_on(date.fromisoformat(str(ev["fig_post_date"])))
scenes = [pre, post]
if mode == "recovery":
    rec, rec_d = clearest_near(post_d.replace(year=RECOVERY_YEAR))
    scenes.append(rec)

# one stretch for all true-colour panels and all three bands, so the colour
# balance is kept: 1st-99.5th percentile over the whole view, with a gamma that
# lifts the mid-tones and a slight desaturation for a softer, natural look
v = np.concatenate([s[b][np.isfinite(s[b])] for s in scenes for b in ("B04", "B03", "B02")])
LO, HI = np.percentile(v, [1, 99.5])
GAMMA, SAT = 0.75, 0.85


def rgb(s):
    x = np.dstack([np.clip((s[b] - LO) / (HI - LO), 0, 1) for b in ("B04", "B03", "B02")]) ** GAMMA
    grey = x.mean(axis=2, keepdims=True)
    return np.nan_to_num(np.clip(grey + SAT * (x - grey), 0, 1), nan=1)


plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5, "text.color": INK})
fig, axs = plt.subplots(1, 3, figsize=(11, 4.3), constrained_layout=True)


def outline(ax, burn_line=True):
    ax.contour(poly.astype(float), levels=[0.5], colors=SITE_EDGE, linewidths=1.0, extent=ext, origin="upper")
    if burn_line:
        ax.contour(burnt.astype(float), levels=[0.5], colors=BURN_EDGE, linewidths=1.6, extent=ext, origin="upper")
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


axs[0].imshow(rgb(pre), extent=ext); outline(axs[0], burn_line=False)
axs[0].set_title(f"Before  ·  {pre_d:%d %B %Y}".replace(" 0", " "), loc="left")
axs[1].imshow(rgb(post), extent=ext); outline(axs[1])
axs[1].set_title(f"After  ·  {post_d:%d %B %Y}".replace(" 0", " "), loc="left")

if mode == "recovery":
    axs[2].imshow(rgb(rec), extent=ext); outline(axs[2])
    axs[2].set_title(f"Six years later  ·  {rec_d:%d %B %Y}".replace(" 0", " ") if rec_d.year - year == 6
                     else f"{rec_d.year - year} years later  ·  {rec_d:%d %B %Y}".replace(" 0", " "), loc="left")
else:
    import mi_aesj_annual_products as aesj
    from aef_loader import AEFIndex, DataSource, VirtualTiffReader

    st = aoi_geobox(site.bounds, crs=burn.MI_CRS, resolution=10.0, bbox_crs=burn.MI_CRS)

    async def load():
        index = AEFIndex(source=DataSource.SOURCE_COOP)
        await index.download()
        index.load()
        async with VirtualTiffReader() as reader:
            return [await aesj._load_year_scores(reader, index, site.bounds, y, st) for y in (AESJ_FROM, AESJ_TO)]

    a0, a1 = asyncio.run(load())
    spoly = rasterio.features.geometry_mask([site], out_shape=st.shape, transform=st.affine, invert=True)
    change = np.where(spoly, a1 - a0, np.nan)
    sext = [st.affine.c, st.affine.c + st.shape[1] * 10, st.affine.f - st.shape[0] * 10, st.affine.f]
    axs[2].set_facecolor("#F2F2EE")
    im = axs[2].imshow(change, cmap="PuOr", norm=TwoSlopeNorm(0, -0.5, 0.5), extent=sext, interpolation="nearest")
    outline(axs[2])
    axs[2].contour(poly.astype(float), levels=[0.5], colors=INK2, linewidths=0.6, extent=ext, origin="upper")
    axs[2].set_title(f"Change in AE-SJ  ·  {AESJ_FROM} to {AESJ_TO}", loc="left")
    cb = fig.colorbar(im, ax=axs[2], shrink=0.8, pad=0.02)
    cb.set_label("AE-SJ change (+ = more open, SJ1-like)", fontsize=8.5, color=INK2)
    cb.outline.set_visible(False)

# scale bar, top right of the first panel
w, h = ext[1] - ext[0], ext[3] - ext[2]
L_m = 200 if w < 1500 else 500 if w < 4000 else 1000
xe, ys = ext[1] - 0.05 * w, ext[3] - 0.08 * h
halo = [pe.withStroke(linewidth=2.5, foreground=INK)]
axs[0].plot([xe - L_m, xe], [ys, ys], color="white", lw=3.5, solid_capstyle="butt",
            path_effects=[pe.Stroke(linewidth=5, foreground=INK), pe.Normal()])
axs[0].text(xe - L_m / 2, ys - 0.025 * h, f"{L_m} m", ha="center", va="top", color="white", weight="bold",
            path_effects=halo)
fig.savefig(out, dpi=180, facecolor="white")
print(out, f"burn {burnt.sum() / 100:.1f} ha, event {eid}")
