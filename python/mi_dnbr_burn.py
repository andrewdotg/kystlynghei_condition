"""mi_dnbr_burn.py - burn detection for MI kystlynghei sites: FIRMS-triggered,
Sentinel-2 dNBR-mapped.

A burn is a patch of Sentinel-2 pixels whose NBR dropped between the same
season of the previous year and the weeks after a FIRMS active-fire
detection, lying within reach of that detection. FIRMS says that and when a
fire happened; dNBR says where. The AE-SJ succession score plays no part in
detection, so its response to burning can be measured on the burns
afterwards (stage 2).

STAGE 1 - MAP BURNS (`python mi_dnbr_burn.py map`)
  1. Events. FIRMS detections (all four products, from the project's
     active-fire input, active_fires.py) linked
     to an MI polygon - within 750 m (VIIRS) or 1,000 m (MODIS) - are grouped into
     fire events: same year, within EVENT_LINK_M and EVENT_GAP_DAYS of each
     other. Each event is mapped once, over an area covering its detections
     (buffered by their reach) and every linked MI polygon.
  2. Pre-fire: per-pixel median NBR and NDVI of clear observations in Y-1,
     within PRE_HALF_WINDOW_DAYS of the event's first detection date (same
     season, so green-up is not read as a burn); widened once if there are
     fewer than MIN_OBS. Post-fire: median from the day after the first
     detection to POST_WINDOW_DAYS later (early, before regrowth); widened
     once if needed. Clear = SCL 4/5, eroded by EROSION_PX. All tiles of an
     acquisition date are mosaicked.
  3. dNBR and dNDVI are corrected by the unburnt mode of the MI heath
     pixels in the area (histogram peak over the less green 80 %, refined by
     the median within +-DNBR_GROW), which removes green-up and differences
     between years. Pixels brighter in B02 than the BRIGHTNESS_PERCENTILE of
     the heath pixels' mean over April Y-1 are excluded.
  4. Candidate pixel: dNBR_corr >= DNBR_GROW, dNDVI_corr >= DNDVI_MIN,
     pre-fire NDVI >= PRE_NDVI_MIN, post-fire NBR <= POST_NBR_MAX.
  5. Patches: a connected region of candidates is kept only if it holds at
     least MIN_CORE_PX pixels with dNBR_corr >= DNBR_CORE (grows weak burn
     margins from a strong core); gaps of CLOSING_PX are closed to group
     patchy burns; patches need MIN_PATCH_HA of flagged pixels.
  6. A patch is a burn when some pixel lies within reach of one of the
     event's detections (750 m VIIRS, 1,000 m MODIS). A site burnt in year Y
     has >= MIN_SITE_BURN_HA of burn pixels inside its polygon.
  7. Fire type: an event with any detection in June-August is a probable
     wildfire; otherwise it falls in the lyngbrenning season and is a
     probable prescribed burn. A site-year takes the wildfire type if any
     event that burnt it was a wildfire.

RE-CLASSIFY (`python mi_dnbr_burn.py reclassify`): steps 4-7 again from the
  saved event layers, after a threshold change (minutes instead of hours).

STAGE 2 - SUCCESSION RESPONSE (`python mi_dnbr_burn.py aesj`)
  For every burnt site, AE-SJ 2017-2025 is scored on the site's own grid.
  For each burn, the burnt pixels are compared with the site's pixels that
  never burnt (in any year): mean AE-SJ and SJ-class composition per year,
  written per burn and year relative to the burn.

OUTPUTS (OUTPUT_DIR)
  MI_firms_detections.csv         linked FIRMS detections with their event
  MI_firms_events.csv             one row per fire event (dates, sensors, result)
  MI_dnbr_patches.csv             one row per dNBR patch (area, dNBR, reach)
  MI_burn_site_events.csv         burnt area per MI site and event
  MI_burns.csv                    one row per burnt site-year (area, type, dates)
  MI_burn_aesj_trajectories.csv   stage 2
  events/<event_id>.tif           per-event layers (dNBR, masks, patches)
  site_events/<fid>_<year>_<event>.tif   each event's burn cropped to a site
  rasters/<fid>_<year>_burn.tif   site burn mask on the site's AE-SJ grid
                                  (read by mi_aesj_annual_products.py)
"""

from __future__ import annotations

import asyncio
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import planetary_computer
import pystac_client
import rasterio
import rasterio.features
import rasterio.warp
import scipy.ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.geometry import box

import active_fires
from aef_geobox import aoi_geobox

# =============================================================================
# 1. SETTINGS
# =============================================================================

REPO = Path(__file__).resolve().parent.parent
MI_SHAPEFILE = REPO / "inputs" / "MI_kystlynghei.shp"
MI_CRS = "EPSG:25833"
OUTPUT_DIR = Path(__file__).parent / "mi_dnbr_burn_results"

RESOLUTION_M = 10.0
PIXEL_AREA_HA = RESOLUTION_M**2 / 10000.0

# Fire events
EVENT_LINK_M = 2000.0
EVENT_GAP_DAYS = 30

# Sentinel-2 (Planetary Computer)
PC_STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
PC_COLLECTION = "sentinel-2-l2a"
MAX_CLOUD_COVER_PCT = 70
S2_HARMONIZE_CUTOFF = datetime(2022, 1, 25, tzinfo=timezone.utc)  # baseline 04.00 -1000 DN offset
S2_HARMONIZE_OFFSET = 1000.0
PRE_HALF_WINDOW_DAYS, PRE_HALF_WINDOW_WIDE_DAYS = 45, 75
POST_WINDOW_DAYS, POST_WINDOW_WIDE_DAYS = 60, 90
MIN_OBS = 2
EROSION_PX = 2
BRIGHTNESS_MONTH = 4
BRIGHTNESS_PERCENTILE = 95
FIGURE_BANDS = ["B02", "B03", "B04", "B8A", "B12"]

# Burn classification, chosen against registered burns
# (validation_registered_burns.py, threshold grid; Methodology Section 11.4,
# Table 9). The thresholds lie below standard burn-severity bounds (0.10 /
# 0.27), because prescribed heath burns are often of low severity. dNDVI only
# rejects pixels that got greener; the vegetation floor sits above bare rock
# (NDVI ~0.05-0.15) but below dormant spring heather.
TREE_NDVI_PERCENTILE = 80
DNBR_GROW = 0.075
DNBR_CORE = 0.20
DNDVI_MIN = 0.0
PRE_NDVI_MIN = 0.15
POST_NBR_MAX = 0.30
MIN_CORE_PX = 4
CLOSING_PX = 1
MIN_PATCH_HA = 0.25
MIN_SITE_BURN_HA = 0.1

# AE-SJ response (stage 2)
AEF_YEARS = list(range(2017, 2026))

N_WORKERS = 2          # Planetary Computer returned 403s at 4
MAX_RETRIES = 3
RETRY_DELAY_SECONDS = 30
READ_WORKERS = 8


# =============================================================================
# 2. SENTINEL-2 ACCESS
# =============================================================================

def open_catalog():
    return pystac_client.Client.open(PC_STAC_URL, modifier=planetary_computer.sign_inplace)


def search_by_date(catalog, bbox_ll, start: date, end: date) -> list[list]:
    """All tiles per acquisition date, least cloudy first. Keeping every tile
    matters: an area crossing an MGRS tile edge is otherwise empty on part of
    every date (found on fid 298)."""
    search = catalog.search(collections=[PC_COLLECTION], bbox=bbox_ll,
                            datetime=f"{start.isoformat()}/{end.isoformat()}",
                            query={"eo:cloud_cover": {"lt": MAX_CLOUD_COVER_PCT}})
    by_date: dict[date, list] = {}
    for it in search.item_collection():
        by_date.setdefault(it.datetime.date(), []).append(it)
    return [sorted(v, key=lambda it: it.properties.get("eo:cloud_cover", 100))
            for _, v in sorted(by_date.items())]


def read_band(item, band: str, target) -> np.ndarray:
    """One band reprojected (nearest) onto `target`; reflectance, or raw SCL
    codes; NaN where the source is nodata. Baseline 04.00 scenes (from 25 Jan
    2022) get the -1000 DN offset removed."""
    with rasterio.open(item.assets[band].href) as src:
        out = np.full(target.shape, np.nan, dtype=np.float64)
        rasterio.warp.reproject(source=rasterio.band(src, 1), destination=out,
                                dst_transform=target.affine, dst_crs=str(target.crs),
                                resampling=rasterio.warp.Resampling.nearest,
                                src_nodata=0, dst_nodata=np.nan)
    if band == "SCL":
        return out
    if item.datetime >= S2_HARMONIZE_CUTOFF:
        f = np.isfinite(out)
        out[f] = np.clip(out[f], S2_HARMONIZE_OFFSET, None) - S2_HARMONIZE_OFFSET
    return out / 10000.0


def read_mosaic(items, band: str, target) -> np.ndarray:
    out = read_band(items[0], band, target)
    for it in items[1:]:
        gap = ~np.isfinite(out)
        if not gap.any():
            break
        out[gap] = read_band(it, band, target)[gap]
    return out


def read_scene(items, target, heath_mask) -> dict | None:
    scl = read_mosaic(items, "SCL", target)
    clear = (scl == 4) | (scl == 5)
    if EROSION_PX:
        clear = ndi.binary_erosion(clear, structure=np.ones((3, 3)), iterations=EROSION_PX)
    if not (clear & heath_mask).any():
        return None
    b = {k: read_mosaic(items, k, target) for k in ("B02", "B04", "B08", "B12")}
    with np.errstate(invalid="ignore", divide="ignore"):
        nbr = (b["B08"] - b["B12"]) / (b["B08"] + b["B12"])
        ndvi = (b["B08"] - b["B04"]) / (b["B08"] + b["B04"])
    ok = clear & np.isfinite(nbr) & np.isfinite(ndvi)
    return {"items": items, "date": items[0].datetime.date(),
            "clear_frac": float((clear & heath_mask).sum() / max(heath_mask.sum(), 1)),
            "nbr": np.where(ok, nbr, np.nan).astype(np.float32),
            "ndvi": np.where(ok, ndvi, np.nan).astype(np.float32),
            "b02": np.where(ok, b["B02"], np.nan).astype(np.float32)}


def read_scenes(catalog, bbox_ll, start, end, target, heath_mask) -> list[dict]:
    groups = search_by_date(catalog, bbox_ll, start, end)
    with ThreadPoolExecutor(READ_WORKERS) as pool:
        scenes = list(pool.map(lambda g: read_scene(g, target, heath_mask), groups))
    return [s for s in scenes if s is not None]


def nanmedian_stack(scenes, key):
    stack = np.stack([s[key] for s in scenes])
    count = np.isfinite(stack).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(stack, axis=0)
    return np.where(count >= MIN_OBS, med, np.nan), count


# =============================================================================
# 3. BURN MAPPING CORE (shared with validation_registered_burns.py)
# =============================================================================

def unburnt_offset(delta, pre_ndvi, mask) -> float:
    """Unburnt mode of a change layer over `mask`: histogram peak over the
    less green pixels, refined by the median within +-DNBR_GROW."""
    sel = mask & np.isfinite(delta) & np.isfinite(pre_ndvi)
    if sel.sum() < 50:
        return 0.0
    sel &= pre_ndvi <= np.nanpercentile(pre_ndvi[sel], TREE_NDVI_PERCENTILE)
    v = delta[sel]
    hist, edges = np.histogram(v, bins=np.arange(-1.0, 1.0001, 0.01))
    k = int(np.argmax(ndi.uniform_filter1d(hist.astype(float), 5)))
    mode = 0.5 * (edges[k] + edges[k + 1])
    near = v[np.abs(v - mode) < DNBR_GROW]
    return float(np.median(near)) if near.size else float(mode)


def burn_patches(candidate, core):
    """Label kept patches 1..k; returns (labels, kept mask)."""
    eight = np.ones((3, 3), dtype=bool)
    lab, n = ndi.label(candidate, structure=eight)
    if n == 0:
        return np.zeros(candidate.shape, np.int32), np.zeros(candidate.shape, bool)
    n_core = ndi.sum(core, lab, index=np.arange(n + 1))
    keep = n_core >= MIN_CORE_PX
    keep[0] = False
    grown = keep[lab]
    pad = CLOSING_PX + 1
    closed = ndi.binary_closing(np.pad(grown, pad), structure=eight, iterations=CLOSING_PX)[pad:-pad, pad:-pad]
    plab, pn = ndi.label(closed, structure=eight)
    flagged = np.where(grown, plab, 0)
    keep_p = np.bincount(flagged.ravel(), minlength=pn + 1) * PIXEL_AREA_HA >= MIN_PATCH_HA
    keep_p[0] = False
    remap = np.zeros(pn + 1, np.int32)
    remap[np.flatnonzero(keep_p)] = np.arange(1, keep_p.sum() + 1)
    final = keep_p[flagged] & grown
    return np.where(final, remap[plab], 0), final


def classify_layers(dnbr_corr, dndvi_corr, pre_ndvi, post_nbr, bright) -> np.ndarray:
    """Burn-patch labels from the corrected change layers (steps 4-5)."""
    valid = np.isfinite(dnbr_corr) & np.isfinite(dndvi_corr)
    cand = (valid & ~bright.astype(bool) & (dnbr_corr >= DNBR_GROW) & (dndvi_corr >= DNDVI_MIN)
            & (pre_ndvi >= PRE_NDVI_MIN) & (post_nbr <= POST_NBR_MAX))
    labels, _ = burn_patches(cand, cand & (dnbr_corr >= DNBR_CORE))
    return labels


def map_burn(catalog, target, heath_mask, event_date: date) -> dict:
    """dNBR burn patches on `target` for a fire on `event_date`. `heath_mask`
    (MI or registered polygons) sets the offset and brightness references."""
    year = event_date.year
    bbox_ll = rasterio.warp.transform_bounds(str(target.crs), "EPSG:4326", *target.boundingbox)
    out: dict = {"status": "ok"}

    bs = read_scenes(catalog, bbox_ll, date(year - 1, BRIGHTNESS_MONTH, 1),
                     date(year - 1, BRIGHTNESS_MONTH + 1, 1), target, heath_mask)
    bright = np.zeros(target.shape, bool)
    if bs:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mb2 = np.nanmean(np.stack([s["b02"] for s in bs]), axis=0)
        v = mb2[heath_mask & np.isfinite(mb2)]
        if v.size:
            bright = np.isfinite(mb2) & (mb2 > np.percentile(v, BRIGHTNESS_PERCENTILE))
    out["n_brightness_scenes"] = len(bs)

    e_prev = date(year - 1, event_date.month, min(event_date.day, 28))
    for half in (PRE_HALF_WINDOW_DAYS, PRE_HALF_WINDOW_WIDE_DAYS):
        pre = read_scenes(catalog, bbox_ll, e_prev - timedelta(days=half), e_prev + timedelta(days=half),
                          target, heath_mask)
        if len(pre) >= MIN_OBS:
            break
    for w in (POST_WINDOW_DAYS, POST_WINDOW_WIDE_DAYS):
        post = read_scenes(catalog, bbox_ll, event_date + timedelta(days=1), event_date + timedelta(days=w),
                           target, heath_mask)
        if len(post) >= MIN_OBS:
            break
    out.update(pre_half_window_days=half, post_window_days=w, n_pre_scenes=len(pre), n_post_scenes=len(post),
               pre_dates=";".join(str(s["date"]) for s in pre), post_dates=";".join(str(s["date"]) for s in post))
    if len(pre) < MIN_OBS or len(post) < MIN_OBS:
        out["status"] = "too_few_scenes"
        return out

    pre_nbr, n_pre = nanmedian_stack(pre, "nbr")
    pre_ndvi, _ = nanmedian_stack(pre, "ndvi")
    post_nbr, n_post = nanmedian_stack(post, "nbr")
    post_ndvi, _ = nanmedian_stack(post, "ndvi")
    dnbr = pre_nbr - post_nbr
    dndvi = pre_ndvi - post_ndvi
    valid = np.isfinite(dnbr) & np.isfinite(dndvi)
    off_nbr = unburnt_offset(dnbr, pre_ndvi, valid & heath_mask)
    off_ndvi = unburnt_offset(dndvi, pre_ndvi, valid & heath_mask)
    dnbr_corr, dndvi_corr = dnbr - off_nbr, dndvi - off_ndvi

    labels = classify_layers(dnbr_corr, dndvi_corr, pre_ndvi, post_nbr, bright)

    out.update(dnbr_offset=off_nbr, dndvi_offset=off_ndvi,
               n_valid_heath_px=int((valid & heath_mask).sum()), n_heath_px=int(heath_mask.sum()),
               pre_nbr=pre_nbr, post_nbr=post_nbr, pre_ndvi=pre_ndvi, dnbr_corr=dnbr_corr,
               dndvi_corr=dndvi_corr, bright=bright, n_pre=n_pre, n_post=n_post, patch_id=labels,
               fig_pre_date=max(pre, key=lambda s: s["clear_frac"])["date"],
               fig_post_date=max(post, key=lambda s: s["clear_frac"])["date"])
    return out


def write_layers(path: Path, target, layers: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=target.shape[0], width=target.shape[1],
                       count=len(layers), dtype="float32", crs=str(target.crs), transform=target.affine,
                       nodata=np.nan, compress="deflate", predictor=3) as dst:
        for i, (k, v) in enumerate(layers.items(), start=1):
            dst.write(np.asarray(v, np.float32), i)
            dst.set_band_description(i, k)


def read_layers(path: Path) -> tuple[dict, object]:
    with rasterio.open(path) as src:
        return {src.descriptions[i]: src.read(i + 1) for i in range(src.count)}, src.transform


# =============================================================================
# 4. MI SITES AND FIRE EVENTS
# =============================================================================

def load_mi() -> gpd.GeoDataFrame:
    mi = gpd.read_file(MI_SHAPEFILE)
    mi = (mi.set_crs(MI_CRS) if mi.crs is None else mi).to_crs(MI_CRS)
    mi = mi[mi.geometry.notna() & ~mi.geometry.is_empty][["fid", "geometry"]].copy()
    mi["fid"] = mi["fid"].astype(int)
    return mi.reset_index(drop=True)


def build_events(mi: gpd.GeoDataFrame, firms: gpd.GeoDataFrame) -> tuple[pd.DataFrame, gpd.GeoDataFrame]:
    """Linked detections grouped into events. Returns (events, detections
    with event_id). Event ids are <year>_<nnn>, numbered by first date."""
    buf = mi.copy()
    buf["geometry"] = buf.geometry.buffer(max(active_fires.REACH_M[s] for s in active_fires.FIRMS_SENSORS))
    j = gpd.sjoin(firms, buf, predicate="within")
    geom = mi.set_index("fid").geometry
    j["dist_m"] = [p.distance(geom.loc[f]) for p, f in zip(j.geometry, j["fid"])]
    j = j[j["dist_m"] <= j["reach_m"]]
    det = firms.loc[sorted(set(j.index))].copy()
    x, y = det.geometry.x.values, det.geometry.y.values
    d = pd.to_datetime(det["date"]).values.astype("datetime64[D]").astype(int)
    yr = det["year"].values
    I, J = [], []
    for i in range(len(det)):
        k = np.where((np.hypot(x - x[i], y - y[i]) <= EVENT_LINK_M)
                     & (np.abs(d - d[i]) <= EVENT_GAP_DAYS) & (yr == yr[i]))[0]
        I += [i] * len(k)
        J += list(k)
    _, lab = connected_components(coo_matrix((np.ones(len(I)), (I, J)), shape=(len(det), len(det))))
    first = pd.Series(pd.to_datetime(det["date"]).values).groupby(lab).min().sort_values()
    count: dict[int, int] = {}
    name = {}
    for l, t in first.items():
        count[t.year] = count.get(t.year, 0) + 1
        name[l] = f"{t.year}_{count[t.year]:03d}"
    det["event_id"] = [name[l] for l in lab]
    links = pd.DataFrame({"fid": j["fid"].values, "event_id": det.loc[j.index, "event_id"].values})
    ev = det.groupby("event_id").agg(
        year=("year", "first"), first_date=("date", "min"), last_date=("date", "max"),
        n_detections=("date", "size"),
        months=("month", lambda s: ";".join(map(str, sorted(set(s))))),
        sensors=("satellite", lambda s: ";".join(sorted(set(s)))),
        summer=("month", lambda s: bool(s.isin([6, 7, 8]).any())),
    ).reset_index()
    lf = links.groupby("event_id")["fid"].apply(lambda s: sorted(set(s)))
    ev["linked_fids"] = ev["event_id"].map(lambda e: ";".join(map(str, lf.get(e, []))))
    ev["n_linked_sites"] = ev["event_id"].map(lambda e: len(lf.get(e, [])))
    return ev, det


# =============================================================================
# 5. STAGE 1: MAP ONE EVENT
# =============================================================================

def process_event(ev: dict, det: gpd.GeoDataFrame, mi: gpd.GeoDataFrame, catalog) -> dict:
    eid, year = ev["event_id"], int(ev["year"])
    pts = det[det["event_id"] == eid]
    linked = mi[mi["fid"].isin([int(f) for f in str(ev["linked_fids"]).split(";") if f])]
    reach = [active_fires.reach_shape(p, r) for p, r in zip(pts.geometry, pts["reach_m"])]
    aoi = box(*gpd.GeoSeries(list(linked.geometry) + reach, crs=MI_CRS).total_bounds).buffer(100).bounds
    target = aoi_geobox(aoi, crs=MI_CRS, resolution=RESOLUTION_M, bbox_crs=MI_CRS)
    heath = rasterio.features.geometry_mask(list(linked.geometry), out_shape=target.shape,
                                            transform=target.affine, invert=True)
    first_date = ev["first_date"] if isinstance(ev["first_date"], date) else date.fromisoformat(str(ev["first_date"]))
    res = map_burn(catalog, target, heath, first_date)
    row = {k: ev[k] for k in ("event_id", "year", "first_date", "last_date", "n_detections", "months",
                              "sensors", "summer", "n_linked_sites", "linked_fids")}
    row.update({k: v for k, v in res.items() if not isinstance(v, np.ndarray)})
    if res["status"] != "ok":
        return {"event": row, "patches": [], "site_events": []}

    out = finish_event(eid, year, pts, linked, mi, target, res, heath)
    row.update(out["row"])
    return {"event": row, "patches": out["patches"], "site_events": out["site_events"]}


def finish_event(eid, year, pts, linked, mi, target, res, heath) -> dict:
    """Steps 6-7 for one mapped event: accept patches within reach of its
    detections, record burnt MI sites, write the event layers and site masks.
    `res` holds the change layers and patch labels (map_burn() output, or the
    saved event layers re-classified by reclassify_events())."""
    labels = res["patch_id"]
    a = target.affine
    rr, cc = np.mgrid[0:target.shape[0], 0:target.shape[1]]
    xs, ys = a.c + (cc + 0.5) * a.a, a.f + (rr + 0.5) * a.e
    in_reach = active_fires.reach_mask(xs, ys, pts)
    patches = []
    burn = np.zeros(target.shape, bool)
    for pid in range(1, int(labels.max()) + 1):
        m = labels == pid
        near = bool((m & in_reach).any())
        d_min = float(min(np.hypot(xs[m] - p.x, ys[m] - p.y).min() for p in pts.geometry))
        patches.append({"event_id": eid, "year": year, "patch_id": pid, "area_ha": m.sum() * PIXEL_AREA_HA,
                        "mean_dnbr_corr": float(np.nanmean(res["dnbr_corr"][m])),
                        "p90_dnbr_corr": float(np.nanpercentile(res["dnbr_corr"][m], 90)),
                        "dist_to_firms_m": d_min, "firms_near": near})
        if near:
            burn |= m
    row = {"burn_area_ha": float(burn.sum() * PIXEL_AREA_HA), "n_patches": len(patches),
           "n_burn_patches": int(sum(p["firms_near"] for p in patches))}

    # every MI polygon in the area, linked or not, can be burnt by a patch
    site_events = []
    linked_set = set(linked["fid"])
    for _, s in mi[mi.geometry.intersects(box(*target.boundingbox))].iterrows():
        smask = rasterio.features.geometry_mask([s.geometry], out_shape=target.shape,
                                                transform=target.affine, invert=True)
        n = int((burn & smask).sum())
        if n * PIXEL_AREA_HA < MIN_SITE_BURN_HA:
            continue
        dist = pts.distance(s.geometry)
        site_events.append({"fid": int(s["fid"]), "year": year, "event_id": eid,
                            "burn_area_ha": n * PIXEL_AREA_HA, "polygon_area_ha": s.geometry.area / 1e4,
                            "n_valid_px": int((smask & np.isfinite(res["dnbr_corr"])).sum()),
                            "n_polygon_px": int(smask.sum()),
                            "n_firms_in_polygon": int((dist == 0).sum()),
                            "min_firms_dist_m": float(dist.min()), "linked": int(s["fid"]) in linked_set})
        write_site_mask(s, year, eid, burn, target)

    layers = {k: res[k] for k in ("pre_nbr", "post_nbr", "pre_ndvi", "dnbr_corr", "dndvi_corr")}
    layers.update(bright=res["bright"], n_pre=res["n_pre"], n_post=res["n_post"], patch_id=labels,
                  burn=burn, heath=heath)
    write_layers(OUTPUT_DIR / "events" / f"{eid}.tif", target, layers)
    return {"row": row, "patches": patches, "site_events": site_events}



def write_site_mask(site, year: int, eid: str, burn: np.ndarray, target) -> None:
    """Crop the event's burn mask onto the site's own AE-SJ grid (both grids
    share the 10 m lattice, so the crop is exact): site_events/<fid>_<year>_
    <event>.tif, merged per site-year by write_burns_table()."""
    st = aoi_geobox(site.geometry.bounds, crs=MI_CRS, resolution=RESOLUTION_M, bbox_crs=MI_CRS)
    col0 = int(round((st.affine.c - target.affine.c) / RESOLUTION_M))
    row0 = int(round((target.affine.f - st.affine.f) / RESOLUTION_M))
    crop = np.zeros(st.shape, np.uint8)
    r0, c0 = max(row0, 0), max(col0, 0)
    r1, c1 = min(row0 + st.shape[0], burn.shape[0]), min(col0 + st.shape[1], burn.shape[1])
    if r1 > r0 and c1 > c0:
        crop[r0 - row0:r1 - row0, c0 - col0:c1 - col0] = burn[r0:r1, c0:c1]
    smask = rasterio.features.geometry_mask([site.geometry], out_shape=st.shape, transform=st.affine, invert=True)
    crop &= smask.astype(np.uint8)
    path = OUTPUT_DIR / "site_events" / f"{int(site['fid'])}_{year}_{eid}.tif"
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=st.shape[0], width=st.shape[1], count=1,
                       dtype="uint8", crs=str(st.crs), transform=st.affine, compress="deflate") as dst:
        dst.write(crop, 1)


_W: dict = {}


def _init_worker():
    _W["mi"] = load_mi()
    _W["events"], _W["det"] = build_events(_W["mi"], active_fires.load(active_fires.FIRMS_SENSORS))
    _W["catalog"] = open_catalog()


def _run_event(ev: dict) -> dict:
    last = None
    for _ in range(MAX_RETRIES):
        try:
            return process_event(ev, _W["det"], _W["mi"], _W["catalog"])
        except Exception as exc:  # noqa: BLE001 - one bad event must not stop the run
            last = exc
            time.sleep(RETRY_DELAY_SECONDS)  # Planetary Computer 403s are throttling
    return {"event": {"event_id": ev["event_id"], "year": ev["year"], "status": f"error: {last}"},
            "patches": [], "site_events": []}


def run_map() -> None:
    import multiprocessing as mp

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mi = load_mi()
    events, det = build_events(mi, active_fires.load(active_fires.FIRMS_SENSORS))
    det.drop(columns="geometry").assign(x=det.geometry.x, y=det.geometry.y).to_csv(
        OUTPUT_DIR / "MI_firms_detections.csv", index=False)
    paths = {k: OUTPUT_DIR / f"{n}.csv" for k, n in
             (("event", "MI_firms_events"), ("patches", "MI_dnbr_patches"), ("site_events", "MI_burn_site_events"))}
    tables = {k: (pd.read_csv(p).to_dict("records") if p.exists() else []) for k, p in paths.items()}
    done = {r["event_id"] for r in tables["event"] if not str(r.get("status", "")).startswith("error")}
    tables["event"] = [r for r in tables["event"] if r["event_id"] in done]
    todo = events[~events["event_id"].isin(done)].to_dict("records")
    print(f"{len(events)} events ({len(det)} linked detections); {len(todo)} to map", flush=True)

    with mp.get_context("spawn").Pool(N_WORKERS, initializer=_init_worker) as pool:
        for i, res in enumerate(pool.imap_unordered(_run_event, todo), start=1):
            e = res["event"]
            print(f"[{i}/{len(todo)}] {e['event_id']} {e.get('status')} "
                  f"burn={e.get('burn_area_ha', float('nan')):.1f} ha sites={len(res['site_events'])}", flush=True)
            tables["event"].append(e)
            tables["patches"].extend(res["patches"])
            tables["site_events"].extend(res["site_events"])
            for k, p in paths.items():
                pd.DataFrame(tables[k]).to_csv(p, index=False)
    write_burns_table()


def reclassify_events() -> None:
    """Re-run steps 4-7 for every mapped event from its saved layers
    (events/<id>.tif) with the current thresholds; the Sentinel-2 composites
    and corrections do not depend on them, so this equals a full re-map."""
    from odc.geo.geobox import GeoBox
    mi = load_mi()
    events = pd.read_csv(OUTPUT_DIR / "MI_firms_events.csv")
    linked_by_event = build_events(mi, active_fires.load(active_fires.FIRMS_SENSORS))[0].set_index("event_id")["linked_fids"]
    events["linked_fids"] = events["event_id"].map(linked_by_event)
    dets = pd.read_csv(OUTPUT_DIR / "MI_firms_detections.csv")
    det = gpd.GeoDataFrame(dets, geometry=gpd.points_from_xy(dets.x, dets.y), crs=MI_CRS)
    tables = {"event": [], "patches": [], "site_events": []}
    for i, ev in events.iterrows():
        row = ev.to_dict()
        if ev["status"] != "ok":
            tables["event"].append(row)
            continue
        L, tf = read_layers(OUTPUT_DIR / "events" / f"{ev.event_id}.tif")
        target = GeoBox(L["dnbr_corr"].shape, tf, MI_CRS)
        L["patch_id"] = classify_layers(L["dnbr_corr"], L["dndvi_corr"], L["pre_ndvi"], L["post_nbr"],
                                        L["bright"] > 0)
        L["bright"] = L["bright"] > 0
        linked = mi[mi["fid"].isin([int(f) for f in str(ev["linked_fids"]).split(";") if f])]
        out = finish_event(ev.event_id, int(ev.year), det[det.event_id == ev.event_id], linked, mi, target, L,
                           L["heath"] > 0)
        row.update(out["row"])
        tables["event"].append(row)
        tables["patches"].extend(out["patches"])
        tables["site_events"].extend(out["site_events"])
        if (i + 1) % 25 == 0:
            print(f"[{i + 1}/{len(events)}] reclassified", flush=True)
    for k, n in (("event", "MI_firms_events"), ("patches", "MI_dnbr_patches"), ("site_events", "MI_burn_site_events")):
        pd.DataFrame(tables[k]).to_csv(OUTPUT_DIR / f"{n}.csv", index=False)
    write_burns_table()


def write_burns_table() -> pd.DataFrame:
    """One row per burnt site-year. Area from the site-year mask, so events
    that overlap in a year are not double counted."""
    se = pd.read_csv(OUTPUT_DIR / "MI_burn_site_events.csv")
    ev = pd.read_csv(OUTPUT_DIR / "MI_firms_events.csv").set_index("event_id")
    rows = []
    (OUTPUT_DIR / "rasters").mkdir(exist_ok=True)
    for (fid, year), g in se.groupby(["fid", "year"]):
        mask, profile = None, None
        for eid in g["event_id"]:
            with rasterio.open(OUTPUT_DIR / "site_events" / f"{fid}_{year}_{eid}.tif") as src:
                m = src.read(1)
                mask = m if mask is None else (mask | m)
                profile = src.profile
        with rasterio.open(OUTPUT_DIR / "rasters" / f"{fid}_{year}_burn.tif", "w", **profile) as dst:
            dst.write(mask, 1)
        n = int(mask.sum())
        e = ev.loc[g["event_id"]]
        rows.append({"fid": fid, "year": year, "burn_area_ha": n * PIXEL_AREA_HA,
                     "polygon_area_ha": g["polygon_area_ha"].iloc[0],
                     "burn_percent": 100 * n / g["n_polygon_px"].iloc[0],
                     "n_events": len(g), "event_ids": ";".join(g["event_id"]),
                     "firms_first_date": e["first_date"].min(), "firms_last_date": e["last_date"].max(),
                     "firms_months": ";".join(sorted({m for s in e["months"].astype(str) for m in s.split(";")}, key=int)),
                     "sensors": ";".join(sorted({x for s in e["sensors"] for x in s.split(";")})),
                     "n_firms_in_polygon": int(g["n_firms_in_polygon"].sum()),
                     "min_firms_dist_m": float(g["min_firms_dist_m"].min()),
                     "probable_fire_type": "wildfire" if e["summer"].astype(bool).any() else "prescribed_burn"})
    out = pd.DataFrame(rows).sort_values(["year", "fid"])
    out.to_csv(OUTPUT_DIR / "MI_burns.csv", index=False)
    print(f"MI_burns.csv: {len(out)} burnt site-years on {out.fid.nunique()} sites, "
          f"{out.burn_area_ha.sum():.1f} ha", flush=True)
    return out


# =============================================================================
# 6. STAGE 2: AE-SJ RESPONSE OF BURNT PIXELS
# =============================================================================

async def _aesj_site(fid, geometry, years_burnt, reader, index, boundaries):
    import mi_aesj_annual_products as aesj
    target = aoi_geobox(geometry.bounds, crs=MI_CRS, resolution=RESOLUTION_M, bbox_crs=MI_CRS)
    poly = rasterio.features.geometry_mask([geometry], out_shape=target.shape, transform=target.affine, invert=True)
    masks = {}
    for y in years_burnt:
        with rasterio.open(OUTPUT_DIR / "rasters" / f"{fid}_{y}_burn.tif") as src:
            masks[y] = src.read(1).astype(bool)
    ever = np.any(np.stack(list(masks.values())), axis=0)
    scores = {}
    for y in AEF_YEARS:
        s = await aesj._load_year_scores(reader, index, geometry.bounds, y, target)
        if s is not None and s.shape == poly.shape:
            scores[y] = s
    rows = []
    for by, bm in masks.items():
        for gname, gm in (("burnt", bm & poly), ("reference", poly & ~ever)):
            for oy, s in scores.items():
                v = s[gm & np.isfinite(s)]
                if v.size == 0:
                    continue
                codes = aesj.classify_sj(v, boundaries)
                rows.append({"fid": fid, "burn_year": by, "obs_year": oy, "years_since_burn": oy - by,
                             "group": gname, "n_px": int(v.size), "mean_AE_SJ": float(v.mean()),
                             "median_AE_SJ": float(np.median(v)), "sd_AE_SJ": float(v.std()),
                             **{f"prop_SJ{k}": float(np.mean(codes == k)) for k in (1, 2, 3, 4)}})
    return rows


async def _run_aesj(burns: pd.DataFrame, mi: gpd.GeoDataFrame) -> list[dict]:
    from aef_loader import AEFIndex, DataSource, VirtualTiffReader
    import mi_aesj_annual_products as aesj
    boundaries = aesj.load_sj_boundaries()
    index = AEFIndex(source=DataSource.SOURCE_COOP)
    await index.download()
    index.load()
    geom = mi.set_index("fid").geometry
    sem = asyncio.Semaphore(8)
    rows: list[dict] = []
    async with VirtualTiffReader() as reader:
        async def one(fid, years):
            async with sem:
                rows.extend(await _aesj_site(fid, geom.loc[fid], years, reader, index, boundaries))
        await asyncio.gather(*[one(f, sorted(g["year"])) for f, g in burns.groupby("fid")])
    return rows


def run_aesj() -> None:
    burns = pd.read_csv(OUTPUT_DIR / "MI_burns.csv")
    rows = asyncio.run(_run_aesj(burns, load_mi()))
    traj = pd.DataFrame(rows).sort_values(["fid", "burn_year", "group", "obs_year"])
    traj.to_csv(OUTPUT_DIR / "MI_burn_aesj_trajectories.csv", index=False)
    print(f"MI_burn_aesj_trajectories.csv: {len(traj)} rows, {traj.fid.nunique()} sites", flush=True)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("map", "all"):
        run_map()
    if stage == "burns":
        write_burns_table()
    if stage == "reclassify":
        reclassify_events()
    if stage in ("aesj", "all"):
        run_aesj()
