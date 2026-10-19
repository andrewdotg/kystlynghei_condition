#!/usr/bin/env python
"""mi_aesj_annual_products.py - national annual AE-SJ scoring of the MI
kystlynghei polygons, and the 2018-2019 -> 2024-2025 pixel transitions.

For every MI polygon and every year 2018-2025, the AlphaEarth embeddings
(source.coop archive, read with aef_loader) are scored with the fixed AE-SJ
model (ae_sj_model.py, generated from R/model.R) on a lattice-snapped 10 m
grid (aef_geobox.aoi_geobox), and summarised per polygon-year.

OUTPUT (OUTPUT_DIR)
  MI_aesj_annual_<year>.csv  per polygon: AE-SJ mean/median/percentiles,
      within-polygon sd, model-uncertainty CI (each pixel carries its SJ
      class's field-validated sd, propagated as sqrt(sum(sd^2))/n), and the
      number and share of SJ1-4-like pixels.
  MI_aesj_pixel_transitions_2018_2019_to_2024_2025.csv  per polygon: the
      4x4 SJ-class transition matrix between the baseline (mean of 2018 and
      2019) and final (mean of 2024 and 2025) score of every pixel, and the
      SJ3/SJ4 -> SJ1/SJ2 transitions split into scrub and tree removal and
      by mechanism:
        * tree removal: SJ4-source pixels and all pixels in the national
          Picea mask (spruce are trees whatever their SJ-like class)
        * scrub removal: SJ3-source pixels outside the Picea mask
        * burn-coincident: scrub pixels inside any burn mapped by
          mi_dnbr_burn.py in 2018-2025, and tree pixels inside a burn of a
          probable wildfire (June-August FIRMS detection)
        * clearance: everything else
      Prescribed burning is applied to open heath and does not remove trees,
      and felling causes the same NBR drop as fire, so tree pixels inside a
      prescribed burn are counted as clearance. Uncontrolled summer fires can
      kill trees, so tree pixels inside a probable-wildfire burn are counted
      as burn-coincident. A burn-coincident pixel may still have been cleared
      rather than burnt.

RUN MODES
  "all"            the full national run (14,330 polygons, ~15 h at 10
                   polygons in parallel).
  "refresh_burnt"  re-processes only the polygons whose burn split can
                   change - polygons with a burn in mi_dnbr_burn_results/
                   MI_burns.csv, and polygons with burn-coincident pixels in
                   the existing output - and merges them into the existing
                   outputs. Everything else is independent of the burn layer.
                   Run this after mi_dnbr_burn.py.

CONCURRENCY: the blocking AEF fetch runs in asyncio.to_thread, so polygons
overlap (benchmark: 10 in parallel, 7.1x faster than serial, identical
output). Batching neighbouring polygons into one combined read was tested
and rejected: slices of a combined read differed by up to 0.16 AE-SJ units
from the same window read on its own.
"""

from __future__ import annotations

import asyncio
import traceback
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.features
import rasterio.windows
from odc.geo.geobox import GeoBox

from aef_loader import AEFIndex, DataSource, VirtualTiffReader
from aef_loader.utils import dequantize_aef, reproject_datatree

from ae_sj_model import score_ae_sj
from aef_geobox import aoi_geobox

# =============================================================================
# 1. SETTINGS
# =============================================================================

REPO = Path(__file__).resolve().parent.parent
MI_SHAPEFILE = REPO / "inputs" / "MI_kystlynghei.shp"
MI_CRS = "EPSG:25833"
# National Picea mask, EPSG:25833 10 m, on the same 0,0-anchored lattice as
# aoi_geobox, so read_picea_mask() is a plain windowed read.
PICEA_MASK_TIF = REPO / "inputs" / "picea_mask_norway_25833_10m.tif"
BURN_DIR = Path(__file__).parent / "mi_dnbr_burn_results"
MODEL_CLASS_REFERENCE_CSV = REPO / "outputs" / "model" / "FINAL_AE_SJ_covadj_class_reference_positions.csv"
OUTPUT_DIR = Path(__file__).parent / "mi_aesj_annual_results"

RUN_MODE = "refresh_burnt"   # "all" or "refresh_burnt"
MAX_CONCURRENT_FEATURES = 10
CHECKPOINT_EVERY_N_FEATURES = 100

ANALYSIS_YEARS = list(range(2018, 2026))
BASELINE_YEARS = (2018, 2019)
FINAL_YEARS = (2024, 2025)
BURN_YEARS = list(range(2018, 2026))

RESOLUTION_M = 10.0
PIXEL_AREA_HA = RESOLUTION_M**2 / 10000.0
TRANSITION_PAIRS = [(a, b) for a in (1, 2, 3, 4) for b in (1, 2, 3, 4)]
TRANSITIONS_CSV = "MI_aesj_pixel_transitions_2018_2019_to_2024_2025.csv"


# =============================================================================
# 2. MODEL REFERENCE VALUES
# =============================================================================

def load_sj_boundaries() -> tuple[float, float, float]:
    """SJ1/2, SJ2/3, SJ3/4 boundaries: midpoints between adjacent field-class
    medians (outputs/model/, written by R/model.R)."""
    m = pd.read_csv(MODEL_CLASS_REFERENCE_CSV).set_index("sj")["median_AE_SJ"]
    return (m[1] + m[2]) / 2.0, (m[2] + m[3]) / 2.0, (m[3] + m[4]) / 2.0


def load_class_sd() -> np.ndarray:
    """Field-validated sd of AE-SJ per SJ class, indexed by classify_sj() code."""
    sd = pd.read_csv(MODEL_CLASS_REFERENCE_CSV).set_index("sj")["sd_AE_SJ"]
    return np.array([np.nan, sd[1], sd[2], sd[3], sd[4]])


def classify_sj(scores: np.ndarray, boundaries: tuple[float, float, float]) -> np.ndarray:
    sj12, sj23, sj34 = boundaries
    codes = np.select([scores >= sj12, scores >= sj23, scores >= sj34], [1, 2, 3], default=4)
    return np.where(np.isfinite(scores), codes, 0).astype(np.int8)


# =============================================================================
# 3. INPUTS
# =============================================================================

def load_mi_features() -> gpd.GeoDataFrame:
    mi = gpd.read_file(MI_SHAPEFILE)
    mi = (mi.set_crs(MI_CRS) if mi.crs is None else mi).to_crs(MI_CRS)
    mi = mi[mi.geometry.notna() & ~mi.geometry.is_empty][["fid", "geometry"]].reset_index(drop=True)
    mi["fid"] = mi["fid"].astype(int)
    return mi


def read_picea_mask(target: GeoBox) -> np.ndarray:
    """Picea mask on `target`'s grid. The window is built from target's own
    snapped extent; the raw polygon bounds round differently and come out
    one pixel short. Opened per call: GDAL handles are not thread-safe."""
    left, top = target.affine.c, target.affine.f
    right, bottom = left + target.shape.x * RESOLUTION_M, top - target.shape.y * RESOLUTION_M
    with rasterio.open(PICEA_MASK_TIF) as src:
        window = rasterio.windows.from_bounds(left, bottom, right, top, transform=src.transform
                                              ).round_offsets().round_lengths()
        arr = src.read(1, window=window, boundless=True, fill_value=255)
    if arr.shape != (target.shape.y, target.shape.x):
        raise ValueError(f"picea window {arr.shape} != target {(target.shape.y, target.shape.x)}")
    return arr == 1


_WILDFIRE_YEARS: dict | None = None


def wildfire_years(fid: int) -> set:
    """Years in which the polygon's burn was a probable wildfire (MI_burns.csv)."""
    global _WILDFIRE_YEARS
    if _WILDFIRE_YEARS is None:
        b = pd.read_csv(BURN_DIR / "MI_burns.csv")
        w = b[b["probable_fire_type"] == "wildfire"]
        _WILDFIRE_YEARS = w.groupby("fid")["year"].apply(set).to_dict()
    return _WILDFIRE_YEARS.get(fid, set())


def read_burn_mask(fid: int, target: GeoBox) -> tuple[np.ndarray, np.ndarray]:
    """Union of the polygon's mapped burns 2018-2025, and of those that were
    probable wildfires (mi_dnbr_burn.py writes them on this same grid)."""
    out = np.zeros(target.shape, bool)
    wild = np.zeros(target.shape, bool)
    wy = wildfire_years(fid)
    for y in BURN_YEARS:
        p = BURN_DIR / "rasters" / f"{fid}_{y}_burn.tif"
        if p.exists():
            with rasterio.open(p) as src:
                arr = src.read(1)
            if arr.shape != out.shape:
                raise ValueError(f"burn mask {p} {arr.shape} != target {out.shape}")
            out |= arr.astype(bool)
            if y in wy:
                wild |= arr.astype(bool)
    return out, wild


# =============================================================================
# 4. SCORING AND STATISTICS
# =============================================================================

async def _load_year_scores(reader, index, bounds, year: int, target: GeoBox) -> np.ndarray | None:
    tiles = await index.query(bbox=bounds, years=year, bbox_crs=MI_CRS)
    if not tiles:
        return None
    tree = await reader.open_tiles_by_zone(tiles)
    combined = reproject_datatree(tree, target, resampling="nearest")
    combined = await asyncio.to_thread(combined.compute)   # blocking S3 read
    emb = dequantize_aef(combined["embeddings"])
    if "time" in emb.dims:
        emb = emb.isel(time=0)
    return score_ae_sj(emb.transpose("y", "x", "band").values)


def _pixel_stats(scores, mask, prefix, boundaries, class_sd) -> dict:
    """Summary of the selected pixels. sd_AE_SJ is the within-polygon spread
    (ecological heterogeneity); ci95 is model uncertainty only, from each
    pixel's class sd, SE = sqrt(sum(sd^2))/n."""
    vals = scores[mask]
    vals = vals[np.isfinite(vals)]
    n = int(vals.size)
    keys = ["mean", "median", "p10", "p25", "p75", "p90", "sd", "model_se", "ci95_lo", "ci95_hi"]
    if n == 0:
        return {f"{prefix}n": 0, **{f"{prefix}{k}_AE_SJ": np.nan for k in keys}}
    mean = float(vals.mean())
    se = float(np.sqrt(np.sum(class_sd[classify_sj(vals, boundaries)] ** 2)) / n)
    p10, p25, p75, p90 = np.percentile(vals, [10, 25, 75, 90])
    return {f"{prefix}n": n, f"{prefix}mean_AE_SJ": mean, f"{prefix}median_AE_SJ": float(np.median(vals)),
            f"{prefix}p10_AE_SJ": float(p10), f"{prefix}p25_AE_SJ": float(p25),
            f"{prefix}p75_AE_SJ": float(p75), f"{prefix}p90_AE_SJ": float(p90),
            f"{prefix}sd_AE_SJ": float(np.std(vals, ddof=1)) if n > 1 else np.nan,
            f"{prefix}model_se_AE_SJ": se, f"{prefix}ci95_lo_AE_SJ": mean - 1.96 * se,
            f"{prefix}ci95_hi_AE_SJ": mean + 1.96 * se}


async def process_feature(reader, index, fid, geometry, boundaries, class_sd) -> tuple[list[dict], dict]:
    target = aoi_geobox(geometry.bounds, crs=MI_CRS, resolution=RESOLUTION_M, bbox_crs=MI_CRS)
    poly = rasterio.features.geometry_mask([geometry], out_shape=target.shape, transform=target.affine, invert=True)
    picea = await asyncio.to_thread(read_picea_mask, target)
    burnt, wild = await asyncio.to_thread(read_burn_mask, fid, target)

    base_sum = np.zeros(target.shape)
    base_n = np.zeros(target.shape, np.int16)
    final_sum = np.zeros(target.shape)
    final_n = np.zeros(target.shape, np.int16)
    rows = []
    for year in ANALYSIS_YEARS:
        scores = await _load_year_scores(reader, index, geometry.bounds, year, target)
        if scores is None:
            rows.append({"fid": fid, "year": year, "status": f"no_aef_tiles_{year}"})
            continue
        if scores.shape != poly.shape:
            rows.append({"fid": fid, "year": year, "status": f"grid_mismatch_{scores.shape}"})
            continue
        valid = poly & np.isfinite(scores)
        row = {"fid": fid, "year": year, "status": "ok", **_pixel_stats(scores, valid, "aesj_", boundaries, class_sd)}
        codes = classify_sj(scores, boundaries)
        nv = int(valid.sum())
        row["n_valid"] = nv
        row["polygon_area_ha"] = int(poly.sum()) * PIXEL_AREA_HA
        for sj in (1, 2, 3, 4):
            row[f"n_SJ{sj}_like"] = int((valid & (codes == sj)).sum())
        for sj in (1, 2, 3, 4):
            row[f"prop_SJ{sj}_like"] = row[f"n_SJ{sj}_like"] / nv if nv else np.nan
        rows.append(row)
        if year in BASELINE_YEARS:
            base_sum[valid] += scores[valid]
            base_n[valid] += 1
        if year in FINAL_YEARS:
            final_sum[valid] += scores[valid]
            final_n[valid] += 1

    tr: dict = {"fid": fid}
    eligible = poly & (base_n == len(BASELINE_YEARS)) & (final_n == len(FINAL_YEARS))
    n_el = int(eligible.sum())
    tr["n_valid"] = n_el
    tr["n_burnt_in_transition_pool"] = int((eligible & burnt).sum())
    if n_el == 0:
        tr["status"] = "no_transition_eligible_pixels"
        return rows, tr
    tr["status"] = "ok"
    with np.errstate(invalid="ignore"):
        b_codes = classify_sj(np.where(eligible, base_sum / np.maximum(base_n, 1), np.nan), boundaries)
        f_codes = classify_sj(np.where(eligible, final_sum / np.maximum(final_n, 1), np.nan), boundaries)
    for a, b in TRANSITION_PAIRS:
        tr[f"SJ{a}_to_SJ{b}"] = int((eligible & (b_codes == a) & (f_codes == b)).sum())

    to12 = np.isin(f_codes, [1, 2])
    sj3 = eligible & (b_codes == 3) & to12
    sj4 = eligible & (b_codes == 4) & to12
    sj34 = sj3 | sj4
    pic = sj34 & picea
    scrub = sj3 & ~picea                 # scrub removal
    tree = sj4 | pic                     # tree removal: SJ4-source or spruce
    scrub_burnt = scrub & burnt          # any mapped burn
    tree_burnt = tree & wild             # probable wildfires only
    burn_coincident = scrub_burnt | tree_burnt
    counts = {"n_SJ3_to_SJ12": sj3, "n_SJ4_to_SJ12": sj4, "n_SJ34_to_SJ12": sj34,
              "n_scrub_to_SJ12": scrub, "n_tree_to_SJ12": tree,
              "n_scrub_to_SJ12_burnt": scrub_burnt, "n_tree_to_SJ12_burnt": tree_burnt,
              "n_SJ34_to_SJ12_burnt": burn_coincident, "n_SJ34_to_SJ12_clearance": sj34 & ~burn_coincident,
              "n_SJ3_to_SJ12_picea": sj3 & picea, "n_SJ4_to_SJ12_picea": sj4 & picea,
              "n_SJ34_to_SJ12_picea": pic, "n_SJ34_to_SJ12_picea_burnt": pic & wild,
              "n_SJ34_to_SJ12_picea_clearance": pic & ~wild}
    for k, m in counts.items():
        tr[k] = int(m.sum())
    tr["n_picea_px"] = int((eligible & picea).sum())
    for k in ("scrub_to_SJ12", "tree_to_SJ12", "scrub_to_SJ12_burnt", "tree_to_SJ12_burnt"):
        tr[f"area_{k}_ha"] = tr[f"n_{k}"] * PIXEL_AREA_HA
    for k in ("SJ34_to_SJ12", "SJ34_to_SJ12_burnt", "SJ34_to_SJ12_clearance", "SJ34_to_SJ12_picea",
              "SJ34_to_SJ12_picea_burnt", "SJ34_to_SJ12_picea_clearance"):
        tr[f"area_{k}_ha"] = tr[f"n_{k}"] * PIXEL_AREA_HA
    tr["prop_SJ34_to_SJ12_picea"] = (tr["n_SJ34_to_SJ12_picea"] / tr["n_SJ34_to_SJ12"]
                                     if tr["n_SJ34_to_SJ12"] else np.nan)
    n1 = sum(tr[f"SJ{a}_to_SJ1"] for a in (2, 3, 4))
    n0 = sum(tr[f"SJ{a}_to_SJ{a}"] for a in (1, 2, 3, 4))
    n4 = sum(tr[f"SJ{a}_to_SJ4"] for a in (1, 2, 3))
    tr.update(n_toward_SJ1=n1, n_same_zone=n0, n_toward_SJ4=n4, prop_toward_SJ1=n1 / n_el,
              prop_same_zone=n0 / n_el, prop_toward_SJ4=n4 / n_el)
    return rows, tr


# =============================================================================
# 5. OUTPUTS
# =============================================================================

# Columns not produced by this script (found in older result files), dropped on write.
LEGACY_PREFIXES = ("firms_", "confirmed_burn", "nbr_", "n_nbr_", "n_same_year", "same_year",
                   "n_preburn", "n_postburn", "mean_ae_sj_gain", "n_confirmed_burn")


def _drop_legacy(df: pd.DataFrame) -> pd.DataFrame:
    return df[[c for c in df.columns if not c.startswith(LEGACY_PREFIXES)]]


def write_outputs(annual: list[dict], transitions: list[dict]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    a = _drop_legacy(pd.DataFrame(annual))
    for y in ANALYSIS_YEARS:
        a[a["year"] == y].sort_values("fid").to_csv(OUTPUT_DIR / f"MI_aesj_annual_{y}.csv", index=False)
    t = _drop_legacy(pd.DataFrame(transitions))
    # Polygons scored before the scrub/tree columns existed have no burns
    # (every burnt polygon is re-scored by refresh_burnt), so their values
    # follow exactly from the SJ3/SJ4 and spruce counts.
    if "n_scrub_to_SJ12" in t.columns:
        miss = t["n_scrub_to_SJ12"].isna() & t["n_SJ3_to_SJ12"].notna()
        t.loc[miss, "n_scrub_to_SJ12"] = t.loc[miss, "n_SJ3_to_SJ12"] - t.loc[miss, "n_SJ3_to_SJ12_picea"]
        t.loc[miss, "n_tree_to_SJ12"] = t.loc[miss, "n_SJ4_to_SJ12"] + t.loc[miss, "n_SJ3_to_SJ12_picea"]
        t.loc[miss, ["n_scrub_to_SJ12_burnt", "n_tree_to_SJ12_burnt"]] = 0
        for k in ("scrub_to_SJ12", "tree_to_SJ12", "scrub_to_SJ12_burnt", "tree_to_SJ12_burnt"):
            t[f"area_{k}_ha"] = t[f"n_{k}"] * PIXEL_AREA_HA
    t.sort_values("fid").to_csv(OUTPUT_DIR / TRANSITIONS_CSV, index=False)


def load_existing() -> tuple[list[dict], list[dict]]:
    annual = []
    for y in ANALYSIS_YEARS:
        p = OUTPUT_DIR / f"MI_aesj_annual_{y}.csv"
        if p.exists():
            annual.extend(pd.read_csv(p).to_dict("records"))
    p = OUTPUT_DIR / TRANSITIONS_CSV
    return annual, (pd.read_csv(p).to_dict("records") if p.exists() else [])


async def main() -> None:
    mi = load_mi_features()
    boundaries, class_sd = load_sj_boundaries(), load_class_sd()
    annual, transitions = load_existing()
    if RUN_MODE == "refresh_burnt":
        burns = pd.read_csv(BURN_DIR / "MI_burns.csv")
        tr = pd.DataFrame(transitions)
        prev = set(tr.loc[tr["n_SJ34_to_SJ12_burnt"].fillna(0) > 0, "fid"].astype(int))
        todo = set(burns["fid"].astype(int)) | prev
        annual = [r for r in annual if int(r["fid"]) not in todo]
        transitions = [r for r in transitions if int(r["fid"]) not in todo]
        mi = mi[mi["fid"].isin(todo)].reset_index(drop=True)
        print(f"refresh_burnt: {len(mi)} polygons ({burns.fid.nunique()} burnt, {len(prev)} with earlier "
              f"burn-coincident pixels)", flush=True)
    elif RUN_MODE == "all":
        done = {int(r["fid"]) for r in transitions}
        mi = mi[~mi["fid"].isin(done)].reset_index(drop=True)
        print(f"all: {len(mi)} polygons to score ({len(done)} already done)", flush=True)
    else:
        raise ValueError(f"RUN_MODE must be 'all' or 'refresh_burnt', got {RUN_MODE!r}")

    index = AEFIndex(source=DataSource.SOURCE_COOP)
    await index.download()
    index.load()
    sem, lock = asyncio.Semaphore(MAX_CONCURRENT_FEATURES), asyncio.Lock()
    done_n = 0

    async def run_one(fid, geometry):
        nonlocal done_n
        async with sem:
            try:
                rows, tr = await process_feature(reader, index, fid, geometry, boundaries, class_sd)
            except Exception as exc:  # noqa: BLE001 - one bad polygon must not stop the run
                traceback.print_exc()
                rows = [{"fid": fid, "year": None, "status": f"error: {exc}"}]
                tr = {"fid": fid, "status": f"error: {exc}"}
        async with lock:
            annual.extend(rows)
            transitions.append(tr)
            done_n += 1
            print(f"[{done_n}/{len(mi)}] fid={fid}", flush=True)
            if done_n % CHECKPOINT_EVERY_N_FEATURES == 0:
                write_outputs(annual, transitions)

    async with VirtualTiffReader() as reader:
        await asyncio.gather(*[run_one(int(r.fid), r.geometry) for r in mi.itertuples()])
    write_outputs(annual, transitions)
    print(f"Wrote {len(transitions)} polygons to {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    asyncio.run(main())
