"""validation_registered_burns.py - the burn detection tested against burning
registered in naturtiltak (Methodology Section 11.4).

The burns are mapped the way mi_dnbr_burn.py maps MI sites - an active-fire
detection gives the date, Sentinel-2 dNBR the burnt area - on polygons where
burning was registered under the grant scheme for threatened nature types.
Two triggers are tested: FIRMS alone (as used for the MI screening) and
FIRMS plus Sentinel-3 SLSTR (active_fires.py).

STAGE 1 - POOL (`python validation_registered_burns.py pool`)
  From the grant register (TILTAK_SRC: polygons with the final-report text),
  every polygon whose measure includes burning and whose registered year is
  2021-2025. The final report is read for each: registrations whose report
  says no heath burning was done (only clearing, grazing, brush or waste
  burning, burning planned for later, or burning stopped by the weather;
  NOT_BURNT) are removed, unless the report says the heath was burnt in
  earlier years and FIRMS has a detection there up to the registered year
  (EARLIER_BURNING). Registrations without a report are kept, flagged by
  whether the grant was paid. Writes the burns-verified pool to
  inputs/naturtiltak/tiltak_burns_verified.gpkg (all register columns, plus
  burn_evidence and reported_burnt_daa).

STAGE 2 - MAP (`python validation_registered_burns.py map`)
  Registrations overlapping by more than DEDUP_OVERLAP (share of the smaller
  polygon) are the same site registered again and are merged into one
  feature (union, latest registered year). The registered year is usually
  the reporting date at the end of a grant period, so every detection within
  reach of a feature from 2018 up to and including its registered year is
  used; detections less than FIRE_GAP_DAYS apart in the same year form one
  fire. This is done per feature: a fire within reach of several
  neighbouring features is a fire at each of them (and one distinct fire in
  the results). Each fire is mapped with mi_dnbr_burn.map_burn() at its first
  detection date and its per-pixel layers are saved, together with a
  control: the same feature and detections moved CONTROL_SHIFT_KM in a
  random direction (centroid on land), mapped at the same date.

STAGE 3 - THRESHOLDS (`python validation_registered_burns.py thresholds`)
  The classification repeated offline over a grid of thresholds (dNBR
  candidate and core, minimum patch, minimum burn in the polygon, pre-fire
  NDVI, post-fire NBR) for the FIRMS-triggered fires: fires detected against
  controls with a "burn".

STAGE 4 - FINAL (`python validation_registered_burns.py final`)
  With the operational thresholds (mi_dnbr_burn.py). A fire is a hit when
  >= MIN_SITE_BURN_HA of burn patches within reach of its detections lie
  inside the feature or within HIT_BUFFER_M of it (registered polygons are
  roughly drawn). A fire has no usable imagery when fewer than
  MIN_VALID_SHARE of the pixels in the feature and its HIT_BUFFER_M zone
  have clear observations before and after the fire - a different failure
  from a burn that is not visible, so results are given for all fires and
  for fires with usable imagery. Also: detection by polygon size and by the
  burnt area given in the final reports, with a logistic fit of the FIRMS
  detection probability against reported burnt area (the smallest burn the
  method can be expected to find).

STAGE 5 - SPATIAL REPRESENTATIVENESS (`python validation_registered_burns.py spatial`)
  Does FIRMS fire more often at registered burning than in the surrounding
  landscape? Every registration is compared with N_COPIES copies moved 5, 10
  and 20 km (centroid on land), counting FIRMS detections within reach per
  year relative to the registered year (Y-3..Y+1) and over 2018..Y, plus the
  density of detections in rings around each polygon.

OUTPUT: python/validation_registered_burns/ (validation_final.json,
validation_fires.csv, validation_features.csv, threshold_grid.csv,
spatial_representativeness.json, pool_summary.json, pool_report_check.csv,
fire_layers/).
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio.features
from shapely.affinity import translate
from shapely.geometry import Point, box

import active_fires
import mi_dnbr_burn as burn
from aef_geobox import aoi_geobox

REPO = Path(__file__).resolve().parent.parent
TILTAK_SRC = REPO.parent.parent / "Data" / "Kystlynghei" / "kystlynghei_naturtiltak_tiltakEN.gpkg"
POOL = REPO / "inputs" / "naturtiltak" / "tiltak_burns_verified.gpkg"
POOL_LAYER = "tiltak_burns_verified"
OUT = Path(__file__).parent / "validation_registered_burns"
LAYER_DIR = OUT / "fire_layers"
YEARS = (2021, 2025)
YEAR_FIX = {350: 2022}          # register year 2099 (placeholder); 2022 as in tiltakPolys_GEE.shp
SIZE_BINS = [0, 5, 20, 100, 1e6]
SIZE_LABELS = ["< 5 ha", "5-20 ha", "20-100 ha", "> 100 ha"]
REPORTED_BINS = [0, 1, 5, 20, 1e6]
REPORTED_LABELS = ["< 1 ha", "1-5 ha", "5-20 ha", "> 20 ha"]
TRIGGERS = {"firms": active_fires.FIRMS_SENSORS, "firms_s3": ("VIIRS", "MODIS", "SLSTR")}

# Final-report reading (stage 1), by register fid. No heath burning in the
# reported period: 105-106 not burnt that year; 109-113, 482-483, 542
# clearing or juniper cutting only; 132, 146-147, 424, 465, 973-974 brush,
# waste or ditch burning only; 330, 439, 540, 618, 636, 760 burning attempted
# or planned, not done; 541, 778, 779, 799-800, 821 weather or a deliberate
# pause; 361, 473 grazing only; 453 firebreaks only.
NOT_BURNT = {105, 106, 109, 110, 111, 112, 113, 132, 146, 147, 330, 361, 424, 439, 453, 465, 473, 482,
             483, 540, 541, 542, 618, 636, 760, 778, 779, 799, 800, 821, 973, 974}
EARLIER_BURNING = {541, 779}    # report: heath burnt in earlier years (kept if FIRMS agrees)
UNCLEAR = {732}                 # "se vedlagt sluttrapport", no burnt area given
OTHER_PART = {157, 158}         # burnt, but east of the cairn, not the area applied for
EVIDENCE_ORDER = ["report: burnt", "report: burnt (other part of site)", "earlier burning reported + FIRMS",
                  "report: unclear", "no report, grant paid", "no report, not paid"]

DEDUP_OVERLAP = 0.9
FIRE_GAP_DAYS = 30
CONTROL_SHIFT_KM = 5
HIT_BUFFER_M = 250
MIN_VALID_SHARE = 0.25
SEED = 42


# =============================================================================
# STAGE 1: BURNS-VERIFIED POOL
# =============================================================================

def run_pool() -> dict:
    src = gpd.read_file(TILTAK_SRC, fid_as_index=True)
    src.index.name = "src_fid"
    src = src.reset_index().to_crs(burn.MI_CRS)
    src["tiltakYear"] = [str(YEAR_FIX.get(f, y)) for f, y in zip(src.src_fid, src.tiltakYear)]
    year = pd.to_numeric(src["tiltakYear"], errors="coerce")
    reg = src[src["tiltakTypes"].fillna("").str.contains("burn") & year.between(*YEARS)].copy()
    firms = active_fires.load(active_fires.FIRMS_SENSORS)
    reg["firms_up_to_year"] = [bool((d["year"] <= int(y)).any()) for d, y in
                               ((active_fires.detections_near(g, firms), y)
                                for g, y in zip(reg.geometry, reg["tiltakYear"]))]
    has_report = reg["tiltakSrc"] != "local_gpkg_specific_fields"
    paid = reg["sum_utbetalt"].fillna(0) > 0

    def evidence(fid, rep, pd_, fr):
        if fid in NOT_BURNT:
            return "earlier burning reported + FIRMS" if fid in EARLIER_BURNING and fr else None
        if fid in UNCLEAR:
            return "report: unclear"
        if fid in OTHER_PART:
            return "report: burnt (other part of site)"
        if rep:
            return "report: burnt"
        return "no report, grant paid" if pd_ else "no report, not paid"
    reg["burn_evidence"] = [evidence(*a) for a in zip(reg.src_fid, has_report, paid, reg.firms_up_to_year)]
    daa = reg["tiltakEN"].fillna("").str.extract(r"burned area: ([0-9]+(?:[,.][0-9]+)?)")[0]
    reg["reported_burnt_daa"] = pd.to_numeric(daa.str.replace(",", "."), errors="coerce")
    reg.drop(columns="geometry").assign(kept=reg.burn_evidence.notna()).to_csv(OUT / "pool_report_check.csv",
                                                                               index=False)
    pool = reg[reg.burn_evidence.notna()]
    POOL.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:          # GeoPackage cannot be written on the network drive
        p = Path(tmp) / POOL.name
        pool.to_file(p, layer=POOL_LAYER, driver="GPKG")
        shutil.copyfile(p, POOL)
    removed = reg[reg.burn_evidence.isna()]
    out = {"registered_burning": len(reg), "kept": len(pool), "removed": len(removed),
           "removed_with_firms": int(removed.firms_up_to_year.sum()),
           "removed_fids": removed.src_fid.tolist(),
           "kept_by_evidence": pool.burn_evidence.value_counts().to_dict(),
           "features_registered": int(len(dedup_features(_as_reg(reg.assign(burn_evidence=None))))),
           "features_pool": int(len(dedup_features(_as_reg(pool))))}
    json.dump(out, open(OUT / "pool_summary.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    return out


def _as_reg(g: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame({"reg_id": g["src_fid"].astype(int).values,
                             "reg_year": pd.to_numeric(g["tiltakYear"]).astype(int).values,
                             "burn_evidence": g["burn_evidence"].values,
                             "reported_burnt_daa": g["reported_burnt_daa"].values},
                            geometry=g.geometry.values, crs=g.crs)


def load_pool() -> gpd.GeoDataFrame:
    return _as_reg(gpd.read_file(POOL, layer=POOL_LAYER).to_crs(burn.MI_CRS))


def dedup_features(reg: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    r = reg.reset_index(drop=True)
    sj = gpd.sjoin(r[["geometry"]], r[["geometry"]], predicate="intersects")
    a, b = sj.index.values, sj.index_right.values
    ov = np.array([r.geometry[i].intersection(r.geometry[j]).area / min(r.geometry[i].area, r.geometry[j].area)
                   for i, j in zip(a, b)])
    m = ov > DEDUP_OVERLAP
    _, lab = connected_components(coo_matrix((np.ones(m.sum()), (a[m], b[m])), shape=(len(r), len(r))))
    r["feature"] = lab
    agg = {"reg_id": lambda s: ";".join(map(str, sorted(s))), "reg_year": "max", "reported_burnt_daa": "max",
           "burn_evidence": lambda s: next((o for o in EVIDENCE_ORDER if o in set(s)), None)}
    f = r.dissolve(by="feature", aggfunc=agg)
    f = f.rename(columns={"reg_id": "reg_ids", "reg_year": "year"}).reset_index()
    f["n_registrations"] = f["reg_ids"].str.count(";") + 1
    f["polygon_ha"] = f.geometry.area / 1e4
    return f


# =============================================================================
# STAGE 2: FIRES AND MAPPING
# =============================================================================

def _split_fires(det: pd.DataFrame) -> list[pd.DataFrame]:
    out = []
    for _, d in det.sort_values("date").groupby("year"):
        dates = pd.to_datetime(d["date"])
        grp = (dates.diff().dt.days.fillna(0) > FIRE_GAP_DAYS).cumsum()
        out.extend(g for _, g in d.groupby(grp.values))
    return out


def _point_ids(pts: pd.DataFrame) -> frozenset:
    return frozenset(f"{p.x:.0f},{p.y:.0f},{d}" for p, d in zip(pts.geometry, pts["date"]))


def feature_fires(feats: gpd.GeoDataFrame) -> tuple[pd.DataFrame, dict]:
    """Fires per feature for each trigger, and the mapping jobs. Identical
    fires under both triggers share one layer key (mapped once)."""
    from scipy.sparse.csgraph import connected_components
    fires_all = active_fires.load()
    rows, jobs = [], {}
    for trig, sensors in TRIGGERS.items():
        det = fires_all[fires_all["sensor"].isin(sensors)]
        for _, f in feats.iterrows():
            d = active_fires.detections_near(f.geometry, det)
            d = d[d["year"] <= f["year"]]
            for k, pts in enumerate(_split_fires(d) if len(d) else [], start=1):
                ids = _point_ids(pts)
                key = (f"{f['feature']}_{min(pts['date'])}_"
                       f"{hashlib.md5('|'.join(sorted(ids)).encode()).hexdigest()[:6]}")
                jobs.setdefault(key, {"key": key, "feature": f["feature"], "fire_date": min(pts["date"]),
                                      "geometry": f.geometry, "points": pts})
                rows.append({"trigger": trig, "feature": f["feature"], "fire": k, "layer_key": key,
                             "year": int(pts["year"].iloc[0]), "fire_date": min(pts["date"]),
                             "last_date": max(pts["date"]), "n_detections": len(pts),
                             "sensors": ";".join(sorted(set(pts["sensor"]))),
                             "satellites": ";".join(sorted(set(pts["satellite"]))),
                             "months": ";".join(map(str, sorted(set(pts["month"])))),
                             "inside": bool((pts["dist_m"] == 0).any()), "_ids": ids})
    fires = pd.DataFrame(rows)
    fires["distinct_fire"] = ""
    for trig in TRIGGERS:                                # fires at different features sharing a detection
        idx = fires.index[fires.trigger == trig]
        ids = fires.loc[idx, "_ids"].tolist()
        adj = np.array([[bool(a & b) for b in ids] for a in ids])
        _, lab = connected_components(adj)
        fires.loc[idx, "distinct_fire"] = [f"{trig}_{x}" for x in lab]
    return fires.drop(columns="_ids"), jobs


_W: dict = {}


def _init():
    _W["catalog"] = burn.open_catalog()


def _save_layers(job: dict) -> dict:
    g, pts = job["geometry"], job["points"]
    dx, dy = job["shift"]
    if dx or dy:
        g = translate(g, dx, dy)
        pts = pts.copy()
        pts["geometry"] = pts.geometry.translate(dx, dy)
    reach = [active_fires.reach_shape(p, s, r) for p, s, r in zip(pts.geometry, pts["sensor"], pts["reach_m"])]
    aoi = box(*gpd.GeoSeries([g] + reach, crs=burn.MI_CRS).total_bounds).buffer(800).bounds
    target = aoi_geobox(aoi, crs=burn.MI_CRS, resolution=burn.RESOLUTION_M, bbox_crs=burn.MI_CRS)
    poly = rasterio.features.geometry_mask([g], out_shape=target.shape, transform=target.affine, invert=True)
    a = target.affine
    rr_, cc_ = np.mgrid[0:target.shape[0], 0:target.shape[1]]
    xs, ys = a.c + (cc_ + 0.5) * a.a, a.f + (rr_ + 0.5) * a.e
    in_reach = active_fires.reach_mask(xs, ys, pts)
    from scipy.ndimage import distance_transform_edt
    dist_poly = distance_transform_edt(~poly) * burn.RESOLUTION_M
    res = None
    for _ in range(burn.MAX_RETRIES):
        try:
            res = burn.map_burn(_W["catalog"], target, poly, job["fire_date"])
            break
        except Exception:  # noqa: BLE001
            time.sleep(burn.RETRY_DELAY_SECONDS)
    kind = "control" if (dx or dy) else "fire"
    row = {"layer_key": job["key"], "kind": kind, "status": "error" if res is None else res["status"]}
    if res is not None and res["status"] == "ok":
        np.savez_compressed(LAYER_DIR / f"{job['key']}_{kind}.npz",
                            dnbr=res["dnbr_corr"].astype(np.float32), dndvi=res["dndvi_corr"].astype(np.float32),
                            pre_ndvi=res["pre_ndvi"].astype(np.float32), post_nbr=res["post_nbr"].astype(np.float32),
                            bright=res["bright"], poly=poly, in_reach=in_reach, dist_poly=dist_poly.astype(np.float32))
    return row


def run_map() -> None:
    import multiprocessing as mp
    feats = dedup_features(load_pool())
    fires, jobs = feature_fires(feats)
    feats.drop(columns="geometry").to_csv(OUT / "validation_features_input.csv", index=False)
    fires.to_csv(OUT / "validation_fires_input.csv", index=False)
    LAYER_DIR.mkdir(parents=True, exist_ok=True)
    ne = gpd.read_file(REPO / "inputs" / "natural_earth" / "ne_10m_admin_0_countries.shp")
    land = ne[ne["ADMIN"] == "Norway"].to_crs(burn.MI_CRS).geometry.union_all()
    st = OUT / "fire_layers_status.csv"
    done = pd.read_csv(st) if st.exists() else pd.DataFrame(columns=["layer_key", "kind", "status"])
    have = set(zip(done.layer_key, done.kind))
    work = []
    for key, j in jobs.items():
        rng = np.random.default_rng([int(hashlib.md5(key.encode()).hexdigest()[:8], 16), SEED])
        for _ in range(200):
            th = rng.uniform(0, 2 * np.pi)
            dx, dy = CONTROL_SHIFT_KM * 1000 * np.cos(th), CONTROL_SHIFT_KM * 1000 * np.sin(th)
            if land.contains(Point(j["geometry"].centroid.x + dx, j["geometry"].centroid.y + dy)):
                break
        for kind, shift in (("fire", (0.0, 0.0)), ("control", (dx, dy))):
            if (key, kind) not in have:
                work.append({**j, "shift": shift})
    print(f"{len(feats)} features, fires per trigger {fires.groupby('trigger').size().to_dict()}, "
          f"{len(jobs)} layer keys, {len(work)} maps to run", flush=True)
    rows = []
    if work:
        with mp.get_context("spawn").Pool(burn.N_WORKERS, initializer=_init) as pool:
            for i, r in enumerate(pool.imap_unordered(_save_layers, work), start=1):
                rows.append(r)
                print(f"[{i}/{len(work)}] {r}", flush=True)
    allst = pd.concat([done, pd.DataFrame(rows, columns=["layer_key", "kind", "status"])])
    allst = allst.drop_duplicates(["layer_key", "kind"], keep="last")
    allst[allst.layer_key.isin(jobs)].to_csv(st, index=False)


def _layers(key: str, kind: str) -> dict:
    z = np.load(LAYER_DIR / f"{key}_{kind}.npz")
    return {k: z[k] for k in z.files}


def classify(L: dict, grow, core, min_patch_ha, pre_ndvi_min, post_nbr_max, closing=1, min_core_px=4):
    """Burn patches within reach of the fire's detections, for given thresholds
    (the rule of mi_dnbr_burn.classify_layers() and burn_patches())."""
    import scipy.ndimage as ndi
    valid = np.isfinite(L["dnbr"]) & np.isfinite(L["dndvi"])
    cand = (valid & ~L["bright"] & (L["dnbr"] >= grow) & (L["dndvi"] >= 0)
            & (L["pre_ndvi"] >= pre_ndvi_min) & (L["post_nbr"] <= post_nbr_max))
    eight = np.ones((3, 3), bool)
    lab, n = ndi.label(cand, structure=eight)
    if n == 0:
        return np.zeros(cand.shape, bool)
    keep = ndi.sum(cand & (L["dnbr"] >= core), lab, index=np.arange(n + 1)) >= min_core_px
    keep[0] = False
    grown = keep[lab]
    pad = closing + 1
    closed = ndi.binary_closing(np.pad(grown, pad), structure=eight, iterations=closing)[pad:-pad, pad:-pad]
    plab, pn = ndi.label(closed, structure=eight)
    flagged = np.where(grown, plab, 0)
    big = np.bincount(flagged.ravel(), minlength=pn + 1) * burn.PIXEL_AREA_HA >= min_patch_ha
    big[0] = False
    near = np.zeros(pn + 1, bool)
    near[np.unique(flagged[L["in_reach"] & (flagged > 0)])] = True
    return big[flagged] & near[flagged] & grown


# =============================================================================
# STAGE 3: THRESHOLD GRID
# =============================================================================

def run_thresholds() -> pd.DataFrame:
    import itertools
    status = pd.read_csv(OUT / "fire_layers_status.csv")
    fires = pd.read_csv(OUT / "validation_fires_input.csv")
    keys = set(fires.loc[fires.trigger == "firms", "layer_key"])
    ok = status[(status.status == "ok") & status.layer_key.isin(keys)]
    layers = {(r.layer_key, r.kind): _layers(r.layer_key, r.kind) for r in ok.itertuples()}
    grid = list(itertools.product((0.05, 0.075, 0.10), (0.15, 0.20, 0.27), (0.1, 0.25, 0.5), (0.1, 0.2),
                                  (0.15, 0.20), (0.30, 0.40)))
    rows = []
    for grow, core, mp_ha, site_ha, ndvi_min, nbr_max in grid:
        if core < grow:
            continue
        fk, ck = [], []
        for (key, kind), L in layers.items():
            b = classify(L, grow, core, mp_ha, ndvi_min, nbr_max)
            (fk if kind == "fire" else ck).append((b & L["poly"]).sum() * burn.PIXEL_AREA_HA >= site_ha)
        rows.append({"dnbr_candidate": grow, "dnbr_core": core, "min_patch_ha": mp_ha, "min_site_ha": site_ha,
                     "pre_ndvi_min": ndvi_min, "post_nbr_max": nbr_max,
                     "fires_detected": int(sum(fk)), "fires": len(fk),
                     "controls_detected": int(sum(ck)), "controls": len(ck)})
    g = pd.DataFrame(rows)
    g["fire_rate"] = g.fires_detected / g.fires
    g["false_rate"] = g.controls_detected / g.controls
    g = g.sort_values(["fires_detected", "controls_detected"], ascending=[False, True])
    g.to_csv(OUT / "threshold_grid.csv", index=False)
    print(g.head(25).to_string(index=False))
    return g


# =============================================================================
# STAGE 4: FINAL VALIDATION
# =============================================================================

def logistic_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Maximum-likelihood logistic regression y ~ a + b*x (Newton-Raphson)."""
    X = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2)
    for _ in range(100):
        p = 1 / (1 + np.exp(-X @ beta))
        W = p * (1 - p)
        step = np.linalg.solve(X.T @ (X * W[:, None]) + 1e-9 * np.eye(2), X.T @ (y - p))
        beta += step
        if np.abs(step).max() < 1e-8:
            break
    return float(beta[0]), float(beta[1])


def _classify_layer(key: str, kind: str) -> dict:
    L = _layers(key, kind)
    b = classify(L, burn.DNBR_GROW, burn.DNBR_CORE, burn.MIN_PATCH_HA, burn.PRE_NDVI_MIN, burn.POST_NBR_MAX)
    dp = L["dist_poly"]
    zone = dp <= HIT_BUFFER_M
    rec = {"inside_ha": (b & L["poly"]).sum() * burn.PIXEL_AREA_HA,
           "adjacent_ha": (b & zone & ~L["poly"]).sum() * burn.PIXEL_AREA_HA,
           "beyond_ha": (b & (dp > HIT_BUFFER_M) & (dp <= 750)).sum() * burn.PIXEL_AREA_HA,
           "valid_share": float(np.isfinite(L["dnbr"])[zone].mean())}
    rec["hit_inside"] = rec["inside_ha"] >= burn.MIN_SITE_BURN_HA
    rec["hit"] = rec["inside_ha"] + rec["adjacent_ha"] >= burn.MIN_SITE_BURN_HA
    rec["usable"] = rec["valid_share"] >= MIN_VALID_SHARE
    return rec


def _block(f: pd.DataFrame, ctrl: pd.DataFrame) -> dict:
    fh = f.groupby("feature")["hit"].any()
    dh = f.groupby("distinct_fire")["hit"].any()
    return {"fires": int(len(f)), "fires_hit": int(f.hit.sum()), "fires_hit_inside": int(f.hit_inside.sum()),
            "fires_hit_pct": 100 * f.hit.mean() if len(f) else None,
            "distinct_fires": int(dh.size), "distinct_fires_hit": int(dh.sum()),
            "features": int(fh.size), "features_hit": int(fh.sum()),
            "features_hit_inside": int(f.groupby("feature")["hit_inside"].any().sum()),
            "controls": int(len(ctrl)), "controls_hit": int(ctrl.hit.sum()),
            "controls_hit_pct": 100 * ctrl.hit.mean() if len(ctrl) else None,
            "controls_hit_inside": int(ctrl.hit_inside.sum()),
            "controls_hit_inside_pct": 100 * ctrl.hit_inside.mean() if len(ctrl) else None}


def run_final() -> dict:
    status = pd.read_csv(OUT / "fire_layers_status.csv")
    fe = pd.read_csv(OUT / "validation_features_input.csv")
    fires = pd.read_csv(OUT / "validation_fires_input.csv")
    lay = pd.DataFrame([{"layer_key": r.layer_key, "kind": r.kind, "status": r.status,
                         **(_classify_layer(r.layer_key, r.kind) if r.status == "ok"
                            else {"hit": False, "hit_inside": False, "usable": False, "valid_share": 0.0,
                                  "beyond_ha": 0.0})}
                        for r in status.itertuples()])
    fi = fires.merge(lay[lay.kind == "fire"].drop(columns="kind"), on="layer_key", how="left")
    for c in ("hit", "hit_inside", "usable"):
        fi[c] = fi[c].fillna(False).astype(bool)
    fi["result"] = np.select([fi.hit_inside, fi.hit, ~fi.usable, fi.beyond_ha.fillna(0) >= burn.MIN_SITE_BURN_HA],
                             ["inside polygon", "within 250 m", "no usable imagery", "250-750 m only"],
                             "no burn scar")
    fi.to_csv(OUT / "validation_fires.csv", index=False)
    ct = lay[(lay.kind == "control") & (lay.status == "ok")]   # failed controls are uninformative

    out = {"thresholds": {"dnbr_candidate": burn.DNBR_GROW, "dnbr_core": burn.DNBR_CORE,
                          "min_patch_ha": burn.MIN_PATCH_HA, "min_site_ha": burn.MIN_SITE_BURN_HA,
                          "pre_ndvi_min": burn.PRE_NDVI_MIN, "post_nbr_max": burn.POST_NBR_MAX,
                          "hit_buffer_m": HIT_BUFFER_M, "min_valid_share": MIN_VALID_SHARE},
           "pool": json.load(open(OUT / "pool_summary.json")),
           "registrations": int(fe["n_registrations"].sum()), "features": int(len(fe))}
    for trig in TRIGGERS:
        f = fi[fi.trigger == trig]
        c = ct[ct.layer_key.isin(set(f.layer_key))]
        cu = ct[ct.layer_key.isin(set(f.loc[f.usable, "layer_key"])) & ct.usable]
        fe[trig] = fe.feature.isin(set(f.feature))
        fe[f"{trig}_hit"] = fe.feature.isin(set(f.loc[f.hit, "feature"]))
        out[trig] = {"features_detected": int(fe[trig].sum()), "features_detected_pct": 100 * fe[trig].mean(),
                     "all": _block(f, c), "usable_imagery": _block(f[f.usable], cu),
                     "fires_no_usable_imagery": int((~f.usable).sum()),
                     "result_counts": f.result.value_counts().to_dict()}
    # FIRMS-only figures at the top level: the trigger used for the MI screening
    out.update({"features_firms": out["firms"]["features_detected"],
                "features_firms_pct": out["firms"]["features_detected_pct"],
                "all": out["firms"]["all"], "usable_imagery": out["firms"]["usable_imagery"],
                "result_counts": out["firms"]["result_counts"]})
    new = fi[(fi.trigger == "firms_s3") & ~fi.layer_key.isin(set(fi.loc[fi.trigger == "firms", "layer_key"]))]
    out["s3_added"] = {"features_detected": int((fe.firms_s3 & ~fe.firms).sum()),
                       "features_hit": int((fe.firms_s3_hit & ~fe.firms_hit).sum()),
                       "fires_new_or_changed": int(len(new)), "fires_new_or_changed_hit": int(new.hit.sum()),
                       "fires_s3_only": int((new.sensors == "SLSTR").sum()),
                       "fires_s3_only_hit": int(new.loc[new.sensors == "SLSTR", "hit"].sum())}
    agg = dict(features=("feature", "size"), firms=("firms", "sum"), dnbr=("firms_hit", "sum"),
               firms_s3=("firms_s3", "sum"), dnbr_s3=("firms_s3_hit", "sum"))
    fe["size"] = pd.cut(fe.polygon_ha, SIZE_BINS, labels=SIZE_LABELS)
    out["by_polygon_size"] = fe.groupby("size", observed=True).agg(**agg).reset_index().astype(
        {"size": str}).to_dict("records")
    out["by_evidence"] = fe.groupby("burn_evidence").agg(**agg).reset_index().to_dict("records")
    # detection against the burnt area given in the final reports (decares;
    # one value per application, so a maximum where several polygons share it)
    rep = fe[fe.burn_evidence.isin(EVIDENCE_ORDER[:2]) & fe.reported_burnt_daa.notna()].copy()
    rep["reported_ha"] = rep.reported_burnt_daa / 10
    rep["rsize"] = pd.cut(rep.reported_ha, REPORTED_BINS, labels=REPORTED_LABELS)
    out["by_reported_burnt_area"] = rep.groupby("rsize", observed=True).agg(**agg).reset_index().astype(
        {"rsize": str}).to_dict("records")
    a, b = logistic_fit(np.log10(rep.reported_ha.values), rep.firms.values.astype(float))
    out["logistic_reported_area"] = {"n": int(len(rep)), "n_detected": int(rep.firms.sum()),
                                     "intercept": a, "slope_per_log10_ha": b,
                                     "area_ha_p50": float(10 ** (-a / b)) if b > 0 else None,
                                     "area_ha_p25": float(10 ** ((np.log(1 / 3) - a) / b)) if b > 0 else None}
    fy = fi[fi.trigger == "firms"]
    out["fire_year_minus_registered"] = (fy["year"] - fy.feature.map(fe.set_index("feature")["year"])
                                         ).value_counts().sort_index().to_dict()
    # MI burns on sites overlapping registered features (the MI method itself)
    feats_g = dedup_features(load_pool())
    mi = burn.load_mi()
    mib = pd.read_csv(burn.OUTPUT_DIR / "MI_burns.csv")
    mi_geo = gpd.GeoDataFrame(mib.merge(mi, on="fid"), geometry="geometry", crs=burn.MI_CRS)
    j = gpd.sjoin(mi_geo, feats_g[["feature", "year", "geometry"]].rename(columns={"year": "reg_year"}),
                  predicate="intersects")
    j = j[j.year <= j.reg_year]
    ov = gpd.overlay(feats_g[["feature", "geometry"]], mi, how="intersection")
    ov = ov[ov.geometry.area >= 5000]
    out["operational"] = {"features_overlapping_mi": int(ov.feature.nunique()),
                          "features_with_mi_burn_up_to_registered_year": int(j.feature.nunique())}
    fe.drop(columns="size").to_csv(OUT / "validation_features.csv", index=False)
    json.dump(out, open(OUT / "validation_final.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "pool"}, indent=1, default=str))
    return out


# =============================================================================
# STAGE 5: SPATIAL REPRESENTATIVENESS (FIRMS)
# =============================================================================

SHIFT_KM = (5, 10, 20)
N_COPIES = 20
REL_YEARS = range(-3, 2)
RINGS_KM = ((0.75, 5), (5, 10), (10, 20))


def run_spatial() -> dict:
    reg = load_pool()
    firms = active_fires.load(active_fires.FIRMS_SENSORS)
    ne = gpd.read_file(REPO / "inputs" / "natural_earth" / "ne_10m_admin_0_countries.shp")
    land = ne[ne["ADMIN"] == "Norway"].to_crs(burn.MI_CRS).geometry.union_all()
    rng = np.random.default_rng(SEED)
    by_year = {y: firms[firms["year"] == y].reset_index(drop=True) for y in range(2018, 2026)}

    def hit(g, f):
        return bool(len(active_fires.detections_near(g, f)))

    def flags(g, y_reg):
        out = {k: (hit(g, by_year[y_reg + k]) if y_reg + k in by_year else None) for k in REL_YEARS}
        out["pre"] = any(hit(g, by_year[y]) for y in range(2018, y_reg + 1))
        return out

    rows = []
    for _, r in reg.iterrows():
        g, y = r.geometry, int(r.reg_year)
        rows.append({"reg_id": r.reg_id, "reg_year": y, "kind": "registered", "shift_km": 0, **flags(g, y)})
        for d in SHIFT_KM:
            n, tries = 0, 0
            while n < N_COPIES and tries < 20 * N_COPIES:
                tries += 1
                th = rng.uniform(0, 2 * np.pi)
                gs = translate(g, xoff=d * 1000 * np.cos(th), yoff=d * 1000 * np.sin(th))
                if not land.contains(gs.centroid):
                    continue
                n += 1
                rows.append({"reg_id": r.reg_id, "reg_year": y, "kind": "displaced", "shift_km": d, **flags(gs, y)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "spatial_representativeness_flags.csv", index=False)

    res: dict = {"n_registered": len(reg), "n_copies": N_COPIES, "shift_km": list(SHIFT_KM)}
    tab = []
    for col in [*REL_YEARS, "pre"]:
        sub = df[df[col].notna()]
        row = {"period": "2018 to Y" if col == "pre" else (f"Y{col:+d}" if col else "Y")}
        obs = sub[sub.kind == "registered"][col].astype(float)
        row["registered"], row["n"], row["k"] = float(obs.mean()), int(obs.size), int(obs.sum())
        z, n_, ph = 1.96, obs.size, obs.mean()             # Wilson 95 % interval
        c = (ph + z * z / (2 * n_)) / (1 + z * z / n_)
        h = z * np.sqrt(ph * (1 - ph) / n_ + z * z / (4 * n_ * n_)) / (1 + z * z / n_)
        row["registered_lo"], row["registered_hi"] = float(c - h), float(c + h)
        for d in SHIFT_KM:     # copies averaged per registration first, so each weighs the same
            x = sub[(sub.kind == "displaced") & (sub.shift_km == d)].groupby("reg_id")[col].apply(
                lambda v: v.astype(float).mean())
            row[f"displaced_{d}km"] = float(x.mean())
        tab.append(row)
    res["share_with_detection"] = tab

    dens = []
    for _, r in reg.iterrows():
        g, y = r.geometry, int(r.reg_year)
        f = gpd.GeoDataFrame(pd.concat([by_year[k] for k in range(2018, y + 1)], ignore_index=True),
                             geometry="geometry", crs=burn.MI_CRS)
        near = f.iloc[f.sindex.query(g.buffer(RINGS_KM[-1][1] * 1000), predicate="intersects")]
        dist = near.distance(g)
        rec = {"reg_id": r.reg_id, "years": y - 2018 + 1,
               "n_reach": int((dist <= 750).sum()), "km2_reach": g.buffer(750).area / 1e6}
        for a, b in RINGS_KM:      # rings clipped to land (coarse coastline), pooled over polygons
            ring = g.buffer(b * 1000).difference(g.buffer(a * 1000)).intersection(land)
            rec[f"n_{a}_{b}km"] = int(((dist > a * 1000) & (dist <= b * 1000)).sum())
            rec[f"km2_{a}_{b}km"] = ring.area / 1e6
        dens.append(rec)
    dd = pd.DataFrame(dens)
    res["density_per_100km2_per_year_pooled"] = {
        key: round(100 * dd[f"n_{key}"].sum() / (dd[f"km2_{key}"] * dd["years"]).sum(), 3)
        for key in ["reach"] + [f"{a}_{b}km" for a, b in RINGS_KM]}
    json.dump(res, open(OUT / "spatial_representativeness.json", "w"), indent=1, default=str)
    print(pd.DataFrame(tab).round(3).to_string())
    print(json.dumps(res["density_per_100km2_per_year_pooled"], indent=1))
    return res


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    stages = {"pool": run_pool, "map": run_map, "thresholds": run_thresholds, "final": run_final,
              "spatial": run_spatial}
    for name in (stages if stage == "all" else [stage]):
        stages[name]()
