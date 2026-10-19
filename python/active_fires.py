"""active_fires.py - the active-fire detections used to trigger burn mapping.

One input file, inputs/active_fires/active_fires_2018_2025.gpkg, holds every
detection used in the project, one row per hot pixel:

    sensor  product   satellite      pixel   period             used by
    VIIRS   J1V-C2    NOAA-20        375 m   2018-04 to 2025    MI screening + validation
    VIIRS   SV-C2     Suomi NPP      375 m   2018 to 2025       MI screening + validation
    VIIRS   J2V-C2    NOAA-21        375 m   2023-03 to 2025    MI screening + validation
    MODIS   M-C61     Terra/Aqua     1 km    2018 to 2025-09    MI screening + validation
    SLSTR   SL_2_FRP  Sentinel-3A/B  1 km    2020-08 to 2025    validation only

FIRMS (NASA, the first four) keeps TYPE 0 (presumed vegetation fire); MODIS
also needs CONFIDENCE >= MODIS_MIN_CONFIDENCE. Sentinel-3 (ESA, s3_frp.py)
keeps hotspot class 1 (vegetation fire; gas flares, volcanic and industrial
hotspots are dropped). Sentinel-3 is not used for the MI screening: it adds
few burns and has many false alarms on summer evenings (Methodology
Section 11.4).

Reach: a detection is linked to a polygon when the polygon lies within
REACH_M of the detection - a circle of two VIIRS pixels (750 m) or one MODIS
pixel (1,000 m) around the pixel centre, or for SLSTR the 3 x 3 pixel square
(+-1,500 m) around it.

Raw inputs (outside the repo, read only by `python active_fires.py build`):
the FIRMS archive shapefiles as downloaded (FIRMS_RAW) and the Sentinel-3
monthly fire lists written by `python s3_frp.py download` (s3_frp.CACHE).
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

REPO = Path(__file__).resolve().parent.parent
BURNING = REPO.parent / "Burning"
FIRMS_RAW = {
    BURNING / "DL_FIRE_J1V-C2_807416" / "fire_archive_J1V-C2_807416.shp": ("VIIRS", "J1V-C2", "NOAA-20"),
    BURNING / "extraFIRMS" / "DL_FIRE_SV-C2_823020" / "fire_archive_SV-C2_823020.shp": ("VIIRS", "SV-C2", "S-NPP"),
    BURNING / "extraFIRMS" / "DL_FIRE_J2V-C2_823019" / "fire_archive_J2V-C2_823019.shp": ("VIIRS", "J2V-C2", "NOAA-21"),
    BURNING / "extraFIRMS" / "DL_FIRE_M-C61_823022" / "fire_archive_M-C61_823022.shp": ("MODIS", "M-C61", "Terra/Aqua"),
}
ACTIVE_FIRES = REPO / "inputs" / "active_fires" / "active_fires_2018_2025.gpkg"
LAYER = "active_fires"
MI_CRS = "EPSG:25833"
FIRST_YEAR, LAST_YEAR = 2018, 2025

FIRMS_TYPE_FILTER = [0]
MODIS_MIN_CONFIDENCE = 30
S3_CLASS_FILTER = [1]
REACH_M = {"VIIRS": 750.0, "MODIS": 1000.0, "SLSTR": 1500.0}
SQUARE_REACH = {"SLSTR"}          # 3 x 3 pixel footprint rather than a circle
FIRMS_SENSORS = ("VIIRS", "MODIS")


# ---------------------------------------------------------------- build
def _firms() -> gpd.GeoDataFrame:
    parts = []
    for path, (sensor, product, satellite) in FIRMS_RAW.items():
        g = gpd.read_file(path)
        g = g[g["TYPE"].isin(FIRMS_TYPE_FILTER)]
        if sensor == "MODIS":
            g = g[pd.to_numeric(g["CONFIDENCE"], errors="coerce") >= MODIS_MIN_CONFIDENCE]
        parts.append(gpd.GeoDataFrame({
            "date": pd.to_datetime(g["ACQ_DATE"]).dt.strftime("%Y-%m-%d"),
            "time_utc": g["ACQ_TIME"].astype(str).str.zfill(4),
            "sensor": sensor, "product": product, "satellite": satellite,
            "daynight": g["DAYNIGHT"], "confidence": g["CONFIDENCE"].astype(str),
            "frp_mw": pd.to_numeric(g["FRP"], errors="coerce"),
        }, geometry=g.geometry, crs=g.crs).to_crs(MI_CRS))
    return pd.concat(parts, ignore_index=True)


def _sentinel3() -> gpd.GeoDataFrame:
    import s3_frp
    s = s3_frp.load_s3(crs=MI_CRS)
    s = s[s["classification"].isin(S3_CLASS_FILTER)]
    t = s["time"].dt.tz_convert("UTC")
    return gpd.GeoDataFrame({
        "date": pd.Series([str(d) for d in s["date"]], index=s.index),   # local (Oslo) date
        "time_utc": t.dt.strftime("%H%M"),
        "sensor": "SLSTR", "product": "SL_2_FRP", "satellite": s["satellite"],
        "daynight": s["daynight"], "confidence": "",                      # SLSTR confidence is always 0
        "frp_mw": s["FRP_MWIR"],
    }, geometry=s.geometry, crs=MI_CRS)


def build() -> None:
    f = pd.concat([_firms(), _sentinel3()], ignore_index=True)
    d = pd.to_datetime(f["date"])
    f["year"], f["month"] = d.dt.year, d.dt.month
    f = f[f["year"].between(FIRST_YEAR, LAST_YEAR)]
    f["reach_m"] = f["sensor"].map(REACH_M)
    # rows stay in archive order (FIRMS products in FIRMS_RAW order, then
    # Sentinel-3): mi_dnbr_burn.py numbers same-day fire events in this order
    f = gpd.GeoDataFrame(f.reset_index(drop=True), geometry="geometry", crs=MI_CRS)
    ACTIVE_FIRES.parent.mkdir(parents=True, exist_ok=True)
    # GeoPackage (SQLite) cannot be written in place on the network drive
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / ACTIVE_FIRES.name
        f.to_file(out, layer=LAYER, driver="GPKG")
        shutil.copyfile(out, ACTIVE_FIRES)
    print(f"{len(f)} detections -> {ACTIVE_FIRES}")
    print(f.groupby(["sensor", "satellite"]).size().to_string())


# ---------------------------------------------------------------- use
def load(sensors=None, crs: str = MI_CRS) -> gpd.GeoDataFrame:
    """Detections, one row per hot pixel. `sensors` limits to e.g.
    FIRMS_SENSORS (MI screening). Columns: date (datetime.date), time_utc,
    year, month, sensor, product, satellite, daynight, confidence, frp_mw,
    reach_m, geometry."""
    f = gpd.read_file(ACTIVE_FIRES, layer=LAYER)
    if sensors is not None:
        f = f[f["sensor"].isin(sensors)]
    f["date"] = pd.to_datetime(f["date"]).dt.date
    return f.to_crs(crs).reset_index(drop=True)


def reach_shape(point, sensor: str, reach_m: float):
    """The area a detection can stand for: a circle, or a square for SLSTR."""
    if sensor in SQUARE_REACH:
        return box(point.x - reach_m, point.y - reach_m, point.x + reach_m, point.y + reach_m)
    return point.buffer(reach_m)


def detections_near(geometry, fires: gpd.GeoDataFrame, year: int | None = None) -> gpd.GeoDataFrame:
    """Detections whose reach touches `geometry`, with the distance to the
    polygon (0 = inside) in dist_m."""
    f = fires if year is None else fires[fires["year"] == year]
    f = f.iloc[f.sindex.query(geometry.buffer(max(REACH_M.values()) * 1.5), predicate="intersects")]
    if f.empty:
        return f.assign(dist_m=pd.Series(dtype=float))
    f = f.assign(dist_m=f.distance(geometry))
    near = [reach_shape(p, s, r).intersects(geometry) for p, s, r in zip(f.geometry, f["sensor"], f["reach_m"])]
    return f[np.array(near, dtype=bool)]


def reach_mask(xs: np.ndarray, ys: np.ndarray, pts: gpd.GeoDataFrame) -> np.ndarray:
    """Pixels (centre coordinates xs, ys) within reach of any detection."""
    m = np.zeros(xs.shape, bool)
    for p, s, r in zip(pts.geometry, pts["sensor"], pts["reach_m"]):
        if s in SQUARE_REACH:
            m |= (np.abs(xs - p.x) <= r) & (np.abs(ys - p.y) <= r)
        else:
            m |= np.hypot(xs - p.x, ys - p.y) <= r
    return m


if __name__ == "__main__":
    if (sys.argv[1] if len(sys.argv) > 1 else "build") == "build":
        build()
