"""active_fires.py - the FIRMS active-fire detections used to trigger burn mapping.

One input file, inputs/active_fires/active_fires_2018_2025.gpkg, holds every
detection used in the project, one row per hot pixel, from four NASA FIRMS
archive products:

    sensor  product   satellite      pixel   period
    VIIRS   J1V-C2    NOAA-20        375 m   2018-04 to 2025
    VIIRS   SV-C2     Suomi NPP      375 m   2018 to 2025
    VIIRS   J2V-C2    NOAA-21        375 m   2023-03 to 2025
    MODIS   M-C61     Terra/Aqua     1 km    2018 to 2025-09

Only TYPE 0 (presumed vegetation fire) is kept; MODIS detections also need
CONFIDENCE >= MODIS_MIN_CONFIDENCE.

Reach: a detection is linked to a polygon when the polygon lies within
REACH_M of the detection's pixel centre - two VIIRS pixels (750 m) or one
MODIS pixel (1,000 m).

Raw inputs (outside the repo, read only by `python active_fires.py build`):
the FIRMS archive shapefiles as downloaded (FIRMS_RAW).
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

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
REACH_M = {"VIIRS": 750.0, "MODIS": 1000.0}
FIRMS_SENSORS = ("VIIRS", "MODIS")


# ---------------------------------------------------------------- build
def build() -> None:
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
    f = pd.concat(parts, ignore_index=True)
    d = pd.to_datetime(f["date"])
    f["year"], f["month"] = d.dt.year, d.dt.month
    f = f[f["year"].between(FIRST_YEAR, LAST_YEAR)]
    f["reach_m"] = f["sensor"].map(REACH_M)
    # rows stay in archive order (products in FIRMS_RAW order):
    # mi_dnbr_burn.py numbers same-day fire events in this order
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
def load(sensors=FIRMS_SENSORS, crs: str = MI_CRS) -> gpd.GeoDataFrame:
    """Detections, one row per hot pixel. Columns: date (datetime.date),
    time_utc, year, month, sensor, product, satellite, daynight, confidence,
    frp_mw, reach_m, geometry."""
    f = gpd.read_file(ACTIVE_FIRES, layer=LAYER)
    f = f[f["sensor"].isin(sensors)]
    f["date"] = pd.to_datetime(f["date"]).dt.date
    return f.to_crs(crs).reset_index(drop=True)


def reach_shape(point, reach_m: float):
    """The area a detection can stand for: a circle of radius reach_m."""
    return point.buffer(reach_m)


def detections_near(geometry, fires: gpd.GeoDataFrame, year: int | None = None) -> gpd.GeoDataFrame:
    """Detections within their reach of `geometry`, with the distance to the
    polygon (0 = inside) in dist_m."""
    f = fires if year is None else fires[fires["year"] == year]
    f = f.iloc[f.sindex.query(geometry.buffer(max(REACH_M.values())), predicate="intersects")]
    if f.empty:
        return f.assign(dist_m=pd.Series(dtype=float))
    f = f.assign(dist_m=f.distance(geometry))
    return f[f["dist_m"] <= f["reach_m"]]


def reach_mask(xs: np.ndarray, ys: np.ndarray, pts: gpd.GeoDataFrame) -> np.ndarray:
    """Pixels (centre coordinates xs, ys) within reach of any detection."""
    m = np.zeros(xs.shape, bool)
    for p, r in zip(pts.geometry, pts["reach_m"]):
        m |= np.hypot(xs - p.x, ys - p.y) <= r
    return m


if __name__ == "__main__":
    if (sys.argv[1] if len(sys.argv) > 1 else "build") == "build":
        build()
