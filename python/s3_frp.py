"""s3_frp.py - Sentinel-3 SLSTR active fires (SL_2_FRP) over the Norwegian coast.

A second active-fire record alongside FIRMS, used for the validation against
registered burns only (active_fires.py). SLSTR on
Sentinel-3A and 3B detects fires in 1 km pixels (MWIR, with a 500 m SWIR
channel for hot fires by day), about four overpasses a day at 60 N between
the two satellites, mostly mid-morning and late evening (local time). The
non-time-critical (NT) archive on the Copernicus Data Space (CDSE) starts in
mid-2020.

Only the fire list in each granule (FRP_in.nc, ~0.5 MB) is downloaded, not
the full product. Fires inside the search area (30 km around every
naturtiltak burning polygon and MI site) are kept, with all quality fields,
as one CSV per month in CACHE (outside the repo). `python active_fires.py
build` then adds the class-1 (vegetation fire) pixels to the project's single
active-fire input.

    python s3_frp.py download   catalogue search + FRP_in.nc fires (resumable)

Login: CDSE username/password from ~/cdse.netrc
(machine identity.dataspace.copernicus.eu login ... password ...).
"""
from __future__ import annotations

import netrc
import sys
import tempfile
import threading
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
import shapely
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
CACHE = REPO.parent / "Burning" / "S3" / "sentinel3_frp_granules"   # one CSV of fires per month + ledger
NETRC = Path.home() / "cdse.netrc"
IDP = "identity.dataspace.copernicus.eu"
TOKEN_URL = f"https://{IDP}/auth/realms/CDSE/protocol/openid-connect/token"
CATALOGUE = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD = "https://download.dataspace.copernicus.eu/odata/v1/Products"
START, END = date(2020, 6, 1), date(2026, 1, 1)
SEARCH_BUFFER_M = 30_000
N_WORKERS = 4                            # CDSE allows 4 concurrent downloads
MI_CRS = "EPSG:25833"
KEEP = ["time", "latitude", "longitude", "FRP_MWIR", "FRP_SWIR", "confidence", "classification",
        "used_channel", "n_cloud", "n_water", "n_SWIR_fire", "IFOV_area"]


# ---------------------------------------------------------------- search area
def search_area():
    import validation_registered_burns as v
    t = gpd.read_file(v.TILTAK_SRC).to_crs(MI_CRS)
    mi = gpd.read_file(REPO / "inputs" / "MI_kystlynghei.shp")
    mi = (mi.set_crs(MI_CRS) if mi.crs is None else mi).to_crs(MI_CRS)
    g = pd.concat([t.geometry, mi.geometry])
    g = g[g.notna() & ~g.is_empty]
    full = gpd.GeoSeries([g.buffer(SEARCH_BUFFER_M).union_all().simplify(1_000)], crs=MI_CRS).to_crs(4326).iloc[0]
    coarse = gpd.GeoSeries([g.buffer(SEARCH_BUFFER_M).union_all().simplify(20_000)], crs=MI_CRS).to_crs(4326).iloc[0]
    return full, coarse


# ---------------------------------------------------------------- auth
class Token:
    def __init__(self):
        a = netrc.netrc(NETRC).authenticators(IDP)
        self._user, self._pw = a[0], a[2]
        self._tok, self._exp = None, 0
        self._lock = threading.Lock()

    def get(self):
        with self._lock:
            if time.time() > self._exp - 60:
                r = requests.post(TOKEN_URL, data={"grant_type": "password", "client_id": "cdse-public",
                                                   "username": self._user, "password": self._pw}, timeout=60)
                r.raise_for_status()
                j = r.json()
                self._tok, self._exp = j["access_token"], time.time() + j["expires_in"]
            return self._tok


# ---------------------------------------------------------------- catalogue
def search_month(aoi_wkt: str, y: int, m: int) -> pd.DataFrame:
    a = date(y, m, 1)
    b = date(y + (m == 12), m % 12 + 1, 1)
    q = ("Collection/Name eq 'SENTINEL-3' and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq "
         "'productType' and att/OData.CSC.StringAttribute/Value eq 'SL_2_FRP___') "
         f"and OData.CSC.Intersects(area=geography'SRID=4326;{aoi_wkt}') "
         f"and ContentDate/Start ge {a}T00:00:00.000Z and ContentDate/Start lt {b}T00:00:00.000Z")
    rows, skip = [], 0
    while True:
        u = CATALOGUE + "?" + urllib.parse.urlencode({"$filter": q, "$top": 1000, "$skip": skip,
                                                      "$orderby": "ContentDate/Start"})
        for attempt in range(5):
            try:
                v = requests.get(u, timeout=180).json()["value"]
                break
            except Exception:  # noqa: BLE001
                time.sleep(15 * (attempt + 1))
        else:
            raise RuntimeError(f"catalogue failed {y}-{m:02d}")
        rows += [{"id": x["Id"], "name": x["Name"], "online": x.get("Online", True)} for x in v]
        if len(v) < 1000:
            break
        skip += 1000
    df = pd.DataFrame(rows, columns=["id", "name", "online"])
    if df.empty:
        return df
    # one granule per platform + sensing start: prefer NT over NR, then the latest processing
    df["platform"] = df.name.str[:3]
    df["start"] = df.name.str[16:31]
    df["created"] = df.name.str[48:63]
    df["timeliness"] = df.name.str.extract(r"_(NT|NR|ST)_\d{3}\.SEN3")[0]
    df["rank"] = df.timeliness.map({"NT": 0, "ST": 1, "NR": 2}).fillna(3)
    df = df.sort_values(["rank", "created"], ascending=[True, False])
    return df.drop_duplicates(["platform", "start"]).sort_values("start").reset_index(drop=True)


# ---------------------------------------------------------------- download
_NC_LOCK = threading.Lock()


def fetch_fires(tok: Token, pid: str, name: str, aoi) -> pd.DataFrame:
    u = f"{DOWNLOAD}({pid})/Nodes({name})/Nodes(FRP_in.nc)/$value"
    for attempt in range(6):
        try:
            s = requests.Session()
            s.headers["Authorization"] = "Bearer " + tok.get()
            r = s.get(u, allow_redirects=False, timeout=180)
            while r.status_code in (301, 302, 303, 307, 308):
                r = s.get(r.headers["Location"], allow_redirects=False, timeout=180)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            break
        except Exception:  # noqa: BLE001
            time.sleep(20 * (attempt + 1))
    else:
        raise RuntimeError(f"download failed {name}")
    with _NC_LOCK, tempfile.NamedTemporaryFile(suffix=".nc") as tmp:   # netCDF/HDF5 is not thread-safe
        tmp.write(r.content)
        tmp.flush()
        with xr.open_dataset(tmp.name) as ds:
            if ds.sizes.get("fires", 0) == 0:
                return pd.DataFrame(columns=KEEP + ["platform", "granule"])
            df = ds[[k for k in KEEP if k in ds and ds[k].dims == ("fires",)]].load().to_dataframe()
    df = df.reset_index(drop=True)
    x0, y0, x1, y1 = aoi.bounds
    df = df[df.longitude.between(x0, x1) & df.latitude.between(y0, y1)]
    if len(df):
        df = df[shapely.contains_xy(aoi, df.longitude.values, df.latitude.values)]
    df["platform"] = name[:3]
    df["granule"] = name
    return df


def run_download():
    CACHE.mkdir(parents=True, exist_ok=True)
    full, coarse = search_area()
    tok = Token()
    y, m = START.year, START.month
    while date(y, m, 1) < END:
        tag = f"{y}-{m:02d}"
        done_f, fires_f = CACHE / f"{tag}_granules.csv", CACHE / f"{tag}_fires.csv"
        if done_f.exists() and fires_f.exists():
            y, m = (y + (m == 12), m % 12 + 1)
            continue
        g = search_month(coarse.wkt, y, m)
        out, ok = [], []
        with ThreadPoolExecutor(N_WORKERS) as ex:
            futs = {ex.submit(fetch_fires, tok, r.id, r.name, full): r for r in g.itertuples()}
            for f in as_completed(futs):
                r = futs[f]
                d = f.result()
                ok.append({"name": r.name, "status": "missing" if d is None else "ok",
                           "n_fires": 0 if d is None else len(d)})
                if d is not None and len(d):
                    out.append(d)
        fires = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=KEEP + ["platform", "granule"])
        fires.to_csv(fires_f, index=False)
        pd.DataFrame(ok).to_csv(done_f, index=False)
        print(f"{tag}: {len(g)} granules, {sum(o['status'] == 'missing' for o in ok)} missing, "
              f"{len(fires)} fires in area", flush=True)
        y, m = (y + (m == 12), m % 12 + 1)


# ---------------------------------------------------------------- products
def solar_elevation(t: pd.Series, lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Approximate solar elevation (deg), enough to tell day from night."""
    doy = t.dt.dayofyear.values
    hour = (t.dt.hour + t.dt.minute / 60).values
    decl = np.radians(23.44) * np.sin(2 * np.pi * (284 + doy) / 365)
    ha = np.radians(15 * (hour + lon / 15 - 12))
    la = np.radians(lat)
    return np.degrees(np.arcsin(np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)))


def load_s3(crs: str = MI_CRS) -> gpd.GeoDataFrame:
    """All cached Sentinel-3 fires, one row per fire pixel, in `crs`."""
    parts = [pd.read_csv(f) for f in sorted(CACHE.glob("*_fires.csv"))]
    parts = [p for p in parts if len(p)]
    df = pd.concat(parts, ignore_index=True)
    df["time"] = pd.to_datetime(df["time"], utc=True)
    df["date"] = df["time"].dt.tz_convert("Europe/Oslo").dt.date
    df["year"] = df["time"].dt.year
    df["month"] = df["time"].dt.month
    df["daynight"] = np.where(solar_elevation(df["time"], df.latitude.values, df.longitude.values) > 0, "D", "N")
    df["satellite"] = df.platform.map({"S3A": "Sentinel-3A", "S3B": "Sentinel-3B"})
    df["sensor"] = "SLSTR"
    g = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs=4326)
    return g.to_crs(crs)


if __name__ == "__main__":
    run_download()
