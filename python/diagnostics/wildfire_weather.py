"""Weather before each probable wildfire (qualitative supporting evidence).

For every fire event that burnt an MI site and was classed as a probable
wildfire (FIRMS detection in June-August), daily precipitation and maximum
temperature are taken from the ERA5 reanalysis (Open-Meteo historical
weather API, ~25 km grid; no key needed) at the mean position of the event's
detections. For the WINDOW_DAYS before the first detection, the precipitation
total, the number of dry days (< 1 mm) and the mean daily maximum
temperature are compared with the same calendar window in 1991-2020, and the
year is ranked among 1991-2025 (rank 1 = driest).

met.no station data (Frost API) need a registered client id; ERA5 is used
instead and should be read as a regional indication, not a station record.

Writes python/mi_dnbr_burn_results/wildfire_weather.csv.
"""
import json
import time
import urllib.request
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1] / "mi_dnbr_burn_results"
WINDOW_DAYS = 30
NORMAL = (1991, 2020)
URL = ("https://archive-api.open-meteo.com/v1/archive?latitude={lat:.3f}&longitude={lon:.3f}"
       "&start_date=1991-01-01&end_date={end}&daily=precipitation_sum,temperature_2m_max&timezone=Europe%2FOslo")


def fetch(lat, lon, end):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(URL.format(lat=lat, lon=lon, end=end), timeout=60) as r:
                d = json.load(r)["daily"]
            return pd.DataFrame({"date": pd.to_datetime(d["time"]), "p": d["precipitation_sum"],
                                 "tmax": d["temperature_2m_max"]})
        except Exception:  # noqa: BLE001
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("Open-Meteo request failed")


def window(df, end_day: date):
    start = end_day - timedelta(days=WINDOW_DAYS)
    w = df[(df.date.dt.date >= start) & (df.date.dt.date < end_day)]
    return w.p.sum(), int((w.p < 1).sum()), w.tmax.mean()


def main():
    burns = pd.read_csv(R / "MI_burns.csv")
    ev = pd.read_csv(R / "MI_firms_events.csv").set_index("event_id")
    det = pd.read_csv(R / "MI_firms_detections.csv")
    wf = burns[burns.probable_fire_type == "wildfire"]
    ids = sorted({e for s in wf.event_ids for e in s.split(";") if ev.loc[e, "summer"]})
    import pyproj
    to_ll = pyproj.Transformer.from_crs("EPSG:25833", "EPSG:4326", always_xy=True)
    rows = []
    for eid in ids:
        e = ev.loc[eid]
        d = det[det.event_id == eid]
        lon, lat = to_ll.transform(d.x.mean(), d.y.mean())
        first = date.fromisoformat(str(e.first_date))
        df = fetch(lat, lon, f"{max(first.year, 2025)}-12-31")
        p, dry, tmax = window(df, first)
        clim = [window(df, first.replace(year=y)) for y in range(NORMAL[0], NORMAL[1] + 1)]
        allyrs = {y: window(df, first.replace(year=y))[0] for y in range(1991, min(first.year, 2025) + 1)}
        rank = 1 + sum(v < p for v in allyrs.values())
        sites = burns[burns.event_ids.str.contains(eid)]
        rows.append({"event_id": eid, "first_date": first.isoformat(), "lat": round(lat, 3), "lon": round(lon, 3),
                     "sites": len(sites), "burnt_ha": round(sites.burn_area_ha.sum(), 1),
                     "fids": ";".join(map(str, sites.fid)),
                     f"precip_{WINDOW_DAYS}d_mm": round(p, 1),
                     "precip_normal_mm": round(np.mean([c[0] for c in clim]), 1),
                     "precip_pct_of_normal": round(100 * p / np.mean([c[0] for c in clim]), 0),
                     "dry_days": dry, "dry_days_normal": round(np.mean([c[1] for c in clim]), 1),
                     "tmax_mean": round(tmax, 1), "tmax_anomaly": round(tmax - np.mean([c[2] for c in clim]), 1),
                     "precip_rank_1991_on": f"{rank} of {len(allyrs)}"})
        print(rows[-1], flush=True)
        time.sleep(2)
    out = pd.DataFrame(rows)
    out.to_csv(R / "wildfire_weather.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
