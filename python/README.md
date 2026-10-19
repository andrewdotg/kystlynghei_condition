# Python: national scoring and burn detection

Environment: `~/miniforge3/envs/openeo/bin/python` (aef-loader, rasterio,
geopandas, scipy, odc-geo, pystac-client, planetary-computer, xarray,
netCDF4, requests).

| Script | What it does | Run time |
|---|---|---|
| `ae_sj_model.py` | AE-SJ scoring function; reads the fitted model from `outputs/model/FINAL_AE_SJ_covadj_model.json`, written by `R/model.R` | - |
| `aef_geobox.py` | 10 m grid snapped to a global 0,0 lattice, so every product aligns pixel for pixel | - |
| `mi_aesj_annual_products.py` | AlphaEarth (source.coop) -> AE-SJ for every MI polygon and year 2018-2025; annual statistics and 2018-19 -> 2024-25 transitions, split into clearance and burn-coincident with the burn layer. `RUN_MODE = "all"` (national) or `"refresh_burnt"` (re-scores only polygons with burns, after `mi_dnbr_burn.py`) | ~15 h / ~15 min |
| `active_fires.py` | `python active_fires.py build` writes the single active-fire input `inputs/active_fires/active_fires_2018_2025.gpkg` from the FIRMS archives in `../Burning/` (VIIRS NOAA-20, Suomi NPP, NOAA-21, MODIS; TYPE 0, MODIS confidence >= 30). Also the loader and the reach rule (750 m VIIRS, 1,000 m MODIS) used by the scripts below | minutes |
| `mi_dnbr_burn.py map` | Groups FIRMS detections linked to MI polygons into fire events and maps each with Sentinel-2 dNBR (Planetary Computer); writes `mi_dnbr_burn_results/` (`MI_burns.csv`, events, patches, rasters). Resumes where it stopped | ~2 h |
| `mi_dnbr_burn.py reclassify` | Re-applies the thresholds to the saved event layers (after a threshold change) | minutes |
| `mi_dnbr_burn.py aesj` | AE-SJ 2017-2025 of burnt vs never-burnt pixels per burn -> `MI_burn_aesj_trajectories.csv` | ~15 min |
| `validation_registered_burns.py` | Burn detection tested against naturtiltak burning: `pool` (burns-verified pool from the grant register and its final reports -> `inputs/naturtiltak/tiltak_burns_verified.gpkg`), `map` (Sentinel-2 layers for every FIRMS fire and a 5 km control), `thresholds`, `final`, `spatial` -> `validation_registered_burns/` | ~15 min |
| `diagnostics/burn_example_figure.py` | Report figure of one burn (`report/fig_burn_example_<fid>_<year>_<mode>.png`) | minutes |
| `diagnostics/methodology_workflow_figure.py` | Methodology Figure 1 | seconds |
| `diagnostics/wildfire_weather.py` | ERA5 weather before each probable wildfire (report discussion) | minutes |

Order: `mi_aesj_annual_products.py` (all) once; `active_fires.py build`
when the FIRMS archives are updated; then
`mi_dnbr_burn.py map`, `mi_aesj_annual_products.py` (refresh_burnt),
`mi_dnbr_burn.py aesj`, `validation_registered_burns.py` (all stages) and
`diagnostics/wildfire_weather.py`.

Planetary Computer returns HTTP 403 when hit too hard: keep
`mi_dnbr_burn.N_WORKERS` at 2 and run only one Sentinel-2 job at a time. The
AlphaEarth jobs read from source.coop and can run alongside.

The burn method is described in `mi_dnbr_burn.py`'s module docstring and in
Methodology Section 7; the validation in `validation_registered_burns.py` and
Methodology Section 11.4.
