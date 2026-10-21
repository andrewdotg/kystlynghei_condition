# Kystlynghei condition

Annual, national assessment of the condition of mapped Norwegian coastal
heathland (*kystlynghei*, 14,330 MI polygons, 2018-2025). A central model
turns each 10 m pixel's annual AlphaEarth Foundations embedding into a
continuous succession score (**AE-SJ**: 0 = median of field class SJ4,
1 = median of SJ1). Modules use the score, and independent satellite data,
to assess trajectories, grazing, burning, scrub and tree clearance, and loss
to development.

**Report:** [Coastal heathland
succession](https://andrewdotg.github.io/kystlynghei_condition/) (version 0.1.0,
draft for internal review). The rendered report is `docs/index.html`, published
with GitHub Pages; its source is
`report/kystlynghei_condition_report.qmd`. After re-rendering, copy
`report/kystlynghei_condition_report.html` to `docs/index.html`.

The methodology is a separate document (Word).

## Pipeline

```
inputs/                     MI polygons, active-fire detections, Picea mask, naturtiltak burns
   |
R/model.R                   fits the AE-SJ axis -> outputs/model/ (model JSON read by python/ae_sj_model.py)
   |
python/mi_aesj_annual_products.py      national AE-SJ scoring 2018-2025, pixel transitions
python/active_fires.py                 active fires: FIRMS (4 products) -> one input
python/mi_dnbr_burn.py                 burns: FIRMS-triggered Sentinel-2 dNBR mapping;
                                       AE-SJ response of burnt pixels
python/validation_registered_burns.py  burn detection tested against naturtiltak burning
   |
R/modules/                  temporal_trajectories.R -> grazing.R -> scrub_clearance.R,
                            development_loss.R, SJscore_validation.R,
                            unified_site_summary.R      -> outputs/<module>/
R/naturTiltak/              AE-SJ response to registered burning and clearance
   |
report/                     quarto render report/kystlynghei_condition_report.qmd
```

Run order and run times are in `run_all.R`.

**Software.** Python: `conda env create -f environment.yml` creates the
`kystlynghei` environment with the package versions used (Python 3.12).
R 4.4: the packages and versions used are listed in `R-packages.txt`.

## Methods in brief

**Succession score.** Covariance-adjusted (shrinkage whitening) AlphaEarth
vectors projected on the SJ4 -> SJ1 prototype axis and rescaled to the field
class medians; pixels classed into SJ-like classes at the midpoints between
class medians (`R/model.R`, Methodology Section 3).

**Trajectories.** Theil-Sen trend of the annual polygon median, after removing
regional (county) year effects, classed with thresholds anchored in the field
classes. Each trend has a 90 % confidence interval, with the residual variance
pooled across similar sites (county x starting SJ class, empirical Bayes), and
a confidence level (`R/modules/temporal_trajectories.R`, Methodology Sections 5
and 12).

**Burns** (`python/mi_dnbr_burn.py`, Methodology Section 7). Independent of
the AE-SJ score:

1. FIRMS detections (VIIRS NOAA-20, Suomi NPP, NOAA-21; MODIS) within 750 m
   (VIIRS) or 1,000 m (MODIS) of an MI polygon, grouped into fire events
   (2 km, 30 days).
2. Per event: per-pixel median NBR/NDVI from Sentinel-2 L2A (Planetary
   Computer), same season of the previous year vs. 1-60 days after the first
   detection; SCL 4/5 clear pixels, eroded; tiles mosaicked per date.
3. dNBR and dNDVI corrected by the unburnt mode of the heath pixels; bright
   pixels masked.
4. Candidate: dNBR >= 0.075, dNDVI >= 0, pre-fire NDVI >= 0.15, post-fire
   NBR <= 0.30; regions with >= 4 pixels of dNBR >= 0.20; 1-pixel closing;
   patches >= 0.25 ha.
5. A patch within reach of the event's detections is a burn; a polygon with
   >= 0.1 ha of burn is burnt that year. Fire type: any June-August detection
   = probable wildfire, otherwise lyngbrenning season (probable prescribed).

Stage 2 scores AE-SJ 2017-2025 for each burnt polygon and compares burnt with
never-burnt pixels (`MI_burn_aesj_trajectories.csv`), which the report uses
for the effect of burning on succession by fire type.

**Scrub and tree clearance.** SJ3/SJ4 -> SJ1/SJ2 pixels between 2018-2019
and 2024-2025 (`python/mi_aesj_annual_products.py`). SJ3, non-spruce pixels
inside a mapped burn are "burn-coincident"; SJ4 and spruce pixels count as
burn-coincident only inside a probable wildfire, otherwise as clearance. Headline population: restoring sites (`R/modules/scrub_clearance.R`).

**Validation.** Cross-validation and field agreement (`R/validation.R`,
`R/modules/SJscore_validation.R`); registered clearance (`R/naturTiltak/`);
registered burning (`python/validation_registered_burns.py`: burns-verified
naturtiltak pool, FIRMS-triggered mapping, 5 km controls).

## Main outputs

| Path | Content |
|---|---|
| `python/mi_aesj_annual_results/` | annual polygon statistics and pixel transitions |
| `python/mi_dnbr_burn_results/MI_burns.csv` | one row per burnt site-year: area, dates, satellites, fire type |
| `python/mi_dnbr_burn_results/MI_firms_events.csv`, `MI_firms_detections.csv` | fire events and their detections |
| `python/mi_dnbr_burn_results/events/`, `rasters/` | per-event dNBR layers; per-polygon burn masks on the AE-SJ grid |
| `python/mi_dnbr_burn_results/MI_burn_aesj_trajectories.csv` | AE-SJ of burnt and never-burnt pixels by year |
| `python/validation_registered_burns/` | `validation_final.json`, fire and feature tables, threshold grid, spatial representativeness |
| `outputs/unified_site_summary/` | one record per MI polygon (CSV, polygon and point layers), with a data dictionary and QGIS metadata (`.qmd`) for the layers |

## Repository layout

- `R/` - model, validation, shared setup (`00_setup.R`), theme (`theme.R`).
- `R/modules/` - MI-wide analyses. `grazing.R` and `scrub_clearance.R` use
  in-memory objects from `temporal_trajectories.R` and must run after it in
  the same session.
- `R/naturTiltak/` - analyses of the naturtiltak register.
- `python/` - national scoring, active fires, burn mapping, burn validation;
  `python/diagnostics/` - figure scripts for the report and methodology.

`R/00_setup.R` resolves `REPO_ROOT` (this folder) and `PROJECT_ROOT` (the
data folder above it; override with `KYSTLYNGHEI_PROJECT_ROOT`).

## Licence

GNU General Public License v3.0 (`LICENSE`).
