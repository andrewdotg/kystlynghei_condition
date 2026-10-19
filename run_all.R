# =============================================================================
# run_all.R - run the whole pipeline end to end
# =============================================================================
#
# Convenience wrapper, not a requirement - every script below can also be
# opened and run/sourced individually. Order matters: each module after
# model.R reads files that model.R (and, for temporal_trajectories.R,
# grazing.R and modules/scrub_clearance.R) writes out.
#
# PYTHON STEPS (run before the R modules; environment: miniforge3 "openeo"):
#   1. python/mi_aesj_annual_products.py   RUN_MODE = "all"
#        national annual AE-SJ scores and pixel transitions (~15 h)
#   2. python/active_fires.py build
#        only when the FIRMS archives are updated: the single active-fire
#        input inputs/active_fires/active_fires_2018_2025.gpkg
#   3. python/mi_dnbr_burn.py map
#        FIRMS-triggered Sentinel-2 dNBR burn mapping (~1.5 h)
#   4. python/mi_aesj_annual_products.py   RUN_MODE = "refresh_burnt"
#        re-scores the burnt polygons so the clearance/burn split uses the
#        burn layer (or run step 1 after step 3 instead)
#   5. python/mi_dnbr_burn.py aesj
#        AE-SJ before and after each burn, burnt vs never-burnt pixels
#   6. python/validation_registered_burns.py   (Methodology Section 11.4)
#        burns-verified naturtiltak pool, FIRMS-triggered mapping
#   7. python/diagnostics/wildfire_weather.py  (report discussion)
#
# Every source() after the first uses file.path(REPO_ROOT, ...), because
# 00_setup.R changes the working directory to PROJECT_ROOT.
#
# R/modules/scrub_clearance.R and grazing.R must run in the same R session
# right after temporal_trajectories.R - they use in-memory objects that
# script builds rather than saves.
# =============================================================================

source("R/00_setup.R")

# 1. Central model - fits the AE-SJ condition axis, writes the model bundle,
#    reference CSVs and outputs/model/FINAL_AE_SJ_covadj_model.json, the
#    fitted model that the Python scoring (python/ae_sj_model.py) reads.
source(file.path(REPO_ROOT, "R", "model.R"))

# 2. Model validation: cross-validation and leave-one-year-out (Methodology
#    Section 11.1, Tables 3-4) and the comparison of model variants. Only
#    needs re-running when the model is refitted, not for annual updates.
# source(file.path(REPO_ROOT, "R", "validation.R"))

# 3. Add-on modules - each analyses the MI polygons, reading model.R's output.
source(file.path(REPO_ROOT, "R", "modules", "temporal_trajectories.R"))
source(file.path(REPO_ROOT, "R", "modules", "grazing.R"))            # depends on objects from the script above
source(file.path(REPO_ROOT, "R", "modules", "scrub_clearance.R"))    # depends on objects from temporal_trajectories.R above
source(file.path(REPO_ROOT, "R", "modules", "development_loss.R"))
source(file.path(REPO_ROOT, "R", "modules", "SJscore_validation.R"))  # field SJ vs. modal AE-SJ, same site-year
source(file.path(REPO_ROOT, "R", "modules", "unified_site_summary.R"))  # one row per site; reads the burn layer

# =============================================================================
# 4. naturTiltak validation scripts (Methodology Sections 7 and 11)
# =============================================================================
#
# Not part of the MI pipeline - these test the AE-SJ response to registered
# management (burning, scrub and tree clearance) in naturTiltak polygons.
# Nothing downstream depends on them.

source(file.path(REPO_ROOT, "R", "naturTiltak", "naturTiltak_burn.R"))     # AE-SJ response to registered burns (Methodology Figure 3)

# source(file.path(REPO_ROOT, "R", "naturTiltak", "naturTiltak_scrub.R"))  # registered clearance (Methodology Table 6); interactive (file.choose()) - run separately

# 5. Report: quarto render report/kystlynghei_condition_report.qmd
