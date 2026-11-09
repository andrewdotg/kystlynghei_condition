# =============================================================================
# export_render_data.R - copy the data the report reads into report/render_data/
# =============================================================================
#
# report/kystlynghei_condition_report.qmd reads only from report/render_data/
# (a path relative to the qmd), so the report renders anywhere the qmd, this
# folder, report/img/, fig_burn_example_2679_2019_recovery.png, references.bib
# and ambio.csl are copied together. The pipeline outputs it comes from are
# ~515 MB and not in git; this folder is a few MB and is.
#
# Run from the repo root after each pipeline run (run_all.R and the Python
# steps):
#   Rscript R/export_render_data.R
# It overwrites the previous export.
#
# Tables are written as gzip-compressed CSV, which readr::read_csv reads
# directly. Large tables are trimmed to the columns the report uses. Values are
# copied as text, so they read back exactly as from the source files.
#
# The MI polygons are stored only in the form the report uses them (the report
# computes no areas from them):
#   - layer "mi_kystlynghei": the 8 attribute fields, with the geometry made
#     valid and simplified with dTolerance = 5 (the web map outlines);
#   - layer "site_points": one point per site, st_point_on_surface() of the
#     valid full geometry (map markers and the E/N coordinates in the popups).
# Field names are stored as the report reads them from the shapefile
# (områdenav, X7RA.SJ, X7JB.BT).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

stopifnot("Run from the repo root" = file.exists("R/theme.R"))

RD <- file.path("report", "render_data")
unlink(RD, recursive = TRUE)
dir.create(RD, recursive = TRUE)

manifest <- list()
note <- function(file, source, columns) {
  manifest[[length(manifest) + 1]] <<- tibble(file = file, source = source, columns = columns)
}

# Copy a file unchanged, gzip-compressed if it is a CSV.
copy_as_is <- function(src) {
  stopifnot(file.exists(src))
  if (grepl("\\.csv$", src)) {
    dest <- file.path(RD, paste0(basename(src), ".gz"))
    con <- gzfile(dest, "wb")
    writeBin(readBin(src, "raw", n = file.size(src)), con)
    close(con)
  } else {
    dest <- file.path(RD, basename(src))
    file.copy(src, dest, overwrite = TRUE)
  }
  note(basename(dest), src, "all (unchanged)")
}

# Read selected columns as text and write them unchanged.
as_text <- function(src, cols) {
  read_csv(src, col_select = all_of(cols), col_types = cols(.default = col_character()),
           na = character(), show_col_types = FALSE)
}
write_trimmed <- function(df, dest, source) {
  write_csv(df, file.path(RD, dest), na = "NA")
  note(dest, source, paste(names(df), collapse = ", "))
}

# ---- site-level results and model / trajectory summaries ---------------------
for (src in c(
  "outputs/unified_site_summary/MI_unified_site_summary.csv",
  "outputs/temporal_trajectories/07_site_level_results_ALL.csv",
  "outputs/temporal_trajectories/02b_regional_year_effects.csv",
  "outputs/temporal_trajectories/02_triage_threshold_derivation.csv",
  "outputs/temporal_trajectories/grazing/24_grazing_slope_median_bootstrap_ci.csv",
  "outputs/temporal_trajectories/grazing/24_grazing_slope_test_result.csv",
  "outputs/temporal_trajectories/grazing/23_no_grazing_vs_grazing_by_tilstand.csv",
  "outputs/scrub_clearance/SJ34_to_SJ12_restorative/15_restorative_sites_SJ34_to_SJ12.csv",
  "outputs/scrub_clearance/SJ34_to_SJ12_restorative/20_SJ34_to_SJ12_transition_pathways.csv",
  "outputs/sj_score_validation/02_sj_score_validation_confusion_matrix.csv",
  "outputs/model/FINAL_AE_SJ_covadj_model_summary.csv",
  "outputs/model/FINAL_AE_SJ_covadj_class_reference_positions.csv",
  # burns (python/mi_dnbr_burn.py, python/diagnostics/wildfire_weather.py)
  "python/mi_dnbr_burn_results/MI_burns.csv",
  "python/mi_dnbr_burn_results/MI_firms_events.csv",
  "python/mi_dnbr_burn_results/MI_firms_detections.csv",
  "python/mi_dnbr_burn_results/MI_burn_aesj_trajectories.csv",
  "python/mi_dnbr_burn_results/wildfire_weather.csv",
  # registered-burn validation (python/validation_registered_burns.py)
  "python/validation_registered_burns/validation_final.json",
  "python/validation_registered_burns/spatial_representativeness.json"
)) copy_as_is(src)

# ---- trimmed tables ----------------------------------------------------------
annual_files <- sprintf("python/mi_aesj_annual_results/MI_aesj_annual_%d.csv", 2018:2025)
purrr::map_dfr(annual_files, as_text, cols = c("fid", "year", "status", "aesj_median_AE_SJ")) |>
  write_trimmed("aesj_annual.csv.gz", "python/mi_aesj_annual_results/MI_aesj_annual_2018.csv ... _2025.csv")

dev_src <- "inputs/MI_DW_persistent_development_annual.csv"
dev_cols <- names(read_csv(dev_src, n_max = 0, show_col_types = FALSE))
as_text(dev_src, c("fid", grep("^dev_ha_", dev_cols, value = TRUE))) |>
  write_trimmed("MI_DW_persistent_development_annual.csv.gz", dev_src)

as_text("outputs/model/FINAL_AE_SJ_covadj_point_scores.csv", c("sj", "AE_SJ")) |>
  write_trimmed("FINAL_AE_SJ_covadj_point_scores.csv.gz", "outputs/model/FINAL_AE_SJ_covadj_point_scores.csv")

# ---- MI polygons -------------------------------------------------------------
sf::sf_use_s2(FALSE)
mi <- sf::st_read(
  "inputs/MI_kystlynghei.shp",
  query = 'SELECT fid, "områdenav", kommunenum, kommuner, year, "7RA-SJ", "7JB-BT", faktaark FROM MI_kystlynghei',
  quiet = TRUE
)
mi_valid <- sf::st_make_valid(mi)
outlines <- sf::st_simplify(mi_valid, dTolerance = 5, preserveTopology = TRUE)
pts <- suppressWarnings(sf::st_point_on_surface(mi_valid |> select(fid)))
# GeoPackages can't be written in place on the network drive: write locally, then copy
gpkg <- file.path(tempdir(), "mi_kystlynghei.gpkg")
unlink(gpkg)
sf::st_write(outlines, gpkg, layer = "mi_kystlynghei", layer_options = "FID=gpkg_fid", quiet = TRUE)
sf::st_write(pts, gpkg, layer = "site_points", layer_options = "FID=gpkg_fid", append = TRUE, quiet = TRUE)
stopifnot(file.copy(gpkg, file.path(RD, "mi_kystlynghei.gpkg"), overwrite = TRUE))
note("mi_kystlynghei.gpkg", "inputs/MI_kystlynghei.shp",
     paste0("layer mi_kystlynghei: ", paste(setdiff(names(mi), attr(mi, "sf_column")), collapse = ", "),
            " + outline (made valid, simplified 5 m); layer site_points: fid + point on surface"))

# ---- theme -------------------------------------------------------------------
writeLines(c("# Generated by R/export_render_data.R from R/theme.R - edit R/theme.R, not this copy.",
             readLines("R/theme.R", encoding = "UTF-8")),
           file.path(RD, "theme.R"), useBytes = TRUE)
note("theme.R", "R/theme.R", "ggplot theme and colours")

# ---- manifest ----------------------------------------------------------------
m <- bind_rows(manifest) |>
  mutate(size_kb = ceiling(file.size(file.path(RD, file)) / 1024))
writeLines(c(
  "# Report data",
  "",
  paste0("Generated by `R/export_render_data.R` on ", format(Sys.Date()),
         ". Do not edit; re-run the script from the repo root after each pipeline run."),
  "",
  "`report/kystlynghei_condition_report.qmd` reads only these files. Tables are gzip-compressed CSV;",
  "trimmed tables keep only the columns the report uses, with values copied unchanged.",
  "",
  "| File | Source | Columns | Size (kB) |",
  "|---|---|---|---:|",
  sprintf("| `%s` | `%s` | %s | %s |", m$file, m$source, m$columns, format(m$size_kb, big.mark = ","))
), file.path(RD, "README.md"), useBytes = TRUE)

sizes <- file.size(list.files(RD, full.names = TRUE))
cat(sprintf("Wrote %d files to %s: %.1f MB in total, largest %s (%.1f MB)\n",
            length(sizes), RD, sum(sizes) / 1e6,
            basename(list.files(RD)[which.max(sizes)]), max(sizes) / 1e6))
