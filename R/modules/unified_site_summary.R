# =============================================================================
# KYSTLYNGHEI MI - UNIFIED SITE SUMMARY (trajectory + scrub/tree + burn +
# picea, one row per MI site)
# =============================================================================
#
# Joins the per-site results into one table, one row per MI polygon. All
# inputs are files, so this module can be re-run on its own:
#   - trajectories and pixel transitions (scrub, tree, spruce, burn split),
#     for every site: outputs/temporal_trajectories/07_site_level_results_ALL.csv
#     (R/modules/temporal_trajectories.R)
#   - burns: python/mi_dnbr_burn_results/MI_burns.csv (one row per burnt
#     site-year: dNBR-mapped burnt area, FIRMS dates, provisional fire type)
#   - spruce cover: inputs/kystlynghei_picea_points.shp (site-wide %)
#
# OUTPUT: one CSV with full column names and two shapefiles with the same
# attributes (abbreviated field names): points for national maps, polygons
# for detailed viewing. A data dictionary (CSV column, shapefile field,
# description, units) is written beside them, and each shapefile gets a QGIS
# metadata file (.qmd, loaded by QGIS with the layer) listing its fields.
#   outputs/unified_site_summary/MI_unified_site_summary.csv
#   outputs/unified_site_summary/MI_unified_site_summary_polygons.shp (+ .qmd)
#   outputs/unified_site_summary/MI_unified_site_summary_points.shp (+ .qmd)
#   outputs/unified_site_summary/MI_unified_site_summary_data_dictionary.csv
# =============================================================================


# =============================================================================
# 0. SHARED SETUP
# =============================================================================

if (!exists("REPO_ROOT")) source("R/00_setup.R")
required_packages(c("readr", "dplyr", "sf", "stringr"))


# =============================================================================
# 1. SETTINGS
# =============================================================================

STATE_TRENDS_CSV <- file.path(REPO_ROOT, "outputs", "temporal_trajectories", "07_site_level_results_ALL.csv")

PICEA_SHAPEFILE <- file.path(REPO_ROOT, "inputs", "kystlynghei_picea_points.shp")

BURNS_CSV <- file.path(REPO_ROOT, "python", "mi_dnbr_burn_results", "MI_burns.csv")

MI_SHAPEFILE <- file.path(REPO_ROOT, "inputs", "MI_kystlynghei.shp")

OUTPUT_DIR <- file.path(OUTPUTS_ROOT, "unified_site_summary")

# Priority flags for targeting measures (report Discussion):
#   "Poor condition, open heath remaining": tilstand Dårlig / Svært redusert,
#     with at least PRIORITY_OPEN_MIN of the site SJ1/SJ2-like in the final
#     period (2024-2025 AE-SJ composition, p1_final + p2_final).
#   "Good/moderate condition, overgrowing": tilstand God / Moderat, with a
#     Moderate or Strong succession trajectory.
PRIORITY_OPEN_MIN <- 0.5
PRIORITY_LABELS <- c(poor_open = "Poor condition, open heath remaining",
                     good_overgrowing = "Good/moderate condition, overgrowing")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 2. TRAJECTORY (ALL SITES)
# =============================================================================

if (!file.exists(STATE_TRENDS_CSV)) {
  stop(
    "Could not find ", STATE_TRENDS_CSV,
    " - run R/modules/temporal_trajectories.R first (it writes this file)."
  )
}

# sj_change is the Theil-Sen trend magnitude (sj_slope * year_span), the
# quantity behind triage_class; sj_base and sj_final are the endpoint values.
trajectory <- readr::read_csv(STATE_TRENDS_CSV, show_col_types = FALSE) %>%
  transmute(
    fid = trimws(as.character(fid)),
    tilstand = as.character(tilstand),
    triage_class = as.character(triage_class),
    sj_base = as.numeric(sj_base),
    sj_final = as.numeric(sj_final),
    sj_change = as.numeric(sj_change),
    # 90 % interval on sj_change and the confidence flag from it
    # (temporal_trajectories.R header items 1 and 6).
    sj_change_se = as.numeric(sj_change_se),
    sj_change_lo = as.numeric(sj_change_lo),
    sj_change_hi = as.numeric(sj_change_hi),
    trend_confidence = as.character(trend_confidence),
    model_se_mean = as.numeric(model_se_mean),
    n_years = as.numeric(n_years),
    open_final = as.numeric(p1_final) + as.numeric(p2_final)
  ) %>%
  mutate(
    pri_poor_open = tilstand %in% c("Dårlig", "Svært redusert") & coalesce(open_final >= PRIORITY_OPEN_MIN, FALSE),
    pri_good_overgrowing = tilstand %in% c("God", "Moderat") &
      triage_class %in% c("Moderate succession", "Strong succession"),
    priority = case_when(
      pri_poor_open ~ PRIORITY_LABELS[["poor_open"]],
      pri_good_overgrowing ~ PRIORITY_LABELS[["good_overgrowing"]],
      TRUE ~ NA_character_
    ),
    # Confidence in the site's trend, shown with the priority flag.
    priority_confidence = if_else(is.na(priority), NA_character_, trend_confidence)
  )


# =============================================================================
# 3. SCRUB/TREE CLEARANCE + PICEA (every site with valid pixel-transition
#    data; sites without are "none detected" after the join in section 6)
# =============================================================================

# Every site is included, whatever its trajectory class: a site can have
# pixel-level clearance or spruce removal without a restoring trend overall.
# (The headline clearance figures in scrub_clearance.R use only restoring
# sites.) transition_qc_ok keeps sites with enough valid transition pixels
# for the counts to be trusted.
state_trends_raw <- readr::read_csv(STATE_TRENDS_CSV, show_col_types = FALSE)

scrub_tree_raw <- state_trends_raw %>%
  filter(coalesce(transition_qc_ok, FALSE))

# Clearance type from the scrub and tree pixel counts.
clear_type_from_counts <- function(n_sj3, n_sj4) {
  case_when(
    n_sj3 > 0 & n_sj4 > 0 ~ "shrub and tree removal",
    n_sj3 > 0 ~ "shrub removal",
    n_sj4 > 0 ~ "tree removal",
    TRUE ~ NA_character_
  )
}

scrub_tree <- scrub_tree_raw %>%
  transmute(
    fid = trimws(as.character(fid)),
    # scrub = SJ3-source, non-spruce; trees = SJ4-source or spruce
    clear_type = clear_type_from_counts(n_scrub_to_SJ12, n_tree_to_SJ12),
    n_sj3_12 = n_SJ3_to_SJ12, n_sj4_12 = n_SJ4_to_SJ12, n_sj34_12 = n_SJ34_to_SJ12,
    # Area by source class (10 m pixels = 0.01 ha).
    ha_sj3_12 = n_SJ3_to_SJ12 * 0.01, ha_sj4_12 = n_SJ4_to_SJ12 * 0.01, ha_sj34_12 = area_SJ34_to_SJ12_ha,
    prop_sj34 = if_else(n_transition_pixels > 0, n_SJ34_to_SJ12 / n_transition_pixels, NA_real_),
    tree_removed = n_tree_to_SJ12 > 0,
    ha_scrub_12 = area_scrub_to_SJ12_ha, ha_tree_12 = area_tree_to_SJ12_ha,

    # Spruce removal: transitioning pixels that lie in the national Picea
    # mask (inputs/picea_mask_norway_25833_10m.tif), whatever the mechanism.
    n_sj34pic = n_SJ34_to_SJ12_picea,
    ha_sj34pic = area_SJ34_to_SJ12_picea_ha,
    prop_pic = prop_SJ34_to_SJ12_picea,

    # Burn-coincident vs clearance split, per pixel, from python/
    # mi_aesj_annual_products.py: scrub (SJ3-source, non-spruce) pixels inside
    # any dNBR-mapped burn (2018-2025) and tree (SJ4-source or spruce) pixels
    # inside a probable-wildfire burn are burn-coincident; the rest is clearance.
    n_sj34_brn = n_SJ34_to_SJ12_burnt, n_sj34_clr = n_SJ34_to_SJ12_clearance,
    ha_sj34_brn = area_SJ34_to_SJ12_burnt_ha, ha_sj34_clr = area_SJ34_to_SJ12_clearance_ha,
    n_sj34pic_brn = n_SJ34_to_SJ12_picea_burnt, n_sj34pic_clr = n_SJ34_to_SJ12_picea_clearance,
    ha_sj34pic_brn = area_SJ34_to_SJ12_picea_burnt_ha, ha_sj34pic_clr = area_SJ34_to_SJ12_picea_clearance_ha
  )


# =============================================================================
# 4. PICEA (one point per fid: joined by ID)
# =============================================================================

if (!file.exists(PICEA_SHAPEFILE)) {
  stop("Could not find ", PICEA_SHAPEFILE)
}

picea_raw <- sf::st_read(PICEA_SHAPEFILE, quiet = TRUE) %>%
  sf::st_drop_geometry()

# fid is stored as a decimal number in this file; it is rounded to match.
picea <- picea_raw %>%
  transmute(
    fid = as.character(as.integer(round(fid))),
    picea_pct = as.numeric(picea_pct)
  )


# =============================================================================
# 5. BURNS (one row per burnt site-year in python/mi_dnbr_burn_results/
#    MI_burns.csv, collapsed to one row per site)
# =============================================================================

if (!file.exists(BURNS_CSV)) {
  stop("Could not find ", BURNS_CSV, " - run python/mi_dnbr_burn.py first.")
}

burns <- readr::read_csv(BURNS_CSV, show_col_types = FALSE) %>%
  mutate(fid = trimws(as.character(fid)))

burn_by_site <- burns %>%
  group_by(fid) %>%
  summarise(
    burn_confirmed = TRUE,
    n_burn_years = n(),
    # A site burnt more than once counts each burn (reburning is a separate event).
    burn_area_ha_total = sum(burn_area_ha, na.rm = TRUE),
    most_recent_burn_year = max(year),
    burn_years = paste(sort(unique(year)), collapse = ";"),
    # A site with any summer (June-August) fire is a probable wildfire site.
    fire_type = if (any(probable_fire_type == "wildfire")) "wildfire" else "prescribed_burn",
    .groups = "drop"
  )


# =============================================================================
# 6. JOIN EVERYTHING ONTO MI GEOMETRY
# =============================================================================

mi_sf <- sf::st_read(MI_SHAPEFILE, quiet = TRUE) %>%
  mutate(fid = trimws(as.character(fid))) %>%
  select(fid)

unified <- mi_sf %>%
  left_join(trajectory, by = "fid") %>%
  left_join(scrub_tree, by = "fid") %>%
  left_join(picea, by = "fid") %>%
  left_join(burn_by_site, by = "fid") %>%
  mutate(
    burn_confirmed = coalesce(burn_confirmed, FALSE),
    tree_removed = coalesce(tree_removed, FALSE),
    n_sj34pic = coalesce(n_sj34pic, 0L),
    ha_sj34_12 = coalesce(ha_sj34_12, 0),
    ha_sj34pic = coalesce(ha_sj34pic, 0),
    burn_area_ha_total = coalesce(burn_area_ha_total, 0),

    # Pixel-exact burn-coincident/clearance split from section 3, coalesced
    # to 0 for sites with no transition-pixel data.
    n_sj34_brn = coalesce(n_sj34_brn, 0L), n_sj34_clr = coalesce(n_sj34_clr, 0L),
    ha_sj34_burnt = coalesce(ha_sj34_brn, 0), ha_sj34_clearance = coalesce(ha_sj34_clr, 0),
    n_sj34pic_brn = coalesce(n_sj34pic_brn, 0L), n_sj34pic_clr = coalesce(n_sj34pic_clr, 0L),
    ha_picea_burnt = coalesce(ha_sj34pic_brn, 0), ha_picea_clearance = coalesce(ha_sj34pic_clr, 0),

    # Spruce removed: transitioning pixels overlap the Picea mask.
    # picea_pct (site-wide spruce cover) is separate information.
    picea_removal = n_sj34pic > 0,

    # Which mechanisms are present at the site (burnt, scrub, tree removal,
    # or combinations), for mapping.
    site_type = case_when(
      burn_confirmed & !is.na(clear_type) ~ paste0("burnt + ", clear_type),
      burn_confirmed ~ "burnt",
      !is.na(clear_type) ~ clear_type,
      TRUE ~ "none detected"
    )
  )

cat("\npriority breakdown:\n")
print(sf::st_drop_geometry(unified) %>% count(priority))

cat("\nsite_type breakdown:\n")
print(sf::st_drop_geometry(unified) %>% count(site_type, sort = TRUE))


# =============================================================================
# 7. WRITE OUTPUTS
# =============================================================================

readr::write_csv(sf::st_drop_geometry(unified), file.path(OUTPUT_DIR, "MI_unified_site_summary.csv"))

# Shapefile (DBF) field names are kept <= 10 characters - see the transmute/
# rename below. The CSV above has the full names if these abbreviations
# aren't self-explanatory enough on their own.
unified_shp_fields <- unified %>%
  transmute(
    fid,
    tilstand,
    triage = triage_class,
    trend_conf = trend_confidence,
    sj_base, sj_final, sj_change,
    chg_se = sj_change_se, chg_lo90 = sj_change_lo, chg_hi90 = sj_change_hi,
    mod_se = model_se_mean,
    n_years,
    # Priority flags - see section 1.
    open_fin = open_final,
    pri_poor = pri_poor_open,
    pri_overgr = pri_good_overgrowing,
    priority,
    pri_conf = priority_confidence,
    clear_type,
    n_sj3_12, n_sj4_12, n_sj34_12,
    ha_sj3_12, ha_sj4_12, ha_sj34_12,
    # Pixel-exact burn-coincident/clearance split - see section 3.
    n_34_brn = n_sj34_brn, n_34_clr = n_sj34_clr,
    ha_34_brn = ha_sj34_burnt, ha_34_clr = ha_sj34_clearance,
    prop_sj34,
    picea_pct,
    n_sj34pic, ha_sj34pic, prop_pic,
    # Same pixel-verified split, restricted to the picea overlap.
    n_picbrn = n_sj34pic_brn, n_picclr = n_sj34pic_clr,
    ha_pic_brn = ha_picea_burnt, ha_pic_clr = ha_picea_clearance,
    picea_rem = picea_removal,
    burn_conf = burn_confirmed,
    burn_yrs = n_burn_years,
    brn_years = burn_years,
    burn_ha = burn_area_ha_total,
    burn_last = most_recent_burn_year,
    fire_type,
    site_type
  )

POLYGON_SHP <- file.path(OUTPUT_DIR, "MI_unified_site_summary_polygons.shp")
sf::st_write(unified_shp_fields, POLYGON_SHP, delete_layer = TRUE, quiet = TRUE, layer_options = "ENCODING=UTF-8")
cat("\nPolygon shapefile written to:\n  ", POLYGON_SHP, "\n", sep = "")

POINT_SHP <- file.path(OUTPUT_DIR, "MI_unified_site_summary_points.shp")
sf::st_write(sf::st_point_on_surface(unified_shp_fields), POINT_SHP, delete_layer = TRUE, quiet = TRUE, layer_options = "ENCODING=UTF-8")
cat("Point shapefile written to:\n  ", POINT_SHP, "\n", sep = "")


# =============================================================================
# 8. DATA DICTIONARY AND QGIS METADATA
# =============================================================================

# One row per column: CSV name, shapefile field (NA = CSV only), description,
# units. AE-SJ score: 1 = SJ1-like (open heath), 0 = SJ4-like (late succession).
data_dictionary <- tibble::tribble(
  ~csv_column, ~shp_field, ~description, ~units,
  "fid", "fid", "MI polygon ID (fid in the MI kystlynghei polygons)", "",
  "tilstand", "tilstand", "Field condition assessment from the MI survey (God, Moderat, Dårlig, Svært redusert)", "",
  "triage_class", "triage", "Trajectory class 2018-2025 from sj_change: Strong/Moderate succession, No strong signal, Moderate/Strong restoration. Empty = excluded (heavy grazing BT5/BT6, development, too few years or implausible annual jumps)", "",
  "trend_confidence", "trend_conf", "Confidence in the trajectory from the 90 % interval on sj_change: High = interval wholly beyond the moderate threshold on the side of the class (at least a moderate trend in that direction), Medium = interval excludes zero (direction clear, size uncertain), Low = interval includes zero. Empty for No strong signal and excluded sites", "",
  "sj_base", "sj_base", "Mean AE-SJ score of the first two years (normally 2018-2019), from the annual site medians", "AE-SJ score",
  "sj_final", "sj_final", "Mean AE-SJ score of the last two years (normally 2024-2025), from the annual site medians", "AE-SJ score",
  "sj_change", "sj_change", "Trend-equivalent change over the period: Theil-Sen slope of the annual site median AE-SJ, corrected for regional (county) year effects, x year span. Negative = towards SJ4 (succession), positive = towards SJ1 (restoration)", "AE-SJ score",
  "sj_change_se", "chg_se", "Standard error of sj_change, from the scatter of the year-corrected annual medians around the Theil-Sen line, pooled with similar sites (county x starting SJ class, empirical Bayes)", "AE-SJ score",
  "sj_change_lo", "chg_lo90", "Lower limit of the 90 % interval on sj_change (t distribution, residual plus prior df)", "AE-SJ score",
  "sj_change_hi", "chg_hi90", "Upper limit of the 90 % interval on sj_change", "AE-SJ score",
  "model_se_mean", "mod_se", "Mean over years of the calibration standard error of the site's mean AE-SJ score (field-calibrated spread of each pixel's SJ class, pixels assumed independent; a lower bound). Uncertainty of the score level, not of the trend", "AE-SJ score",
  "n_years", "n_years", "Number of years with a valid annual score", "years",
  "open_final", "open_fin", "Share of the site SJ1- or SJ2-like in 2024-2025", "proportion 0-1",
  "pri_poor_open", "pri_poor", "Priority flag: condition Dårlig or Svært redusert with at least half the site still SJ1/SJ2-like in 2024-2025", "TRUE/FALSE",
  "pri_good_overgrowing", "pri_overgr", "Priority flag: condition God or Moderat with a Moderate or Strong succession trajectory", "TRUE/FALSE",
  "priority", "priority", "Priority group for targeting measures (one of the two flags above; poor/open takes precedence)", "",
  "priority_confidence", "pri_conf", "trend_confidence of priority sites (empty for other sites)", "",
  "clear_type", "clear_type", "Clearance from pixel transitions 2018-2019 to 2024-2025: shrub removal (SJ3-source, not spruce), tree removal (SJ4-source or spruce), or both", "",
  "n_sj3_12", "n_sj3_12", "Pixels SJ3-like in 2018-2019 and SJ1/SJ2-like in 2024-2025", "10 m pixels",
  "n_sj4_12", "n_sj4_12", "Pixels SJ4-like in 2018-2019 and SJ1/SJ2-like in 2024-2025", "10 m pixels",
  "n_sj34_12", "n_sj34_12", "Pixels SJ3/SJ4-like in 2018-2019 and SJ1/SJ2-like in 2024-2025", "10 m pixels",
  "ha_sj3_12", "ha_sj3_12", "Area of n_sj3_12", "ha",
  "ha_sj4_12", "ha_sj4_12", "Area of n_sj4_12", "ha",
  "ha_sj34_12", "ha_sj34_12", "Area of n_sj34_12", "ha",
  "prop_sj34", "prop_sj34", "Share of the site's transition pixels that moved from SJ3/SJ4-like to SJ1/SJ2-like", "proportion 0-1",
  "tree_removed", NA, "Any tree removal (SJ4-source or spruce pixels to SJ1/SJ2)", "TRUE/FALSE",
  "ha_scrub_12", NA, "Scrub removal area (SJ3-source, not spruce)", "ha",
  "ha_tree_12", NA, "Tree removal area (SJ4-source or spruce)", "ha",
  "n_sj34pic", "n_sj34pic", "SJ3/SJ4 to SJ1/SJ2 pixels inside the national Picea mask (spruce removal)", "10 m pixels",
  "ha_sj34pic", "ha_sj34pic", "Area of n_sj34pic", "ha",
  "prop_pic", "prop_pic", "Share of the SJ3/SJ4 to SJ1/SJ2 pixels inside the Picea mask", "proportion 0-1",
  "n_sj34_brn", "n_34_brn", "SJ3/SJ4 to SJ1/SJ2 pixels coinciding with a mapped burn (scrub inside any burn, trees inside a probable-wildfire burn)", "10 m pixels",
  "n_sj34_clr", "n_34_clr", "SJ3/SJ4 to SJ1/SJ2 pixels attributed to clearance (not burn-coincident)", "10 m pixels",
  "ha_sj34_brn", NA, "Area of n_sj34_brn (sites with transition data only)", "ha",
  "ha_sj34_clr", NA, "Area of n_sj34_clr (sites with transition data only)", "ha",
  "n_sj34pic_brn", "n_picbrn", "Spruce-mask transition pixels coinciding with a mapped burn", "10 m pixels",
  "n_sj34pic_clr", "n_picclr", "Spruce-mask transition pixels attributed to clearance", "10 m pixels",
  "ha_sj34pic_brn", NA, "Area of n_sj34pic_brn (sites with transition data only)", "ha",
  "ha_sj34pic_clr", NA, "Area of n_sj34pic_clr (sites with transition data only)", "ha",
  "picea_pct", "picea_pct", "Site-wide spruce cover from the national Picea mask", "%",
  "burn_confirmed", "burn_conf", "At least one burn mapped at the site 2018-2025 (FIRMS-triggered Sentinel-2 dNBR)", "TRUE/FALSE",
  "n_burn_years", "burn_yrs", "Number of years with a mapped burn", "years",
  "burn_area_ha_total", "burn_ha", "Mapped burnt area, summed over burns", "ha",
  "most_recent_burn_year", "burn_last", "Year of the most recent mapped burn", "year",
  "burn_years", "brn_years", "Years with a mapped burn, separated by ;", "",
  "fire_type", "fire_type", "prescribed_burn (all detections outside June-August) or wildfire (any detection in June-August); provisional", "",
  "ha_sj34_burnt", "ha_34_brn", "Area of n_sj34_brn, 0 where no transition data", "ha",
  "ha_sj34_clearance", "ha_34_clr", "Area of n_sj34_clr, 0 where no transition data", "ha",
  "ha_picea_burnt", "ha_pic_brn", "Area of n_sj34pic_brn, 0 where no transition data", "ha",
  "ha_picea_clearance", "ha_pic_clr", "Area of n_sj34pic_clr, 0 where no transition data", "ha",
  "picea_removal", "picea_rem", "Spruce removed: any SJ3/SJ4 to SJ1/SJ2 pixel inside the Picea mask", "TRUE/FALSE",
  "site_type", "site_type", "Mechanisms present at the site: burnt, shrub/tree removal, combinations, or none detected", ""
)

csv_cols <- names(sf::st_drop_geometry(unified))
shp_cols <- setdiff(names(unified_shp_fields), attr(unified_shp_fields, "sf_column"))
stopifnot(
  setequal(data_dictionary$csv_column, csv_cols),
  setequal(na.omit(data_dictionary$shp_field), shp_cols)
)
data_dictionary <- data_dictionary[match(csv_cols, data_dictionary$csv_column), ]
readr::write_csv(data_dictionary, file.path(OUTPUT_DIR, "MI_unified_site_summary_data_dictionary.csv"), na = "")

# QGIS layer metadata (.qmd beside the .shp; QGIS loads it with the layer and
# shows it under Layer Properties > Metadata). The field list is in the abstract.
xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

write_qmd <- function(shp, geometry_kind) {
  dd <- data_dictionary[!is.na(data_dictionary$shp_field), ]
  dd <- dd[match(shp_cols, dd$shp_field), ]
  fields <- paste0(dd$shp_field, " (", dd$csv_column, ")",
                   ifelse(dd$units == "", "", paste0(" [", dd$units, "]")),
                   ": ", dd$description, collapse = "\n")
  abstract <- paste0(
    "AE-SJ succession monitoring of mapped kystlynghei (MI) sites, 2018-2025, one ", geometry_kind,
    " per MI site. AE-SJ score: 1 = SJ1-like (open heath), 0 = SJ4-like (late succession). ",
    "Full column names are in MI_unified_site_summary.csv; field descriptions are also in ",
    "MI_unified_site_summary_data_dictionary.csv.\n\nFields:\n", fields
  )
  qmd <- c(
    "<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>",
    "<qgis version=\"3.28.0\">",
    paste0("  <identifier>", tools::file_path_sans_ext(basename(shp)), "</identifier>"),
    "  <parentidentifier></parentidentifier>",
    "  <language>ENG</language>",
    "  <type>dataset</type>",
    paste0("  <title>Kystlynghei site summary (", geometry_kind, "s)</title>"),
    paste0("  <abstract>", xml_escape(abstract), "</abstract>"),
    "  <keywords vocabulary=\"gmd:topicCategory\">",
    "    <keyword>environment</keyword>",
    "  </keywords>",
    "  <contact>",
    "    <name>Andrew Møller Gray</name>",
    "    <organization>Norwegian Institute for Nature Research (NINA)</organization>",
    "    <position></position><voice></voice><fax></fax>",
    "    <email>andrew.gray@nina.no</email>",
    "    <role>author</role>",
    "  </contact>",
    "  <links/>",
    paste0("  <history>Written by R/modules/unified_site_summary.R on ", Sys.Date(), "</history>"),
    "  <fees></fees>",
    "  <encoding>UTF-8</encoding>",
    "  <crs>",
    "    <spatialrefsys><authid>EPSG:25833</authid></spatialrefsys>",
    "  </crs>",
    "</qgis>"
  )
  writeLines(enc2utf8(qmd), sub("\\.shp$", ".qmd", shp), useBytes = TRUE)
}

write_qmd(POLYGON_SHP, "polygon")
write_qmd(POINT_SHP, "point")
cat("Data dictionary and QGIS metadata (.qmd) written to:\n  ", OUTPUT_DIR, "\n", sep = "")
