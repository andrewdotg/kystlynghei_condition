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
# for detailed viewing.
#   outputs/unified_site_summary/MI_unified_site_summary.csv
#   outputs/unified_site_summary/MI_unified_site_summary_polygons.shp
#   outputs/unified_site_summary/MI_unified_site_summary_points.shp
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
    )
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
    sj_base, sj_final, sj_change,
    n_years,
    # Priority flags - see section 1.
    open_fin = open_final,
    pri_poor = pri_poor_open,
    pri_overgr = pri_good_overgrowing,
    priority,
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
