# =============================================================================
# KYSTLYNGHEI AE-SJ TEMPORAL TRIAGE ANALYSIS
# =============================================================================
#
# Classifies the AE-SJ trajectory of every MI site, 2018-2025, and describes
# how the change is distributed within each site.
#
# 1. Trend: the Theil-Sen slope of the annual site median AE-SJ, expressed as
#    the change over the period (sj_change = slope x year span).
#
# 2. Triage classes, from two thresholds anchored in the field classes
#    (TRIAGE_* settings below):
#       moderate threshold = TRIAGE_MODERATE_FRACTION x |SJ1 - SJ2| gap
#       large    threshold = TRIAGE_LARGE_FRACTION    x |SJ3 - SJ4| gap
#    giving Strong succession | Moderate succession | No strong signal |
#    Moderate restoration | Strong restoration. "Succession" = towards
#    SJ4-like structure, "restoration" = towards SJ1-like structure. "No
#    strong signal" means the moderate threshold was not crossed in either
#    direction; it does not mean the site is unchanged.
#
# 3. Within-site mechanism, separately for succeeding and restoring sites:
#    whether the change is widespread across the site or concentrated in part
#    of it (consistent with a localised disturbance or targeted work), from
#    the SJ-class composition and pixel transitions.
#
# 4. BT5/BT6 (heavy grazing) sites are excluded: heavy grazing drives a
#    grassland trajectory that is not on the succession axis of the score.
#
# 5. Kendall's tau and BH-FDR p-values are written as diagnostics only. With
#    8 annual observations per site, significance testing has little power
#    and is not used to screen sites.
#
# Plot styling comes from R/theme.R. R/modules/grazing.R and
# R/modules/scrub_clearance.R use objects built here and must be sourced in
# the same R session after this script (run_all.R does this).
# =============================================================================


# =============================================================================
# 0. SHARED SETUP
# =============================================================================
#
# 00_setup.R resolves PROJECT_ROOT (or KYSTLYNGHEI_PROJECT_ROOT if set) and
# sets the working directory there, so relative paths below work wherever this
# is run from. theme.R provides the shared theme and the mechanism-block
# plotting helpers (make_composition_change_plot, make_mechanism_scatter,
# etc.).
if (!exists("REPO_ROOT")) source("R/00_setup.R")
source(file.path(REPO_ROOT, "R/theme.R"))

required_packages(c(
  "readr", "dplyr", "tidyr", "purrr",
  "stringr", "tibble", "ggplot2", "sf", "ggrepel"
))


# =============================================================================
# 1. PATHS
# =============================================================================

# python/mi_aesj_annual_products.py writes one CSV per year with both the
# zonal statistics and the SJ-class composition; sections 7 and 19 below
# each read the columns they need from the same files.
ANNUAL_INPUT_DIR <- file.path(REPO_ROOT, "python", "mi_aesj_annual_results")

ANNUAL_REGEX <-
  "^MI_aesj_annual_[0-9]{4}\\.csv$"

COMPOSITION_SEARCH_ROOT <- ANNUAL_INPUT_DIR

COMPOSITION_REGEX <- ANNUAL_REGEX

TRANSITION_SEARCH_ROOT <- ANNUAL_INPUT_DIR

TRANSITION_REGEX <-
  "^MI_aesj_pixel_transitions_2018_2019_to_2024_2025\\.csv$"

# Field-class reference positions from R/model.R.
MODEL_REFERENCE_CSV <- file.path(
  REPO_ROOT,
  "outputs",
  "model",
  "FINAL_AE_SJ_covadj_class_reference_positions.csv"
)

# Built-up development per site (Dynamic World), one row per fid.
DEVELOPMENT_CSV <- file.path(REPO_ROOT, "inputs", "MI_DW_persistent_development_annual.csv")

# MI polygons (the same file the Python pipeline uses).
MI_SHAPEFILE <- file.path(REPO_ROOT, "inputs", "MI_kystlynghei.shp")

# grazing.R (sourced after this script) writes to the same OUTPUT_DIR.
OUTPUT_DIR <- file.path(OUTPUTS_ROOT, "temporal_trajectories")


# =============================================================================
# 2. ANALYSIS SETTINGS
# =============================================================================

MIN_YEARS <- 6

# ---- Triage thresholds ------------------------------------------------------
# Both thresholds are defined as a fraction of the empirical AE-SJ gap
# between a chosen pair of adjacent SJ reference medians. Pair indices refer
# to sj = 1..4 (SJ1 Intakt ... SJ4 Sein gjenvekstsuksesjonsfase).
#
# Defaults reproduce the original screening values: half the SJ1-SJ2 gap for
# "moderate", half the SJ3-SJ4 gap for "large".
TRIAGE_MODERATE_PAIR <- c(1L, 2L)
TRIAGE_LARGE_PAIR <- c(3L, 4L)
TRIAGE_MODERATE_FRACTION <- 0.5
TRIAGE_LARGE_FRACTION <- 0.5

# Minimum valid pixels.
MIN_PIXELS_STATE <- 1
MIN_COMPOSITION_PIXELS <- 5
MIN_TRANSITION_PIXELS <- 5

# Descriptive thresholds used only to turn continuous pixel-transition
# metrics into plain-language pattern classes. They are NOT statistical
# thresholds.
WIDESPREAD_CHANGED_THRESHOLD <- 0.50
DIRECTIONAL_COHERENCE_THRESHOLD <- 0.75

# Number of strongest sites labelled in each ranked figure. Full site lists
# are always written to CSV and shapefile.
SITE_GRAPH_TOP_N_PER_DIRECTION <- 30

# Exclude BT5 and BT6 from ALL temporal analyses: overgrazing drives a
# grassland-bound trajectory, not the forest-succession axis this score
# measures.
EXCLUDE_BT <- c(5L, 6L)

# The supplied development script defines detected development as
# dev_ha_total > 0. Setting this to 0 reproduces that conservative rule.
DEVELOPMENT_EXCLUDE_MIN_HA <- 0

# Maximum plausible absolute year-to-year change in polygon median AE-SJ.
MAX_ABS_ANNUAL_STEP <- 0.50

# Scores outside 0-1 are valid because the final AE-SJ model is unclipped.
REFERENCE_SCORE_MIN <- 0.0
REFERENCE_SCORE_MAX <- 1.0
EXTREME_SCORE_MIN <- -0.5
EXTREME_SCORE_MAX <- 1.5


# =============================================================================
# 3. OUTPUT DIRECTORY
# =============================================================================
# Packages and the shared theme come from R/00_setup.R and R/theme.R.

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 4. GENERAL HELPERS
# =============================================================================

rename_if_needed <- function(dat, old, new) {
  if (old %in% names(dat) && !new %in% names(dat)) {
    names(dat)[names(dat) == old] <- new
  }
  dat
}

find_id_field <- function(x) {
  candidates <- c("fid", "FID", "mi_id", "id", "lokalId", "lokalid")
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) {
    stop(
      "Could not identify a polygon ID field. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

normalise_tilstand <- function(x) {
  z <- trimws(as.character(x))
  z_ascii <- iconv(z, from = "", to = "ASCII//TRANSLIT")
  z_key <- tolower(gsub("[^a-zA-Z0-9]+", "", z_ascii))
  dplyr::case_when(
    z_key == "god" ~ "God",
    z_key == "moderat" ~ "Moderat",
    z_key %in% c("darlig", "daarlig") ~ "Dårlig",
    z_key %in% c("svartredusert", "svaertredusert") ~ "Svært redusert",
    TRUE ~ NA_character_
  )
}

find_bt_field <- function(x) {
  nm <- names(x)
  if ("X7JB.BT" %in% nm) return("X7JB.BT")
  key <- tolower(gsub("[^a-zA-Z0-9]+", "", nm))
  preferred_keys <- c("x7jbbt", "7jbbt", "jbbt", "bt", "grazingpressure")
  idx <- which(key %in% preferred_keys)
  if (length(idx) == 0) return(NA_character_)
  nm[idx[1]]
}

parse_bt <- function(x) {
  z <- trimws(as.character(x))
  out <- suppressWarnings(as.integer(z))
  missing <- is.na(out)
  if (any(missing)) {
    extracted <- stringr::str_extract(z[missing], "(?<![0-9])[1-6](?![0-9])")
    out[missing] <- suppressWarnings(as.integer(extracted))
  }
  out
}

read_annual_file <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  year_from_file <- suppressWarnings(as.integer(stringr::str_extract(basename(path), "20[0-9]{2}")))
  if (!"year" %in% names(x)) {
    x$year <- year_from_file
  } else {
    x$year <- suppressWarnings(as.integer(x$year))
    x$year[!is.finite(x$year)] <- year_from_file
  }
  x$source_file <- basename(path)
  x
}

safe_kendall <- function(y, x, min_n = MIN_YEARS) {
  keep <- is.finite(y) & is.finite(x)
  if (sum(keep) < min_n) return(c(tau = NA_real_, p = NA_real_))
  tst <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "kendall", exact = FALSE))
  c(tau = unname(tst$estimate), p = unname(tst$p.value))
}

safe_theilsen_slope <- function(y, x, min_n = MIN_YEARS) {
  keep <- is.finite(y) & is.finite(x)
  y <- y[keep]; x <- x[keep]
  if (length(y) < min_n) return(NA_real_)
  ij <- utils::combn(seq_along(y), 2)
  dx <- x[ij[2, ]] - x[ij[1, ]]
  dy <- y[ij[2, ]] - y[ij[1, ]]
  valid <- is.finite(dx) & is.finite(dy) & dx != 0
  if (!any(valid)) return(NA_real_)
  stats::median(dy[valid] / dx[valid], na.rm = TRUE)
}

trajectory_step_qc <- function(dat) {
  dat %>%
    arrange(fid, year) %>%
    group_by(fid) %>%
    summarise(
      n_years_step = n(),
      max_abs_step = {
        d <- diff(AE_SJ_p50)
        if (length(d) == 0 || all(!is.finite(d))) NA_real_ else max(abs(d), na.rm = TRUE)
      },
      prop_nonnegative_steps = {
        d <- diff(AE_SJ_p50); d <- d[is.finite(d)]
        if (length(d) == 0) NA_real_ else mean(d >= 0)
      },
      prop_nonpositive_steps = {
        d <- diff(AE_SJ_p50); d <- d[is.finite(d)]
        if (length(d) == 0) NA_real_ else mean(d <= 0)
      },
      min_score = min(AE_SJ_p50, na.rm = TRUE),
      max_score = max(AE_SJ_p50, na.rm = TRUE),
      any_outside_reference = any(AE_SJ_p50 < REFERENCE_SCORE_MIN | AE_SJ_p50 > REFERENCE_SCORE_MAX, na.rm = TRUE),
      any_strongly_extreme = any(AE_SJ_p50 < EXTREME_SCORE_MIN | AE_SJ_p50 > EXTREME_SCORE_MAX, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(implausible_step = is.finite(max_abs_step) & max_abs_step > MAX_ABS_ANNUAL_STEP)
}

normalised_entropy4 <- function(p1, p2, p3, p4) {
  p <- cbind(p1, p2, p3, p4)
  out <- rep(NA_real_, nrow(p))
  ok <- rowSums(is.finite(p)) == 4 & rowSums(p, na.rm = TRUE) > 0
  if (!any(ok)) return(out)
  pp <- p[ok, , drop = FALSE]
  pp <- pp / rowSums(pp)
  terms <- ifelse(pp > 0, pp * log(pp), 0)
  out[ok] <- -rowSums(terms) / log(4)
  out
}

dominant_sj_zone <- function(p1, p2, p3, p4) {
  p <- cbind(p1, p2, p3, p4)
  out <- rep(NA_character_, nrow(p))
  ok <- rowSums(is.finite(p)) == 4
  if (!any(ok)) return(out)
  pp <- p[ok, , drop = FALSE]
  max_p <- apply(pp, 1, max)
  tied <- rowSums(abs(pp - max_p) < 1e-12) > 1
  idx <- max.col(pp, ties.method = "first")
  ok_index <- which(ok)
  out[ok_index[!tied]] <- paste0("SJ", idx[!tied], "-like")
  out
}

summarise_composition_period <- function(dat, years, suffix) {
  out <- dat %>%
    filter(year %in% years) %>%
    group_by(fid) %>%
    summarise(
      n_comp_years = n_distinct(year),
      n_pixels_min = min(n_valid, na.rm = TRUE),
      p1 = mean(prop_SJ1_like, na.rm = TRUE),
      p2 = mean(prop_SJ2_like, na.rm = TRUE),
      p3 = mean(prop_SJ3_like, na.rm = TRUE),
      p4 = mean(prop_SJ4_like, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      complete_period = n_comp_years == length(years),
      across(c(p1, p2, p3, p4), ~ if_else(complete_period, .x, NA_real_)),
      entropy = normalised_entropy4(p1, p2, p3, p4),
      dominant = dominant_sj_zone(p1, p2, p3, p4)
    )
  names(out)[names(out) != "fid"] <- paste0(names(out)[names(out) != "fid"], "_", suffix)
  out
}


# =============================================================================
# 6. CHECK INPUT PATHS
# =============================================================================

required_paths <- c(ANNUAL_INPUT_DIR, MODEL_REFERENCE_CSV, DEVELOPMENT_CSV, MI_SHAPEFILE)
missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0) {
  stop(
    "One or more required input paths do not exist:\n\n",
    paste0("  ", missing_paths, collapse = "\n"),
    "\n\nCheck the file-location settings at the top of the script."
  )
}


# =============================================================================
# 7. READ ANNUAL COVARIANCE-ADJUSTED AE-SJ DATA
# =============================================================================

input_files <- sort(list.files(ANNUAL_INPUT_DIR, pattern = ANNUAL_REGEX, full.names = TRUE))

if (length(input_files) == 0) {
  stop("No annual covariance-adjusted AE-SJ CSVs matched:\n", ANNUAL_REGEX, "\nin:\n", ANNUAL_INPUT_DIR)
}

cat("\nAnnual files found:\n")
cat(paste0("  ", basename(input_files), collapse = "\n"), "\n")

ae_mi <- purrr::map_dfr(input_files, read_annual_file)
annual_id_field <- find_id_field(ae_mi)
ae_mi$fid <- trimws(as.character(ae_mi[[annual_id_field]]))

# Left column = python/mi_aesj_annual_products.py's _pixel_stats()
# "aesj_" prefix, renamed to the field names used below.
for (pair in list(
  c("aesj_mean_AE_SJ", "AE_SJ_mean"), c("aesj_sd_AE_SJ", "AE_SJ_stdDev"),
  c("aesj_p10_AE_SJ", "AE_SJ_p10"), c("aesj_p25_AE_SJ", "AE_SJ_p25"),
  c("aesj_median_AE_SJ", "AE_SJ_p50"), c("aesj_p75_AE_SJ", "AE_SJ_p75"),
  c("aesj_p90_AE_SJ", "AE_SJ_p90"), c("aesj_n", "AE_SJ_count")
)) {
  ae_mi <- rename_if_needed(ae_mi, pair[1], pair[2])
}

required_fields <- c(
  "AE_SJ_mean", "AE_SJ_stdDev", "AE_SJ_p10", "AE_SJ_p25",
  "AE_SJ_p50", "AE_SJ_p75", "AE_SJ_p90", "AE_SJ_count"
)

missing_fields <- setdiff(required_fields, names(ae_mi))
if (length(missing_fields) > 0) {
  stop("Missing expected annual AE-SJ fields: ", paste(missing_fields, collapse = ", "))
}

ae_mi <- ae_mi %>%
  mutate(year = as.integer(year), across(all_of(required_fields), as.numeric))


# =============================================================================
# 8. READ FINAL EMPIRICAL SJ REFERENCE POSITIONS AND DERIVE TRIAGE THRESHOLDS
# =============================================================================

sj_reference_table <- readr::read_csv(MODEL_REFERENCE_CSV, show_col_types = FALSE) %>%
  transmute(sj = as.integer(sj), reference = as.numeric(median_AE_SJ)) %>%
  arrange(sj)

if (!identical(sj_reference_table$sj, 1:4)) {
  stop("Model reference CSV must contain SJ1-SJ4.")
}

ref_dat <- tibble::tibble(
  sj = 1:4,
  reference = sj_reference_table$reference,
  label = c(
    "SJ1  Intakt",
    "SJ2  Brakkleggingsfase",
    "SJ3  Tidlig gjenvekstsuksesjonsfase",
    "SJ4  Sein gjenvekstsuksesjonsfase"
  )
)

adjacent_gap <- function(ref_dat, pair) {
  abs(ref_dat$reference[pair[1]] - ref_dat$reference[pair[2]])
}

moderate_gap <- adjacent_gap(ref_dat, TRIAGE_MODERATE_PAIR)
large_gap    <- adjacent_gap(ref_dat, TRIAGE_LARGE_PAIR)

moderate_threshold <- TRIAGE_MODERATE_FRACTION * moderate_gap
large_threshold    <- TRIAGE_LARGE_FRACTION * large_gap

if (!is.finite(moderate_threshold) || !is.finite(large_threshold)) {
  stop("Could not compute triage thresholds from the SJ reference table.")
}

if (large_threshold <= moderate_threshold) {
  stop(
    "The 'large' triage threshold (", signif(large_threshold, 4),
    ") is not greater than the 'moderate' threshold (",
    signif(moderate_threshold, 4), "). Check TRIAGE_MODERATE_PAIR, ",
    "TRIAGE_LARGE_PAIR and their fractions."
  )
}

triage_threshold_summary <- tibble::tibble(
  tier = c("moderate", "large"),
  sj_pair = c(
    paste0("SJ", TRIAGE_MODERATE_PAIR[1], "-SJ", TRIAGE_MODERATE_PAIR[2]),
    paste0("SJ", TRIAGE_LARGE_PAIR[1], "-SJ", TRIAGE_LARGE_PAIR[2])
  ),
  reference_gap = c(moderate_gap, large_gap),
  fraction_of_gap = c(TRIAGE_MODERATE_FRACTION, TRIAGE_LARGE_FRACTION),
  threshold_AE_SJ_units = c(moderate_threshold, large_threshold),
  description = c(
    paste0(
      "A site is flagged 'Moderate' if its Theil-Sen trend-equivalent change over the ",
      "observed period is at least ", signif(TRIAGE_MODERATE_FRACTION, 2),
      " of the empirical AE-SJ gap between SJ", TRIAGE_MODERATE_PAIR[1],
      " and SJ", TRIAGE_MODERATE_PAIR[2], "."
    ),
    paste0(
      "A site is flagged 'Large' if its Theil-Sen trend-equivalent change over the ",
      "observed period is at least ", signif(TRIAGE_LARGE_FRACTION, 2),
      " of the empirical AE-SJ gap between SJ", TRIAGE_LARGE_PAIR[1],
      " and SJ", TRIAGE_LARGE_PAIR[2], "."
    )
  )
)

readr::write_csv(triage_threshold_summary, file.path(OUTPUT_DIR, "02_triage_threshold_derivation.csv"))

cat(
  "\nTriage thresholds derived from SJ reference spacing:\n",
  "  Moderate: ", signif(moderate_threshold, 4), " AE-SJ units (",
  signif(TRIAGE_MODERATE_FRACTION, 2), " x SJ", TRIAGE_MODERATE_PAIR[1], "-SJ",
  TRIAGE_MODERATE_PAIR[2], " gap of ", signif(moderate_gap, 4), ")\n",
  "  Large:    ", signif(large_threshold, 4), " AE-SJ units (",
  signif(TRIAGE_LARGE_FRACTION, 2), " x SJ", TRIAGE_LARGE_PAIR[1], "-SJ",
  TRIAGE_LARGE_PAIR[2], " gap of ", signif(large_gap, 4), ")\n\n",
  sep = ""
)


# =============================================================================
# 9. READ MI GEOMETRY, TILSTAND AND BT
# =============================================================================

mi_sf <- sf::st_read(MI_SHAPEFILE, quiet = TRUE)
mi_id_field <- find_id_field(mi_sf)
mi_sf$fid <- trimws(as.character(mi_sf[[mi_id_field]]))

if ("tilstand" %in% names(mi_sf)) {
  tilstand_lookup <- mi_sf %>%
    sf::st_drop_geometry() %>%
    transmute(fid, tilstand = normalise_tilstand(tilstand)) %>%
    distinct(fid, .keep_all = TRUE)
} else if ("tilstand" %in% names(ae_mi)) {
  tilstand_lookup <- ae_mi %>%
    transmute(fid, tilstand = normalise_tilstand(tilstand)) %>%
    distinct(fid, .keep_all = TRUE)
} else {
  stop("No tilstand field found in the MI shapefile or annual CSVs.")
}

bt_field_mi <- find_bt_field(mi_sf)
bt_field_annual <- find_bt_field(ae_mi)

if (!is.na(bt_field_mi)) {
  bt_lookup <- mi_sf %>%
    sf::st_drop_geometry() %>%
    transmute(fid, bt = parse_bt(.data[[bt_field_mi]])) %>%
    distinct(fid, .keep_all = TRUE)
} else if (!is.na(bt_field_annual)) {
  bt_lookup <- ae_mi %>%
    transmute(fid, bt = parse_bt(.data[[bt_field_annual]])) %>%
    distinct(fid, .keep_all = TRUE)
} else {
  stop(
    "Could not find a BT field in the MI shapefile or annual CSVs. ",
    "BT5/6 must be excluded, so the analysis stops rather than guessing."
  )
}

ae_mi <- ae_mi %>%
  select(-any_of(c("tilstand", "bt"))) %>%
  left_join(tilstand_lookup, by = "fid") %>%
  left_join(bt_lookup, by = "fid") %>%
  mutate(tilstand = factor(tilstand, levels = TILSTAND_LEVELS))


# =============================================================================
# 10. READ DEVELOPMENT RESULTS AND DEFINE EXCLUDED FIDS
# =============================================================================

development <- readr::read_csv(DEVELOPMENT_CSV, show_col_types = FALSE)
development_id_field <- find_id_field(development)

development <- development %>%
  mutate(
    fid = trimws(as.character(.data[[development_id_field]])),
    dev_ha_total = dplyr::coalesce(suppressWarnings(as.numeric(dev_ha_total)), 0)
  )

development_excluded <- development %>%
  filter(dev_ha_total > DEVELOPMENT_EXCLUDE_MIN_HA) %>%
  distinct(fid) %>%
  pull(fid)


# =============================================================================
# 11. APPLY OUT-OF-SCOPE FILTERS
# =============================================================================

all_fids_before <- n_distinct(ae_mi$fid)

bt_excluded <- ae_mi %>%
  filter(bt %in% EXCLUDE_BT) %>%
  distinct(fid) %>%
  pull(fid)

excluded_fids <- union(bt_excluded, development_excluded)

ae_analysis <- ae_mi %>% filter(!fid %in% excluded_fids)

filter_summary <- tibble::tibble(
  filter = c(
    "All annual AE-SJ FIDs", "Excluded BT5/BT6", "Excluded persistent development",
    "Unique FIDs excluded by either rule", "FIDs retained"
  ),
  n_fids = c(
    all_fids_before, length(unique(bt_excluded)), length(unique(development_excluded)),
    length(unique(excluded_fids)), n_distinct(ae_analysis$fid)
  )
)

cat("\nOut-of-scope filtering:\n")
print(filter_summary)

readr::write_csv(filter_summary, file.path(OUTPUT_DIR, "00_filter_summary.csv"))


# =============================================================================
# 12. ANALYSIS DATASET
# =============================================================================

state_dat <- ae_analysis %>%
  filter(
    !is.na(fid), fid != "",
    is.finite(AE_SJ_p50), is.finite(AE_SJ_count),
    AE_SJ_count >= MIN_PIXELS_STATE
  )


# =============================================================================
# 13. TRAJECTORY QC: ANNUAL JUMPS AND OUT-OF-REFERENCE SCORES
# =============================================================================

trajectory_qc <- trajectory_step_qc(state_dat)

readr::write_csv(trajectory_qc, file.path(OUTPUT_DIR, "01_trajectory_QC.csv"))

implausible_step_fids <- trajectory_qc %>% filter(implausible_step) %>% pull(fid)

state_dat <- state_dat %>% filter(!fid %in% implausible_step_fids)

cat(
  "\nTrajectory QC:\n",
  "  FIDs excluded for max |annual step| > ", MAX_ABS_ANNUAL_STEP, ": ",
  length(unique(implausible_step_fids)), "\n",
  "  FIDs with at least one AE-SJ score outside 0-1: ",
  sum(trajectory_qc$any_outside_reference, na.rm = TRUE), "\n",
  "  FIDs with at least one AE-SJ score outside ", EXTREME_SCORE_MIN, " to ",
  EXTREME_SCORE_MAX, ": ", sum(trajectory_qc$any_strongly_extreme, na.rm = TRUE), "\n",
  sep = ""
)


# =============================================================================
# 14. FIGURE: ANNUAL TRAJECTORIES BY TILSTAND (context figure)
# =============================================================================

tilstand_trajectory <- state_dat %>%
  filter(!is.na(tilstand)) %>%
  group_by(tilstand, year) %>%
  summarise(
    n_fids = n_distinct(fid),
    p25_P50 = quantile(AE_SJ_p50, 0.25, na.rm = TRUE),
    median_P50 = median(AE_SJ_p50, na.rm = TRUE),
    p75_P50 = quantile(AE_SJ_p50, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

first_plot_year <- min(tilstand_trajectory$year, na.rm = TRUE)
last_plot_year <- max(tilstand_trajectory$year, na.rm = TRUE)

tilstand_end_labels <- tilstand_trajectory %>%
  filter(year == last_plot_year)

p_tilstand <- ggplot(
  tilstand_trajectory,
  aes(x = year, y = median_P50, colour = tilstand, fill = tilstand, group = tilstand)
) +
  geom_ribbon(aes(ymin = p25_P50, ymax = p75_P50), alpha = 0.16, colour = NA) +
  sj_reference_layers(ref_dat, first_plot_year) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.1) +
  geom_text(
    data = tilstand_end_labels,
    aes(label = tilstand),
    hjust = -0.12, size = gg_text_size("row_label"), fontface = "bold", show.legend = FALSE
  ) +
  scale_colour_manual(values = TILSTAND_COLOURS, guide = "none") +
  scale_fill_manual(values = TILSTAND_COLOURS, guide = "none") +
  scale_x_continuous(
    breaks = sort(unique(tilstand_trajectory$year)),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Annual succession trajectories by condition class",
    subtitle = "Median AE-SJ score each year, by mapped condition class",
    x = "Year",
    y = "Site median AE-SJ score",
    caption = "Each site receives equal weight; ribbons show the between-site range (25th-75th percentile). Lower AE-SJ means more successional structure."
  ) +
  theme_report() +
  theme(plot.margin = margin(18, 62, 14, 18))


# =============================================================================
# 15. POLYGON-LEVEL ROBUST TEMPORAL CHANGE
# =============================================================================
#
# The temporal-change statistic is the Theil-Sen slope of annual site median
# AE-SJ. Kendall's tau and BH-adjusted p-values are retained as supplementary
# diagnostics only - see note at the top of the script on statistical power
# at n = 8.
# =============================================================================

state_trends <- state_dat %>%
  group_by(fid) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    n_years = n(),
    first_year = min(year),
    last_year = max(year),
    sj_base = median(head(AE_SJ_p50, 2), na.rm = TRUE),
    sj_final = median(tail(AE_SJ_p50, 2), na.rm = TRUE),
    sj_mean = mean(AE_SJ_p50, na.rm = TRUE),
    sj_tau = safe_kendall(AE_SJ_p50, year)[["tau"]],
    sj_p = safe_kendall(AE_SJ_p50, year)[["p"]],
    sj_slope = safe_theilsen_slope(AE_SJ_p50, year),
    tilstand = dplyr::first(as.character(tilstand[!is.na(tilstand)]), default = NA_character_),
    bt = dplyr::first(bt[!is.na(bt)], default = NA_integer_),
    .groups = "drop"
  ) %>%
  left_join(
    trajectory_qc %>%
      select(fid, max_abs_step, prop_nonnegative_steps, prop_nonpositive_steps,
             min_score, max_score, any_outside_reference, any_strongly_extreme, implausible_step),
    by = "fid"
  ) %>%
  mutate(
    sj_padj = p.adjust(sj_p, method = "BH"),
    year_span = last_year - first_year,
    # sj_change is the Theil-Sen trend magnitude, the quantity that drives
    # triage_class. sj_base/sj_final are the raw endpoint values.
    sj_change = sj_slope * year_span
  )


# =============================================================================
# 16. TRIAGE CLASSIFICATION
# =============================================================================

state_trends <- state_trends %>%
  mutate(
    direction_eligible = n_years >= MIN_YEARS & is.finite(sj_slope) & is.finite(sj_change),
    
    triage_direction = case_when(
      !direction_eligible ~ NA_character_,
      sj_change > 0 ~ DIRECTION_RESTORATIVE,
      sj_change < 0 ~ DIRECTION_SUCCESSIONAL,
      TRUE ~ DIRECTION_NO_CHANGE
    ),
    
    triage_tier = case_when(
      !direction_eligible ~ NA_character_,
      abs(sj_change) >= large_threshold ~ "Strong",
      abs(sj_change) >= moderate_threshold ~ "Moderate",
      TRUE ~ "No strong signal"
    ),
    
    triage_class = case_when(
      !direction_eligible ~ NA_character_,
      triage_tier == "No strong signal" ~ "No strong signal",
      triage_tier == "Strong" & triage_direction == DIRECTION_RESTORATIVE ~ "Strong restoration",
      triage_tier == "Strong" & triage_direction == DIRECTION_SUCCESSIONAL ~ "Strong succession",
      triage_tier == "Moderate" & triage_direction == DIRECTION_RESTORATIVE ~ "Moderate restoration",
      triage_tier == "Moderate" & triage_direction == DIRECTION_SUCCESSIONAL ~ "Moderate succession",
      TRUE ~ NA_character_
    ),
    
    triage_class = factor(triage_class, levels = TRIAGE_LEVELS),
    tilstand = factor(tilstand, levels = TILSTAND_LEVELS)
  )



# =============================================================================
# 17. OVERVIEW FIGURES
# =============================================================================

# ---- 17a. How the thresholds split the dataset ------------------------------
#
# Four dashed lines in total: an inner pair (the moderate threshold, one on
# each side of zero) and an outer pair (the strong threshold, one on each
# side of zero). Each line is labelled directly on the plot so the meaning
# of all four does not depend on the caption alone.

threshold_label_dat <- tibble::tibble(
  x = c(-large_threshold, -moderate_threshold, moderate_threshold, large_threshold),
  label = c(
    "Strong\nsuccession", "Moderate\nsuccession",
    "Moderate\nrestoration", "Strong\nrestoration"
  )
)

p_threshold_distribution <- ggplot(
  state_trends %>% filter(direction_eligible),
  aes(x = sj_change, fill = triage_class)
) +
  geom_histogram(bins = 70, colour = SKOG_THEME$surface, linewidth = 0.15) +
  geom_vline(
    xintercept = threshold_label_dat$x,
    linetype = "dashed", linewidth = 0.5, colour = SKOG_THEME$rule2
  ) +
  geom_text(
    data = threshold_label_dat,
    aes(x = x, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.2, size = gg_text_size("eyebrow"), colour = SKOG_THEME$ink2, lineheight = 0.85
  ) +
  scale_fill_manual(values = TRIAGE_COLOURS, drop = FALSE, name = NULL) +
  labs(
    title = "How the triage thresholds split the dataset",
    subtitle = "Estimated AE-SJ change per site, compared with thresholds set from SJ reference spacing",
    x = "Estimated AE-SJ change over the monitoring period",
    y = "Number of sites",
    caption = wrap_caption(paste0(
      "Each site is compared with a moderate threshold (\u00b1", signif(moderate_threshold, 3),
      ") and a strong threshold (\u00b1", signif(large_threshold, 3), ") in AE-SJ units, in either direction."
    ))
  ) +
  theme_report() +
  theme(legend.position = "bottom")

# ---- 17b. Tilstand x triage matrix ------------------------------------------
#
# Tile colour marks the triage category (as elsewhere). Tile darkness scales
# with the number of sites in that cell, RELATIVE TO ITS OWN COLUMN - not
# the whole grid - so a column with small absolute counts (e.g. Strong
# restoration) still shows a readable contrast between its own tilstand
# rows, instead of being washed out by comparison to a much larger column
# like "No strong signal".

matrix_dat <- state_trends %>%
  filter(direction_eligible, !is.na(tilstand), !is.na(triage_class)) %>%
  count(tilstand, triage_class, name = "n_fids") %>%
  tidyr::complete(
    tilstand = factor(TILSTAND_LEVELS, levels = TILSTAND_LEVELS),
    triage_class = factor(TRIAGE_LEVELS, levels = TRIAGE_LEVELS),
    fill = list(n_fids = 0L)
  ) %>%
  group_by(tilstand) %>%
  mutate(
    row_total = sum(n_fids),
    row_pct = if_else(row_total > 0, 100 * n_fids / row_total, NA_real_),
    tile_label = if_else(
      row_total > 0,
      sprintf("n = %s\n(%.1f%%)", skog_num(n_fids), row_pct),
      paste0("n = ", skog_num(n_fids))
    )
  ) %>%
  ungroup() %>%
  group_by(triage_class) %>%
  mutate(
    col_max = max(n_fids),
    tile_alpha = if_else(col_max > 0, 0.15 + 0.85 * (n_fids / col_max), 0.15)
  ) %>%
  ungroup()


p_triage_matrix <- ggplot(matrix_dat, aes(x = triage_class, y = tilstand, fill = triage_class, alpha = tile_alpha)) +
  geom_tile(colour = SKOG_THEME$surface, linewidth = 1.0) +
  geom_text(aes(label = tile_label), colour = SKOG_THEME$ink, size = gg_text_size("value"), lineheight = 0.95, alpha = 1) +
  scale_fill_manual(values = TRIAGE_COLOURS, guide = "none") +
  scale_alpha_identity() +
  scale_y_discrete(limits = rev(TILSTAND_LEVELS)) +
  labs(
    title = "Where are successional and restorative sites concentrated?",
    subtitle = "Number of sites in each triage category, by mapped condition class (tilstand)",
    x = NULL,
    y = "Mapped condition (tilstand)",
    caption = "Percentages are of the row total (within each tilstand class). Within each column, darker tiles hold more sites."
  ) +
  theme_report() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


# =============================================================================
# 18. READ ANNUAL SJ-LIKE PIXEL COMPOSITION PRODUCT
# =============================================================================

if (!dir.exists(COMPOSITION_SEARCH_ROOT)) {
  stop("Composition search root does not exist:\n", COMPOSITION_SEARCH_ROOT)
}

composition_files <- sort(list.files(COMPOSITION_SEARCH_ROOT, pattern = COMPOSITION_REGEX, recursive = TRUE, full.names = TRUE))

if (length(composition_files) == 0) {
  stop("No annual composition CSVs matched:\n", COMPOSITION_REGEX, "\nunder:\n", COMPOSITION_SEARCH_ROOT)
}

composition_file_years <- suppressWarnings(as.integer(stringr::str_extract(basename(composition_files), "20[0-9]{2}")))

if (anyDuplicated(composition_file_years)) {
  duplicate_years <- unique(composition_file_years[duplicated(composition_file_years)])
  stop(
    "Multiple composition files were found for year(s): ", paste(duplicate_years, collapse = ", "),
    "\nNarrow COMPOSITION_SEARCH_ROOT so there is one export per year.\n\n",
    paste(composition_files, collapse = "\n")
  )
}

missing_composition_years <- setdiff(2018:2025, composition_file_years)
if (length(missing_composition_years) > 0) {
  stop("Composition exports are missing year(s): ", paste(missing_composition_years, collapse = ", "))
}

cat("\nComposition files found:\n")
cat(paste0("  ", composition_files, collapse = "\n"), "\n")

composition_raw <- purrr::map_dfr(composition_files, read_annual_file)
composition_id_field <- find_id_field(composition_raw)
composition_raw$fid <- trimws(as.character(composition_raw[[composition_id_field]]))

for (pair in list(
  c("n_valid_sum", "n_valid"), c("n_SJ1_like_sum", "n_SJ1_like"),
  c("n_SJ2_like_sum", "n_SJ2_like"), c("n_SJ3_like_sum", "n_SJ3_like"), c("n_SJ4_like_sum", "n_SJ4_like")
)) {
  composition_raw <- rename_if_needed(composition_raw, pair[1], pair[2])
}

COMPOSITION_PROP_FIELDS <- c("prop_SJ1_like", "prop_SJ2_like", "prop_SJ3_like", "prop_SJ4_like")
required_composition_fields <- c("fid", "year", "n_valid", COMPOSITION_PROP_FIELDS)
missing_composition_fields <- setdiff(required_composition_fields, names(composition_raw))

if (length(missing_composition_fields) > 0) {
  stop("Composition CSVs are missing required fields: ", paste(missing_composition_fields, collapse = ", "))
}

composition_raw <- composition_raw %>%
  mutate(
    year = as.integer(year),
    n_valid = suppressWarnings(as.numeric(n_valid)),
    across(all_of(COMPOSITION_PROP_FIELDS), ~ suppressWarnings(as.numeric(.x)))
  )

composition_duplicates <- composition_raw %>% count(fid, year, name = "n_rows") %>% filter(n_rows > 1)
if (nrow(composition_duplicates) > 0) {
  stop(
    "Composition input contains duplicate fid-year rows. First duplicates:\n",
    paste(utils::capture.output(print(utils::head(composition_duplicates, 20))), collapse = "\n")
  )
}

composition_raw <- composition_raw %>%
  mutate(
    prop_sum_calc = prop_SJ1_like + prop_SJ2_like + prop_SJ3_like + prop_SJ4_like,
    composition_qc_ok =
      is.finite(n_valid) & n_valid >= MIN_COMPOSITION_PIXELS &
      if_all(all_of(COMPOSITION_PROP_FIELDS), ~ is.finite(.x) & .x >= -0.001 & .x <= 1.001) &
      is.finite(prop_sum_calc) & abs(prop_sum_calc - 1) <= 0.02,
    entropy = normalised_entropy4(prop_SJ1_like, prop_SJ2_like, prop_SJ3_like, prop_SJ4_like),
    dominant_zone = dominant_sj_zone(prop_SJ1_like, prop_SJ2_like, prop_SJ3_like, prop_SJ4_like)
  )

composition_qc_summary <- tibble::tibble(
  metric = c("Composition fid-year rows", "Rows with >= minimum valid pixels", "Rows passing composition QC", "Rows failing composition QC"),
  value = c(
    nrow(composition_raw),
    sum(is.finite(composition_raw$n_valid) & composition_raw$n_valid >= MIN_COMPOSITION_PIXELS, na.rm = TRUE),
    sum(composition_raw$composition_qc_ok, na.rm = TRUE),
    sum(!composition_raw$composition_qc_ok, na.rm = TRUE)
  )
)

readr::write_csv(composition_qc_summary, file.path(OUTPUT_DIR, "05_composition_QC_summary.csv"))

composition_annual <- composition_raw %>% filter(composition_qc_ok, fid %in% state_trends$fid)


# =============================================================================
# 19. BASELINE AND FINAL COMPOSITION FOR EACH POLYGON
# =============================================================================

composition_base <- summarise_composition_period(composition_annual, c(2018, 2019), "base")
composition_final <- summarise_composition_period(composition_annual, c(2024, 2025), "final")

composition_periods <- composition_base %>%
  full_join(composition_final, by = "fid") %>%
  mutate(
    dp1 = p1_final - p1_base, dp2 = p2_final - p2_base,
    dp3 = p3_final - p3_base, dp4 = p4_final - p4_base,
    entropy_change = entropy_final - entropy_base
  )


# =============================================================================
# 20. READ PIXEL-TRANSITION PRODUCT
# =============================================================================

if (!dir.exists(TRANSITION_SEARCH_ROOT)) {
  stop("Transition search root does not exist:\n", TRANSITION_SEARCH_ROOT)
}

transition_files <- sort(list.files(TRANSITION_SEARCH_ROOT, pattern = TRANSITION_REGEX, recursive = TRUE, full.names = TRUE))

if (length(transition_files) != 1) {
  stop(
    "Expected exactly one pixel-transition CSV but found ", length(transition_files),
    ".\nPattern:\n", TRANSITION_REGEX, "\nSearch root:\n", TRANSITION_SEARCH_ROOT,
    if (length(transition_files) > 0) paste0("\nMatches:\n", paste(transition_files, collapse = "\n")) else ""
  )
}

TRANSITION_CSV <- transition_files[1]
cat("\nPixel-transition file:\n  ", TRANSITION_CSV, "\n", sep = "")

transition_raw <- readr::read_csv(TRANSITION_CSV, show_col_types = FALSE, progress = FALSE)
transition_id_field <- find_id_field(transition_raw)
transition_raw$fid <- trimws(as.character(transition_raw[[transition_id_field]]))

for (pair in list(
  c("n_valid_sum", "n_valid"), c("n_toward_SJ1_sum", "n_toward_SJ1"),
  c("n_same_zone_sum", "n_same_zone"), c("n_toward_SJ4_sum", "n_toward_SJ4")
)) {
  transition_raw <- rename_if_needed(transition_raw, pair[1], pair[2])
}

TRANSITION_MATRIX_FIELDS <- c(
  "SJ1_to_SJ1", "SJ1_to_SJ2", "SJ1_to_SJ3", "SJ1_to_SJ4",
  "SJ2_to_SJ1", "SJ2_to_SJ2", "SJ2_to_SJ3", "SJ2_to_SJ4",
  "SJ3_to_SJ1", "SJ3_to_SJ2", "SJ3_to_SJ3", "SJ3_to_SJ4",
  "SJ4_to_SJ1", "SJ4_to_SJ2", "SJ4_to_SJ3", "SJ4_to_SJ4"
)

required_transition_fields <- c("fid", "n_valid", "prop_toward_SJ1", "prop_same_zone", "prop_toward_SJ4", TRANSITION_MATRIX_FIELDS)
missing_transition_fields <- setdiff(required_transition_fields, names(transition_raw))

if (length(missing_transition_fields) > 0) {
  stop("Pixel-transition CSV is missing required fields: ", paste(missing_transition_fields, collapse = ", "))
}

transition_numeric_fields <- intersect(
  c("n_valid", "n_toward_SJ1", "n_same_zone", "n_toward_SJ4", "prop_toward_SJ1", "prop_same_zone", "prop_toward_SJ4", TRANSITION_MATRIX_FIELDS),
  names(transition_raw)
)

transition_raw <- transition_raw %>%
  mutate(across(all_of(transition_numeric_fields), ~ suppressWarnings(as.numeric(.x))))

transition_duplicates <- transition_raw %>% count(fid, name = "n_rows") %>% filter(n_rows > 1)
if (nrow(transition_duplicates) > 0) {
  stop(
    "Pixel-transition input contains duplicate FIDs. First duplicates:\n",
    paste(utils::capture.output(print(utils::head(transition_duplicates, 20))), collapse = "\n")
  )
}


# =============================================================================
# 21. PIXEL-TRANSITION EXTENT AND COHERENCE METRICS
# =============================================================================

transition_analysis <- transition_raw %>%
  filter(fid %in% state_trends$fid) %>%
  mutate(
    # Completeness check: the 16 SJa_to_SJb transition counts must sum to
    # n_valid (every eligible pixel falls in exactly one pair).
    transition_prop_sum = rowSums(across(all_of(TRANSITION_MATRIX_FIELDS))) / n_valid,
    transition_qc_ok =
      is.finite(n_valid) & n_valid >= MIN_TRANSITION_PIXELS &
      is.finite(transition_prop_sum) & abs(transition_prop_sum - 1) <= 0.02,
    mean_signed_zone_shift = (
      1 * SJ1_to_SJ2 + 2 * SJ1_to_SJ3 + 3 * SJ1_to_SJ4 -
        1 * SJ2_to_SJ1 + 1 * SJ2_to_SJ3 + 2 * SJ2_to_SJ4 -
        2 * SJ3_to_SJ1 - 1 * SJ3_to_SJ2 + 1 * SJ3_to_SJ4 -
        3 * SJ4_to_SJ1 - 2 * SJ4_to_SJ2 - 1 * SJ4_to_SJ3
    ) / n_valid,
    mean_abs_zone_shift = (
      1 * SJ1_to_SJ2 + 2 * SJ1_to_SJ3 + 3 * SJ1_to_SJ4 +
        1 * SJ2_to_SJ1 + 1 * SJ2_to_SJ3 + 2 * SJ2_to_SJ4 +
        2 * SJ3_to_SJ1 + 1 * SJ3_to_SJ2 + 1 * SJ3_to_SJ4 +
        3 * SJ4_to_SJ1 + 2 * SJ4_to_SJ2 + 1 * SJ4_to_SJ3
    ) / n_valid,
    long_jump_fraction = (
      SJ1_to_SJ3 + SJ1_to_SJ4 + SJ2_to_SJ4 + SJ3_to_SJ1 + SJ4_to_SJ1 + SJ4_to_SJ2
    ) / n_valid
  )

transition_qc_summary <- tibble::tibble(
  metric = c("Transition FIDs", "FIDs with >= minimum matched pixels", "FIDs passing transition QC", "FIDs failing transition QC"),
  value = c(
    nrow(transition_analysis),
    sum(is.finite(transition_analysis$n_valid) & transition_analysis$n_valid >= MIN_TRANSITION_PIXELS, na.rm = TRUE),
    sum(transition_analysis$transition_qc_ok, na.rm = TRUE),
    sum(!transition_analysis$transition_qc_ok, na.rm = TRUE)
  )
)

readr::write_csv(transition_qc_summary, file.path(OUTPUT_DIR, "06_transition_QC_summary.csv"))


# =============================================================================
# 22. MASTER SITE-LEVEL RESULTS TABLE
# =============================================================================


# Scrub/tree clearance counts (SJ3/4 -> SJ1/2) from
# python/mi_aesj_annual_products.py, split into burn-coincident and clearance
# with its pixel-level burn mask, and the same split for spruce pixels.
# any_of() keeps this working with transition files that lack these columns.
PICEA_BURN_SPLIT_FIELDS <- c(
  "n_SJ3_to_SJ12", "n_SJ4_to_SJ12", "n_SJ34_to_SJ12",
  "n_SJ34_to_SJ12_burnt", "n_SJ34_to_SJ12_clearance",
  "area_SJ34_to_SJ12_ha", "area_SJ34_to_SJ12_burnt_ha", "area_SJ34_to_SJ12_clearance_ha",
  "n_picea_px", "n_SJ3_to_SJ12_picea", "n_SJ4_to_SJ12_picea", "n_SJ34_to_SJ12_picea",
  "n_SJ34_to_SJ12_picea_burnt", "n_SJ34_to_SJ12_picea_clearance",
  "area_SJ34_to_SJ12_picea_ha", "area_SJ34_to_SJ12_picea_burnt_ha", "area_SJ34_to_SJ12_picea_clearance_ha",
  "prop_SJ34_to_SJ12_picea",
  # scrub removal = SJ3-source, non-spruce; tree removal = SJ4-source or spruce
  "n_scrub_to_SJ12", "n_tree_to_SJ12", "n_scrub_to_SJ12_burnt", "n_tree_to_SJ12_burnt",
  "area_scrub_to_SJ12_ha", "area_tree_to_SJ12_ha", "area_scrub_to_SJ12_burnt_ha", "area_tree_to_SJ12_burnt_ha"
)

site_results <- state_trends %>%
  left_join(composition_periods, by = "fid") %>%
  left_join(
    transition_analysis %>%
      select(fid, n_valid, transition_qc_ok, prop_toward_SJ1, prop_same_zone, prop_toward_SJ4,
             mean_signed_zone_shift, mean_abs_zone_shift, long_jump_fraction, all_of(TRANSITION_MATRIX_FIELDS),
             any_of(PICEA_BURN_SPLIT_FIELDS)) %>%
      rename(n_transition_pixels = n_valid),
    by = "fid"
  ) %>%
  mutate(
    prop_changed_any = prop_toward_SJ1 + prop_toward_SJ4,
    prop_concordant = case_when(
      triage_direction == DIRECTION_RESTORATIVE ~ prop_toward_SJ1,
      triage_direction == DIRECTION_SUCCESSIONAL ~ prop_toward_SJ4,
      TRUE ~ NA_real_
    ),
    prop_opposite = case_when(
      triage_direction == DIRECTION_RESTORATIVE ~ prop_toward_SJ4,
      triage_direction == DIRECTION_SUCCESSIONAL ~ prop_toward_SJ1,
      TRUE ~ NA_real_
    ),
    directional_coherence_changed = case_when(
      triage_tier %in% c("Moderate", "Strong") & is.finite(prop_changed_any) & prop_changed_any > 0 ~
        prop_concordant / prop_changed_any,
      TRUE ~ NA_real_
    ),
    change_pattern = case_when(
      !direction_eligible ~ NA_character_,
      triage_tier == "No strong signal" ~ "Below triage threshold",
      !(is.finite(n_transition_pixels) & dplyr::coalesce(transition_qc_ok, FALSE)) ~ NA_character_,
      is.finite(prop_changed_any) & prop_changed_any <= 1e-12 ~ "Within-class shift",
      is.finite(directional_coherence_changed) & directional_coherence_changed < DIRECTIONAL_COHERENCE_THRESHOLD ~ "Mixed-direction change",
      is.finite(prop_changed_any) & prop_changed_any >= WIDESPREAD_CHANGED_THRESHOLD ~ "Widespread change",
      is.finite(prop_changed_any) ~ "Partial change",
      TRUE ~ NA_character_
    ),
    # Separate from triage_direction (which is based on the continuous
    # AE-SJ trend): this compares only the site's dominant SJ Class in the
    # baseline vs final period, i.e. whether that single label moved at
    # all - see section 27 for the figure this also feeds.
    dominant_direction = case_when(
      is.na(dominant_base) | is.na(dominant_final) ~ NA_character_,
      match(dominant_final, SJ_ZONE_LEVELS) == match(dominant_base, SJ_ZONE_LEVELS) ~ "No change",
      match(dominant_final, SJ_ZONE_LEVELS) > match(dominant_base, SJ_ZONE_LEVELS) ~ "Succession",
      TRUE ~ "Restoration"
    )
  )

readr::write_csv(site_results, file.path(OUTPUT_DIR, "07_site_level_results_ALL.csv"))


# =============================================================================
# 23. SPLIT INTO SUCCESSIONAL / RESTORATIVE / NO-STRONG-SIGNAL POPULATIONS
# =============================================================================

decline_sites <- site_results %>%
  filter(triage_class %in% c("Moderate succession", "Strong succession")) %>%
  arrange(sj_change)

recovery_sites <- site_results %>%
  filter(triage_class %in% c("Moderate restoration", "Strong restoration")) %>%
  arrange(desc(sj_change))

no_signal_sites <- site_results %>%
  filter(triage_class == "No strong signal") %>%
  arrange(abs(sj_change))

readr::write_csv(decline_sites, file.path(OUTPUT_DIR, "08_successional_sites.csv"))
readr::write_csv(recovery_sites, file.path(OUTPUT_DIR, "09_restorative_sites.csv"))
readr::write_csv(no_signal_sites, file.path(OUTPUT_DIR, "10_no_strong_signal_sites.csv"))




# =============================================================================
# 25. SUCCESSIONAL (DECLINE) MECHANISM BLOCK
# =============================================================================
#
# Sites moving toward SJ4-like structure. The key question this block is
# built to answer: is the succession spread evenly across the site
# (consistent with a site-wide pressure and little visible localised
# management), or concentrated in part of the site (consistent with a
# localised disturbance or an untreated pocket within an otherwise stable
# site)?
# =============================================================================

decline_composition_summary <- build_period_composition_summary(decline_sites, group_vars = "tilstand")

p_decline_composition <- make_composition_change_plot(
  decline_composition_summary,
  title = "Successional sites: how the structure shifts",
  subtitle = "Mean SJ-like composition, 2018-2019 vs 2024-2025, by condition class"
)

decline_mechanism_dat <- decline_sites %>%
  filter(!is.na(tilstand), transition_qc_ok, is.finite(prop_changed_any), is.finite(directional_coherence_changed))

p_decline_mechanism <- make_mechanism_scatter(
  decline_mechanism_dat,
  tier_colours = DECLINE_TIER_COLOURS,
  title = "How is succession happening within each site?",
  subtitle = "Sites moving toward SJ4-like structure: how much of the site changed, and how consistently"
)

p_decline_pattern <- make_pattern_summary_plot(
  decline_sites,
  palette = DECLINE_PATTERN_COLOURS,
  title = "Succession pattern by condition class",
  subtitle = "Widespread change points to a site-wide pressure; partial change points to a localised one"
)

decline_pattern_counts <- decline_sites %>%
  filter(!is.na(change_pattern), change_pattern != "Below triage threshold") %>%
  count(tilstand, change_pattern, name = "n_fids")


p_decline_trajectories <- make_trajectory_plot(
  state_dat, decline_sites, ref_dat, line_colour = SKOG_THEME$p3,
  title = "Strongest successional sites",
  subtitle = paste0("Annual AE-SJ score, top ", SITE_GRAPH_TOP_N_PER_DIRECTION, " sites by estimated change; full list in CSV/shapefile"),
  top_n = SITE_GRAPH_TOP_N_PER_DIRECTION
)

p_decline_entropy <- make_diversity_trajectory_plot(
  composition_annual, decline_sites, colour_palette = TILSTAND_COLOURS,
  title = "Structural diversity of successional sites",
  subtitle = "Median mix of SJ Classes within each site, by year"
)


# =============================================================================
# 26. RESTORATIVE (RECOVERY) MECHANISM BLOCK
# =============================================================================
#
# Sites moving toward SJ1-like structure. The key question this block is
# built to answer: is the restoration spread evenly across the site
# (consistent with a site-wide change, e.g. reduced pressure across the
# whole site), or concentrated in part of the site (consistent with
# targeted restoration work in that part of the site)?
# =============================================================================

recovery_composition_summary <- build_period_composition_summary(recovery_sites, group_vars = "tilstand")

p_recovery_composition <- make_composition_change_plot(
  recovery_composition_summary,
  title = "Restorative sites: how the structure shifts",
  subtitle = "Mean SJ-like composition, 2018-2019 vs 2024-2025, by condition class"
)

recovery_mechanism_dat <- recovery_sites %>%
  filter(!is.na(tilstand), transition_qc_ok, is.finite(prop_changed_any), is.finite(directional_coherence_changed))

p_recovery_mechanism <- make_mechanism_scatter(
  recovery_mechanism_dat,
  tier_colours = RECOVERY_TIER_COLOURS,
  title = "How is restoration happening within each site?",
  subtitle = "Sites moving toward SJ1-like structure: how much of the site changed, and how consistently"
)

p_recovery_pattern <- make_pattern_summary_plot(
  recovery_sites,
  palette = RECOVERY_PATTERN_COLOURS,
  title = "Restoration pattern by condition class",
  subtitle = "Widespread change points to a site-wide driver; partial change points to targeted intervention"
)

recovery_pattern_counts <- recovery_sites %>%
  filter(!is.na(change_pattern), change_pattern != "Below triage threshold") %>%
  count(tilstand, change_pattern, name = "n_fids")


p_recovery_trajectories <- make_trajectory_plot(
  state_dat, recovery_sites, ref_dat, line_colour = SKOG_THEME$g3,
  title = "Strongest restorative sites",
  subtitle = paste0("Annual AE-SJ score, top ", SITE_GRAPH_TOP_N_PER_DIRECTION, " sites by estimated change; full list in CSV/shapefile"),
  top_n = SITE_GRAPH_TOP_N_PER_DIRECTION
)

p_recovery_entropy <- make_diversity_trajectory_plot(
  composition_annual, recovery_sites, colour_palette = TILSTAND_COLOURS,
  title = "Structural diversity of restorative sites",
  subtitle = "Median mix of SJ Classes within each site, by year"
)


# =============================================================================
#
# The full 2018 x 2025 transition detail (exactly which class each site
# moved to) is still written to CSV. The plotted figure summarises it
# instead: for each 2018 starting class, the share of sites that kept the
# same dominant class, drifted toward succession, or drifted toward
# restoration - using the same stacked-bar grammar as the succession/
# restoration pattern charts (05/10), rather than a 16-cell matrix. That's
# a deliberate simplicity-for-detail trade-off: the 3-way split answers
# "did it hold, and which way did it move" at a glance; the CSV is there
# for anyone who needs to know exactly which class it moved to.

dominant_matrix <- site_results %>%
  filter(direction_eligible, !is.na(dominant_base), !is.na(dominant_final)) %>%
  mutate(
    dominant_base = factor(dominant_base, levels = SJ_ZONE_LEVELS),
    dominant_final = factor(dominant_final, levels = SJ_ZONE_LEVELS)
  ) %>%
  count(dominant_base, dominant_final, name = "n_fids") %>%
  tidyr::complete(
    dominant_base = factor(SJ_ZONE_LEVELS, levels = SJ_ZONE_LEVELS),
    dominant_final = factor(SJ_ZONE_LEVELS, levels = SJ_ZONE_LEVELS),
    fill = list(n_fids = 0L)
  ) %>%
  group_by(dominant_base) %>%
  mutate(
    row_total = sum(n_fids),
    row_pct = if_else(row_total > 0, 100 * n_fids / row_total, NA_real_),
    tile_label = if_else(
      row_total > 0,
      sprintf("n = %s\n(%.1f%%)", skog_num(n_fids), row_pct),
      paste0("n = ", skog_num(n_fids))
    )
  ) %>%
  ungroup() %>%
  mutate(
    direction_group = case_when(
      as.integer(dominant_final) == as.integer(dominant_base) ~ "No change",
      as.integer(dominant_final) > as.integer(dominant_base) ~ "Succession",
      TRUE ~ "Restoration"
    ),
    direction_group = factor(direction_group, levels = c("Restoration", "No change", "Succession"))
  )

readr::write_csv(
  dominant_matrix %>% mutate(dominant_base = as.character(dominant_base), dominant_final = as.character(dominant_final)),
  file.path(OUTPUT_DIR, "13_dominant_SJ_class_baseline_final_matrix.csv")
)

dominant_summary <- dominant_matrix %>%
  group_by(dominant_base, direction_group) %>%
  summarise(n_fids = sum(n_fids), .groups = "drop") %>%
  group_by(dominant_base) %>%
  mutate(pct = 100 * n_fids / sum(n_fids)) %>%
  ungroup()

readr::write_csv(dominant_summary, file.path(OUTPUT_DIR, "14_dominant_SJ_class_summary_by_direction.csv"))

dominant_summary_labels <- dominant_summary %>% filter(pct >= 5)

p_dominant_matrix <- ggplot(dominant_summary, aes(x = dominant_base, y = pct, fill = direction_group)) +
  geom_col(width = 0.65, position = "stack") +
  geom_text(
    data = dominant_summary_labels,
    aes(label = sprintf("%.0f%%", pct)),
    position = position_stack(vjust = 0.5),
    colour = SKOG_THEME$surface, fontface = "bold", size = gg_text_size("footnote")
  ) +
  scale_fill_manual(values = DOMINANT_DIRECTION_COLOURS, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(
    title = "Dominant SJ Class: what happened to each 2018 starting class by 2025",
    subtitle = "Share of sites that kept the same dominant SJ Class, or drifted toward succession or restoration",
    x = "2018 dominant SJ Class",
    y = "Share of sites",
    caption = wrap_caption(paste0(
      "Grey = same dominant SJ Class in 2025. Red = dominant class moved toward SJ4 (succession). Blue = dominant ",
      "class moved toward SJ1 (restoration). Exactly which class each site moved to is in the accompanying CSV."
    ))
  ) +
  theme_report() +
  theme(legend.position = "top")


# =============================================================================
# 28. SAVE FIGURES
# =============================================================================
# The figures written to OUTPUT_DIR. The other ggplot objects built above
# (sections 17, 25, 26, 27) are available in the session but not saved.
figure_paths <- list(
  annual_tilstand         = file.path(OUTPUT_DIR, "00_annual_trajectories_by_tilstand.png"),
  decline_composition     = file.path(OUTPUT_DIR, "03_succession_composition_change.png"),
  decline_trajectories    = file.path(OUTPUT_DIR, "06_succession_trajectories.png"),
  recovery_composition    = file.path(OUTPUT_DIR, "08_restoration_composition_change.png")
)

plots_to_save <- list(
  list(p_tilstand, figure_paths$annual_tilstand, 10, 7),
  list(p_decline_composition, figure_paths$decline_composition, 12, 6),
  list(p_decline_trajectories, figure_paths$decline_trajectories, 9, 7),
  list(p_recovery_composition, figure_paths$recovery_composition, 12, 6)
)

for (x in plots_to_save) {
  if (!is.null(x[[1]])) {
    ggsave(filename = x[[2]], plot = x[[1]], width = x[[3]], height = x[[4]], dpi = 350, bg = "white")
  }
}


# =============================================================================
# 29. METHOD NOTES
# =============================================================================

method_notes <- c(
  "KYSTLYNGHEI AE-SJ TEMPORAL TRIAGE ANALYSIS",
  "",
  "Triage classification:",
  paste0(
    "  Moderate threshold: ", signif(moderate_threshold, 6), " AE-SJ units ",
    "(", signif(TRIAGE_MODERATE_FRACTION, 2), " x the SJ", TRIAGE_MODERATE_PAIR[1],
    "-SJ", TRIAGE_MODERATE_PAIR[2], " reference gap)."
  ),
  paste0(
    "  Strong threshold:   ", signif(large_threshold, 6), " AE-SJ units ",
    "(", signif(TRIAGE_LARGE_FRACTION, 2), " x the SJ", TRIAGE_LARGE_PAIR[1],
    "-SJ", TRIAGE_LARGE_PAIR[2], " reference gap)."
  ),
  "  Every direction-eligible site is classified by comparing its Theil-Sen",
  "  trend-equivalent change (slope x observed year span) against these two",
  "  absolute thresholds, in either direction. This is an exhaustive,",
  "  ecologically anchored triage rule - not a percentile screen - so the",
  "  number of sites in each tier reflects the actual population, not a",
  "  fixed quota.",
  "",
  "  'No strong signal' means the site did not cross the moderate threshold",
  "  in either direction over the monitoring period. It does not mean the",
  "  site is literally unchanged.",
  "",
  "Kendall's tau and BH-adjusted p-values are retained in the site-level",
  "table as diagnostic-only columns. They are not used to flag or filter any",
  "site: with 8 annual observations per site, formal significance testing",
  "after multiple-testing correction has very little statistical power.",
  "",
  "BT (grazing pressure) is read and used only to exclude BT5/BT6, which",
  "follow a different (grassland-bound) trajectory driven by overgrazing",
  "and are not on the AE-SJ forest-succession axis this score measures.",
  "BT1-4 is not currently used as a covariate elsewhere in the analysis.",
  "",
  "Successional and restorative mechanism blocks:",
  "  Every changing site (Moderate or Strong tier) is analysed in one of two",
  "  parallel blocks, split by direction:",
  "    Successional block - sites moving toward SJ4-like structure.",
  "                          Widespread, directionally coherent change across",
  "                          most of the site is consistent with a site-wide",
  "                          pressure and little visible localised",
  "                          management; change confined to part of the site",
  "                          is consistent with a localised disturbance or an",
  "                          untreated pocket.",
  "    Restorative block  - sites moving toward SJ1-like structure.",
  "                          Widespread change is consistent with a",
  "                          site-wide shift (e.g. reduced grazing pressure",
  "                          across the whole site); change confined to part",
  "                          of the site is consistent with targeted",
  "                          restoration work in that part of the site.",
  "  These patterns are descriptive spatial evidence, not proof of",
  "  management cause - they indicate whether the change is uniform or",
  "  localised, which is a necessary (but not sufficient) signature of",
  "  targeted intervention.",
  "",
  paste0(
    "  Widespread change: >= ", 100 * WIDESPREAD_CHANGED_THRESHOLD,
    "% of matched pixels cross an SJ Class boundary, with the coherence",
    " criterion below also met."
  ),
  paste0(
    "  Coherent: >= ", 100 * DIRECTIONAL_COHERENCE_THRESHOLD,
    "% of boundary-crossing pixels agree with the site-level direction."
  ),
  "  These are descriptive thresholds, not statistical cutoffs, and do not",
  "  measure spatial clustering (adjacency) of changed pixels.",
  "",
  "Composition and structural diversity (entropy):",
  "  Baseline = pixelwise mean AE-SJ across separately scored 2018 and 2019",
  "  images. Final = pixelwise mean AE-SJ across separately scored 2024 and",
  "  2025 images. Structural diversity (normalised Shannon entropy of the",
  "  four SJ-like proportions): 0 = composition concentrated in one SJ",
  "  Class, 1 = equal proportions across all four SJ Classes. Plotted as",
  "  an annual median by tilstand (2018-2025), not just a before/after pair.",
  ""
)

writeLines(method_notes, file.path(OUTPUT_DIR, "METHOD_NOTES.txt"))


# =============================================================================
# ADD-ON FIGURE: ALL SITE TRAJECTORIES BY TRIAGE CLASS
# =============================================================================
#
# Every direction-eligible site's annual AE-SJ trajectory.
#
# Colours:
#   Strong succession      = dark red
#   Moderate succession    = light red
#   No strong signal       = grey
#   Moderate restoration   = light blue
#   Strong restoration     = dark blue
#
# Lines are intentionally very thin and highly transparent because this
# figure is intended to show the full population distribution rather than
# individual sites.
# =============================================================================


all_trajectory_dat <- state_dat %>%
  inner_join(
    state_trends %>%
      filter(
        direction_eligible,
        !is.na(triage_class)
      ) %>%
      select(
        fid,
        triage_class
      ),
    by = "fid"
  ) %>%
  mutate(
    triage_class = factor(
      triage_class,
      levels = TRIAGE_LEVELS
    )
  )


# Number of sites in each category for the legend.
trajectory_class_n <- state_trends %>%
  filter(
    direction_eligible,
    !is.na(triage_class)
  ) %>%
  count(
    triage_class,
    name = "n_fids"
  ) %>%
  tidyr::complete(
    triage_class = factor(
      TRIAGE_LEVELS,
      levels = TRIAGE_LEVELS
    ),
    fill = list(
      n_fids = 0L
    )
  )


trajectory_legend_labels <- setNames(
  paste0(
    as.character(trajectory_class_n$triage_class),
    "  (n = ",
    skog_num(trajectory_class_n$n_fids),
    ")"
  ),
  as.character(trajectory_class_n$triage_class)
)


first_all_trajectory_year <- min(
  all_trajectory_dat$year,
  na.rm = TRUE
)


p_all_trajectories <- ggplot(
  all_trajectory_dat,
  aes(
    x = year,
    y = AE_SJ_p50,
    group = fid,
    colour = triage_class
  )
) +
  
  # Empirical SJ reference positions.
  sj_reference_layers(
    ref_dat,
    first_all_trajectory_year
  ) +
  
  # Draw no-signal sites first so coloured trajectories remain readable.
  geom_line(
    data = all_trajectory_dat %>%
      filter(
        triage_class == "No strong signal"
      ),
    linewidth = 0.18,
    alpha = 0.055
  ) +
  
  # Moderate succession.
  geom_line(
    data = all_trajectory_dat %>%
      filter(
        triage_class == "Moderate succession"
      ),
    linewidth = 0.20,
    alpha = 0.11
  ) +
  
  # Moderate restoration.
  geom_line(
    data = all_trajectory_dat %>%
      filter(
        triage_class == "Moderate restoration"
      ),
    linewidth = 0.20,
    alpha = 0.11
  ) +
  
  # Strong succession.
  geom_line(
    data = all_trajectory_dat %>%
      filter(
        triage_class == "Strong succession"
      ),
    linewidth = 0.24,
    alpha = 0.18
  ) +
  
  # Strong restoration.
  geom_line(
    data = all_trajectory_dat %>%
      filter(
        triage_class == "Strong restoration"
      ),
    linewidth = 0.24,
    alpha = 0.18
  ) +
  
  scale_colour_manual(
    values = TRIAGE_COLOURS,
    breaks = TRIAGE_LEVELS,
    labels = trajectory_legend_labels,
    drop = FALSE,
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = sort(
      unique(
        all_trajectory_dat$year
      )
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.02
      )
    )
  ) +
  
  labs(
    title = "Annual AE-SJ trajectories across all sites",
    subtitle = paste0(
      "Every direction-eligible kystlynghei site, coloured by temporal triage class"
    ),
    x = "Year",
    y = "Site median AE-SJ score",
    caption = wrap_caption(
      paste0(
        "Thin semi-transparent lines show individual site trajectories. ",
        "Red = succession toward SJ4-like structure; blue = restoration toward SJ1-like structure; ",
        "grey = no strong signal. The SJ reference lines show the empirical median positions of SJ1-SJ4."
      )
    )
  ) +
  
  theme_report() +
  
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "vertical"
  )





# =============================================================================
# 30. CONSOLE SUMMARY
# =============================================================================

cat(
  "\n============================================================\n",
  "AE-SJ TEMPORAL TRIAGE ANALYSIS COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

triage_summary_console <- state_trends %>%
  filter(direction_eligible) %>%
  count(triage_class, name = "n_fids")

cat("Direction-eligible sites: ", sum(state_trends$direction_eligible, na.rm = TRUE), "\n\n", sep = "")
print(triage_summary_console, n = Inf)

cat(
  "\nComposition rows passing QC: ", sum(composition_raw$composition_qc_ok, na.rm = TRUE),
  " / ", nrow(composition_raw), "\n",
  "Transition FIDs passing QC: ", sum(transition_analysis$transition_qc_ok, na.rm = TRUE),
  " / ", nrow(transition_analysis), "\n\n",
  sep = ""
)

cat("Succession change-pattern counts:\n")
print(decline_pattern_counts, n = Inf)

cat("\nRestoration change-pattern counts:\n")
print(recovery_pattern_counts, n = Inf)

cat("\nOutput directory:\n", OUTPUT_DIR, "\n", sep = "")

# ---- Trend distribution by tilstand ------------------------------------------


tilstand_dist_dat <- state_trends %>%
  filter(direction_eligible, !is.na(tilstand))

# One consistent bin width across all four panels, so bars line up and are
# comparable - not re-derived per panel.
trend_binwidth <- diff(range(tilstand_dist_dat$sj_change, na.rm = TRUE)) / 22

# Bins and percentages are computed by hand rather than left to
# geom_histogram()'s stat_bin: stat_bin's density is normalised WITHIN each
# fill group (each triage_class), not within each panel, so a group with
# very few sites (e.g. a tilstand's "Strong restoration", sometimes n=1)
# gets treated as 100% of itself and spikes - the percentages then don't
# mean "% of this tilstand" at all. Binning first and dividing by each
# tilstand's true total avoids that.
tilstand_dist_summary <- tilstand_dist_dat %>%
  mutate(
    bin_index = floor(sj_change / trend_binwidth),
    bin_mid = (bin_index + 0.5) * trend_binwidth
  ) %>%
  count(tilstand, bin_mid, triage_class, name = "n_in_bin") %>%
  group_by(tilstand) %>%
  mutate(pct_of_tilstand = 100 * n_in_bin / sum(n_in_bin)) %>%
  ungroup()

p_tilstand_distribution <- ggplot(tilstand_dist_summary, aes(x = bin_mid, y = pct_of_tilstand, fill = triage_class)) +
  geom_col(width = trend_binwidth, colour = SKOG_THEME$surface, linewidth = 0.15) +
  geom_vline(
    xintercept = c(-large_threshold, -moderate_threshold, moderate_threshold, large_threshold),
    linetype = "dashed", linewidth = 0.4, colour = SKOG_THEME$rule2
  ) + xlim(-0.5, 0.5) +
  facet_wrap(~ tilstand, ncol = 1, strip.position = "left") +
  scale_fill_manual(values = TRIAGE_COLOURS, drop = FALSE, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = "Where are successional and restorative sites concentrated?",
    subtitle = "Estimated AE-SJ change per site, by mapped condition class (tilstand)",
    x = "Estimated AE-SJ change over the monitoring period",
    y = "Share of sites in that condition class",
    caption = "Each panel sums to 100% of that tilstand class. Dashed lines mark the same moderate/strong change thresholds."
  ) +
  theme_report() +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  )


