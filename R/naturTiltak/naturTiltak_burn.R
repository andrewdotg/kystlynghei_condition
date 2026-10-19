# =============================================================================
# NATURTILTAK — POST-BURN MANAGEMENT-INTERVENTION TEST, USING FINAL
# COVARIANCE-ADJUSTED AE-SJ MODEL
# =============================================================================
#
# A management-intervention TEST, not the deployed MI-wide burn screen: uses
# naturTiltak polygons (independently logged management interventions) to
# check whether the AE-SJ method detects a real response in burnt cells
# relative to unburnt reference cells. Burn detection on the MI sites is
# python/mi_dnbr_burn.py (FIRMS-triggered Sentinel-2 dNBR), and its validation
# against registered burning is python/validation_registered_burns.py.
#
# It reads the final model bundle written by R/model.R.
#
# Final score orientation:
#   0 = median observed SJ4 reference
#   1 = median observed SJ1 reference
#   higher = more SJ1-like / more intact structural state
#   lower  = more SJ4-like / later-regrowth structural state
#   score is NOT clipped
#
# Temporal design:
#   -2 = two calendar years before burn
#   -1 = primary pre-burn baseline
#    0 = burn year (excluded from main analysis)
#    1+ = complete post-burn years
#
# Baseline-adjusted burnt-reference effect:
#   positive = burnt cells shifted toward SJ1 relative to reference cells
#   negative = burnt cells shifted toward SJ4 relative to reference cells
# =============================================================================


# =============================================================================
# 0. SHARED SETUP
# =============================================================================
#
# Sourced before SETTINGS below, since two of those settings (MODEL_BUNDLE_
# PATH, OUTPUT_DIR) now need REPO_ROOT, which this defines.
if (!exists("REPO_ROOT")) source("R/00_setup.R")


# =============================================================================
# 1. SETTINGS
# =============================================================================

PROJECT_DIR <-
  "/data/P-Prosjekter2/156004_08_effektovervaking_2026_2030_opsjon_2_fjernmaling"

AE_CSV_PATH <-
  "Burning/burntEqualAreaConditionMetrics/kystlynghei_v3_AE_SJ_annual_means_all_years.csv"

# Fitted model bundle from R/model.R.
MODEL_BUNDLE_PATH <-
  file.path(REPO_ROOT, "outputs", "model", "FINAL_AE_SJ_covadj_model_bundle.rds")

# Outputs go to outputs/naturTiltak_burn/.
OUTPUT_DIR <- file.path(REPO_ROOT, "outputs", "naturTiltak_burn")

FID_EXAMPLE <- "201"


# =============================================================================
# 2. PACKAGES
# =============================================================================

required_packages(c("readr", "dplyr", "tidyr", "ggplot2"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 3. PATH RESOLUTION
# =============================================================================

resolve_path <- function(path, project_dir = PROJECT_DIR) {
  candidates <- c(
    path,
    file.path(project_dir, path),
    file.path(project_dir, "Andrew_kystlynghei", path)
  )

  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) return(normalizePath(hit[1]))

  stop(
    "File not found. Checked:\n",
    paste(candidates, collapse = "\n")
  )
}

ae_csv_path <- resolve_path(AE_CSV_PATH)
model_bundle_path <- resolve_path(MODEL_BUNDLE_PATH)


# =============================================================================
# 4. READ FINAL MODEL BUNDLE
# =============================================================================

model <- readRDS(model_bundle_path)

required_model_fields <- c(
  "adjustment_matrix",
  "direction",
  "SJ1_RAW_ANCHOR",
  "SJ4_RAW_ANCHOR",
  "class_reference"
)

missing_model <- setdiff(required_model_fields, names(model))
if (length(missing_model) > 0) {
  stop("Model bundle is missing: ", paste(missing_model, collapse = ", "))
}

W <- as.matrix(model$adjustment_matrix)
direction <- as.numeric(model$direction)
SJ1_RAW_ANCHOR <- as.numeric(model$SJ1_RAW_ANCHOR)
SJ4_RAW_ANCHOR <- as.numeric(model$SJ4_RAW_ANCHOR)
class_reference <- model$class_reference

if (!all(dim(W) == c(64, 64))) stop("Adjustment matrix must be 64x64.")
if (length(direction) != 64) stop("Model direction must have 64 coefficients.")


# =============================================================================
# 5. READ ANNUAL CELL EMBEDDINGS
# =============================================================================

ae_raw <- readr::read_csv(
  ae_csv_path,
  guess_max = 100000,
  show_col_types = FALSE
)

# Simple clean-names equivalent for the fields used below.
names(ae_raw) <- names(ae_raw) %>%
  tolower() %>%
  gsub("[^a-z0-9]+", "_", .) %>%
  gsub("^_|_$", "", .)

ae_dat <- ae_raw
rm(ae_raw)

required_id_cols <- c(
  "sample_id",
  "parent_fid",
  "sample_class",
  "intervention_year",
  "observation_year",
  "years_since_intervention"
)

missing_id_cols <- setdiff(required_id_cols, names(ae_dat))
if (length(missing_id_cols) > 0) {
  stop("Missing identifier columns: ", paste(missing_id_cols, collapse = ", "))
}

ae_dat <- ae_dat %>%
  mutate(
    sample_id = as.character(sample_id),
    parent_fid = as.character(parent_fid),
    sample_class = as.character(sample_class),
    intervention_year = as.integer(intervention_year),
    observation_year = as.integer(observation_year),
    years_since_burn = as.integer(years_since_intervention)
  )

# Burn year is intentionally excluded.
ae_dat_analysis <- ae_dat %>%
  filter(years_since_burn %in% c(-2, -1) | years_since_burn >= 1)


# =============================================================================
# 6. APPLY FINAL COVARIANCE-ADJUSTED AE-SJ MODEL
# =============================================================================

# After name cleaning, AE_A00_mean becomes ae_a00_mean.
ae_mean_vars <- sprintf("ae_a%02d_mean", 0:63)
missing_ae_vars <- setdiff(ae_mean_vars, names(ae_dat_analysis))

if (length(missing_ae_vars) > 0) {
  stop("Missing AlphaEarth mean variables: ", paste(missing_ae_vars, collapse = ", "))
}

X <- ae_dat_analysis %>%
  select(all_of(ae_mean_vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

complete_vector <- apply(X, 1, function(z) all(is.finite(z)))

X_adj <- matrix(NA_real_, nrow = nrow(X), ncol = 64)
X_adj[complete_vector, ] <- X[complete_vector, , drop = FALSE] %*% W

adj_norm <- sqrt(rowSums(X_adj^2))
valid_vector <- complete_vector & is.finite(adj_norm) & adj_norm > 1e-12

X_unit <- matrix(NA_real_, nrow = nrow(X), ncol = 64)
X_unit[valid_vector, ] <-
  X_adj[valid_vector, , drop = FALSE] /
  adj_norm[valid_vector]

raw_score <- rep(NA_real_, nrow(X))
raw_score[valid_vector] <- as.numeric(
  X_unit[valid_vector, , drop = FALSE] %*% direction
)

# 0 = SJ4 median, 1 = SJ1 median. Unclipped by design.
AE_SJ <-
  (raw_score - SJ4_RAW_ANCHOR) /
  (SJ1_RAW_ANCHOR - SJ4_RAW_ANCHOR)

sj_cells <- ae_dat_analysis %>%
  mutate(
    AE_SJ_raw = raw_score,
    AE_SJ = AE_SJ,
    AE_SJ_clipped01 = pmin(1, pmax(0, AE_SJ)),
    treatment = case_when(
      sample_class == "burnt" ~ "Burnt",
      sample_class %in% c("unburnt", "unburnt_control") ~ "Reference",
      TRUE ~ NA_character_
    ),
    treatment = factor(treatment, levels = c("Burnt", "Reference"))
  ) %>%
  filter(!is.na(treatment), is.finite(AE_SJ))

cat("\nFinal AE-SJ score summary:\n")
print(summary(sj_cells$AE_SJ))
cat("Below 0: ", sum(sj_cells$AE_SJ < 0), "\n", sep = "")
cat("Above 1: ", sum(sj_cells$AE_SJ > 1), "\n", sep = "")


# =============================================================================
# 7. CELL -> FID x YEAR x TREATMENT
# =============================================================================

sj_fid_year <- sj_cells %>%
  group_by(
    parent_fid,
    intervention_year,
    observation_year,
    years_since_burn,
    treatment
  ) %>%
  summarise(
    n_cells = n_distinct(sample_id),
    mean_AE_SJ = mean(AE_SJ, na.rm = TRUE),
    median_AE_SJ = median(AE_SJ, na.rm = TRUE),
    sd_AE_SJ = sd(AE_SJ, na.rm = TRUE),
    .groups = "drop"
  )

sj_overall <- sj_fid_year %>%
  group_by(years_since_burn, treatment) %>%
  summarise(
    n_fids = n_distinct(parent_fid),
    overall_mean = mean(mean_AE_SJ, na.rm = TRUE),
    sd_between_fids = sd(mean_AE_SJ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se = sd_between_fids / sqrt(n_fids),
    ci_lower = overall_mean - 1.96 * se,
    ci_upper = overall_mean + 1.96 * se
  )

cell_counts <- sj_cells %>%
  group_by(years_since_burn, treatment) %>%
  summarise(n_cells = n_distinct(sample_id), .groups = "drop")

sj_overall <- sj_overall %>%
  left_join(cell_counts, by = c("years_since_burn", "treatment"))


# =============================================================================
# 8. BURNT - REFERENCE CONTRAST AND BASELINE ADJUSTMENT
# =============================================================================

sj_contrast <- sj_fid_year %>%
  select(
    parent_fid,
    intervention_year,
    observation_year,
    years_since_burn,
    treatment,
    mean_AE_SJ
  ) %>%
  tidyr::pivot_wider(names_from = treatment, values_from = mean_AE_SJ) %>%
  mutate(burnt_minus_reference = Burnt - Reference)

baseline_difference <- sj_contrast %>%
  filter(years_since_burn == -1) %>%
  select(
    parent_fid,
    baseline_burnt_minus_reference = burnt_minus_reference
  )

sj_contrast <- sj_contrast %>%
  left_join(baseline_difference, by = "parent_fid") %>%
  mutate(
    baseline_adjusted_difference =
      burnt_minus_reference - baseline_burnt_minus_reference
  )

sj_effect_overall <- sj_contrast %>%
  group_by(years_since_burn) %>%
  summarise(
    n_fids = sum(is.finite(baseline_adjusted_difference)),
    mean_effect = mean(baseline_adjusted_difference, na.rm = TRUE),
    sd_effect = sd(baseline_adjusted_difference, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se = sd_effect / sqrt(n_fids),
    ci_lower = mean_effect - 1.96 * se,
    ci_upper = mean_effect + 1.96 * se
  )


# =============================================================================
# 9. PRE-BURN STABILITY (-2 VS -1)
# =============================================================================

preburn_comparison <- sj_fid_year %>%
  filter(years_since_burn %in% c(-2, -1)) %>%
  select(parent_fid, treatment, years_since_burn, mean_AE_SJ) %>%
  mutate(year_label = paste0("year_", years_since_burn)) %>%
  select(-years_since_burn) %>%
  pivot_wider(names_from = year_label, values_from = mean_AE_SJ)

if (all(c("year_-2", "year_-1") %in% names(preburn_comparison))) {
  preburn_comparison <- preburn_comparison %>%
    mutate(preburn_change = `year_-1` - `year_-2`)
}


# =============================================================================
# 10. REFERENCE POSITIONS AND PLOT HELPERS
# =============================================================================

ref_tbl <- class_reference %>%
  transmute(
    sj = as.integer(sj),
    SJ_class = as.character(SJ_class),
    reference = as.numeric(median_AE_SJ)
  ) %>%
  arrange(sj)

ref_values <- ref_tbl$reference

relative_years <- sort(unique(sj_cells$years_since_burn))
x_breaks <- sort(unique(c(relative_years, 0)))
x_labels <- ifelse(x_breaks == 0, "Burn", as.character(x_breaks))

theme_sj <- theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.caption = element_text(size = 9, colour = "grey35", hjust = 0)
  )

SJ_CLASS_DESCRIPTIONS <- c(
  SJ1 = "SJ1  Intakt",
  SJ2 = "SJ2  Brakkleggingsfase",
  SJ3 = "SJ3  Tidlig gjenvekstsuksesjonsfase",
  SJ4 = "SJ4  Sein gjenvekstsuksesjonsfase"
)

ref_label_tbl <- ref_tbl %>%
  mutate(
    sj_key = paste0("SJ", sj),
    label = unname(SJ_CLASS_DESCRIPTIONS[sj_key])
  )

# =============================================================================
# 11. FIGURE 1 — OVERALL BURNT VS REFERENCE TRAJECTORY
# =============================================================================

p_overall <- ggplot(
  sj_overall,
  aes(
    x = years_since_burn,
    y = overall_mean,
    colour = treatment,
    fill = treatment,
    group = treatment
  )
) +
  geom_hline(
    yintercept = ref_values,
    linetype = "dotted",
    linewidth = 0.45,
    colour = "grey45"
  ) +
  
  # SJ class labels
  geom_text(
    data = ref_label_tbl,
    aes(
      x = min(x_breaks),
      y = reference,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -0.35,
    size = 3.2,
    colour = "grey30"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  
  geom_ribbon(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    alpha = 0.18,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 1.1
  ) +
  
  geom_point(
    size = 2.6
  ) +
  
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  
  labs(
    title = "Covariance-adjusted AlphaEarth succession trajectories following heathland burning",
    x = "Years relative to burn",
    y = "AlphaEarth SJ score",
    caption = "Higher values are more SJ1-like. Ribbons show 95% confidence intervals across intervention FIDs."
  ) +
  
  theme_sj

# =============================================================================
# 12. FIGURE 2 — BASELINE-ADJUSTED BURNT - REFERENCE EFFECT
# =============================================================================

p_effect <- ggplot(
  sj_effect_overall,
  aes(x = years_since_burn, y = mean_effect)
) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, colour = "grey35") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  labs(
    title = "Burn-associated shift in the AlphaEarth succession gradient",
    subtitle = "Burnt-reference difference after subtracting each FID's year -1 pre-burn difference",
    x = "Years relative to burn",
    y = "Baseline-adjusted burnt - reference difference",
    caption = "Positive values indicate a shift of burnt cells toward the SJ1 end relative to local reference cells."
  ) +
  theme_sj


# =============================================================================
# 13. FIGURE 3 — ONE EXAMPLE FID
# =============================================================================

plot_sj_fid <- function(fid) {
  fid <- as.character(fid)
  cells <- sj_cells %>% filter(parent_fid == fid)
  if (nrow(cells) == 0) stop("FID not found: ", fid)

  means <- cells %>%
    group_by(years_since_burn, treatment) %>%
    summarise(
      n_cells = n_distinct(sample_id),
      mean_AE_SJ = mean(AE_SJ, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, colour = "grey35") +
    geom_hline(yintercept = ref_values, linetype = "dotted", linewidth = 0.35, colour = "grey55") +
    geom_line(
      data = cells,
      aes(x = years_since_burn, y = AE_SJ, colour = treatment, group = sample_id),
      linewidth = 0.45,
      alpha = 0.18
    ) +
    geom_point(
      data = cells,
      aes(x = years_since_burn, y = AE_SJ, colour = treatment),
      size = 0.9,
      alpha = 0.22
    ) +
    geom_line(
      data = means,
      aes(x = years_since_burn, y = mean_AE_SJ, colour = treatment, group = treatment),
      linewidth = 1.3
    ) +
    geom_point(
      data = means,
      aes(x = years_since_burn, y = mean_AE_SJ, colour = treatment),
      size = 2.7
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    labs(
      title = paste0("FID ", fid, ": covariance-adjusted AlphaEarth succession trajectory"),
      subtitle = "Thin lines = individual fixed 0.25-ha cells; thick lines = treatment means",
      x = "Years relative to burn",
      y = "AlphaEarth SJ score",
      caption = "SJ4 median = 0; SJ1 median = 1; burn calendar year omitted."
    ) +
    theme_sj
}

p_fid_example <- plot_sj_fid(FID_EXAMPLE)


# =============================================================================
# 14. FIGURE 4 — ALL BURNT FID TRAJECTORIES
# =============================================================================

burnt_fid <- sj_fid_year %>% filter(treatment == "Burnt")
overall_burnt <- sj_overall %>% filter(treatment == "Burnt")

p_all_burnt <- ggplot(
  burnt_fid,
  aes(x = years_since_burn, y = mean_AE_SJ, group = parent_fid)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, colour = "grey35") +
  geom_hline(yintercept = ref_values, linetype = "dotted", linewidth = 0.35, colour = "grey55") +
  geom_line(alpha = 0.18, linewidth = 0.45) +
  geom_line(
    data = overall_burnt,
    aes(x = years_since_burn, y = overall_mean, group = 1),
    inherit.aes = FALSE,
    linewidth = 1.4
  ) +
  geom_point(
    data = overall_burnt,
    aes(x = years_since_burn, y = overall_mean),
    inherit.aes = FALSE,
    size = 2.6
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  labs(
    title = "Burnt-area covariance-adjusted AlphaEarth succession trajectories",
    subtitle = "Thin lines = intervention FIDs; thick line = mean across FIDs",
    x = "Years relative to burn",
    y = "AlphaEarth SJ score"
  ) +
  theme_sj


# =============================================================================
# 15. FIGURE 5 — PRE-BURN STABILITY
# =============================================================================

p_preburn <- NULL
if (all(c("year_-2", "year_-1") %in% names(preburn_comparison))) {
  p_preburn <- ggplot(
    preburn_comparison,
    aes(x = `year_-2`, y = `year_-1`, colour = treatment)
  ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.6, colour = "grey40") +
    geom_point(alpha = 0.55, size = 1.8) +
    coord_equal() +
    labs(
      title = "Pre-burn stability of the covariance-adjusted AE-SJ signal",
      x = "AE-SJ score at year -2",
      y = "AE-SJ score at year -1"
    ) +
    theme_sj
}

# Methodology Figure 3 (p_overall). p_effect, p_fid_example, p_all_burnt and
# p_preburn are built above but not printed or saved.
print(p_overall)


# =============================================================================
# 16. SAVE TABLES
# =============================================================================

readr::write_csv(sj_cells, file.path(OUTPUT_DIR, "01_AE_SJ_covadj_cell_trajectories.csv"))
readr::write_csv(sj_fid_year, file.path(OUTPUT_DIR, "02_AE_SJ_covadj_fid_year_summary.csv"))
readr::write_csv(sj_overall, file.path(OUTPUT_DIR, "03_AE_SJ_covadj_overall_burnt_reference.csv"))
readr::write_csv(sj_contrast, file.path(OUTPUT_DIR, "04_AE_SJ_covadj_fid_burnt_reference_contrast.csv"))
readr::write_csv(sj_effect_overall, file.path(OUTPUT_DIR, "05_AE_SJ_covadj_baseline_adjusted_effect.csv"))
readr::write_csv(preburn_comparison, file.path(OUTPUT_DIR, "06_AE_SJ_covadj_preburn_minus2_minus1.csv"))


# =============================================================================
# 17. SAVE FIGURES
# =============================================================================

ggsave(file.path(OUTPUT_DIR, "01_AE_SJ_covadj_overall_burnt_reference.png"), p_overall, width = 8.5, height = 6, dpi = 400, bg = "white")


# =============================================================================
# 18. SAVE MASTER OBJECT
# =============================================================================

saveRDS(
  list(
    model_bundle_path = model_bundle_path,
    cell_scores = sj_cells,
    fid_year = sj_fid_year,
    overall = sj_overall,
    fid_contrast = sj_contrast,
    overall_effect = sj_effect_overall,
    preburn = preburn_comparison
  ),
  file.path(OUTPUT_DIR, "AE_SJ_covadj_postburn_analysis_bundle.rds")
)


# =============================================================================
# 19. FINISHED
# =============================================================================

cat("\n============================================\n")
cat("COVARIANCE-ADJUSTED AE-SJ POST-BURN ANALYSIS COMPLETE\n")
cat("============================================\n")
cat("Cells: ", dplyr::n_distinct(sj_cells$sample_id), "\n", sep = "")
cat("FIDs: ", dplyr::n_distinct(sj_cells$parent_fid), "\n", sep = "")
cat("Outputs saved to:\n", OUTPUT_DIR, "\n", sep = "")


# =============================================================================
# PRESCRIBED-BURN MASK SUMMARY
# =============================================================================
#
# Reads the exported burn masks:
#
#   1   = detected prescribed burn
#   0   = valid non-burn
#   255 = NoData / insufficient imagery
#
# Images containing ZERO burn pixels are discarded.
#
# Outputs:
#
#   1. Per-image burn statistics
#   2. Overall summary:
#        - number of retained burn images
#        - total detected burn area
#        - mean burn area per image
#        - mean burn percentage per image
#
# Burn percentage:
#
#   burn pixels / valid pixels (0 + 1) * 100
#
# =============================================================================


# -----------------------------------------------------------------------------
# LOCATION OF DOWNLOADED BURN MASKS
# -----------------------------------------------------------------------------

BURN_MASK_DIR <- file.path("deploy/BurningAnalysis/GEE_prescribed_burn_masks/"
)


if (!dir.exists(BURN_MASK_DIR)) {
  
  stop(
    "Burn-mask directory not found:\n",
    BURN_MASK_DIR
  )
}


# -----------------------------------------------------------------------------
# FIND MASK TIFFS
# -----------------------------------------------------------------------------

burn_mask_files <- list.files(
  BURN_MASK_DIR,
  pattern = "^prescribed_burn_mask_FID_.*\\.tif$",
  full.names = TRUE,
  ignore.case = TRUE
)


if (length(burn_mask_files) == 0) {
  
  stop(
    "No prescribed-burn mask TIFFs found in:\n",
    BURN_MASK_DIR
  )
}


cat(
  "\nBurn-mask TIFFs found: ",
  length(burn_mask_files),
  "\n",
  sep = ""
)


# =============================================================================
# PROCESS ONE MASK
# =============================================================================

summarise_burn_mask <- function(path) {
  
  r <- terra::rast(
    path
  )
  
  
  if (terra::nlyr(r) != 1) {
    
    stop(
      "Expected one raster band in:\n",
      path
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Extract FID and year from filename:
  #
  # prescribed_burn_mask_FID_201_2021.tif
  # ---------------------------------------------------------------------------
  
  filename <- basename(
    path
  )
  
  
  filename_match <- stringr::str_match(
    filename,
    "prescribed_burn_mask_FID_([^_]+)_([0-9]{4})\\.tif$"
  )
  
  
  fid <- filename_match[
    1,
    2
  ]
  
  
  year <- suppressWarnings(
    as.integer(
      filename_match[
        1,
        3
      ]
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Count burn pixels
  # ---------------------------------------------------------------------------
  
  n_burn <- terra::global(
    r == 1,
    fun = "sum",
    na.rm = TRUE
  )[
    1,
    1
  ]
  
  
  # ---------------------------------------------------------------------------
  # Count valid non-burn pixels
  # ---------------------------------------------------------------------------
  
  n_nonburn <- terra::global(
    r == 0,
    fun = "sum",
    na.rm = TRUE
  )[
    1,
    1
  ]
  
  
  n_burn <- as.numeric(
    n_burn
  )
  
  
  n_nonburn <- as.numeric(
    n_nonburn
  )
  
  
  if (!is.finite(n_burn)) {
    n_burn <- 0
  }
  
  
  if (!is.finite(n_nonburn)) {
    n_nonburn <- 0
  }
  
  
  # ---------------------------------------------------------------------------
  # Valid analysis area
  # ---------------------------------------------------------------------------
  
  n_valid <-
    n_burn +
    n_nonburn
  
  
  # ---------------------------------------------------------------------------
  # Pixel area
  #
  # These masks should be 20 m, but calculate it from the raster rather
  # than hard-coding 0.04 ha.
  # ---------------------------------------------------------------------------
  
  raster_resolution <-
    terra::res(
      r
    )
  
  
  pixel_area_m2 <-
    abs(
      raster_resolution[1] *
        raster_resolution[2]
    )
  
  
  pixel_area_ha <-
    pixel_area_m2 /
    10000
  
  
  # ---------------------------------------------------------------------------
  # Areas
  # ---------------------------------------------------------------------------
  
  burn_area_ha <-
    n_burn *
    pixel_area_ha
  
  
  valid_area_ha <-
    n_valid *
    pixel_area_ha
  
  
  # ---------------------------------------------------------------------------
  # Burn percentage
  # ---------------------------------------------------------------------------
  
  burn_percent <-
    if (
      n_valid > 0
    ) {
      
      100 *
        n_burn /
        n_valid
      
    } else {
      
      NA_real_
      
    }
  
  
  tibble::tibble(
    
    file =
      filename,
    
    fid =
      fid,
    
    intervention_year =
      year,
    
    pixel_size_m =
      raster_resolution[1],
    
    pixel_area_ha =
      pixel_area_ha,
    
    n_burn_pixels =
      n_burn,
    
    n_nonburn_pixels =
      n_nonburn,
    
    n_valid_pixels =
      n_valid,
    
    burn_area_ha =
      burn_area_ha,
    
    valid_area_ha =
      valid_area_ha,
    
    burn_percent =
      burn_percent
    
  )
}


# =============================================================================
# PROCESS ALL MASKS
# =============================================================================

burn_mask_all <- purrr::map_dfr(
  burn_mask_files,
  summarise_burn_mask
)


cat(
  "\nAll mask rasters processed: ",
  nrow(burn_mask_all),
  "\n",
  sep = ""
)


# =============================================================================
# DISCARD IMAGES WITH NO DETECTED BURN
# =============================================================================

burn_mask_detected <- burn_mask_all %>%
  
  dplyr::filter(
    n_burn_pixels > 0
  )


n_discarded <-
  nrow(burn_mask_all) -
  nrow(burn_mask_detected)


cat(
  "Masks containing burn: ",
  nrow(burn_mask_detected),
  "\n",
  sep = ""
)


cat(
  "Masks discarded because they contain no value-1 pixels: ",
  n_discarded,
  "\n",
  sep = ""
)


if (nrow(burn_mask_detected) == 0) {
  
  stop(
    "None of the burn masks contained any detected burn pixels."
  )
}


# =============================================================================
# OVERALL SUMMARY
# =============================================================================

prescribed_burn_summary <- burn_mask_detected %>%
  
  summarise(
    
    n_burn_images =
      n(),
    
    total_burn_area_ha =
      sum(
        burn_area_ha,
        na.rm = TRUE
      ),
    
    mean_burn_area_ha_per_image =
      mean(
        burn_area_ha,
        na.rm = TRUE
      ),
    
    median_burn_area_ha_per_image =
      median(
        burn_area_ha,
        na.rm = TRUE
      ),
    
    min_burn_area_ha =
      min(
        burn_area_ha,
        na.rm = TRUE
      ),
    
    max_burn_area_ha =
      max(
        burn_area_ha,
        na.rm = TRUE
      ),
    
    mean_burn_percent_per_image =
      mean(
        burn_percent,
        na.rm = TRUE
      ),
    
    median_burn_percent_per_image =
      median(
        burn_percent,
        na.rm = TRUE
      )
    
  )


# =============================================================================
# PRINT RESULTS
# =============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "PRESCRIBED-BURN MASK SUMMARY\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Downloaded masks examined: ",
  nrow(burn_mask_all),
  "\n",
  sep = ""
)


cat(
  "Masks with detected burn: ",
  prescribed_burn_summary$n_burn_images,
  "\n",
  sep = ""
)


cat(
  "Masks with no detected burn discarded: ",
  n_discarded,
  "\n\n",
  sep = ""
)


cat(
  "TOTAL DETECTED PRESCRIBED-BURN AREA: ",
  round(
    prescribed_burn_summary$total_burn_area_ha,
    2
  ),
  " ha\n",
  sep = ""
)


cat(
  "AVERAGE BURN AREA PER RETAINED IMAGE: ",
  round(
    prescribed_burn_summary$mean_burn_area_ha_per_image,
    2
  ),
  " ha\n",
  sep = ""
)


cat(
  "AVERAGE BURN PERCENT PER RETAINED IMAGE: ",
  round(
    prescribed_burn_summary$mean_burn_percent_per_image,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "\n============================================================\n"
)


# =============================================================================
# SHOW PER-IMAGE RESULTS
# =============================================================================

burn_mask_detected <- burn_mask_detected %>%
  
  arrange(
    fid,
    intervention_year
  )


print(
  burn_mask_detected %>%
    
    select(
      fid,
      intervention_year,
      burn_area_ha,
      valid_area_ha,
      burn_percent
    ),
  
  n = Inf
)


# =============================================================================
# SAVE RESULTS
# =============================================================================

readr::write_csv(
  
  burn_mask_detected,
  
  file.path(
    OUTPUT_DIR,
    "08_prescribed_burn_mask_per_image_summary.csv"
  )
  
)


readr::write_csv(
  
  prescribed_burn_summary,
  
  file.path(
    OUTPUT_DIR,
    "09_prescribed_burn_mask_overall_summary.csv"
  )
  
)


# =============================================================================
# PRESCRIBED BURNS
# YEAR BEFORE VS YEAR AFTER SJ-LIKE COMPOSITION
# =============================================================================
#
# Input:
#
#   prescribed_burn_mean_SJ_composition_minus1_plus1.csv
#
# Expected columns include:
#
#   period
#   years_since_burn
#   SJ_type
#   mean_fraction
#   mean_percentage
#   n_interventions
#   n_unique_fids
#
# The export already gives each FID/year intervention equal weight.
# Therefore this block DOES NOT average the values again.
# =============================================================================


# -----------------------------------------------------------------------------
# CSV PATH
# -----------------------------------------------------------------------------
#
# Path to the exported composition CSV.
# -----------------------------------------------------------------------------

PRESCRIBED_SJ_COMPOSITION_CSV <- file.path("deploy/BurningAnalysis/prescribed_burn_mean_SJ_composition_minus1_plus1.csv"
)


if (!file.exists(
  PRESCRIBED_SJ_COMPOSITION_CSV
)) {
  
  stop(
    "Prescribed-burn SJ composition CSV not found:\n",
    PRESCRIBED_SJ_COMPOSITION_CSV
  )
}


# =============================================================================
# READ CSV
# =============================================================================

prescribed_sj_comp <- readr::read_csv(
  PRESCRIBED_SJ_COMPOSITION_CSV,
  show_col_types = FALSE
)


required_fields <- c(
  "period",
  "SJ_type",
  "mean_fraction"
)


missing_fields <- setdiff(
  required_fields,
  names(
    prescribed_sj_comp
  )
)


if (length(
  missing_fields
) > 0) {
  
  stop(
    "Prescribed SJ CSV is missing: ",
    paste(
      missing_fields,
      collapse = ", "
    )
  )
}


# =============================================================================
# PREPARE GRAPH DATA
# =============================================================================

prescribed_sj_plot_dat <-
  prescribed_sj_comp |>
  
  dplyr::mutate(
    
    period =
      dplyr::recode(
        period,
        "Year before" = "Before",
        "Year after" = "After"
      ),
    
    period =
      factor(
        period,
        levels = c(
          "Before",
          "After"
        )
      ),
    
    SJ_type =
      factor(
        SJ_type,
        levels = c(
          "SJ1",
          "SJ2",
          "SJ3",
          "SJ4"
        )
      ),
    
    mean_fraction =
      as.numeric(
        mean_fraction
      ),
    
    percentage =
      100 *
      mean_fraction
    
  ) |>
  
  dplyr::arrange(
    period,
    SJ_type
  )


cat(
  "\nPrescribed-burn SJ composition used for graph:\n"
)


print(
  prescribed_sj_plot_dat,
  n = Inf
)


# =============================================================================
# CHECK THAT EACH BAR SUMS TO ~100%
# =============================================================================

prescribed_sj_sum_check <-
  prescribed_sj_plot_dat |>
  
  dplyr::group_by(
    period
  ) |>
  
  dplyr::summarise(
    
    total_percentage =
      sum(
        percentage,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


cat(
  "\nComposition sum check:\n"
)


print(
  prescribed_sj_sum_check
)


# =============================================================================
# SJ COLOURS
# =============================================================================

SJ_COLOURS <- c(
  
  "SJ1" =
    "#2a78d6",
  
  "SJ2" =
    "#86b6ef",
  
  "SJ3" =
    "#e09b83",
  
  "SJ4" =
    "#b4442d"
  
)


# =============================================================================
# SAMPLE SIZE FOR SUBTITLE
# =============================================================================

n_interventions_plot <-
  if (
    "n_interventions" %in%
    names(
      prescribed_sj_plot_dat
    )
  ) {
    
    unique(
      prescribed_sj_plot_dat$n_interventions
    )[1]
    
  } else {
    
    NA_integer_
    
  }


n_fids_plot <-
  if (
    "n_unique_fids" %in%
    names(
      prescribed_sj_plot_dat
    )
  ) {
    
    unique(
      prescribed_sj_plot_dat$n_unique_fids
    )[1]
    
  } else {
    
    NA_integer_
    
  }


plot_subtitle <-
  if (
    is.finite(
      n_interventions_plot
    ) &&
    is.finite(
      n_fids_plot
    )
  ) {
    
    paste0(
      "Composition change over 64 kystlynghei sites where prescribed burns occured"
    )
    
  } else {
    
    paste0(
      "Detected prescribed-burn pixels; ",
      "each intervention contributes equally"
    )
    
  }


# =============================================================================
# STACKED BEFORE / AFTER COMPOSITION PLOT
# =============================================================================

p_prescribed_sj <-
  ggplot2::ggplot(
    
    prescribed_sj_plot_dat,
    
    ggplot2::aes(
      x = period,
      y = percentage,
      fill = SJ_type
    )
    
  ) +
  
  ggplot2::geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.25
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      
      label =
        ifelse(
          
          percentage >= 3,
          
          paste0(
            round(
              percentage,
              1
            ),
            "%"
          ),
          
          ""
          
        )
      
    ),
    
    position =
      ggplot2::position_stack(
        vjust = 0.5
      ),
    
    size = 3.4
    
  ) +
  
  ggplot2::scale_fill_manual(
    values = SJ_COLOURS,
    drop = FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    
    limits = c(
      0,
      100
    ),
    
    breaks = seq(
      0,
      100,
      20
    ),
    
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.02
        )
      )
    
  ) +
  
  ggplot2::labs(
    
    title =
      "Prescribed burns: SJ composition before and after burning",
    
    subtitle =
      plot_subtitle,
    
    x =
      NULL,
    
    y =
      "Mean fraction of detected burn pixels",
    
    fill =
      "SJ type",
  
    
  ) +
  
  theme_sj


# Composition of the burnt cells before and after burning (prescribed burns).
print(
  p_prescribed_sj
)


# =============================================================================
# SAVE FIGURE
# =============================================================================

ggplot2::ggsave(

  filename =
    file.path(OUTPUT_DIR, "prescribed_burn_SJ_composition_before_after.png"
    ),
  
  plot =
    p_prescribed_sj,
  
  width =
    8,
  
  height =
    7,
  
  dpi =
    400,
  
  bg =
    "white"
  
)