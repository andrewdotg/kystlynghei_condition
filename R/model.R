# =============================================================================
# KYSTLYNGHEI ALPHAEARTH SJ — FINAL COVARIANCE-ADJUSTED AXIS
# =============================================================================
#
# Purpose
# -------
# Fit the final production AE-SJ model selected after method comparison.
#
# Validated configuration reproduced here:
#   * field calibration uses AlphaEarth from the year of assessment
#   * point weighted: every retained 0.25-ha sample contributes once
#   * 64-D AlphaEarth vectors receive a shrunk covariance adjustment
#   * shrinkage = 0.20
#   * covariance transform is estimated without using SJ labels
#   * BT5 and BT6 do NOT contribute to SJ endpoint calibration
#   * SJ1 and SJ4 alone define the structural succession axis
#   * final score is endpoint calibrated:
#         0 = median observed SJ4 reference
#         1 = median observed SJ1 reference
#   * score is NOT clipped
#
# IMPORTANT ABOUT BT5/BT6
# -----------------------
# The validation run that selected this model fitted the UNSUPERVISED covariance
# adjustment on all complete AE rows, including BT5/BT6, but excluded BT5/BT6
# from the supervised SJ endpoint calibration. This script deliberately
# reproduces that validated model exactly. BT5/BT6 are also excluded from the
# final point-score analysis CSV.
#
# If you instead want BT5/BT6 to have zero influence even on the unsupervised
# covariance estimate, change COVARIANCE_FIT_SCOPE to "analysis_eligible" and
# revalidate that altered model before treating it as the final production fit.
#
# Outputs
# -------
#   FINAL_AE_SJ_covadj_point_scores.csv
#   FINAL_AE_SJ_covadj_class_reference_positions.csv
#   FINAL_AE_SJ_covadj_model_summary.csv
#   FINAL_AE_SJ_covadj_adjustment_matrix_long.csv
#   FINAL_AE_SJ_covadj_adjustment_matrix_wide.csv
#   FINAL_AE_SJ_covadj_direction_64D.csv
#   FINAL_AE_SJ_covadj_combined_numerator_64D.csv
#   FINAL_AE_SJ_covadj_scale_parameters.csv
#   FINAL_AE_SJ_covadj_model_bundle.rds
#   FINAL_AE_SJ_covadj_model.json      the fitted model for python/ae_sj_model.py
#   FINAL_AE_SJ_covadj_score_by_SJ.png
#
# Scoring needs BOTH the 64x64 covariance-adjustment matrix and the 64-D axis
# direction. A single 64-D direction is not sufficient because the adjusted
# vector is L2-normalised after the matrix transform.
# =============================================================================


# =============================================================================
# 1. SETTINGS
# =============================================================================

INPUT_DIR <- "equalArea2/GEEpointStats/"
INPUT_REGEX <- "^kystlynghei_condition_points_0p25ha_surveyYear.*\\.csv$"

SHRINKAGE <- 0.20
EIGEN_FLOOR_REL <- 1e-6

EXCLUDE_BT_FROM_CALIBRATION <- c(5, 6)
ALLOW_MISSING_BT <- TRUE

# "all_complete" reproduces the validated best model.
# "analysis_eligible" also excludes BT5/BT6 from the covariance estimate.
COVARIANCE_FIT_SCOPE <- "all_complete"

MIN_ENDPOINT_POINTS <- 30
PLOT_MAX_POINTS <- 40000
RANDOM_SEED <- 42


# =============================================================================
# 2. PACKAGES
# =============================================================================

if (!exists("REPO_ROOT")) source("R/00_setup.R")
source(file.path(REPO_ROOT, "R/ae_sj_core.R"))

required_packages(c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "ggplot2",
  "jsonlite"
))
set.seed(RANDOM_SEED)


# =============================================================================
# 3. HELPERS
# =============================================================================
#
# clean_character, clean_mi_character, code_to_numeric, get_col_or_na,
# normalize_vector and clamp_cosine (below) now live in R/ae_sj_core.R,
# shared with validation.R. fit_endpoint_axis and score_axis stay here
# unchanged - validation.R has its own differently-behaved versions (CV
# train/test splitting, a different min-points default), so they were kept
# as two separate copies rather than merged. See ae_sj_core.R's header.

row_l2_normalize <- function(m) {
  nrm <- sqrt(rowSums(m^2))
  out <- matrix(
    NA_real_,
    nrow = nrow(m),
    ncol = ncol(m),
    dimnames = dimnames(m)
  )
  good <- is.finite(nrm) & nrm > 1e-12
  out[good, ] <- m[good, , drop = FALSE] / nrm[good]
  list(unit = out, norm = nrm, good = good)
}

fit_covariance_adjustment <- function(X, shrink = SHRINKAGE,
                                      eigen_floor_rel = EIGEN_FLOOR_REL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (nrow(X) < 100) stop("Too few rows to estimate covariance adjustment.")

  # cov() centres internally when estimating Sigma. We intentionally do NOT
  # subtract the mean when applying the transform. This reproduces the model
  # that was validated in the comparison analysis.
  S <- stats::cov(X)
  if (any(!is.finite(S))) stop("Non-finite covariance matrix.")

  p <- ncol(S)
  target_scale <- mean(diag(S))
  if (!is.finite(target_scale) || target_scale <= 0) {
    stop("Invalid covariance scale.")
  }

  S_adj <- (1 - shrink) * S + shrink * target_scale * diag(p)

  ee <- eigen(S_adj, symmetric = TRUE)
  eigen_floor <- max(ee$values) * eigen_floor_rel
  vals <- pmax(ee$values, eigen_floor)

  W <- ee$vectors %*%
    diag(1 / sqrt(vals), nrow = p) %*%
    t(ee$vectors)

  # Numerical symmetry check.
  W <- (W + t(W)) / 2

  list(
    W = W,
    shrink = shrink,
    raw_covariance = S,
    adjusted_covariance = S_adj,
    eigenvalues = ee$values,
    eigen_floor = eigen_floor,
    n_fit = nrow(X)
  )
}

apply_covariance_adjustment <- function(X, fit) {
  X %*% fit$W
}

fit_endpoint_axis <- function(unit_matrix, eligible, sj) {
  i1 <- which(eligible & sj == 1)
  i4 <- which(eligible & sj == 4)

  if (length(i1) < MIN_ENDPOINT_POINTS || length(i4) < MIN_ENDPOINT_POINTS) {
    stop(
      "Too few endpoint observations after filtering: SJ1=", length(i1),
      ", SJ4=", length(i4)
    )
  }

  p1 <- normalize_vector(colMeans(unit_matrix[i1, , drop = FALSE]))
  p4 <- normalize_vector(colMeans(unit_matrix[i4, , drop = FALSE]))
  direction <- normalize_vector(p1 - p4)

  raw1 <- as.numeric(unit_matrix[i1, , drop = FALSE] %*% direction)
  raw4 <- as.numeric(unit_matrix[i4, , drop = FALSE] %*% direction)

  a1 <- median(raw1, na.rm = TRUE)
  a4 <- median(raw4, na.rm = TRUE)

  if (!is.finite(a1) || !is.finite(a4) || a1 <= a4) {
    stop("Unexpected endpoint orientation: SJ1 anchor should exceed SJ4 anchor.")
  }

  list(
    prototype_SJ1 = p1,
    prototype_SJ4 = p4,
    direction = direction,
    anchor_SJ1 = a1,
    anchor_SJ4 = a4,
    n_SJ1 = length(i1),
    n_SJ4 = length(i4)
  )
}

score_axis <- function(unit_matrix, idx, fit) {
  raw <- as.numeric(unit_matrix[idx, , drop = FALSE] %*% fit$direction)

  score01 <-
    (raw - fit$anchor_SJ4) /
    (fit$anchor_SJ1 - fit$anchor_SJ4)

  cos1 <- clamp_cosine(as.numeric(
    unit_matrix[idx, , drop = FALSE] %*% fit$prototype_SJ1
  ))
  cos4 <- clamp_cosine(as.numeric(
    unit_matrix[idx, , drop = FALSE] %*% fit$prototype_SJ4
  ))

  tibble(
    AE_SJ_raw_projection = raw,
    AE_SJ = score01,
    cosine_to_SJ1_prototype = cos1,
    cosine_to_SJ4_prototype = cos4,
    prototype_cosine_margin = cos1 - cos4
  )
}


# =============================================================================
# 4. INPUT FILES
# =============================================================================

if (is.null(INPUT_DIR)) {
  chosen_file <- file.choose()
  input_dir <- dirname(chosen_file)
} else {
  input_dir <- normalizePath(INPUT_DIR, mustWork = TRUE)
}

input_files <- sort(list.files(
  input_dir,
  pattern = INPUT_REGEX,
  full.names = TRUE,
  ignore.case = TRUE
))

if (length(input_files) == 0) {
  stop("No matching calibration CSVs found in: ", input_dir)
}

# Outputs go to outputs/model/ inside the repo.
output_dir <- file.path(OUTPUTS_ROOT, "model")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nInput files: ", length(input_files), "\n", sep = "")
cat(paste0("  ", basename(input_files), collapse = "\n"), "\n")

# read_one_csv() is defined in R/ae_sj_core.R.
raw <- dplyr::bind_rows(lapply(input_files, read_one_csv))

required_fields <- c("parentID", "_survey_year_num", "7RA-SJ")
missing_required <- setdiff(required_fields, names(raw))
if (length(missing_required) > 0) {
  stop("Missing required fields: ", paste(missing_required, collapse = ", "))
}


# =============================================================================
# 5. STANDARDISE FIELD VARIABLES AND DEDUPLICATE
# =============================================================================

raw_year <- code_to_numeric(raw[["_survey_year_num"]])
year_std <- ifelse(
  is.finite(raw_year),
  raw_year,
  raw$year_from_filename
)

parent_id <- clean_character(raw[["parentID"]])
sj <- code_to_numeric(raw[["7RA-SJ"]])
bt <- code_to_numeric(get_col_or_na(raw, "7JB-BT"))

survey_parent_key <- paste(year_std, parent_id, sep = "__")

x_25833 <- code_to_numeric(get_col_or_na(raw, "sample_center_x_25833"))
y_25833 <- code_to_numeric(get_col_or_na(raw, "sample_center_y_25833"))
sample_id <- clean_character(get_col_or_na(raw, "sample_id"))
system_index <- clean_character(get_col_or_na(raw, "system:index"))

coord_ok <- is.finite(x_25833) & is.finite(y_25833)
sample_key <- rep(NA_character_, nrow(raw))

sample_key[coord_ok] <- paste(
  year_std[coord_ok],
  parent_id[coord_ok],
  sprintf("%.2f", x_25833[coord_ok]),
  sprintf("%.2f", y_25833[coord_ok]),
  sep = "__"
)

no_key <- is.na(sample_key)
has_sample_id <- no_key & !is.na(sample_id)
sample_key[has_sample_id] <- paste(
  year_std[has_sample_id], parent_id[has_sample_id], sample_id[has_sample_id],
  sep = "__"
)

no_key <- is.na(sample_key)
has_system_index <- no_key & !is.na(system_index)
sample_key[has_system_index] <- paste(
  year_std[has_system_index], parent_id[has_system_index],
  system_index[has_system_index], sep = "__"
)

no_key <- is.na(sample_key)
sample_key[no_key] <- paste(
  year_std[no_key], raw$source_file[no_key], which(no_key), sep = "__"
)

dat <- raw %>%
  mutate(
    year_std = year_std,
    parent_id = parent_id,
    survey_parent_key = survey_parent_key,
    sj = sj,
    bt = bt,
    x_25833 = x_25833,
    y_25833 = y_25833,
    sample_key = sample_key
  ) %>%
  arrange(source_file) %>%
  distinct(sample_key, .keep_all = TRUE)

rm(raw)
gc()

# Confirm assessment-year pairing where the filename contains a year.
year_mismatch <-
  is.finite(dat$year_std) &
  is.finite(dat$year_from_filename) &
  dat$year_std != dat$year_from_filename

cat("\nAssessment-year mismatches: ", sum(year_mismatch, na.rm = TRUE), "\n", sep = "")

if (any(year_mismatch, na.rm = TRUE)) {
  stop(
    "Assessment year disagrees with filename year for at least one row. ",
    "Resolve this before fitting the final model."
  )
}


# =============================================================================
# 6. BUILD 64-D ALPHAEARTH MATRIX
# =============================================================================

ae_columns <- sort(grep(
  "^AE_A[0-9]{2}_mean$",
  names(dat),
  value = TRUE
))

if (length(ae_columns) != 64) {
  stop("Expected 64 AlphaEarth bands; found ", length(ae_columns))
}

ae_raw <- sapply(ae_columns, function(v) code_to_numeric(dat[[v]]))
ae_raw <- as.matrix(ae_raw)
colnames(ae_raw) <- ae_columns

complete_ae <- apply(ae_raw, 1, function(z) all(is.finite(z)))

# Analysis/calibration rows: valid SJ and BT5/BT6 removed.
analysis_eligible <-
  complete_ae &
  is.finite(dat$sj) &
  dat$sj %in% 1:4

if (ALLOW_MISSING_BT) {
  analysis_eligible <- analysis_eligible &
    (is.na(dat$bt) | !(dat$bt %in% EXCLUDE_BT_FROM_CALIBRATION))
} else {
  analysis_eligible <- analysis_eligible &
    is.finite(dat$bt) &
    !(dat$bt %in% EXCLUDE_BT_FROM_CALIBRATION)
}

cat("Complete AE rows: ", format(sum(complete_ae), big.mark = ","), "\n", sep = "")
cat(
  "Final analysis rows after SJ + BT filtering: ",
  format(sum(analysis_eligible), big.mark = ","), "\n", sep = ""
)


# =============================================================================
# 7. FIT COVARIANCE ADJUSTMENT
# =============================================================================

if (COVARIANCE_FIT_SCOPE == "all_complete") {
  cov_fit_idx <- which(complete_ae)
} else if (COVARIANCE_FIT_SCOPE == "analysis_eligible") {
  cov_fit_idx <- which(analysis_eligible)
} else {
  stop("COVARIANCE_FIT_SCOPE must be 'all_complete' or 'analysis_eligible'.")
}

covadj_fit <- fit_covariance_adjustment(
  ae_raw[cov_fit_idx, , drop = FALSE],
  shrink = SHRINKAGE
)

covadj_raw <- matrix(
  NA_real_,
  nrow = nrow(ae_raw),
  ncol = ncol(ae_raw),
  dimnames = list(NULL, ae_columns)
)

covadj_raw[complete_ae, ] <- apply_covariance_adjustment(
  ae_raw[complete_ae, , drop = FALSE],
  covadj_fit
)

covadj_norm_obj <- row_l2_normalize(covadj_raw)
covadj_unit <- covadj_norm_obj$unit


# =============================================================================
# 8. FIT FINAL SJ1 -> SJ4 ENDPOINT GEOMETRY
# =============================================================================

final_fit <- fit_endpoint_axis(
  unit_matrix = covadj_unit,
  eligible = analysis_eligible,
  sj = dat$sj
)

analysis_idx <- which(analysis_eligible)
score_tbl <- score_axis(covadj_unit, analysis_idx, final_fit)

# Raw vector norms for diagnostics.
raw_norm <- sqrt(rowSums(ae_raw^2))

score_tbl <- score_tbl %>%
  mutate(
    AE_raw_vector_norm = raw_norm[analysis_idx],
    AE_covadj_vector_norm = covadj_norm_obj$norm[analysis_idx]
  )


# =============================================================================
# 9. CLASS REFERENCE POSITIONS ON FINAL 0 -> 1 SCALE
# =============================================================================

class_reference <- tibble(
  sj = dat$sj[analysis_idx],
  AE_SJ = score_tbl$AE_SJ
) %>%
  filter(is.finite(sj), sj %in% 1:4, is.finite(AE_SJ)) %>%
  group_by(sj) %>%
  summarise(
    n = n(),
    mean_AE_SJ = mean(AE_SJ),
    sd_AE_SJ = sd(AE_SJ),
    q05_AE_SJ = quantile(AE_SJ, 0.05, na.rm = TRUE),
    q25_AE_SJ = quantile(AE_SJ, 0.25, na.rm = TRUE),
    median_AE_SJ = median(AE_SJ),
    q75_AE_SJ = quantile(AE_SJ, 0.75, na.rm = TRUE),
    q95_AE_SJ = quantile(AE_SJ, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SJ_class = paste0("SJ", sj),
    description = c(
      "Intact",
      "Fallow",
      "Early regrowth",
      "Late regrowth"
    )
  ) %>%
  select(SJ_class, sj, description, everything())

cat("\nFinal class reference positions:\n")
print(class_reference, n = Inf, width = Inf)


# =============================================================================
# 10. SAVE TIDY POINT-SCORE CSV
# =============================================================================

# Keep all non-AE metadata from the source files, then append final model scores.
# This avoids duplicating the 64 calibration bands in the analysis CSV while
# preserving the field variables for later ecological analysis.
metadata_columns <- setdiff(names(dat), ae_columns)

point_scores <- bind_cols(
  dat[analysis_idx, metadata_columns, drop = FALSE],
  score_tbl
)

readr::write_csv(
  point_scores,
  file.path(output_dir, "FINAL_AE_SJ_covadj_point_scores.csv")
)

readr::write_csv(
  class_reference,
  file.path(output_dir, "FINAL_AE_SJ_covadj_class_reference_positions.csv")
)


# =============================================================================
# 11. MODEL SUMMARY
# =============================================================================

model_summary <- tibble(
  model = "Covariance-adjusted point-weighted AE-SJ endpoint axis",
  shrinkage = SHRINKAGE,
  covariance_fit_scope = COVARIANCE_FIT_SCOPE,
  covariance_fit_n = covadj_fit$n_fit,
  analysis_n = length(analysis_idx),
  endpoint_SJ1_n = final_fit$n_SJ1,
  endpoint_SJ4_n = final_fit$n_SJ4,
  excluded_BT = paste(EXCLUDE_BT_FROM_CALIBRATION, collapse = ","),
  allow_missing_BT = ALLOW_MISSING_BT,
  SJ4_mapped_value = 0,
  SJ1_mapped_value = 1,
  SJ4_raw_anchor = final_fit$anchor_SJ4,
  SJ1_raw_anchor = final_fit$anchor_SJ1,
  direction_norm = sqrt(sum(final_fit$direction^2)),
  covariance_condition_number =
    max(covadj_fit$eigenvalues) /
    max(min(covadj_fit$eigenvalues), .Machine$double.eps),
  covariance_eigen_floor = covadj_fit$eigen_floor
)

readr::write_csv(
  model_summary,
  file.path(output_dir, "FINAL_AE_SJ_covadj_model_summary.csv")
)


# =============================================================================
# 12. EXPORT MODEL PARAMETERS
# =============================================================================

# R uses row observations: X_adjusted = X %*% W. The exports use column
# vectors (x_adjusted = W_cols %*% x, one row per adjusted dimension).
# W is mathematically symmetric, but the explicit orientation is exported.
W_cols <- t(covadj_fit$W)

ae_bands <- sub("^AE_", "", ae_columns)
ae_bands <- sub("_mean$", "", ae_bands)

if (!identical(ae_bands, sprintf("A%02d", 0:63))) {
  stop("Unexpected AlphaEarth band order. Expected A00 ... A63.")
}

adj_dims <- sprintf("C%02d", 0:63)

# Long adjustment matrix: one row per output-dimension x input-band coefficient.
adjustment_long <- expand.grid(
  output_dimension = adj_dims,
  input_band = ae_bands,
  stringsAsFactors = FALSE
) %>%
  arrange(output_dimension, input_band)

adjustment_long$coefficient <- as.numeric(t(W_cols))

# Safer explicit construction to guarantee row-major mapping.
adjustment_long <- bind_rows(lapply(seq_along(adj_dims), function(j) {
  tibble(
    output_dimension = adj_dims[j],
    input_band = ae_bands,
    coefficient = as.numeric(W_cols[j, ])
  )
}))

readr::write_csv(
  adjustment_long,
  file.path(output_dir, "FINAL_AE_SJ_covadj_adjustment_matrix_long.csv")
)

adjustment_wide <- as.data.frame(W_cols, check.names = FALSE)
names(adjustment_wide) <- ae_bands
adjustment_wide <- tibble::as_tibble(adjustment_wide) %>%
  mutate(output_dimension = adj_dims, .before = 1)

readr::write_csv(
  adjustment_wide,
  file.path(output_dir, "FINAL_AE_SJ_covadj_adjustment_matrix_wide.csv")
)

direction_export <- tibble(
  adjusted_dimension = adj_dims,
  coefficient = as.numeric(final_fit$direction)
)

readr::write_csv(
  direction_export,
  file.path(output_dir, "FINAL_AE_SJ_covadj_direction_64D.csv")
)

# Numerator coefficients in the ORIGINAL AE bands:
# raw numerator = x %*% W %*% direction.
# These are useful for checking an implementation, but they are not sufficient
# by themselves because the denominator is ||x %*% W||.
combined_numerator <- as.numeric(covadj_fit$W %*% final_fit$direction)

readr::write_csv(
  tibble(
    ae_band = ae_bands,
    coefficient = combined_numerator
  ),
  file.path(output_dir, "FINAL_AE_SJ_covadj_combined_numerator_64D.csv")
)

scale_parameters <- tibble(
  parameter = c(
    "SJ4_RAW_ANCHOR",
    "SJ1_RAW_ANCHOR",
    "SJ4_MAPPED_VALUE",
    "SJ1_MAPPED_VALUE",
    "COVARIANCE_SHRINKAGE"
  ),
  value = c(
    final_fit$anchor_SJ4,
    final_fit$anchor_SJ1,
    0,
    1,
    SHRINKAGE
  )
)

readr::write_csv(
  scale_parameters,
  file.path(output_dir, "FINAL_AE_SJ_covadj_scale_parameters.csv")
)


# =============================================================================
# 13. MODEL FILE FOR PYTHON SCORING
# =============================================================================
#
# python/ae_sj_model.py reads this file, so the national Python scoring
# applies exactly the model fitted here (17 significant digits: full double precision).

jsonlite::write_json(
  list(
    description = "Covariance-adjusted AlphaEarth 7RA-SJ succession axis (AE-SJ)",
    scoring = paste(
      "x_adj = adjustment_matrix %*% x (x = 64 AlphaEarth bands in band order);",
      "x_unit = x_adj / ||x_adj||; raw = sum(x_unit * direction);",
      "AE_SJ = (raw - sj4_raw_anchor) / (sj1_raw_anchor - sj4_raw_anchor);",
      "0 = median SJ4 reference, 1 = median SJ1 reference, unclipped"
    ),
    bands = ae_bands,
    adjustment_matrix = unname(W_cols),
    direction = as.numeric(final_fit$direction),
    sj4_raw_anchor = final_fit$anchor_SJ4,
    sj1_raw_anchor = final_fit$anchor_SJ1,
    covariance_shrinkage = SHRINKAGE,
    class_reference = class_reference
  ),
  file.path(output_dir, "FINAL_AE_SJ_covadj_model.json"),
  digits = I(17), auto_unbox = TRUE, pretty = TRUE
)


# =============================================================================
# 14. SAVE MODEL BUNDLE FOR R ANALYSES
# =============================================================================

model_bundle <- list(
  description = "Validated covariance-adjusted point-weighted AlphaEarth 7RA-SJ succession axis",
  score_description = "0 = median observed SJ4 reference; 1 = median observed SJ1 reference; unclipped",
  shrinkage = SHRINKAGE,
  covariance_fit_scope = COVARIANCE_FIT_SCOPE,
  exclude_BT_from_calibration = EXCLUDE_BT_FROM_CALIBRATION,
  allow_missing_BT = ALLOW_MISSING_BT,
  r_bands = ae_columns,
  ae_bands = ae_bands,
  adjusted_dimensions = adj_dims,
  adjustment_matrix = covadj_fit$W,
  adjustment_matrix_cols = W_cols,
  direction = final_fit$direction,
  prototype_SJ1 = final_fit$prototype_SJ1,
  prototype_SJ4 = final_fit$prototype_SJ4,
  SJ1_RAW_ANCHOR = final_fit$anchor_SJ1,
  SJ4_RAW_ANCHOR = final_fit$anchor_SJ4,
  class_reference = class_reference,
  model_summary = model_summary
)

saveRDS(
  model_bundle,
  file.path(output_dir, "FINAL_AE_SJ_covadj_model_bundle.rds")
)


# =============================================================================
# 15. DESCRIPTIVE FINAL-MODEL FIGURE
# =============================================================================

plot_dat <- tibble(
  sj = factor(
    dat$sj[analysis_idx],
    levels = 1:4,
    labels = c(
      "SJ1\nIntact",
      "SJ2\nFallow",
      "SJ3\nEarly regrowth",
      "SJ4\nLate regrowth"
    )
  ),
  AE_SJ = score_tbl$AE_SJ
) %>%
  filter(is.finite(AE_SJ))

if (nrow(plot_dat) > PLOT_MAX_POINTS) {
  set.seed(RANDOM_SEED)
  plot_points <- plot_dat %>% slice_sample(n = PLOT_MAX_POINTS)
} else {
  plot_points <- plot_dat
}

p_final <- ggplot(plot_dat, aes(x = sj, y = AE_SJ)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(
    data = plot_points,
    width = 0.12,
    height = 0,
    alpha = 0.04,
    size = 0.35
  ) +
  geom_hline(yintercept = c(0, 1), linetype = "dotted", linewidth = 0.4) +
  theme_classic(base_size = 12) +
  labs(
    x = "Observed 7RA-SJ",
    y = "Covariance-adjusted AlphaEarth SJ score",
    title = "Final AlphaEarth coastal-heath succession axis",
    subtitle = "Point weighted; BT5/BT6 excluded from endpoint calibration; SJ4 median = 0, SJ1 median = 1"
  )



# NDVItest ----------------------------------------------------------------

# =============================================================================
# FINAL MODEL FIGURE COLOURED BY NDVI
# =============================================================================

plot_dat_ndvi <- tibble(
  sj = factor(
    dat$sj[analysis_idx],
    levels = 1:4,
    labels = c(
      "SJ1\nIntact",
      "SJ2\nFallow",
      "SJ3\nEarly regrowth",
      "SJ4\nLate regrowth"
    )
  ),
  
  AE_SJ = score_tbl$AE_SJ,
  
  NDVI = code_to_numeric(
    dat$MTCI_season_mean_mean[analysis_idx]
  )
) %>%
  filter(
    is.finite(AE_SJ)
  )


# -------------------------------------------------------------------------
# Sample points for display only.
# Boxplots below still use ALL observations.
# -------------------------------------------------------------------------

plot_points_ndvi <- plot_dat_ndvi %>%
  filter(
    is.finite(NDVI)
  )

if (nrow(plot_points_ndvi) > PLOT_MAX_POINTS) {
  
  set.seed(RANDOM_SEED)
  
  plot_points_ndvi <-
    plot_points_ndvi %>%
    slice_sample(
      n = PLOT_MAX_POINTS
    )
}


# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

p_final_ndvi <-
  ggplot(
    plot_dat_ndvi,
    aes(
      x = sj,
      y = AE_SJ
    )
  ) +
  
  # Full-data distributions
  geom_boxplot(
    outlier.shape = NA,
    width = 0.5,
    colour = "black",
    fill = NA
  ) +
  
  # Display sample coloured by NDVI
  geom_jitter(
    data = plot_points_ndvi,
    aes(
      colour = NDVI
    ),
    width = 0.12,
    height = 0,
    alpha = 0.25,
    size = 0.55
  ) +
  
  # Endpoint reference positions
  geom_hline(
    yintercept = c(0, 1),
    linetype = "dotted",
    linewidth = 0.4
  ) +
  
  scale_colour_viridis_c(
    option = "D",
    name = "MTCI mean"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    x = "Observed 7RA-SJ",
    y = "Covariance-adjusted AlphaEarth SJ score",
    
    title =
      "AlphaEarth coastal-heath succession axis and vegetation greenness",
    
    subtitle =
      "Points coloured by seasonal mean NDVI; boxplots show all observations"
  ) +
  
  theme(
    legend.position = "right"
  )


# =============================================================================
# 16. FINISHED
# =============================================================================

cat("\n============================================================\n")
cat("FINAL COVARIANCE-ADJUSTED AE-SJ MODEL COMPLETE\n")
cat("============================================================\n")
cat("Output directory:\n", output_dir, "\n\n", sep = "")
cat("SJ4 raw anchor: ", final_fit$anchor_SJ4, " -> 0\n", sep = "")
cat("SJ1 raw anchor: ", final_fit$anchor_SJ1, " -> 1\n", sep = "")
cat("Direction norm: ", sqrt(sum(final_fit$direction^2)), "\n", sep = "")


###

# =============================================================================
# DESCRIPTIVE FINAL-MODEL FIGURE (standalone - reads the point-scores CSV)
# =============================================================================
library(readr)
library(dplyr)
library(ggplot2)
library(ggforce)   # install.packages("ggforce") if you don't have it - used for rounded box corners

# Set this if `output_dir` isn't already in your session from running the
# full script (it's deploy/AE_SJ_FINAL_COVARIANCE_ADJUSTED under INPUT_DIR).
# output_dir <- "path/to/deploy/AE_SJ_FINAL_COVARIANCE_ADJUSTED"

point_scores <- readr::read_csv(
  file.path(output_dir, "FINAL_AE_SJ_covadj_point_scores.csv"),
  show_col_types = FALSE
)

# ---- Style tokens (subset of the shared skog design system) ----------------
# Matches R/theme.R's 2026-09-22 redesign (kept as a local copy rather than
# sourcing theme.R here, same as before - see that file's own note on why
# this keeps an independent copy).
skog_ink          <- "#202822"
skog_ink3         <- "#8B927E"
skog_rule         <- "#DCDDCC"
skog_surface      <- "#FFFFFF"
skog_g3           <- "#4F7A52"  # "good" (moss)
skog_neutral      <- "#D3D6C4"  # "undecided"
skog_p3           <- "#B14A2A"  # "poor" (ember)
skog_p4           <- "#7A2E17"  # "poor", most severe
skog_accent_text  <- "#7A5289"  # explanatory-caption tone (heather)
skog_font_display <- "Fraunces"
skog_font_sans    <- "IBM Plex Sans"

SJ_CLASS_LEVELS <- c("SJ4", "SJ3", "SJ2", "SJ1")  # bottom -> top
SJ_CLASS_COLOURS <- c(SJ4 = skog_p4, SJ3 = skog_p3, SJ2 = skog_neutral, SJ1 = skog_g3)
SJ_CLASS_DESCRIPTIONS <- c(
  SJ1 = "SJ1  Intakt", SJ2 = "SJ2  Brakkleggingsfase",
  SJ3 = "SJ3  Tidlig gjenvekstsuksesjonsfase", SJ4 = "SJ4  Sein gjenvekstsuksesjonsfase"
)

plot_dat <- point_scores %>%
  transmute(sj_code = paste0("SJ", sj), AE_SJ = AE_SJ) %>%
  filter(is.finite(AE_SJ), sj_code %in% SJ_CLASS_LEVELS) %>%
  mutate(sj_code = factor(sj_code, levels = SJ_CLASS_LEVELS))

# Standard Tukey box-plot statistics (same convention geom_boxplot() uses),
# computed by hand so the box and whiskers can be drawn as separate,
# independently styled layers - rounded/stroke-free box, plain grey
# whiskers with caps.
box_stat <- function(x) {
  x <- x[is.finite(x)]
  q1 <- as.numeric(quantile(x, 0.25, na.rm = TRUE))
  q3 <- as.numeric(quantile(x, 0.75, na.rm = TRUE))
  iqr <- q3 - q1
  in_range <- x[x >= q1 - 1.5 * iqr & x <= q3 + 1.5 * iqr]
  tibble(
    q1 = q1,
    q3 = q3,
    median_AE_SJ = median(x, na.rm = TRUE),
    whisker_min = if (length(in_range) > 0) min(in_range) else q1,
    whisker_max = if (length(in_range) > 0) max(in_range) else q3
  )
}

sj_summary <- plot_dat %>%
  group_by(sj_code) %>%
  summarise(n = n(), box_stat(AE_SJ), .groups = "drop") %>%
  mutate(
    y_pos = as.numeric(sj_code),
    class_label = SJ_CLASS_DESCRIPTIONS[as.character(sj_code)],
    stat_label = sprintf("n = %s   median %.2f", format(n, big.mark = " "), median_AE_SJ)
  )

box_half_height <- 0.22

box_corners <- bind_rows(lapply(seq_len(nrow(sj_summary)), function(i) {
  row <- sj_summary[i, ]
  tibble(
    sj_code = row$sj_code,
    x = c(row$q1, row$q3, row$q3, row$q1),
    y = c(row$y_pos - box_half_height, row$y_pos - box_half_height,
          row$y_pos + box_half_height, row$y_pos + box_half_height)
  )
}))

whisker_caps <- bind_rows(
  transmute(sj_summary, sj_code, x_cap = whisker_min, y_pos),
  transmute(sj_summary, sj_code, x_cap = whisker_max, y_pos)
)

p_final <- ggplot() +
  # Whisker line (min to max) and end caps, drawn first so the box covers
  # the middle portion, leaving only the whisker ends visible beyond it.
  geom_segment(
    data = sj_summary,
    aes(x = whisker_min, xend = whisker_max, y = y_pos, yend = y_pos),
    colour = skog_ink3, linewidth = 0.4
  ) +
  geom_segment(
    data = whisker_caps,
    aes(x = x_cap, xend = x_cap, y = y_pos - 0.08, yend = y_pos + 0.08),
    colour = skog_ink3, linewidth = 0.4
  ) +
  # Box: rounded corners, no stroke.
  ggforce::geom_shape(
    data = box_corners,
    aes(x = x, y = y, group = sj_code, fill = sj_code),
    radius = grid::unit(2, "mm"), colour = NA
  ) +
  # Median mark, drawn on top of the box.
  geom_segment(
    data = sj_summary,
    aes(x = median_AE_SJ, xend = median_AE_SJ,
        y = y_pos - box_half_height * 0.9, yend = y_pos + box_half_height * 0.9),
    colour = skog_surface, linewidth = 1.1, lineend = "round"
  ) +
  geom_vline(xintercept = c(0, 1), linetype = "dotted", linewidth = 0.4, colour = skog_ink3) +
  geom_text(
    data = sj_summary,
    aes(x = -Inf, y = y_pos + 0.14, label = class_label),
    hjust = 1, fontface = "bold", family = skog_font_sans, size = 3.6, colour = skog_ink
  ) +
  geom_text(
    data = sj_summary,
    aes(x = -Inf, y = y_pos - 0.16, label = stat_label),
    hjust = 1, colour = skog_ink3, family = "IBM Plex Mono", size = 3.1
  ) +
  scale_fill_manual(values = SJ_CLASS_COLOURS, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = NULL,
    title = "AlphaEarth coastal-heath succession axis",
    caption = paste0(
      "AE-SJ score 0 = MI SJ4 (late regrowth); ",
      "1 = MI SJ1 (intact)"
    )
  ) +
  theme_minimal(base_size = 12, base_family = skog_font_sans) +
  theme(
    plot.background = element_rect(fill = skog_surface, colour = NA),
    panel.background = element_rect(fill = skog_surface, colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = skog_rule, linewidth = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = skog_ink3),
    axis.line.x = element_line(colour = skog_ink3),
    axis.ticks.x = element_line(colour = skog_ink3),
    plot.title = element_text(colour = skog_ink, family = skog_font_display, face = "bold", size = 18),
    plot.caption = element_text(colour = skog_accent_text, hjust = 0, size = 10),
    plot.caption.position = "plot",
    plot.title.position = "plot",
    plot.margin = margin(18, 30, 14, 320)
  )

print(p_final)

ggsave(
  file.path(output_dir, "FINAL_AE_SJ_covadj_score_by_SJ_styled.png"),
  p_final,
  width = 12,
  height = 5.5,
  dpi = 350,
  bg = skog_surface
)
