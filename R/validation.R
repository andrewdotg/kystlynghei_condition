# =============================================================================
# KYSTLYNGHEI ALPHAEARTH SJ METHOD COMPARISON
# Whitening + reference-cloud cosine distances + BT5/6 exclusion
# =============================================================================
#
# PURPOSE
# -------
# Compare the existing AlphaEarth SJ-axis approach with several amendments:
#
#   M0  Native AE axis, original-style fit using all BT classes
#       (comparison baseline only)
#
#   M1  Native AE axis, but BT5/BT6 excluded from axis construction
#       (candidate baseline)
#
#   M2  Shrunk-whitened AE axis, BT5/BT6 excluded from axis construction
#
#   M3  Native-space endpoint reference-cloud score:
#         mean cosine distance to k nearest SJ1 references
#         mean cosine distance to k nearest SJ4 references
#         cloud score = d_SJ4 / (d_SJ1 + d_SJ4)
#       so SJ4-like -> 0 and SJ1-like -> 1
#
#   M4  Same reference-cloud score in shrunk-whitened AE space
#
# The script also exports diagnostics designed to identify points that project
# as SJ1-like on the linear axis but are nevertheless unlike the actual SJ1
# reference cloud. This is especially relevant to BT5/BT6 grass/herb-dominated
# sites.
#
# IMPORTANT
# ---------
# * Every retained 0.25-ha point remains point weighted for the LINEAR AXES.
# * BT5 and BT6 are excluded from M1-M4 endpoint construction, but are still
#   SCORED and reported as ecological diagnostics.
# * M0 retains BT5/BT6 to show how much the BT filter changes the axis.
# * Reference-cloud kNN pools are capped per parent and globally for computation
#   and to stop one very large parent from dominating the local neighbour pool.
# * All AlphaEarth observations are assumed to be from the survey/assessment
#   year, as in the point extraction. A filename/year mismatch
#
# WHITENING
# ----------
# Shrunk whitening: shrink = 0.2 by default; shrink = 1 means plain cosine
# similarity (no whitening).
#
# Here we use a transparent covariance-shrink whitening metric consistent with
# that behaviour:
#
#   Sigma_shrunk = (1-shrink) * Sigma + shrink * mean(diag(Sigma)) * I
#   W            = Sigma_shrunk^(-1/2)
#   X_white      = X %*% W
#
# At shrink=1, W is a scalar multiple of I, so cosine geometry is unchanged.
# At shrink=0, full covariance whitening is applied (with a small eigenvalue
# floor for numerical stability).
#
# The whitener is fitted WITHOUT SJ labels, by default on all complete
# training AE rows (unlabelled target-ecosystem pixels). Set
# WHITEN_FIT_SCOPE <- "axis_eligible" to exclude BT5/6 from this fit as well.
#
# SCALE
# -----
# Every linear axis is calibrated as:
#
#   median observed SJ4 reference -> 0
#   median observed SJ1 reference -> 1
#
# Values are NOT clipped. A score >1 or <0 therefore remains informative.
#
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

INPUT_DIR <- "equalArea2/GEEpointStats/"
INPUT_REGEX <- "^kystlynghei_condition_points_0p25ha_surveyYear.*\\.csv$"


RANDOM_SEED <- 42

# BT classes excluded from candidate endpoint construction.
EXCLUDE_BT_FROM_AXIS <- c(5, 6)
ALLOW_MISSING_BT_IN_AXIS <- TRUE

# Whitening. 0 = full whitening; 1 = no whitening under cosine geometry.
WHITEN_SHRINK <- 0.20
WHITEN_EIGEN_FLOOR_REL <- 1e-6

# "all_complete" fits the whitening on all complete rows (unsupervised).
# "axis_eligible" also excludes BT5/6 from the whitening fit.
WHITEN_FIT_SCOPE <- "all_complete"

# Reference-cloud settings: k = 25 neighbours, at most 2,000 references per
# endpoint.
K_NN <- 25
MAX_REF_PER_PARENT <- 60
MAX_REF_PER_ENDPOINT <- 2000
KNN_CHUNK_SIZE <- 500

# Validation settings. kNN scoring is considerably slower than a linear axis.
N_CV_FOLDS <- 5
N_CV_REPEATS <- 2
RUN_LEAVE_YEAR_OUT <- TRUE
RUN_SPATIAL_BLOCK_CV <- TRUE
SPATIAL_BLOCK_SIZE_M <- 100000
N_SPATIAL_FOLDS <- 5

# Full-data whitening sensitivity for the LINEAR axis only.
WHITEN_SHRINK_GRID <- c(1.0, 0.5, 0.2, 0.0)

# Plotting only.
PLOT_MAX_POINTS <- 40000


# =============================================================================
# 2. PACKAGES
# =============================================================================

if (!exists("REPO_ROOT")) source("R/00_setup.R")
source(file.path(REPO_ROOT, "R/ae_sj_core.R"))

required_packages(c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "tibble",
  "purrr"
))
set.seed(RANDOM_SEED)


# =============================================================================
# 3. HELPERS
# =============================================================================
#
# clean_character, clean_mi_character, code_to_numeric, get_col_or_na,
# normalize_vector and clamp_cosine now live in R/ae_sj_core.R, shared with
# model.R. row_l2_normalize, fit_endpoint_axis and score_axis (below) stay
# here unchanged - model.R has its own differently-behaved versions (no
# train/test split concept, a different near-zero-norm threshold), so they
# were kept as two separate copies rather than merged. See ae_sj_core.R's
# header for why.

row_l2_normalize <- function(m) {
  norms <- sqrt(rowSums(m^2))
  out <- matrix(
    NA_real_,
    nrow = nrow(m),
    ncol = ncol(m),
    dimnames = dimnames(m)
  )
  good <- is.finite(norms) & norms > 0
  out[good, ] <- m[good, , drop = FALSE] / norms[good]
  list(unit = out, norm = norms, good = good)
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = method))
}

auc_two_groups <- function(score_positive, score_negative) {
  x <- score_positive[is.finite(score_positive)]
  y <- score_negative[is.finite(score_negative)]
  nx <- length(x)
  ny <- length(y)
  if (nx < 2 || ny < 2) return(NA_real_)
  r <- rank(c(x, y), ties.method = "average")
  (sum(r[seq_len(nx)]) - nx * (nx + 1) / 2) / (nx * ny)
}

sample_for_plot <- function(dat, n_max = PLOT_MAX_POINTS) {
  if (nrow(dat) <= n_max) return(dat)
  dat %>% slice_sample(n = n_max)
}

mode_character <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tt <- table(x)
  names(tt)[which.max(tt)]
}

# -----------------------------------------------------------------------------
# Shrunk covariance whitener.
# -----------------------------------------------------------------------------

fit_shrunk_whitener <- function(X, shrink = WHITEN_SHRINK,
                                eigen_floor_rel = WHITEN_EIGEN_FLOOR_REL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (nrow(X) < 100) stop("Too few rows to fit the whitening covariance.")

  # cov() centres internally when estimating Sigma. We intentionally do NOT
  # subtract the mean when applying W below: at shrink=1 the transform is a
  # scalar multiple of identity and therefore exactly preserves cosine geometry.
  S <- stats::cov(X)

  if (any(!is.finite(S))) {
    stop("Non-finite covariance while fitting whitener.")
  }

  p <- ncol(S)
  target_scale <- mean(diag(S))

  if (!is.finite(target_scale) || target_scale <= 0) {
    stop("Invalid covariance scale while fitting whitener.")
  }

  S_shrunk <-
    (1 - shrink) * S +
    shrink * target_scale * diag(p)

  ee <- eigen(S_shrunk, symmetric = TRUE)
  floor_value <- max(ee$values) * eigen_floor_rel
  vals <- pmax(ee$values, floor_value)

  W <- ee$vectors %*%
    diag(1 / sqrt(vals), nrow = p) %*%
    t(ee$vectors)

  list(
    W = W,
    shrink = shrink,
    covariance = S,
    shrunk_covariance = S_shrunk,
    eigenvalues = ee$values,
    eigen_floor = floor_value,
    n_fit = nrow(X)
  )
}

apply_whitener <- function(X, whitener) {
  X %*% whitener$W
}

# -----------------------------------------------------------------------------
# Axis fitting and endpoint calibration.
# -----------------------------------------------------------------------------

fit_endpoint_axis <- function(unit_matrix, train_idx, sj, eligible,
                              min_points = 30) {
  use <- train_idx[
    eligible[train_idx] &
      is.finite(sj[train_idx]) &
      sj[train_idx] %in% c(1, 4)
  ]

  i1 <- use[sj[use] == 1]
  i4 <- use[sj[use] == 4]

  if (length(i1) < min_points || length(i4) < min_points) {
    stop(
      "Too few eligible endpoint points: SJ1=", length(i1),
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
    stop("Unexpected endpoint orientation while fitting axis.")
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

  # 0 = median SJ4, 1 = median SJ1. Unclipped by design.
  score01 <-
    (raw - fit$anchor_SJ4) /
    (fit$anchor_SJ1 - fit$anchor_SJ4)

  cos1 <- as.numeric(unit_matrix[idx, , drop = FALSE] %*% fit$prototype_SJ1)
  cos4 <- as.numeric(unit_matrix[idx, , drop = FALSE] %*% fit$prototype_SJ4)

  tibble(
    idx = idx,
    raw = raw,
    score01 = score01,
    cosine_to_SJ1_prototype = clamp_cosine(cos1),
    cosine_to_SJ4_prototype = clamp_cosine(cos4),
    prototype_cosine_margin = clamp_cosine(cos1) - clamp_cosine(cos4)
  )
}

# -----------------------------------------------------------------------------
# Reference-cloud construction.
# -----------------------------------------------------------------------------

sample_reference_cloud <- function(idx, sj, parent_id, endpoint,
                                   max_per_parent = MAX_REF_PER_PARENT,
                                   max_total = MAX_REF_PER_ENDPOINT,
                                   seed = RANDOM_SEED) {
  ii <- idx[sj[idx] == endpoint]
  if (length(ii) == 0) return(integer())

  set.seed(seed + endpoint * 1000 + length(ii))

  ref_tbl <- tibble(
    idx = ii,
    parent_id = parent_id[ii]
  ) %>%
    mutate(
      parent_id = if_else(
        is.na(parent_id) | parent_id == "",
        paste0("missing_parent__", idx),
        parent_id
      ),
      .random_order = runif(n())
    ) %>%
    group_by(parent_id) %>%
    arrange(.random_order, .by_group = TRUE) %>%
    slice_head(n = max_per_parent) %>%
    ungroup() %>%
    select(-.random_order)

  if (nrow(ref_tbl) > max_total) {
    ref_tbl <- ref_tbl %>% slice_sample(n = max_total)
  }

  ref_tbl$idx
}

# -----------------------------------------------------------------------------
# Chunked kNN cosine cloud score.
#
# For each query row, returns:
#   - mean distance to k nearest references
#   - nearest reference id
#   - nearest reference distance
#
# Self-match is removed when query_ids and ref_ids are the same sample key.
# -----------------------------------------------------------------------------

knn_cloud_distance <- function(query_unit, ref_unit,
                               query_ids = NULL, ref_ids = NULL,
                               k = K_NN,
                               chunk_size = KNN_CHUNK_SIZE) {
  nq <- nrow(query_unit)
  nr <- nrow(ref_unit)

  if (nr < k) {
    return(tibble(
      mean_knn_cosine_distance = rep(NA_real_, nq),
      nearest_cosine_distance = rep(NA_real_, nq),
      nearest_ref_id = rep(NA_character_, nq)
    ))
  }

  mean_d <- rep(NA_real_, nq)
  nearest_d <- rep(NA_real_, nq)
  nearest_id <- rep(NA_character_, nq)

  starts <- seq(1, nq, by = chunk_size)

  for (s in starts) {
    e <- min(nq, s + chunk_size - 1)
    qidx <- s:e

    sim <- query_unit[qidx, , drop = FALSE] %*% t(ref_unit)
    sim <- pmin(pmax(sim, -1), 1)

    if (!is.null(query_ids) && !is.null(ref_ids)) {
      match_col <- match(query_ids[qidx], ref_ids)
      has_self <- which(!is.na(match_col))
      if (length(has_self) > 0) {
        sim[cbind(has_self, match_col[has_self])] <- -Inf
      }
    }

    # Nearest reference.
    j_near <- max.col(sim, ties.method = "first")
    best_sim <- sim[cbind(seq_len(nrow(sim)), j_near)]
    nearest_d[qidx] <- 1 - best_sim

    if (!is.null(ref_ids)) {
      nearest_id[qidx] <- ref_ids[j_near]
    }

    # Mean of k smallest cosine distances = mean of k largest cosine sims.
    dmat <- 1 - sim
    mean_d[qidx] <- apply(
      dmat,
      1,
      function(z) {
        z <- z[is.finite(z)]
        if (length(z) < k) return(NA_real_)
        zz <- sort.int(z, partial = seq_len(k))
        mean(zz[seq_len(k)])
      }
    )

    rm(sim, dmat)
  }

  tibble(
    mean_knn_cosine_distance = mean_d,
    nearest_cosine_distance = nearest_d,
    nearest_ref_id = nearest_id
  )
}

score_reference_clouds <- function(unit_matrix, query_idx, ref1_idx, ref4_idx,
                                   sample_key, k = K_NN) {
  q <- unit_matrix[query_idx, , drop = FALSE]
  r1 <- unit_matrix[ref1_idx, , drop = FALSE]
  r4 <- unit_matrix[ref4_idx, , drop = FALSE]

  d1 <- knn_cloud_distance(
    query_unit = q,
    ref_unit = r1,
    query_ids = sample_key[query_idx],
    ref_ids = sample_key[ref1_idx],
    k = k
  )

  d4 <- knn_cloud_distance(
    query_unit = q,
    ref_unit = r4,
    query_ids = sample_key[query_idx],
    ref_ids = sample_key[ref4_idx],
    k = k
  )

  m1 <- d1$mean_knn_cosine_distance
  m4 <- d4$mean_knn_cosine_distance

  # Distance-ratio score with SJ1 as the good endpoint:
  # close to SJ1 => m1 small, m4 large => score approaches 1.
  score01 <- m4 / (m1 + m4 + 1e-12)

  tibble(
    idx = query_idx,
    d_SJ1_knn = m1,
    d_SJ4_knn = m4,
    cloud_score01 = score01,
    nearest_endpoint = if_else(m1 <= m4, "SJ1", "SJ4"),
    reference_novelty = pmin(m1, m4),
    SJ1_nearest_dist = d1$nearest_cosine_distance,
    SJ1_nearest_ref = d1$nearest_ref_id,
    SJ4_nearest_dist = d4$nearest_cosine_distance,
    SJ4_nearest_ref = d4$nearest_ref_id
  )
}

# -----------------------------------------------------------------------------
# Metric calculation.
# score higher = more SJ1-like / better along the structural succession axis.
# -----------------------------------------------------------------------------

metrics_from_score <- function(score, observed_sj, bt, model,
                               validation, repeat_id = NA_integer_,
                               fold_id = NA_integer_, heldout_year = NA_real_) {
  keep <- is.finite(score) & is.finite(observed_sj) & observed_sj %in% 1:4

  score <- score[keep]
  sj <- observed_sj[keep]
  bt2 <- bt[keep]

  if (length(score) < 10) return(tibble())

  calc_one <- function(k, eval_set) {
    if (sum(k) < 10) return(tibble())

    sc <- score[k]
    ss <- sj[k]
    goodness <- 5 - ss

    meds <- tapply(sc, ss, median, na.rm = TRUE)
    med_vec <- rep(NA_real_, 4)
    if (length(meds) > 0) {
      med_vec[as.integer(names(meds))] <- as.numeric(meds)
    }

    tibble(
      model = model,
      validation = validation,
      cv_repeat = repeat_id,
      fold = fold_id,
      heldout_year = heldout_year,
      evaluation_set = eval_set,
      n = length(sc),
      spearman_vs_goodness = safe_cor(sc, goodness, "spearman"),
      kendall_vs_goodness = safe_cor(sc, goodness, "kendall"),
      SJ1_vs_SJ4_AUC = auc_two_groups(sc[ss == 1], sc[ss == 4]),
      median_SJ1 = med_vec[1],
      median_SJ2 = med_vec[2],
      median_SJ3 = med_vec[3],
      median_SJ4 = med_vec[4],
      endpoint_separation = med_vec[1] - med_vec[4],
      medians_strictly_ordered = all(diff(med_vec) < 0, na.rm = FALSE)
    )
  }

  bind_rows(
    calc_one(rep(TRUE, length(score)), "all_BT"),
    calc_one(!(bt2 %in% c(5, 6)), "BT_not_5_6"),
    calc_one(bt2 %in% c(5, 6), "BT5_6_only")
  )
}

# -----------------------------------------------------------------------------
# Grouped fold creation.
# -----------------------------------------------------------------------------

make_grouped_folds <- function(idx, groups, sj, k, seed) {
  set.seed(seed)

  gt <- tibble(
    idx = idx,
    group = groups[idx],
    sj = as.character(sj[idx])
  ) %>%
    mutate(
      group = if_else(
        is.na(group) | group == "",
        paste0("missing_group__", idx),
        group
      )
    ) %>%
    group_by(group) %>%
    summarise(
      stratum = mode_character(sj),
      n_points = n(),
      .groups = "drop"
    )

  gt$fold <- NA_integer_

  for (st in sort(unique(gt$stratum))) {
    jj <- which(gt$stratum == st)
    jj <- sample(jj)
    gt$fold[jj] <- rep(seq_len(k), length.out = length(jj))
  }

  jj_na <- which(is.na(gt$fold))
  if (length(jj_na) > 0) {
    jj_na <- sample(jj_na)
    gt$fold[jj_na] <- rep(seq_len(k), length.out = length(jj_na))
  }

  mp <- setNames(gt$fold, gt$group)
  as.integer(mp[
    ifelse(is.na(groups[idx]) | groups[idx] == "",
           paste0("missing_group__", idx), groups[idx])
  ])
}

# -----------------------------------------------------------------------------
# Fit and score every comparison model for a given train/test split.
# -----------------------------------------------------------------------------

fit_score_split <- function(train_idx, test_idx,
                            ae_raw, ae_native_unit,
                            sj, bt, parent_id, sample_key,
                            complete_ae, axis_eligible,
                            seed = RANDOM_SEED) {

  out <- list()

  # ---------------------------------------------------------------------------
  # M0: original-style native axis including BT5/6 in endpoint construction.
  # ---------------------------------------------------------------------------

  all_bt_eligible <- complete_ae

  fit0 <- fit_endpoint_axis(
    unit_matrix = ae_native_unit,
    train_idx = train_idx,
    sj = sj,
    eligible = all_bt_eligible
  )

  out$native_axis_allBT <- score_axis(ae_native_unit, test_idx, fit0)$score01

  # ---------------------------------------------------------------------------
  # M1: native axis excluding BT5/6.
  # ---------------------------------------------------------------------------

  fit1 <- fit_endpoint_axis(
    unit_matrix = ae_native_unit,
    train_idx = train_idx,
    sj = sj,
    eligible = axis_eligible
  )

  out$native_axis_exclBT56 <- score_axis(ae_native_unit, test_idx, fit1)$score01

  # ---------------------------------------------------------------------------
  # Whitening fit (unsupervised).
  # ---------------------------------------------------------------------------

  if (WHITEN_FIT_SCOPE == "all_complete") {
    wfit_idx <- train_idx[complete_ae[train_idx]]
  } else if (WHITEN_FIT_SCOPE == "axis_eligible") {
    wfit_idx <- train_idx[complete_ae[train_idx] & axis_eligible[train_idx]]
  } else {
    stop("WHITEN_FIT_SCOPE must be 'all_complete' or 'axis_eligible'.")
  }

  whit <- fit_shrunk_whitener(
    ae_raw[wfit_idx, , drop = FALSE],
    shrink = WHITEN_SHRINK
  )

  white_train_test_idx <- sort(unique(c(train_idx, test_idx)))
  white_mat <- matrix(
    NA_real_,
    nrow = nrow(ae_raw),
    ncol = ncol(ae_raw)
  )

  white_mat[white_train_test_idx, ] <- apply_whitener(
    ae_raw[white_train_test_idx, , drop = FALSE],
    whit
  )

  white_norm <- row_l2_normalize(
    white_mat[white_train_test_idx, , drop = FALSE]
  )

  white_unit <- matrix(
    NA_real_,
    nrow = nrow(ae_raw),
    ncol = ncol(ae_raw)
  )
  white_unit[white_train_test_idx, ] <- white_norm$unit

  # ---------------------------------------------------------------------------
  # M2: whitened linear axis excluding BT5/6.
  # ---------------------------------------------------------------------------

  fit2 <- fit_endpoint_axis(
    unit_matrix = white_unit,
    train_idx = train_idx,
    sj = sj,
    eligible = axis_eligible
  )

  out$white_axis_exclBT56 <- score_axis(white_unit, test_idx, fit2)$score01

  # ---------------------------------------------------------------------------
  # Reference clouds for M3/M4.
  # Endpoint refs are eligible training observations only.
  # ---------------------------------------------------------------------------

  ref_base <- train_idx[
    axis_eligible[train_idx] &
      complete_ae[train_idx] &
      is.finite(sj[train_idx]) &
      sj[train_idx] %in% c(1, 4)
  ]

  ref1 <- sample_reference_cloud(
    ref_base, sj, parent_id, endpoint = 1,
    seed = seed + 101
  )

  ref4 <- sample_reference_cloud(
    ref_base, sj, parent_id, endpoint = 4,
    seed = seed + 404
  )

  # ---------------------------------------------------------------------------
  # M3: native cloud distance ratio.
  # ---------------------------------------------------------------------------

  cloud_native <- score_reference_clouds(
    unit_matrix = ae_native_unit,
    query_idx = test_idx,
    ref1_idx = ref1,
    ref4_idx = ref4,
    sample_key = sample_key,
    k = K_NN
  )

  out$native_cloud_exclBT56 <- cloud_native$cloud_score01

  # ---------------------------------------------------------------------------
  # M4: whitened cloud distance ratio.
  # ---------------------------------------------------------------------------

  cloud_white <- score_reference_clouds(
    unit_matrix = white_unit,
    query_idx = test_idx,
    ref1_idx = ref1,
    ref4_idx = ref4,
    sample_key = sample_key,
    k = K_NN
  )

  out$white_cloud_exclBT56 <- cloud_white$cloud_score01

  out
}


# =============================================================================
# 4. INPUT FILES
# =============================================================================

input_dir <- normalizePath(INPUT_DIR, mustWork = TRUE)

input_files <- list.files(
  input_dir,
  pattern = INPUT_REGEX,
  full.names = TRUE,
  ignore.case = TRUE
)

input_files <- sort(input_files)

if (length(input_files) == 0) {
  stop("No matching input CSV files found in: ", input_dir)
}

# Outputs go to the repo's outputs/ folder.
output_dir <- file.path(OUTPUTS_ROOT, "validation")
qc_dir <- file.path(output_dir, "01_QC")
full_dir <- file.path(output_dir, "02_full_data_scores")
cv_dir <- file.path(output_dir, "03_grouped_CV")
year_dir <- file.path(output_dir, "04_leave_year_out")
spatial_dir <- file.path(output_dir, "05_spatial_block_CV")
sensitivity_dir <- file.path(output_dir, "06_whitening_sensitivity")

for (d in c(output_dir, qc_dir, full_dir, cv_dir, year_dir,
            spatial_dir, sensitivity_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("\nMatched ", length(input_files), " input files.\n", sep = "")


# =============================================================================
# 5. READ AND STANDARDISE DATA
# =============================================================================

# read_one_csv() is defined in R/ae_sj_core.R.
raw <- bind_rows(lapply(input_files, read_one_csv))

required_fields <- c("parentID", "_survey_year_num", "7RA-SJ")
missing_required <- setdiff(required_fields, names(raw))
if (length(missing_required) > 0) {
  stop("Missing required fields: ", paste(missing_required, collapse = ", "))
}

raw_year <- code_to_numeric(raw[["_survey_year_num"]])
year_std <- ifelse(is.finite(raw_year), raw_year, raw$year_from_filename)

parent_id <- clean_character(raw[["parentID"]])
sj <- code_to_numeric(raw[["7RA-SJ"]])
bt <- code_to_numeric(get_col_or_na(raw, "7JB-BT"))

survey_parent_key <- paste(year_std, parent_id, sep = "__")
cv_group <- parent_id
missing_group <- is.na(cv_group) | cv_group == ""
cv_group[missing_group] <- survey_parent_key[missing_group]

x_25833 <- code_to_numeric(get_col_or_na(raw, "sample_center_x_25833"))
y_25833 <- code_to_numeric(get_col_or_na(raw, "sample_center_y_25833"))
sample_id <- clean_character(get_col_or_na(raw, "sample_id"))
system_index <- clean_character(get_col_or_na(raw, "system:index"))

coord_ok <- is.finite(x_25833) & is.finite(y_25833)
sample_key <- rep(NA_character_, nrow(raw))

sample_key[coord_ok] <- paste(
  year_std[coord_ok], parent_id[coord_ok],
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
    cv_group = cv_group,
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

# Year-of-assessment QC.
year_qc <- dat %>%
  transmute(
    source_file,
    sample_key,
    year_std,
    year_from_filename,
    mismatch = is.finite(year_std) &
      is.finite(year_from_filename) &
      year_std != year_from_filename
  )

write_csv(year_qc, file.path(qc_dir, "01_assessment_year_QC.csv"))

cat("\nAssessment-year mismatches: ", sum(year_qc$mismatch, na.rm = TRUE), "\n", sep = "")


# =============================================================================
# 6. ALPHAEARTH MATRIX
# =============================================================================

ae_columns <- grep("^AE_A[0-9]{2}_mean$", names(dat), value = TRUE)
ae_columns <- sort(ae_columns)

if (length(ae_columns) != 64) {
  stop("Expected exactly 64 AlphaEarth mean bands; found ", length(ae_columns))
}

ae_raw <- sapply(ae_columns, function(v) code_to_numeric(dat[[v]]))
ae_raw <- as.matrix(ae_raw)
colnames(ae_raw) <- ae_columns

complete_ae <- apply(ae_raw, 1, function(z) all(is.finite(z)))

native_norm <- row_l2_normalize(ae_raw)
ae_native_unit <- native_norm$unit
dat$AE_vector_norm <- native_norm$norm
dat$AE_complete <- complete_ae & native_norm$good
complete_ae <- dat$AE_complete

# Candidate axis eligibility. Missing BT is retained by default.
axis_eligible <-
  complete_ae &
  is.finite(dat$sj) &
  dat$sj %in% 1:4

if (ALLOW_MISSING_BT_IN_AXIS) {
  axis_eligible <- axis_eligible &
    (is.na(dat$bt) | !(dat$bt %in% EXCLUDE_BT_FROM_AXIS))
} else {
  axis_eligible <- axis_eligible &
    is.finite(dat$bt) &
    !(dat$bt %in% EXCLUDE_BT_FROM_AXIS)
}

dat$axis_eligible_excl_BT56 <- axis_eligible
dat$tilstand_clean <- clean_character(get_col_or_na(dat, "tilstand"))
dat$local_condition_clean <- clean_character(get_col_or_na(dat, "local_condition"))

support <- dat %>%
  filter(AE_complete, is.finite(sj), sj %in% 1:4) %>%
  mutate(
    BT_group = case_when(
      bt %in% c(5, 6) ~ "BT5_6",
      is.na(bt) ~ "BT_missing",
      TRUE ~ "BT_other"
    )
  ) %>%
  count(sj, BT_group, name = "n_points")

write_csv(support, file.path(qc_dir, "02_support_by_SJ_BTgroup.csv"))
print(support, n = Inf)

valid_sj_idx <- which(complete_ae & is.finite(dat$sj) & dat$sj %in% 1:4)


# =============================================================================
# 7. FULL-DATA M0/M1 AXES: EFFECT OF EXCLUDING BT5/BT6
# =============================================================================

fit_native_all <- fit_endpoint_axis(
  ae_native_unit,
  train_idx = valid_sj_idx,
  sj = dat$sj,
  eligible = complete_ae
)

fit_native_excl <- fit_endpoint_axis(
  ae_native_unit,
  train_idx = valid_sj_idx,
  sj = dat$sj,
  eligible = axis_eligible
)

axis_cos <- clamp_cosine(sum(
  fit_native_all$direction * fit_native_excl$direction
))

axis_filter_comparison <- tibble(
  comparison = "native_allBT_vs_native_excluding_BT5_6",
  cosine_similarity = axis_cos,
  angular_difference_deg = acos(axis_cos) * 180 / pi,
  allBT_n_SJ1 = fit_native_all$n_SJ1,
  allBT_n_SJ4 = fit_native_all$n_SJ4,
  exclBT56_n_SJ1 = fit_native_excl$n_SJ1,
  exclBT56_n_SJ4 = fit_native_excl$n_SJ4
)

write_csv(
  axis_filter_comparison,
  file.path(qc_dir, "03_axis_change_after_excluding_BT5_6.csv")
)

cat("\nEffect of excluding BT5/6 from native endpoint construction:\n")
print(axis_filter_comparison, width = Inf)


# =============================================================================
# 8. FULL-DATA WHITENER AND WHITENED AXIS
# =============================================================================

if (WHITEN_FIT_SCOPE == "all_complete") {
  whiten_fit_idx <- which(complete_ae)
} else if (WHITEN_FIT_SCOPE == "axis_eligible") {
  whiten_fit_idx <- which(axis_eligible)
} else {
  stop("WHITEN_FIT_SCOPE must be 'all_complete' or 'axis_eligible'.")
}

full_whitener <- fit_shrunk_whitener(
  ae_raw[whiten_fit_idx, , drop = FALSE],
  shrink = WHITEN_SHRINK
)

ae_white_raw <- apply_whitener(ae_raw, full_whitener)
white_norm <- row_l2_normalize(ae_white_raw)
ae_white_unit <- white_norm$unit

fit_white_excl <- fit_endpoint_axis(
  ae_white_unit,
  train_idx = valid_sj_idx,
  sj = dat$sj,
  eligible = axis_eligible
)

whitener_summary <- tibble(
  shrink = WHITEN_SHRINK,
  fit_scope = WHITEN_FIT_SCOPE,
  n_fit = full_whitener$n_fit,
  raw_cov_condition_number =
    max(full_whitener$eigenvalues) /
    max(min(full_whitener$eigenvalues), .Machine$double.eps),
  eigen_floor = full_whitener$eigen_floor
)

write_csv(whitener_summary, file.path(qc_dir, "04_whitener_summary.csv"))


# =============================================================================
# 9. FULL-DATA REFERENCE CLOUDS
# =============================================================================

ref_base_full <- valid_sj_idx[axis_eligible[valid_sj_idx]]

ref1_full <- sample_reference_cloud(
  ref_base_full, dat$sj, dat$parent_id,
  endpoint = 1, seed = RANDOM_SEED + 101
)

ref4_full <- sample_reference_cloud(
  ref_base_full, dat$sj, dat$parent_id,
  endpoint = 4, seed = RANDOM_SEED + 404
)

ref_support <- tibble(
  endpoint = c("SJ1", "SJ4"),
  n_reference_points = c(length(ref1_full), length(ref4_full)),
  n_reference_parents = c(
    n_distinct(dat$parent_id[ref1_full]),
    n_distinct(dat$parent_id[ref4_full])
  )
)

write_csv(ref_support, file.path(qc_dir, "05_reference_cloud_support.csv"))
print(ref_support)


# =============================================================================
# 10. SCORE ALL VALID OBSERVATIONS WITH ALL FULL-DATA METHODS
# =============================================================================

score0 <- score_axis(ae_native_unit, valid_sj_idx, fit_native_all)
score1 <- score_axis(ae_native_unit, valid_sj_idx, fit_native_excl)
score2 <- score_axis(ae_white_unit, valid_sj_idx, fit_white_excl)

cat("\nScoring native reference clouds...\n")
cloud_native_full <- score_reference_clouds(
  ae_native_unit,
  query_idx = valid_sj_idx,
  ref1_idx = ref1_full,
  ref4_idx = ref4_full,
  sample_key = dat$sample_key,
  k = K_NN
)

cat("\nScoring whitened reference clouds...\n")
cloud_white_full <- score_reference_clouds(
  ae_white_unit,
  query_idx = valid_sj_idx,
  ref1_idx = ref1_full,
  ref4_idx = ref4_full,
  sample_key = dat$sample_key,
  k = K_NN
)

full_scores <- dat[valid_sj_idx, , drop = FALSE] %>%
  transmute(
    source_file,
    year_std,
    parent_id,
    survey_parent_key,
    cv_group,
    sample_key,
    x_25833,
    y_25833,
    sj,
    bt,
    tilstand = tilstand_clean,
    local_condition = local_condition_clean,
    axis_eligible_excl_BT56,
    AE_vector_norm
  ) %>%
  mutate(
    native_axis01_allBT = score0$score01,
    native_axis01_exclBT56 = score1$score01,
    white_axis01_exclBT56 = score2$score01,

    native_cos_SJ1_prototype = score1$cosine_to_SJ1_prototype,
    native_cos_SJ4_prototype = score1$cosine_to_SJ4_prototype,
    native_prototype_margin = score1$prototype_cosine_margin,

    white_cos_SJ1_prototype = score2$cosine_to_SJ1_prototype,
    white_cos_SJ4_prototype = score2$cosine_to_SJ4_prototype,
    white_prototype_margin = score2$prototype_cosine_margin,

    native_d_SJ1_knn = cloud_native_full$d_SJ1_knn,
    native_d_SJ4_knn = cloud_native_full$d_SJ4_knn,
    native_cloud01 = cloud_native_full$cloud_score01,
    native_nearest_endpoint = cloud_native_full$nearest_endpoint,
    native_reference_novelty = cloud_native_full$reference_novelty,
    native_SJ1_nearest_dist = cloud_native_full$SJ1_nearest_dist,
    native_SJ1_nearest_ref = cloud_native_full$SJ1_nearest_ref,
    native_SJ4_nearest_dist = cloud_native_full$SJ4_nearest_dist,
    native_SJ4_nearest_ref = cloud_native_full$SJ4_nearest_ref,

    white_d_SJ1_knn = cloud_white_full$d_SJ1_knn,
    white_d_SJ4_knn = cloud_white_full$d_SJ4_knn,
    white_cloud01 = cloud_white_full$cloud_score01,
    white_nearest_endpoint = cloud_white_full$nearest_endpoint,
    white_reference_novelty = cloud_white_full$reference_novelty,
    white_SJ1_nearest_dist = cloud_white_full$SJ1_nearest_dist,
    white_SJ1_nearest_ref = cloud_white_full$SJ1_nearest_ref,
    white_SJ4_nearest_dist = cloud_white_full$SJ4_nearest_dist,
    white_SJ4_nearest_ref = cloud_white_full$SJ4_nearest_ref
  )

write_csv(
  full_scores,
  file.path(full_dir, "01_full_data_AE_SJ_method_scores.csv")
)


# =============================================================================
# 11. REFERENCE-NOVELTY THRESHOLDS
# =============================================================================
#
# Endpoint references are scored against their own reference clouds with the
# exact self-match removed. The 95th percentile of within-endpoint novelty is
# used as a diagnostic threshold for "unlike either reference cloud".
# This is NOT a condition threshold; it is an out-of-reference diagnostic.
# =============================================================================

ref_rows <- which(full_scores$sample_key %in% dat$sample_key[c(ref1_full, ref4_full)])

native_ref_novelty <- full_scores$native_reference_novelty[ref_rows]
white_ref_novelty <- full_scores$white_reference_novelty[ref_rows]

native_novelty_q95 <- quantile(native_ref_novelty, 0.95, na.rm = TRUE)
white_novelty_q95 <- quantile(white_ref_novelty, 0.95, na.rm = TRUE)

full_scores <- full_scores %>%
  mutate(
    native_outside_endpoint_clouds =
      native_reference_novelty > native_novelty_q95,
    white_outside_endpoint_clouds =
      white_reference_novelty > white_novelty_q95
  )

write_csv(
  full_scores,
  file.path(full_dir, "01_full_data_AE_SJ_method_scores.csv")
)

novelty_thresholds <- tibble(
  space = c("native", "whitened"),
  reference_novelty_q95 = c(native_novelty_q95, white_novelty_q95)
)

write_csv(
  novelty_thresholds,
  file.path(full_dir, "02_reference_novelty_thresholds.csv")
)


# =============================================================================
# 12. FULL-DATA DESCRIPTIVE COMPARISON
# =============================================================================

full_metrics <- bind_rows(
  metrics_from_score(
    full_scores$native_axis01_allBT,
    full_scores$sj, full_scores$bt,
    "M0_native_axis_allBT", "full_data"
  ),
  metrics_from_score(
    full_scores$native_axis01_exclBT56,
    full_scores$sj, full_scores$bt,
    "M1_native_axis_exclBT56", "full_data"
  ),
  metrics_from_score(
    full_scores$white_axis01_exclBT56,
    full_scores$sj, full_scores$bt,
    "M2_white_axis_exclBT56", "full_data"
  ),
  metrics_from_score(
    full_scores$native_cloud01,
    full_scores$sj, full_scores$bt,
    "M3_native_cloud_exclBT56", "full_data"
  ),
  metrics_from_score(
    full_scores$white_cloud01,
    full_scores$sj, full_scores$bt,
    "M4_white_cloud_exclBT56", "full_data"
  )
)

write_csv(full_metrics, file.path(full_dir, "03_full_data_method_metrics.csv"))
print(full_metrics, n = Inf, width = Inf)

method_correlations <- full_scores %>%
  select(
    native_axis01_allBT,
    native_axis01_exclBT56,
    white_axis01_exclBT56,
    native_cloud01,
    white_cloud01
  ) %>%
  cor(use = "pairwise.complete.obs", method = "spearman") %>%
  as.data.frame() %>%
  rownames_to_column("method")

write_csv(
  method_correlations,
  file.path(full_dir, "04_method_score_Spearman_correlations.csv")
)


# =============================================================================
# 13. BT5/BT6 DIAGNOSTICS
# =============================================================================

bt56_summary <- full_scores %>%
  filter(bt %in% c(5, 6)) %>%
  group_by(bt, sj) %>%
  summarise(
    n = n(),
    median_native_axis01 = median(native_axis01_exclBT56, na.rm = TRUE),
    median_white_axis01 = median(white_axis01_exclBT56, na.rm = TRUE),
    median_native_cloud01 = median(native_cloud01, na.rm = TRUE),
    median_white_cloud01 = median(white_cloud01, na.rm = TRUE),
    median_native_novelty = median(native_reference_novelty, na.rm = TRUE),
    median_white_novelty = median(white_reference_novelty, na.rm = TRUE),
    prop_native_outside = mean(native_outside_endpoint_clouds, na.rm = TRUE),
    prop_white_outside = mean(white_outside_endpoint_clouds, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(bt56_summary, file.path(full_dir, "05_BT5_BT6_diagnostics.csv"))
print(bt56_summary, n = Inf)

# Key diagnostic: axis says "SJ1-like" but cloud distance says "unusual".
p_bt56 <- full_scores %>%
  filter(bt %in% c(5, 6)) %>%
  ggplot(
    aes(
      x = white_axis01_exclBT56,
      y = white_reference_novelty,
      shape = factor(bt),
      colour = factor(sj)
    )
  ) +
  geom_point(alpha = 0.55, size = 1.7) +
  geom_hline(yintercept = white_novelty_q95, linetype = "dashed") +
  geom_vline(xintercept = c(0, 1), linetype = "dotted") +
  theme_bw(base_size = 12) +
  labs(
    title = "BT5/BT6: projected succession position versus reference-cloud novelty",
    subtitle = paste0(
      "Whitened AE; score 0 = SJ4 median, 1 = SJ1 median. ",
      "Dashed line = 95th percentile of endpoint-reference novelty"
    ),
    x = "Whitened AE-SJ axis score",
    y = "Distance to nearest endpoint reference cloud",
    shape = "BT",
    colour = "Observed SJ"
  )

print(p_bt56)
ggsave(
  file.path(full_dir, "06_BT5_BT6_axis_vs_reference_novelty.png"),
  p_bt56, width = 9, height = 6.5, dpi = 300
)


# =============================================================================
# 14. FULL-DATA SCORE DISTRIBUTION FIGURE
# =============================================================================

plot_long <- full_scores %>%
  select(
    sj, bt,
    native_axis01_allBT,
    native_axis01_exclBT56,
    white_axis01_exclBT56,
    native_cloud01,
    white_cloud01
  ) %>%
  pivot_longer(
    cols = -c(sj, bt),
    names_to = "method",
    values_to = "score01"
  ) %>%
  filter(is.finite(score01)) %>%
  sample_for_plot()

p_methods <- ggplot(
  plot_long,
  aes(x = factor(sj), y = score01)
) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.025, size = 0.35) +
  geom_hline(yintercept = c(0, 1), linetype = "dotted", linewidth = 0.35) +
  facet_wrap(~ method, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Comparison of AlphaEarth SJ scoring methods",
    subtitle = "Linear axes are endpoint-calibrated: SJ4 median = 0, SJ1 median = 1",
    x = "Observed 7RA-SJ",
    y = "SJ1-likeness score"
  )

print(p_methods)
ggsave(
  file.path(full_dir, "07_method_score_distributions_by_SJ.png"),
  p_methods, width = 12, height = 8, dpi = 300
)


# =============================================================================
# 15. GROUPED CROSS-VALIDATION METHOD COMPARISON
# =============================================================================

cv_metrics_list <- list()
cv_counter <- 0L

for (r in seq_len(N_CV_REPEATS)) {
  cat("\nGrouped CV repeat ", r, " / ", N_CV_REPEATS, "\n", sep = "")

  foldid <- make_grouped_folds(
    valid_sj_idx,
    dat$cv_group,
    dat$sj,
    k = N_CV_FOLDS,
    seed = RANDOM_SEED + r
  )

  for (f in seq_len(N_CV_FOLDS)) {
    cat("  fold ", f, " / ", N_CV_FOLDS, "\n", sep = "")

    train_idx <- valid_sj_idx[foldid != f]
    test_idx <- valid_sj_idx[foldid == f]

    scores <- fit_score_split(
      train_idx = train_idx,
      test_idx = test_idx,
      ae_raw = ae_raw,
      ae_native_unit = ae_native_unit,
      sj = dat$sj,
      bt = dat$bt,
      parent_id = dat$parent_id,
      sample_key = dat$sample_key,
      complete_ae = complete_ae,
      axis_eligible = axis_eligible,
      seed = RANDOM_SEED + r * 100 + f
    )

    for (nm in names(scores)) {
      cv_counter <- cv_counter + 1L
      cv_metrics_list[[cv_counter]] <- metrics_from_score(
        scores[[nm]],
        dat$sj[test_idx],
        dat$bt[test_idx],
        model = nm,
        validation = "parent_grouped_CV",
        repeat_id = r,
        fold_id = f
      )
    }
  }
}

cv_metrics <- bind_rows(cv_metrics_list)

write_csv(
  cv_metrics,
  file.path(cv_dir, "01_grouped_CV_method_metrics.csv")
)

cv_summary <- cv_metrics %>%
  filter(evaluation_set %in% c("all_BT", "BT_not_5_6")) %>%
  group_by(model, evaluation_set) %>%
  summarise(
    n_folds = n(),
    mean_spearman = mean(spearman_vs_goodness, na.rm = TRUE),
    sd_spearman = sd(spearman_vs_goodness, na.rm = TRUE),
    mean_kendall = mean(kendall_vs_goodness, na.rm = TRUE),
    mean_AUC = mean(SJ1_vs_SJ4_AUC, na.rm = TRUE),
    sd_AUC = sd(SJ1_vs_SJ4_AUC, na.rm = TRUE),
    mean_endpoint_separation = mean(endpoint_separation, na.rm = TRUE),
    prop_ordered = mean(medians_strictly_ordered, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cv_summary,
  file.path(cv_dir, "02_grouped_CV_method_summary.csv")
)

cat("\nGrouped-CV method summary:\n")
print(cv_summary, n = Inf, width = Inf)

p_cv <- cv_metrics %>%
  filter(evaluation_set == "BT_not_5_6") %>%
  ggplot(aes(x = model, y = spearman_vs_goodness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 1.2) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Does whitening or reference-cloud scoring improve held-out SJ ordering?",
    subtitle = "Parent-grouped CV; evaluation shown for BT1-4",
    x = NULL,
    y = "Held-out Spearman: score vs SJ goodness"
  )

print(p_cv)
ggsave(
  file.path(cv_dir, "03_grouped_CV_Spearman_method_comparison.png"),
  p_cv, width = 9, height = 6.5, dpi = 300
)


# =============================================================================
# 16. LEAVE-ONE-YEAR-OUT COMPARISON
# =============================================================================

if (RUN_LEAVE_YEAR_OUT) {
  years <- sort(unique(dat$year_std[valid_sj_idx]))
  years <- years[is.finite(years)]

  year_metrics_list <- list()
  counter <- 0L

  for (yy in years) {
    cat("\nLeave-year-out: ", yy, "\n", sep = "")

    train_idx <- valid_sj_idx[dat$year_std[valid_sj_idx] != yy]
    test_idx <- valid_sj_idx[dat$year_std[valid_sj_idx] == yy]

    if (length(train_idx) == 0 || length(test_idx) == 0) next

    scores <- fit_score_split(
      train_idx = train_idx,
      test_idx = test_idx,
      ae_raw = ae_raw,
      ae_native_unit = ae_native_unit,
      sj = dat$sj,
      bt = dat$bt,
      parent_id = dat$parent_id,
      sample_key = dat$sample_key,
      complete_ae = complete_ae,
      axis_eligible = axis_eligible,
      seed = RANDOM_SEED + yy
    )

    for (nm in names(scores)) {
      counter <- counter + 1L
      year_metrics_list[[counter]] <- metrics_from_score(
        scores[[nm]],
        dat$sj[test_idx],
        dat$bt[test_idx],
        model = nm,
        validation = "leave_one_year_out",
        heldout_year = yy
      )
    }
  }

  year_metrics <- bind_rows(year_metrics_list)

  write_csv(
    year_metrics,
    file.path(year_dir, "01_leave_year_out_method_metrics.csv")
  )

  p_year <- year_metrics %>%
    filter(evaluation_set == "BT_not_5_6") %>%
    ggplot(
      aes(
        x = heldout_year,
        y = spearman_vs_goodness,
        colour = model,
        group = model
      )
    ) +
    geom_line() +
    geom_point(size = 1.8) +
    theme_bw(base_size = 11) +
    labs(
      title = "Temporal transfer of alternative AE-SJ scoring methods",
      subtitle = "Each field assessment is paired with AlphaEarth from its assessment year",
      x = "Held-out assessment year",
      y = "Spearman: score vs SJ goodness",
      colour = "Method"
    )

  print(p_year)
  ggsave(
    file.path(year_dir, "02_leave_year_out_method_comparison.png"),
    p_year, width = 10, height = 6.5, dpi = 300
  )
}


# =============================================================================
# 17. SPATIAL-BLOCK CV COMPARISON
# =============================================================================

if (
  RUN_SPATIAL_BLOCK_CV &&
  any(is.finite(dat$x_25833[valid_sj_idx])) &&
  any(is.finite(dat$y_25833[valid_sj_idx]))
) {
  parent_xy <- dat %>%
    slice(valid_sj_idx) %>%
    filter(is.finite(x_25833), is.finite(y_25833)) %>%
    group_by(cv_group) %>%
    summarise(
      x = median(x_25833),
      y = median(y_25833),
      .groups = "drop"
    ) %>%
    mutate(
      block_x = floor(x / SPATIAL_BLOCK_SIZE_M),
      block_y = floor(y / SPATIAL_BLOCK_SIZE_M),
      spatial_block = paste(block_x, block_y, sep = "__")
    )

  block_map <- setNames(parent_xy$spatial_block, parent_xy$cv_group)
  row_block <- unname(block_map[dat$cv_group])

  spatial_idx <- valid_sj_idx[!is.na(row_block[valid_sj_idx])]
  unique_blocks <- unique(row_block[spatial_idx])

  if (length(unique_blocks) >= N_SPATIAL_FOLDS) {
    set.seed(RANDOM_SEED + 50000)
    unique_blocks <- sample(unique_blocks)
    block_fold <- setNames(
      rep(seq_len(N_SPATIAL_FOLDS), length.out = length(unique_blocks)),
      unique_blocks
    )
    foldid <- unname(block_fold[row_block[spatial_idx]])

    spatial_metrics_list <- list()
    counter <- 0L

    for (f in seq_len(N_SPATIAL_FOLDS)) {
      cat("\nSpatial CV fold ", f, " / ", N_SPATIAL_FOLDS, "\n", sep = "")

      train_idx <- spatial_idx[foldid != f]
      test_idx <- spatial_idx[foldid == f]

      train_classes <- sort(unique(dat$sj[train_idx]))
      if (!all(c(1, 4) %in% train_classes)) next

      scores <- fit_score_split(
        train_idx = train_idx,
        test_idx = test_idx,
        ae_raw = ae_raw,
        ae_native_unit = ae_native_unit,
        sj = dat$sj,
        bt = dat$bt,
        parent_id = dat$parent_id,
        sample_key = dat$sample_key,
        complete_ae = complete_ae,
        axis_eligible = axis_eligible,
        seed = RANDOM_SEED + 50000 + f
      )

      for (nm in names(scores)) {
        counter <- counter + 1L
        spatial_metrics_list[[counter]] <- metrics_from_score(
          scores[[nm]],
          dat$sj[test_idx],
          dat$bt[test_idx],
          model = nm,
          validation = "spatial_block_CV",
          fold_id = f
        )
      }
    }

    spatial_metrics <- bind_rows(spatial_metrics_list)

    write_csv(
      parent_xy,
      file.path(spatial_dir, "01_parent_spatial_blocks.csv")
    )

    write_csv(
      spatial_metrics,
      file.path(spatial_dir, "02_spatial_block_method_metrics.csv")
    )

    spatial_summary <- spatial_metrics %>%
      filter(evaluation_set %in% c("all_BT", "BT_not_5_6")) %>%
      group_by(model, evaluation_set) %>%
      summarise(
        n_folds = n(),
        mean_spearman = mean(spearman_vs_goodness, na.rm = TRUE),
        mean_AUC = mean(SJ1_vs_SJ4_AUC, na.rm = TRUE),
        mean_endpoint_separation = mean(endpoint_separation, na.rm = TRUE),
        prop_ordered = mean(medians_strictly_ordered, na.rm = TRUE),
        .groups = "drop"
      )

    write_csv(
      spatial_summary,
      file.path(spatial_dir, "03_spatial_block_method_summary.csv")
    )

    cat("\nSpatial-CV method summary:\n")
    print(spatial_summary, n = Inf, width = Inf)
  } else {
    warning("Too few spatial blocks for spatial CV.")
  }
}


# =============================================================================
# 18. WHITENING-SHRINK SENSITIVITY — FULL-DATA LINEAR AXIS
# =============================================================================
#
# This is intentionally cheap: it compares the linear candidate axis at several
# shrink levels without rerunning kNN or CV. Use it to see whether the result is
# specific to shrink=0.2.
# =============================================================================

shrink_results <- vector("list", length(WHITEN_SHRINK_GRID))

for (ii in seq_along(WHITEN_SHRINK_GRID)) {
  sh <- WHITEN_SHRINK_GRID[ii]
  cat("\nWhitening sensitivity shrink = ", sh, "\n", sep = "")

  ww <- fit_shrunk_whitener(
    ae_raw[whiten_fit_idx, , drop = FALSE],
    shrink = sh
  )

  xw <- apply_whitener(ae_raw, ww)
  xwu <- row_l2_normalize(xw)$unit

  ff <- fit_endpoint_axis(
    xwu,
    train_idx = valid_sj_idx,
    sj = dat$sj,
    eligible = axis_eligible
  )

  ss <- score_axis(xwu, valid_sj_idx, ff)$score01

  mm <- metrics_from_score(
    ss,
    dat$sj[valid_sj_idx],
    dat$bt[valid_sj_idx],
    model = paste0("whiten_shrink_", sh),
    validation = "full_data_shrink_sensitivity"
  ) %>%
    mutate(shrink = sh, .before = 1)

  shrink_results[[ii]] <- mm
}

shrink_results <- bind_rows(shrink_results)

write_csv(
  shrink_results,
  file.path(sensitivity_dir, "01_whitening_shrink_sensitivity.csv")
)


# =============================================================================
# 19. SIMPLE METHOD-RANKING TABLE
# =============================================================================
#
# This does not declare a winner automatically. It puts the key held-out metrics
# side by side. Prefer improvements that are repeated in parent-grouped CV,
# leave-year-out and spatial-block CV rather than a full-data gain alone.
# =============================================================================

ranking <- cv_summary %>%
  filter(evaluation_set == "BT_not_5_6") %>%
  arrange(desc(mean_spearman), desc(mean_AUC))

write_csv(
  ranking,
  file.path(output_dir, "METHOD_COMPARISON_RANKING_groupedCV_BT1_4.csv")
)


# =============================================================================
# 20. FINISHED
# =============================================================================

cat("\n\n============================================================\n")
cat("AE-SJ METHOD COMPARISON COMPLETE\n")
cat("============================================================\n\n")

cat("Outputs written to:\n", output_dir, "\n\n", sep = "")

cat("Most useful files to inspect first:\n")
cat("  1. 03_grouped_CV/02_grouped_CV_method_summary.csv\n")
cat("  2. 05_spatial_block_CV/03_spatial_block_method_summary.csv\n")
cat("  3. 02_full_data_scores/05_BT5_BT6_diagnostics.csv\n")
cat("  4. 02_full_data_scores/06_BT5_BT6_axis_vs_reference_novelty.png\n")
cat("  5. 01_QC/03_axis_change_after_excluding_BT5_6.csv\n")
cat("  6. 06_whitening_sensitivity/01_whitening_shrink_sensitivity.csv\n\n")

cat("Interpretation reminder:\n")
cat("  Linear axis score: 0 = median SJ4, 1 = median SJ1; values are unclipped.\n")
cat("  Cloud score: close to SJ4 -> 0; close to SJ1 -> 1.\n")
cat("  High reference_novelty = unlike BOTH endpoint reference clouds.\n")
