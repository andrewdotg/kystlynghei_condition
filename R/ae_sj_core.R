# =============================================================================
# ae_sj_core.R - shared AE-SJ data-cleaning helpers
# =============================================================================
#
# model.R (the deployed model) and validation.R (the method-comparison /
# cross-validation script that selected it) both started from the same
# analysis boilerplate and had drifted into two near-identical copies of a
# handful of small cleaning helpers. This file is the single copy.
#
# Helpers shared by model.R and validation.R. Functions with the same name in
# both scripts (fit_endpoint_axis, score_axis, row_l2_normalize) differ in
# behaviour and are kept separately in each.
# =============================================================================

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "<NA>")] <- NA_character_
  x
}

clean_mi_character <- function(x) {
  x <- clean_character(x)
  x[toupper(x) == "X"] <- NA_character_
  x
}

code_to_numeric <- function(x) {
  suppressWarnings(readr::parse_number(clean_mi_character(x)))
}

get_col_or_na <- function(dat, nm) {
  if (nm %in% names(dat)) return(dat[[nm]])
  rep(NA_character_, nrow(dat))
}

normalize_vector <- function(v) {
  n <- sqrt(sum(v^2))
  if (!is.finite(n) || n <= 0) return(rep(NA_real_, length(v)))
  v / n
}

clamp_cosine <- function(x) pmin(pmax(x, -1), 1)

# Reads one raw calibration-point CSV, tagging it with the source filename and
# the year parsed from that filename (used as a fallback/cross-check against
# the in-data survey-year field).
read_one_csv <- function(path) {
  x <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN", "null", "NULL"),
    show_col_types = FALSE,
    progress = FALSE
  )

  filename <- basename(path)
  year_match <- stringr::str_match(filename, "(20[0-9]{2})(?=\\.csv$)")[, 2]

  x$source_file <- filename
  x$year_from_filename <- suppressWarnings(as.integer(year_match))
  x
}
