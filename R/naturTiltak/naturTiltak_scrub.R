# =============================================================================
# NATURTILTAK - SCRUB / TREE CLEARANCE VALIDATION ANALYSIS
# =============================================================================
#
# Tests the AE-SJ clearance response against scrub and tree clearance
# registered in naturtiltak (Methodology Section 11.3, Table 6). The MI-wide
# clearance analysis is the separate R/modules/scrub_clearance.R.
#
# PURPOSE
# -------
# Use the tiltak clearance dataset as a VALIDATION dataset for the AE-SJ method:
#
#   1. Join the exported AE-SJ clearance-response CSV to the tiltak DBF.
#   2. Exclude records that explicitly state:
#        "The report states that the measure was not completed as planned."
#   3. Standardise the many tiltakType strings into ecologically useful
#      clearance groups:
#        - Scrub clearance
#        - Tree clearance
#        - Scrub + tree clearance
#        - Unspecified scrub/tree clearance
#   4. Report how many eligible sites contain AE-SJ clearance-response pixels.
#   5. Test whether the PRE-treatment SJ-like source zones correspond to the
#      reported management:
#        scrub clearance -> expected signal from SJ3-like pixels
#        tree clearance  -> expected signal from SJ4-like pixels
#        both            -> expected signal from both SJ3-like and SJ4-like
#   6. Summarise the full response-pixel SJ transition matrix by clearance type.
#   7. Produce the pre/post SJ-composition figures, faceted by
#      standardised clearance type.
#
# IMPORTANT INTERPRETATION
# ------------------------
# The "response pixels" were selected because:
#
#   AE_SJ_POST - AE_SJ_PRE >= 0.075
#
# Therefore the AE-SJ gain inside those selected pixels is a detection
# diagnostic, NOT an unbiased management-effect estimate. Do not compare
# resp_AE_SJ_gain_mean directly with the burn effect size as though the two
# samples had been selected independently.
#
# Whole-polygon AE_SJ_gain_mean is not conditioned on the +0.075 response
# threshold, but it can be diluted because the tiltak polygons are imprecise.
#
# This script deliberately stops BEFORE analysing the MI kystlynghei dataset.
# =============================================================================


# =============================================================================
# 1. REQUIRED PACKAGES
# =============================================================================

if (!exists("REPO_ROOT")) source("R/00_setup.R")

required_packages(c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "purrr",
  "ggplot2",
  "foreign"
))


# =============================================================================
# 2. INPUTS
# =============================================================================
#
# The two file.choose() calls make this script interactive (e.g. RStudio);
# it is run separately from run_all.R.

cat(
  "\nSelect the AE-SJ clearance-response CSV.\n"
)

INPUT_CSV <- file.choose()


cat(
  "\nSelect tiltakPolys_GEE.dbf (or the corresponding tiltakPolys DBF).\n"
)

TILTAK_DBF <- file.choose()


# Outputs go to outputs/naturTiltak_scrub/.
OUTPUT_DIR <- file.path(REPO_ROOT, "outputs", "naturTiltak_scrub")

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# Minimum number of valid matched pixel-equivalents required for a site to
# enter the validation analysis.
MIN_VALID_PIXELS <- 5


# The response threshold applied when the CSV was exported; checked here.
EXPECTED_GAIN_THRESHOLD <- 0.075


# A site is called "detected" when it passes site QC and contains any positive
# area/pixel weight above the +0.075 response threshold.
#
# Because the exported regional summaries can contain fractional edge-pixel
# weights, the script ALSO reports a >= 1 pixel-equivalent sensitivity count.
MIN_RESPONSE_PIXEL_EQUIVALENTS_SENSITIVITY <- 1


# =============================================================================
# 3. SJ REFERENCE VALUES AND PLOT STYLE
# =============================================================================
#
# R/model.R's fitted SJ1-SJ4 reference positions (Methodology Table 1),
# copied here; update them if the model is refitted.

SJ1_REF <- 1.000000
SJ2_REF <- 0.854205
SJ3_REF <- 0.533042
SJ4_REF <- 0.000000

SJ12_BOUNDARY <- mean(c(SJ1_REF, SJ2_REF))
SJ23_BOUNDARY <- mean(c(SJ2_REF, SJ3_REF))
SJ34_BOUNDARY <- mean(c(SJ3_REF, SJ4_REF))


SJ_LEVELS <- c(
  "SJ1-like",
  "SJ2-like",
  "SJ3-like",
  "SJ4-like"
)


SJ_COLOURS <- c(
  "SJ1-like" = "#2a78d6",
  "SJ2-like" = "#86b6ef",
  "SJ3-like" = "#e09b83",
  "SJ4-like" = "#b4442d"
)


CLEARANCE_TYPE_LEVELS <- c(
  "Scrub clearance",
  "Tree clearance",
  "Scrub + tree clearance",
  "Unspecified scrub/tree clearance"
)


theme_sj <-
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 14,
      margin = ggplot2::margin(
        b = 5
      )
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5,
      margin = ggplot2::margin(
        b = 9
      )
    ),
    plot.caption = ggplot2::element_text(
      size = 9,
      colour = "grey35",
      hjust = 0,
      margin = ggplot2::margin(
        t = 8
      )
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    strip.background = ggplot2::element_blank()
  )


# =============================================================================
# 4. GENERAL HELPERS
# =============================================================================

num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}


safe_div <- function(a, b) {
  
  a <- as.numeric(a)
  b <- as.numeric(b)
  
  out <- rep(
    NA_real_,
    max(
      length(a),
      length(b)
    )
  )
  
  a <- rep(
    a,
    length.out = length(out)
  )
  
  b <- rep(
    b,
    length.out = length(out)
  )
  
  ok <-
    is.finite(a) &
    is.finite(b) &
    b > 0
  
  out[ok] <-
    a[ok] /
    b[ok]
  
  out
}


safe_mean <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  mean(x)
}


safe_median <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  stats::median(x)
}


fmt_pct <- function(x, digits = 1) {
  
  ifelse(
    is.finite(x),
    paste0(
      formatC(
        100 * x,
        format = "f",
        digits = digits
      ),
      "%"
    ),
    "NA"
  )
}


save_plot <- function(
    p,
    filename,
    width = 12,
    height = 7
) {
  
  ggplot2::ggsave(
    filename = file.path(
      OUTPUT_DIR,
      filename
    ),
    plot = p,
    width = width,
    height = height,
    dpi = 350,
    bg = "white"
  )
}


normalise_fid <- function(x) {
  
  z <- trimws(
    as.character(x)
  )
  
  # Handles DBF numeric IDs that may appear as "1040.0".
  z <- sub(
    "\\.0+$",
    "",
    z
  )
  
  z
}


# =============================================================================
# 5. ACTION-TYPE NORMALISATION
# =============================================================================
#
# The source contains multiple generations of terminology.
#
# We intentionally distinguish SPECIFIC scrub and tree terms from ambiguous
# generic clearance terms. This is important for a validation dataset:
# "scrub_or_tree_clearance" does not tell us which vegetation type was removed,
# and "clearing cutting" is also not specific enough to assign confidently.
#
# "invasive removal" is NOT automatically treated as tree removal.
# =============================================================================


SCRUB_SPECIFIC_TOKENS <- tolower(
  c(
    "scrub/woody vegetation clearing"
  )
)


TREE_SPECIFIC_TOKENS <- tolower(
  c(
    "Sitka spruce removal",
    "sitka_removal",
    "sitka_spruce_removal",
    "tree_removal"
  )
)


GENERIC_CLEARANCE_TOKENS <- tolower(
  c(
    "scrub_or_tree_clearance",
    "tree_or_scrub_removal",
    "clearing cutting"
  )
)


ALL_CLEARANCE_TOKENS <- unique(
  c(
    SCRUB_SPECIFIC_TOKENS,
    TREE_SPECIFIC_TOKENS,
    GENERIC_CLEARANCE_TOKENS
  )
)


split_action_tokens <- function(x) {
  
  if (
    length(x) == 0 ||
    is.na(x) ||
    !nzchar(trimws(x))
  ) {
    return(character(0))
  }
  
  out <-
    stringr::str_split(
      as.character(x),
      ";",
      simplify = FALSE
    )[[1]]
  
  out <-
    stringr::str_trim(
      out
    )
  
  out <-
    tolower(
      out
    )
  
  unique(
    out[
      nzchar(out)
    ]
  )
}


has_any_token <- function(x, token_set) {
  
  xx <- split_action_tokens(x)
  
  any(
    xx %in%
      token_set
  )
}


non_clearance_token_string <- function(x) {
  
  xx <- split_action_tokens(x)
  
  other <-
    setdiff(
      xx,
      ALL_CLEARANCE_TOKENS
    )
  
  if (length(other) == 0) {
    return("")
  }
  
  paste(
    other,
    collapse = "; "
  )
}


n_non_clearance_tokens <- function(x) {
  
  xx <- split_action_tokens(x)
  
  length(
    setdiff(
      xx,
      ALL_CLEARANCE_TOKENS
    )
  )
}


classify_clearance_type <- function(x) {
  
  has_scrub <-
    has_any_token(
      x,
      SCRUB_SPECIFIC_TOKENS
    )
  
  has_tree <-
    has_any_token(
      x,
      TREE_SPECIFIC_TOKENS
    )
  
  has_generic <-
    has_any_token(
      x,
      GENERIC_CLEARANCE_TOKENS
    )
  
  if (has_scrub && has_tree) {
    return(
      "Scrub + tree clearance"
    )
  }
  
  if (has_scrub) {
    return(
      "Scrub clearance"
    )
  }
  
  if (has_tree) {
    return(
      "Tree clearance"
    )
  }
  
  if (has_generic) {
    return(
      "Unspecified scrub/tree clearance"
    )
  }
  
  NA_character_
}


# =============================================================================
# 6. READ THE CLEARANCE-RESPONSE CSV
# =============================================================================

dat_raw <-
  readr::read_csv(
    INPUT_CSV,
    col_types = readr::cols(
      fid = readr::col_character(),
      .default = readr::col_guess()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )


required_cols <- c(
  "fid",
  "tiltakYear",
  "tiltakType",
  "analysis_pre_year",
  "analysis_post_year",
  "AE_SJ_gain_threshold",
  "n_valid_pixels",
  "n_response_pixels",
  "prop_response_pixels",
  "response_area_ha",
  "AE_SJ_pre_mean",
  "AE_SJ_post_mean",
  "AE_SJ_gain_mean",
  "resp_AE_SJ_pre_mean",
  "resp_AE_SJ_post_mean",
  "resp_AE_SJ_gain_mean"
)


missing_cols <-
  setdiff(
    required_cols,
    names(dat_raw)
  )


if (length(missing_cols) > 0) {
  
  stop(
    "The response CSV is missing required columns:\n  ",
    paste(
      missing_cols,
      collapse = "\n  "
    )
  )
}


transition_cols <-
  as.vector(
    outer(
      paste0(
        "SJ",
        1:4
      ),
      paste0(
        "SJ",
        1:4
      ),
      FUN = function(a, b) {
        paste0(
          "resp_",
          a,
          "_to_",
          b
        )
      }
    )
  )


missing_transition_cols <-
  setdiff(
    transition_cols,
    names(dat_raw)
  )


if (length(missing_transition_cols) > 0) {
  
  stop(
    "The response CSV is missing response-pixel transition columns:\n  ",
    paste(
      missing_transition_cols,
      collapse = "\n  "
    )
  )
}


composition_cols <- c(
  paste0(
    "pre_prop_SJ",
    1:4
  ),
  paste0(
    "post_prop_SJ",
    1:4
  ),
  paste0(
    "resp_pre_prop_SJ",
    1:4
  ),
  paste0(
    "resp_post_prop_SJ",
    1:4
  )
)


missing_composition_cols <-
  setdiff(
    composition_cols,
    names(dat_raw)
  )


if (length(missing_composition_cols) > 0) {
  
  stop(
    "The response CSV is missing SJ-composition columns:\n  ",
    paste(
      missing_composition_cols,
      collapse = "\n  "
    )
  )
}


numeric_cols <-
  intersect(
    c(
      "tiltakYear",
      "analysis_pre_year",
      "analysis_post_year",
      "AE_SJ_gain_threshold",
      "n_valid_pixels",
      "n_response_pixels",
      "prop_response_pixels",
      "response_area_ha",
      "AE_SJ_pre_mean",
      "AE_SJ_post_mean",
      "AE_SJ_gain_mean",
      "resp_AE_SJ_pre_mean",
      "resp_AE_SJ_post_mean",
      "resp_AE_SJ_gain_mean",
      paste0(
        "pre_n_SJ",
        1:4
      ),
      paste0(
        "post_n_SJ",
        1:4
      ),
      paste0(
        "resp_pre_n_SJ",
        1:4
      ),
      paste0(
        "resp_post_n_SJ",
        1:4
      ),
      composition_cols,
      transition_cols
    ),
    names(dat_raw)
  )


dat <-
  dat_raw |>
  dplyr::mutate(
    fid =
      normalise_fid(
        fid
      ),
    dplyr::across(
      dplyr::all_of(
        numeric_cols
      ),
      num
    ),
    tiltakYear =
      as.integer(
        tiltakYear
      ),
    analysis_pre_year =
      as.integer(
        analysis_pre_year
      ),
    analysis_post_year =
      as.integer(
        analysis_post_year
      )
  )


if (anyDuplicated(dat$fid)) {
  
  duplicate_fids <-
    dat |>
    dplyr::count(
      fid,
      name = "n_rows"
    ) |>
    dplyr::filter(
      n_rows > 1
    )
  
  stop(
    "The response CSV contains duplicate FIDs. First duplicates:\n",
    paste(
      utils::capture.output(
        print(
          utils::head(
            duplicate_fids,
            20
          )
        )
      ),
      collapse = "\n"
    )
  )
}


# Check the response threshold.

observed_thresholds <-
  sort(
    unique(
      stats::na.omit(
        dat$AE_SJ_gain_threshold
      )
    )
  )


if (
  length(observed_thresholds) > 0 &&
  any(
    abs(
      observed_thresholds -
      EXPECTED_GAIN_THRESHOLD
    ) > 1e-9
  )
) {
  
  warning(
    "Export threshold differs from EXPECTED_GAIN_THRESHOLD. Observed: ",
    paste(
      observed_thresholds,
      collapse = ", "
    )
  )
}


# =============================================================================
# 7. READ ORIGINAL TILTAK DBF AND EXTRACT IMPLEMENTATION STATUS
# =============================================================================
#
# Inspection of the supplied DBF shows that implementation information is
# embedded in the English tiltakEN text, including the exact phrases:
#
#   "The report states that the measure was completed as planned."
#   "The report states that the measure was not completed as planned."
#
# The negative phrase is tested FIRST because it contains the positive phrase
# "completed as planned".
#
# Only records explicitly stating "not completed as planned" are excluded.
# Records with no explicit completion statement remain, labelled, so that
# sensitivity analyses can restrict to explicitly completed records.
# =============================================================================


tiltak_dbf_raw <-
  foreign::read.dbf(
    TILTAK_DBF,
    as.is = TRUE
  )


dbf_required <- c(
  "fid",
  "tiltakEN",
  "tiltakType"
)


missing_dbf <-
  setdiff(
    dbf_required,
    names(tiltak_dbf_raw)
  )


if (length(missing_dbf) > 0) {
  
  stop(
    "The tiltak DBF is missing required fields:\n  ",
    paste(
      missing_dbf,
      collapse = "\n  "
    )
  )
}


tiltak_status <-
  tibble::as_tibble(
    tiltak_dbf_raw
  ) |>
  dplyr::transmute(
    fid =
      normalise_fid(
        .data$fid
      ),
    tiltakEN_dbf =
      as.character(
        .data$tiltakEN
      ),
    tiltakType_dbf =
      as.character(
        .data$tiltakType
      )
  ) |>
  dplyr::mutate(
    implementation_status =
      dplyr::case_when(
        stringr::str_detect(
          stringr::str_to_lower(
            dplyr::coalesce(
              tiltakEN_dbf,
              ""
            )
          ),
          stringr::fixed(
            "the report states that the measure was not completed as planned."
          )
        ) ~
          "Not completed as planned",
        stringr::str_detect(
          stringr::str_to_lower(
            dplyr::coalesce(
              tiltakEN_dbf,
              ""
            )
          ),
          stringr::fixed(
            "the report states that the measure was completed as planned."
          )
        ) ~
          "Completed as planned",
        TRUE ~
          "No explicit completion statement"
      )
  )


if (anyDuplicated(tiltak_status$fid)) {
  
  stop(
    "FID is not unique in the tiltak DBF. Resolve duplicate IDs before analysis."
  )
}


# =============================================================================
# 8. JOIN STATUS AND STANDARDISE CLEARANCE TYPE
# =============================================================================

dat <-
  dat |>
  dplyr::left_join(
    tiltak_status,
    by = "fid"
  )


if (any(is.na(dat$implementation_status))) {
  
  missing_status_fids <-
    dat |>
    dplyr::filter(
      is.na(
        implementation_status
      )
    ) |>
    dplyr::pull(
      fid
    )
  
  stop(
    "Some response-CSV FIDs were not found in the tiltak DBF:\n  ",
    paste(
      missing_status_fids,
      collapse = ", "
    )
  )
}


# Prefer the tiltakType carried in the response CSV, because that is the
# field used to create the clearance sample. The DBF copy is retained for QC.

dat <-
  dat |>
  dplyr::mutate(
    reported_scrub_specific =
      vapply(
        tiltakType,
        has_any_token,
        logical(1),
        token_set = SCRUB_SPECIFIC_TOKENS
      ),
    reported_tree_specific =
      vapply(
        tiltakType,
        has_any_token,
        logical(1),
        token_set = TREE_SPECIFIC_TOKENS
      ),
    reported_generic_clearance =
      vapply(
        tiltakType,
        has_any_token,
        logical(1),
        token_set = GENERIC_CLEARANCE_TOKENS
      ),
    clearance_type =
      vapply(
        tiltakType,
        classify_clearance_type,
        character(1)
      ),
    clearance_type =
      factor(
        clearance_type,
        levels = CLEARANCE_TYPE_LEVELS
      ),
    non_clearance_actions =
      vapply(
        tiltakType,
        non_clearance_token_string,
        character(1)
      ),
    n_non_clearance_actions =
      vapply(
        tiltakType,
        n_non_clearance_tokens,
        integer(1)
      ),
    has_co_intervention =
      n_non_clearance_actions > 0
  )


# Transparent lookup showing exactly how every original tiltakType string
# maps into the validation groups.

type_lookup <-
  dat |>
  dplyr::count(
    tiltakType,
    clearance_type,
    reported_scrub_specific,
    reported_tree_specific,
    reported_generic_clearance,
    non_clearance_actions,
    has_co_intervention,
    name = "n_fids"
  ) |>
  dplyr::arrange(
    clearance_type,
    dplyr::desc(
      n_fids
    ),
    tiltakType
  )


# =============================================================================
# 9. APPLY IMPLEMENTATION FILTER AND SITE QC
# =============================================================================

site_results <-
  dat |>
  dplyr::mutate(
    excluded_not_completed =
      implementation_status ==
      "Not completed as planned",
    qc_enough_pixels =
      is.finite(
        n_valid_pixels
      ) &
      n_valid_pixels >=
      MIN_VALID_PIXELS,
    analysis_eligible =
      !excluded_not_completed &
      qc_enough_pixels,
    detected_clearance =
      analysis_eligible &
      is.finite(
        n_response_pixels
      ) &
      n_response_pixels > 0,
    detected_clearance_1px_sensitivity =
      analysis_eligible &
      is.finite(
        n_response_pixels
      ) &
      n_response_pixels >=
      MIN_RESPONSE_PIXEL_EQUIVALENTS_SENSITIVITY,
    analysis_status =
      dplyr::case_when(
        excluded_not_completed ~
          "Excluded: not completed as planned",
        !qc_enough_pixels ~
          paste0(
            "Excluded: <",
            MIN_VALID_PIXELS,
            " valid pixel-equivalents"
          ),
        detected_clearance ~
          "Detected AE-SJ clearance response",
        TRUE ~
          "Eligible but no threshold response"
      )
  )


# =============================================================================
# 10. RESPONSE-PIXEL SOURCE-ZONE AND TRANSITION METRICS
# =============================================================================
#
# These transition counts are only for pixels selected by the response threshold
# threshold.
#
# "from SJ3" and "from SJ4" below describe the PRE-treatment SJ-like zone.
#
# Two validation ideas are retained:
#
#   SOURCE SUPPORT:
#     Was there at least one detected response pixel whose PRE zone matched
#     the expected vegetation structure?
#
#   CROSSING SUPPORT:
#     Did at least one such pixel actually cross an SJ midpoint boundary toward
#     SJ1 after treatment?
#
# Crossing support is the stricter criterion. A pixel may move by >0.075 AE-SJ
# yet remain within the same SJ-like score zone.
# =============================================================================


site_results <-
  site_results |>
  dplyr::mutate(
    n_resp_from_SJ1 =
      resp_SJ1_to_SJ1 +
      resp_SJ1_to_SJ2 +
      resp_SJ1_to_SJ3 +
      resp_SJ1_to_SJ4,
    n_resp_from_SJ2 =
      resp_SJ2_to_SJ1 +
      resp_SJ2_to_SJ2 +
      resp_SJ2_to_SJ3 +
      resp_SJ2_to_SJ4,
    n_resp_from_SJ3 =
      resp_SJ3_to_SJ1 +
      resp_SJ3_to_SJ2 +
      resp_SJ3_to_SJ3 +
      resp_SJ3_to_SJ4,
    n_resp_from_SJ4 =
      resp_SJ4_to_SJ1 +
      resp_SJ4_to_SJ2 +
      resp_SJ4_to_SJ3 +
      resp_SJ4_to_SJ4,
    
    prop_resp_from_SJ1 =
      safe_div(
        n_resp_from_SJ1,
        n_response_pixels
      ),
    prop_resp_from_SJ2 =
      safe_div(
        n_resp_from_SJ2,
        n_response_pixels
      ),
    prop_resp_from_SJ3 =
      safe_div(
        n_resp_from_SJ3,
        n_response_pixels
      ),
    prop_resp_from_SJ4 =
      safe_div(
        n_resp_from_SJ4,
        n_response_pixels
      ),
    prop_resp_from_SJ3_or_SJ4 =
      safe_div(
        n_resp_from_SJ3 +
          n_resp_from_SJ4,
        n_response_pixels
      ),
    
    # Strict movement from SJ3 toward an earlier SJ-like zone.
    n_SJ3_cross_toward_SJ1 =
      resp_SJ3_to_SJ1 +
      resp_SJ3_to_SJ2,
    
    # Strict movement from SJ4 toward an earlier SJ-like zone.
    n_SJ4_cross_toward_SJ1 =
      resp_SJ4_to_SJ1 +
      resp_SJ4_to_SJ2 +
      resp_SJ4_to_SJ3,
    
    # Future MI-compatible summary:
    # late-successional source pixels ending specifically in SJ1 or SJ2.
    n_SJ34_to_SJ12 =
      resp_SJ3_to_SJ1 +
      resp_SJ3_to_SJ2 +
      resp_SJ4_to_SJ1 +
      resp_SJ4_to_SJ2,
    
    prop_SJ34_to_SJ12_of_response =
      safe_div(
        n_SJ34_to_SJ12,
        n_response_pixels
      ),
    
    has_response_from_SJ3 =
      detected_clearance &
      n_resp_from_SJ3 > 0,
    
    has_response_from_SJ4 =
      detected_clearance &
      n_resp_from_SJ4 > 0,
    
    has_SJ3_cross_toward_SJ1 =
      detected_clearance &
      n_SJ3_cross_toward_SJ1 > 0,
    
    has_SJ4_cross_toward_SJ1 =
      detected_clearance &
      n_SJ4_cross_toward_SJ1 > 0
  ) |>
  dplyr::mutate(
    expected_source_support =
      dplyr::case_when(
        !detected_clearance ~
          NA,
        clearance_type ==
          "Scrub clearance" ~
          has_response_from_SJ3,
        clearance_type ==
          "Tree clearance" ~
          has_response_from_SJ4,
        clearance_type ==
          "Scrub + tree clearance" ~
          has_response_from_SJ3 &
          has_response_from_SJ4,
        clearance_type ==
          "Unspecified scrub/tree clearance" ~
          has_response_from_SJ3 |
          has_response_from_SJ4,
        TRUE ~
          NA
      ),
    expected_crossing_support =
      dplyr::case_when(
        !detected_clearance ~
          NA,
        clearance_type ==
          "Scrub clearance" ~
          has_SJ3_cross_toward_SJ1,
        clearance_type ==
          "Tree clearance" ~
          has_SJ4_cross_toward_SJ1,
        clearance_type ==
          "Scrub + tree clearance" ~
          has_SJ3_cross_toward_SJ1 &
          has_SJ4_cross_toward_SJ1,
        clearance_type ==
          "Unspecified scrub/tree clearance" ~
          has_SJ3_cross_toward_SJ1 |
          has_SJ4_cross_toward_SJ1,
        TRUE ~
          NA
      )
  )


# =============================================================================
# 11. FILTER / DETECTION REPORT
# =============================================================================

filter_summary <-
  tibble::tibble(
    step = c(
      "Response-CSV FIDs",
      "Explicitly not completed as planned",
      "Retained after implementation filter",
      paste0(
        "Retained with >= ",
        MIN_VALID_PIXELS,
        " valid pixel-equivalents"
      ),
      "Detected: any +0.075 response pixel weight",
      paste0(
        "Detected sensitivity: >= ",
        MIN_RESPONSE_PIXEL_EQUIVALENTS_SENSITIVITY,
        " response pixel-equivalent"
      ),
      "Eligible but no +0.075 response",
      paste0(
        "Excluded for < ",
        MIN_VALID_PIXELS,
        " valid pixel-equivalents"
      )
    ),
    n_fids = c(
      nrow(
        site_results
      ),
      sum(
        site_results$excluded_not_completed,
        na.rm = TRUE
      ),
      sum(
        !site_results$excluded_not_completed,
        na.rm = TRUE
      ),
      sum(
        site_results$analysis_eligible,
        na.rm = TRUE
      ),
      sum(
        site_results$detected_clearance,
        na.rm = TRUE
      ),
      sum(
        site_results$detected_clearance_1px_sensitivity,
        na.rm = TRUE
      ),
      sum(
        site_results$analysis_eligible &
          !site_results$detected_clearance,
        na.rm = TRUE
      ),
      sum(
        !site_results$excluded_not_completed &
          !site_results$qc_enough_pixels,
        na.rm = TRUE
      )
    )
  )


detection_summary_by_type <-
  site_results |>
  dplyr::filter(
    !excluded_not_completed
  ) |>
  dplyr::group_by(
    clearance_type
  ) |>
  dplyr::summarise(
    n_after_implementation_filter =
      dplyr::n(),
    n_qc_eligible =
      sum(
        analysis_eligible,
        na.rm = TRUE
      ),
    n_detected =
      sum(
        detected_clearance,
        na.rm = TRUE
      ),
    n_detected_1px_sensitivity =
      sum(
        detected_clearance_1px_sensitivity,
        na.rm = TRUE
      ),
    detection_rate =
      safe_div(
        n_detected,
        n_qc_eligible
      ),
    median_response_extent =
      safe_median(
        prop_response_pixels[
          analysis_eligible
        ]
      ),
    median_whole_polygon_AE_SJ_gain =
      safe_median(
        AE_SJ_gain_mean[
          analysis_eligible
        ]
      ),
    median_detected_pixel_AE_SJ_gain =
      safe_median(
        resp_AE_SJ_gain_mean[
          detected_clearance
        ]
      ),
    .groups = "drop"
  )


# =============================================================================
# 12. VALIDATION AGAINST REPORTED SCRUB / TREE MANAGEMENT
# =============================================================================

expected_source_validation <-
  site_results |>
  dplyr::filter(
    detected_clearance
  ) |>
  dplyr::group_by(
    clearance_type
  ) |>
  dplyr::summarise(
    n_detected_sites =
      dplyr::n(),
    
    n_sites_with_response_from_SJ3 =
      sum(
        has_response_from_SJ3,
        na.rm = TRUE
      ),
    
    prop_sites_with_response_from_SJ3 =
      safe_div(
        n_sites_with_response_from_SJ3,
        n_detected_sites
      ),
    
    n_sites_with_response_from_SJ4 =
      sum(
        has_response_from_SJ4,
        na.rm = TRUE
      ),
    
    prop_sites_with_response_from_SJ4 =
      safe_div(
        n_sites_with_response_from_SJ4,
        n_detected_sites
      ),
    
    n_sites_with_both_SJ3_and_SJ4_sources =
      sum(
        has_response_from_SJ3 &
          has_response_from_SJ4,
        na.rm = TRUE
      ),
    
    prop_sites_with_both_SJ3_and_SJ4_sources =
      safe_div(
        n_sites_with_both_SJ3_and_SJ4_sources,
        n_detected_sites
      ),
    
    median_prop_response_pixels_from_SJ3 =
      safe_median(
        prop_resp_from_SJ3
      ),
    
    median_prop_response_pixels_from_SJ4 =
      safe_median(
        prop_resp_from_SJ4
      ),
    
    median_prop_response_pixels_from_SJ3_or_SJ4 =
      safe_median(
        prop_resp_from_SJ3_or_SJ4
      ),
    
    n_sites_matching_expected_source =
      sum(
        expected_source_support,
        na.rm = TRUE
      ),
    
    prop_sites_matching_expected_source =
      safe_div(
        n_sites_matching_expected_source,
        n_detected_sites
      ),
    
    n_sites_matching_expected_boundary_crossing =
      sum(
        expected_crossing_support,
        na.rm = TRUE
      ),
    
    prop_sites_matching_expected_boundary_crossing =
      safe_div(
        n_sites_matching_expected_boundary_crossing,
        n_detected_sites
      ),
    
    n_sites_with_SJ34_to_SJ12 =
      sum(
        n_SJ34_to_SJ12 > 0,
        na.rm = TRUE
      ),
    
    prop_sites_with_SJ34_to_SJ12 =
      safe_div(
        n_sites_with_SJ34_to_SJ12,
        n_detected_sites
      ),
    
    median_prop_SJ34_to_SJ12_of_response =
      safe_median(
        prop_SJ34_to_SJ12_of_response
      ),
    
    .groups = "drop"
  )


# =============================================================================
# 13. RESPONSE-PIXEL PRE-TREATMENT SOURCE COMPOSITION
# =============================================================================
#
# Equal site weighting: first calculate each detected site's fraction of
# response pixels originating in each SJ-like zone, then summarise those site
# fractions by clearance type.
# =============================================================================


response_source_long <-
  site_results |>
  dplyr::filter(
    detected_clearance
  ) |>
  dplyr::select(
    fid,
    clearance_type,
    prop_resp_from_SJ1,
    prop_resp_from_SJ2,
    prop_resp_from_SJ3,
    prop_resp_from_SJ4
  ) |>
  tidyr::pivot_longer(
    cols =
      dplyr::starts_with(
        "prop_resp_from_SJ"
      ),
    names_to =
      "source_zone",
    values_to =
      "proportion"
  ) |>
  dplyr::mutate(
    source_zone =
      dplyr::recode(
        source_zone,
        prop_resp_from_SJ1 =
          "SJ1-like",
        prop_resp_from_SJ2 =
          "SJ2-like",
        prop_resp_from_SJ3 =
          "SJ3-like",
        prop_resp_from_SJ4 =
          "SJ4-like"
      ),
    source_zone =
      factor(
        source_zone,
        levels = SJ_LEVELS
      )
  )


response_source_summary <-
  response_source_long |>
  dplyr::group_by(
    clearance_type,
    source_zone
  ) |>
  dplyr::summarise(
    n_sites =
      dplyr::n_distinct(
        fid
      ),
    mean_site_proportion =
      safe_mean(
        proportion
      ),
    median_site_proportion =
      safe_median(
        proportion
      ),
    .groups = "drop"
  )


# =============================================================================
# 14. FULL 4 x 4 RESPONSE-PIXEL TRANSITION SUMMARY BY CLEARANCE TYPE
# =============================================================================

transition_long <-
  site_results |>
  dplyr::filter(
    detected_clearance
  ) |>
  dplyr::select(
    fid,
    clearance_type,
    dplyr::all_of(
      transition_cols
    )
  ) |>
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        transition_cols
      ),
    names_to =
      "transition",
    values_to =
      "n_pixel_equivalents"
  ) |>
  tidyr::extract(
    col =
      transition,
    into = c(
      "from",
      "to"
    ),
    regex =
      "resp_(SJ[1-4])_to_(SJ[1-4])",
    remove =
      TRUE
  ) |>
  dplyr::mutate(
    from =
      factor(
        from,
        levels = paste0(
          "SJ",
          1:4
        )
      ),
    to =
      factor(
        to,
        levels = paste0(
          "SJ",
          1:4
        )
      ),
    n_pixel_equivalents =
      as.numeric(
        n_pixel_equivalents
      )
  )


# Pooled pixel-equivalent transition totals, useful for seeing the actual
# transition pathways. These are AREA/PIXEL weighted and can therefore be
# dominated by larger polygons. Keep alongside the equal-site summaries above.

transition_summary_by_type <-
  transition_long |>
  dplyr::group_by(
    clearance_type,
    from,
    to
  ) |>
  dplyr::summarise(
    n_pixel_equivalents =
      sum(
        n_pixel_equivalents,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>
  dplyr::group_by(
    clearance_type,
    from
  ) |>
  dplyr::mutate(
    from_zone_total =
      sum(
        n_pixel_equivalents,
        na.rm = TRUE
      ),
    prop_of_source_zone =
      safe_div(
        n_pixel_equivalents,
        from_zone_total
      )
  ) |>
  dplyr::ungroup()


# =============================================================================
# 15. PRE/POST COMPOSITION TABLES
# =============================================================================

whole_composition_long <-
  purrr::map_dfr(
    1:4,
    function(k) {
      
      tibble::tibble(
        fid =
          site_results$fid,
        clearance_type =
          site_results$clearance_type,
        analysis_eligible =
          site_results$analysis_eligible,
        SJ_zone =
          paste0(
            "SJ",
            k,
            "-like"
          ),
        Before =
          site_results[[
            paste0(
              "pre_prop_SJ",
              k
            )
          ]],
        After =
          site_results[[
            paste0(
              "post_prop_SJ",
              k
            )
          ]]
      )
    }
  ) |>
  dplyr::filter(
    analysis_eligible
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Before,
      After
    ),
    names_to =
      "period",
    values_to =
      "proportion"
  ) |>
  dplyr::mutate(
    period =
      factor(
        period,
        levels = c(
          "Before",
          "After"
        )
      ),
    SJ_zone =
      factor(
        SJ_zone,
        levels = SJ_LEVELS
      )
  )


response_composition_long <-
  purrr::map_dfr(
    1:4,
    function(k) {
      
      tibble::tibble(
        fid =
          site_results$fid,
        clearance_type =
          site_results$clearance_type,
        detected_clearance =
          site_results$detected_clearance,
        SJ_zone =
          paste0(
            "SJ",
            k,
            "-like"
          ),
        Before =
          site_results[[
            paste0(
              "resp_pre_prop_SJ",
              k
            )
          ]],
        After =
          site_results[[
            paste0(
              "resp_post_prop_SJ",
              k
            )
          ]]
      )
    }
  ) |>
  dplyr::filter(
    detected_clearance
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Before,
      After
    ),
    names_to =
      "period",
    values_to =
      "proportion"
  ) |>
  dplyr::mutate(
    period =
      factor(
        period,
        levels = c(
          "Before",
          "After"
        )
      ),
    SJ_zone =
      factor(
        SJ_zone,
        levels = SJ_LEVELS
      )
  )


whole_composition_summary <-
  whole_composition_long |>
  dplyr::group_by(
    clearance_type,
    period,
    SJ_zone
  ) |>
  dplyr::summarise(
    n_sites =
      dplyr::n_distinct(
        fid
      ),
    mean_proportion =
      safe_mean(
        proportion
      ),
    median_proportion =
      safe_median(
        proportion
      ),
    .groups = "drop"
  )


response_composition_summary <-
  response_composition_long |>
  dplyr::group_by(
    clearance_type,
    period,
    SJ_zone
  ) |>
  dplyr::summarise(
    n_sites =
      dplyr::n_distinct(
        fid
      ),
    mean_proportion =
      safe_mean(
        proportion
      ),
    median_proportion =
      safe_median(
        proportion
      ),
    .groups = "drop"
  )


# =============================================================================
# 16. DESCRIPTIVE AE-SJ EFFECT-SIZE TABLE
# =============================================================================
#
# For future comparison with burning:
#
#   whole_polygon_* = not conditioned on threshold response, but may be diluted
#                     by poorly delineated tiltak polygons.
#
#   response_pixel_* = conditioned on +0.075, so DO NOT use it as an unbiased
#                      effect-size comparison with burning.
#
# Co-intervention counts are included because many tiltak records contain
# burning, grazing, mowing, restoration, etc. in addition to clearance.
# =============================================================================

effect_size_descriptives <-
  site_results |>
  dplyr::filter(
    analysis_eligible
  ) |>
  dplyr::group_by(
    clearance_type
  ) |>
  dplyr::summarise(
    n_sites =
      dplyr::n(),
    n_clearance_only_records =
      sum(
        !has_co_intervention,
        na.rm = TRUE
      ),
    n_with_co_intervention =
      sum(
        has_co_intervention,
        na.rm = TRUE
      ),
    mean_whole_polygon_AE_SJ_gain =
      safe_mean(
        AE_SJ_gain_mean
      ),
    median_whole_polygon_AE_SJ_gain =
      safe_median(
        AE_SJ_gain_mean
      ),
    n_detected =
      sum(
        detected_clearance,
        na.rm = TRUE
      ),
    mean_response_extent =
      safe_mean(
        prop_response_pixels
      ),
    median_response_extent =
      safe_median(
        prop_response_pixels
      ),
    mean_response_pixel_AE_SJ_gain =
      safe_mean(
        resp_AE_SJ_gain_mean[
          detected_clearance
        ]
      ),
    median_response_pixel_AE_SJ_gain =
      safe_median(
        resp_AE_SJ_gain_mean[
          detected_clearance
        ]
      ),
    .groups = "drop"
  )


# =============================================================================
# 17. WRITE TABLES
# =============================================================================

readr::write_csv(
  filter_summary,
  file.path(
    OUTPUT_DIR,
    "00_filter_summary.csv"
  )
)


readr::write_csv(
  site_results |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        )
    ),
  file.path(
    OUTPUT_DIR,
    "01_clearance_validation_site_results.csv"
  )
)


readr::write_csv(
  type_lookup |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        )
    ),
  file.path(
    OUTPUT_DIR,
    "02_tiltakType_unification_lookup.csv"
  )
)


readr::write_csv(
  detection_summary_by_type |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        )
    ),
  file.path(
    OUTPUT_DIR,
    "03_detection_summary_by_clearance_type.csv"
  )
)


readr::write_csv(
  expected_source_validation |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        )
    ),
  file.path(
    OUTPUT_DIR,
    "04_expected_source_validation_by_clearance_type.csv"
  )
)


readr::write_csv(
  response_source_summary |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        ),
      source_zone =
        as.character(
          source_zone
        )
    ),
  file.path(
    OUTPUT_DIR,
    "05_response_source_composition_by_clearance_type.csv"
  )
)


readr::write_csv(
  transition_summary_by_type |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        ),
      from =
        as.character(
          from
        ),
      to =
        as.character(
          to
        )
    ),
  file.path(
    OUTPUT_DIR,
    "06_response_transition_by_clearance_type.csv"
  )
)


readr::write_csv(
  whole_composition_summary |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        ),
      period =
        as.character(
          period
        ),
      SJ_zone =
        as.character(
          SJ_zone
        )
    ),
  file.path(
    OUTPUT_DIR,
    "07_whole_polygon_SJ_composition_by_clearance_type.csv"
  )
)


readr::write_csv(
  response_composition_summary |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        ),
      period =
        as.character(
          period
        ),
      SJ_zone =
        as.character(
          SJ_zone
        )
    ),
  file.path(
    OUTPUT_DIR,
    "08_response_pixel_SJ_composition_by_clearance_type.csv"
  )
)


readr::write_csv(
  effect_size_descriptives |>
    dplyr::mutate(
      clearance_type =
        as.character(
          clearance_type
        )
    ),
  file.path(
    OUTPUT_DIR,
    "09_effect_size_descriptives_by_clearance_type.csv"
  )
)


# =============================================================================
# 18. PLOT FUNCTION - PRE/POST SJ COMPOSITION FACETED BY CLEARANCE TYPE
# =============================================================================

make_composition_plot <- function(
    summary_dat,
    title,
   # subtitle,
    y_title
) {
  
  plot_dat <-
    summary_dat |>
    dplyr::mutate(
      percentage =
        100 *
        mean_proportion,
      clearance_type =
        factor(
          clearance_type,
          levels = CLEARANCE_TYPE_LEVELS
        ),
      period =
        factor(
          period,
          levels = c(
            "Before",
            "After"
          )
        ),
      SJ_zone =
        factor(
          SJ_zone,
          levels = SJ_LEVELS
        )
    )
  
  
  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = period,
      y = percentage,
      fill = SJ_zone
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
      size = 3.1
    ) +
    
    ggplot2::facet_wrap(
      ~ clearance_type,
      nrow = 1,
      drop = FALSE
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
      labels =
        function(x) {
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
        title,
      
      x =
        NULL,
      y =
        y_title,
      fill =
        "SJ type",
      caption =
        paste0(
          "Each FID contributes equally to the group mean. ",
          "SJ1-SJ4 are AE-SJ score zones, not field-observed SJ classes."
        )
    ) +
    
    theme_sj +
    
    ggplot2::theme(
      legend.position =
        "top"
    )
}


# =============================================================================
# 19. FIGURE 1 - WHOLE VALID POLYGON
# =============================================================================

p_whole <-
  make_composition_plot(
    whole_composition_summary,
    title =
      "Reported scrub/tree clearance: pre- and post-treatment SJ composition",
    subtitle =
      paste0(
        "Whole valid masked tiltak polygon; PRE = tiltakYear - 1, ",
        "POST = tiltakYear + 1"
      ),
    y_title =
      "Mean fraction of valid site pixels"
  )




# =============================================================================
# 20. FIGURE 2 - DETECTED RESPONSE PIXELS
# =============================================================================

p_response <-
  make_composition_plot(
    response_composition_summary,
    title =
      "Natur tiltak specifying clearance: pre- and post-treatment SJ composition",
    
    y_title =
      "Mean percentge of SJ type"
  )


print(
  p_response
)


save_plot(
  p_response,
  "02_response_pixel_pre_post_SJ_composition_by_clearance_type.png",
  width = 15,
  height = 6.5
)


# =============================================================================
# 21. TEXT VALIDATION REPORT
# =============================================================================

n_input <-
  nrow(
    site_results
  )

n_not_completed <-
  sum(
    site_results$excluded_not_completed,
    na.rm = TRUE
  )

n_after_status <-
  sum(
    !site_results$excluded_not_completed,
    na.rm = TRUE
  )

n_eligible <-
  sum(
    site_results$analysis_eligible,
    na.rm = TRUE
  )

n_detected <-
  sum(
    site_results$detected_clearance,
    na.rm = TRUE
  )

n_detected_1px <-
  sum(
    site_results$detected_clearance_1px_sensitivity,
    na.rm = TRUE
  )

detection_rate_overall <-
  safe_div(
    n_detected,
    n_eligible
  )


report_lines <- c(
  "KYSTLYNGHEI AE-SJ SCRUB / TREE CLEARANCE VALIDATION",
  "====================================================",
  "",
  "Purpose:",
  "Use independently reported tiltak clearance actions as a validation dataset",
  "for the AE-SJ clearance-response method.",
  "",
  "Implementation-status filter:",
  paste0(
    "  Response-CSV FIDs: ",
    n_input
  ),
  paste0(
    "  Excluded because tiltakEN explicitly says the measure was not completed as planned: ",
    n_not_completed
  ),
  paste0(
    "  Remaining after implementation filter: ",
    n_after_status
  ),
  paste0(
    "  QC-eligible sites (>= ",
    MIN_VALID_PIXELS,
    " valid pixel-equivalents): ",
    n_eligible
  ),
  "",
  "Detection:",
  paste0(
    "  Sites with any +",
    EXPECTED_GAIN_THRESHOLD,
    " AE-SJ response pixel weight: ",
    n_detected,
    " / ",
    n_eligible,
    " (",
    fmt_pct(
      detection_rate_overall
    ),
    ")"
  ),
  paste0(
    "  Sensitivity count requiring >= ",
    MIN_RESPONSE_PIXEL_EQUIVALENTS_SENSITIVITY,
    " response pixel-equivalent: ",
    n_detected_1px
  ),
  "",
  "Important:",
  "  Response-pixel AE-SJ gain is conditioned on the +0.075 threshold and is",
  "  therefore a detection diagnostic, not an unbiased treatment-effect estimate.",
  "  Many tiltak records also contain other management actions; see the",
  "  has_co_intervention and non_clearance_actions fields before causal/effect-size",
  "  comparison with burning.",
  "",
  "Clearance-type detection summary:"
)


for (
  i in seq_len(
    nrow(
      detection_summary_by_type
    )
  )
) {
  
  rr <-
    detection_summary_by_type[
      i,
      ,
      drop = FALSE
    ]
  
  report_lines <-
    c(
      report_lines,
      paste0(
        "  ",
        as.character(
          rr$clearance_type
        ),
        ": ",
        rr$n_detected,
        " detected / ",
        rr$n_qc_eligible,
        " QC-eligible (",
        fmt_pct(
          rr$detection_rate
        ),
        ")"
      )
    )
}


report_lines <-
  c(
    report_lines,
    "",
    "Reported management vs PRE-treatment source-zone validation:"
  )


for (
  i in seq_len(
    nrow(
      expected_source_validation
    )
  )
) {
  
  rr <-
    expected_source_validation[
      i,
      ,
      drop = FALSE
    ]
  
  report_lines <-
    c(
      report_lines,
      paste0(
        "  ",
        as.character(
          rr$clearance_type
        ),
        ": expected source present in ",
        rr$n_sites_matching_expected_source,
        " / ",
        rr$n_detected_sites,
        " detected sites (",
        fmt_pct(
          rr$prop_sites_matching_expected_source
        ),
        "); expected source crossed an SJ boundary toward SJ1 in ",
        rr$n_sites_matching_expected_boundary_crossing,
        " / ",
        rr$n_detected_sites,
        " (",
        fmt_pct(
          rr$prop_sites_matching_expected_boundary_crossing
        ),
        ")."
      )
    )
}


report_lines <-
  c(
    report_lines,
    "",
    "Expected-source rules:",
    "  Scrub clearance: SJ3-like PRE-treatment response pixels.",
    "  Tree clearance: SJ4-like PRE-treatment response pixels.",
    "  Scrub + tree clearance: both SJ3-like and SJ4-like PRE-treatment response pixels.",
    "  Unspecified clearance: either SJ3-like or SJ4-like PRE-treatment response pixels.",
    "",
    "The stricter boundary-crossing version requires those expected pixels to move",
    "into an earlier SJ-like zone. Continuous +0.075 improvement can occur without",
    "crossing a midpoint boundary, so source support and crossing support answer",
    "slightly different validation questions.",
    "",
    paste0(
      "Outputs written to: ",
      OUTPUT_DIR
    )
  )


writeLines(
  report_lines,
  file.path(
    OUTPUT_DIR,
    "VALIDATION_REPORT.txt"
  )
)


cat(
  "\n",
  paste(
    report_lines,
    collapse = "\n"
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 22. CONSOLE TABLES
# =============================================================================

cat(
  "\n\nFILTER SUMMARY\n"
)

print(
  filter_summary,
  n = Inf
)


cat(
  "\n\nDETECTION SUMMARY BY CLEARANCE TYPE\n"
)

print(
  detection_summary_by_type,
  n = Inf
)


cat(
  "\n\nEXPECTED SOURCE-ZONE VALIDATION\n"
)

print(
  expected_source_validation,
  n = Inf
)


cat(
  "\n\nEFFECT-SIZE DESCRIPTIVES (interpret response-pixel gain cautiously)\n"
)

print(
  effect_size_descriptives,
  n = Inf
)


cat(
  "\n\nValidation analysis complete.\nOutput directory:\n",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
