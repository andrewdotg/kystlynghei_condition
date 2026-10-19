# =============================================================================
# KYSTLYNGHEI MI - SJ SCORE VALIDATION: FIELD SJ vs. MODAL AE-SJ CLASS, SAME YEAR
# =============================================================================
#
# A direct, independent check that's different from R/model.R/R/validation.R's
# own validation (which compares the continuous AE_SJ score against field sj
# at the POINT level, pooled across all years): this compares, for every MI
# site with a field survey, the site's FIELD-ASSESSED sj class against the
# MODAL AE-SJ pixel class this pipeline assigns that SAME site IN THE SAME
# YEAR the survey happened - i.e. does the annual per-pixel classification
# actually recover the known field-assessed condition, site by site, year by
# year, not just in the aggregate the original model validation already
# checked.
#
# FIELD SIDE: field_sj and survey_year come from the MI shapefile's own
# X7RA.SJ and year fields, valid for 13,853 of 14,330 sites (96.7 %). The
# model's point-scores CSV covers only the calibration subset and is not
# used here.
#
# MODAL AE-SJ CLASS: python/mi_aesj_annual_products.py's annual
# CSVs already carry n_SJ1_like..n_SJ4_like (pixel counts per class, per
# site, per year) - the modal class is just argmax across those four, no new
# Python computation needed.
#
# =============================================================================


# =============================================================================
# 0. SHARED SETUP
# =============================================================================

if (!exists("REPO_ROOT")) source("R/00_setup.R")
source(file.path(REPO_ROOT, "R/theme.R"))
required_packages(c("readr", "dplyr", "sf", "purrr", "tidyr", "ggplot2", "stringr", "scales"))


# =============================================================================
# 1. SETTINGS
# =============================================================================

MI_SHAPEFILE <- file.path(REPO_ROOT, "inputs", "MI_kystlynghei.shp")
ANNUAL_INPUT_DIR <- file.path(REPO_ROOT, "python", "mi_aesj_annual_results")

OUTPUT_DIR <- file.path(OUTPUTS_ROOT, "sj_score_validation")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SJ_LABELS <- c("SJ1", "SJ2", "SJ3", "SJ4")


# =============================================================================
# 2. FIELD SIDE: (fid, survey_year, field_sj)
# =============================================================================

if (!file.exists(MI_SHAPEFILE)) {
  stop("Could not find ", MI_SHAPEFILE)
}

mi_raw <- sf::st_read(MI_SHAPEFILE, quiet = TRUE) %>% sf::st_drop_geometry()

# X7RA.SJ: the shapefile's own per-polygon field-assessed SJ class (R
# mangles the leading "7RA.SJ" -> "X7RA.SJ" since column names can't start
# with a digit). Values are "1".."4"; a handful are "X" (unclear/no-data,
# excluded below) or blank/NA (475 sites never field-assessed at all).
# `year` is that assessment's survey year, already restricted to 2018-2025
# for every row that also has a valid X7RA.SJ (checked - 0 mismatches).
field_sj <- mi_raw %>%
  transmute(
    fid = trimws(as.character(fid)),
    survey_year = suppressWarnings(as.integer(year)),
    field_sj = suppressWarnings(as.integer(X7RA.SJ))
  ) %>%
  filter(!is.na(field_sj), field_sj %in% 1:4, !is.na(survey_year))

cat("\nField-surveyed (site, year) pairs found: ", nrow(field_sj),
    " across ", dplyr::n_distinct(field_sj$fid), " MI sites (",
    round(100 * nrow(field_sj) / nrow(mi_raw), 1), "% of all ", nrow(mi_raw), " sites)\n", sep = "")


# =============================================================================
# 3. MODEL SIDE: modal AE-SJ class per (fid, year), for the survey years only
# =============================================================================

survey_years <- sort(unique(field_sj$survey_year))

read_modal_year <- function(year) {
  path <- file.path(ANNUAL_INPUT_DIR, sprintf("MI_aesj_annual_%d.csv", year))
  if (!file.exists(path)) {
    cat("NOTE: ", path, " not found - skipping survey_year ", year, ".\n", sep = "")
    return(NULL)
  }
  dat <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  needed <- c("fid", "n_SJ1_like", "n_SJ2_like", "n_SJ3_like", "n_SJ4_like")
  if (!all(needed %in% names(dat))) {
    cat("NOTE: ", path, " is missing one of ", paste(needed, collapse = ", "),
        " - skipping.\n", sep = "")
    return(NULL)
  }
  dat %>%
    filter(status == "ok") %>%
    transmute(
      fid = trimws(as.character(fid)),
      survey_year = year,
      n1 = suppressWarnings(as.numeric(n_SJ1_like)),
      n2 = suppressWarnings(as.numeric(n_SJ2_like)),
      n3 = suppressWarnings(as.numeric(n_SJ3_like)),
      n4 = suppressWarnings(as.numeric(n_SJ4_like))
    ) %>%
    rowwise() %>%
    mutate(
      n_total = sum(n1, n2, n3, n4, na.rm = TRUE),
      modal_sj = if (is.finite(n_total) && n_total > 0) which.max(c(n1, n2, n3, n4)) else NA_integer_
    ) %>%
    ungroup() %>%
    filter(!is.na(modal_sj)) %>%
    select(fid, survey_year, modal_sj, n_total)
}

modal_by_year <- purrr::map_dfr(survey_years, read_modal_year)

cat("Modal AE-SJ class computed for ", nrow(modal_by_year), " (site, year) pairs, years: ",
    paste(survey_years, collapse = ", "), "\n", sep = "")


# =============================================================================
# 4. JOIN AND CONFUSION MATRIX
# =============================================================================

compare <- field_sj %>%
  inner_join(modal_by_year, by = c("fid", "survey_year")) %>%
  mutate(
    field_sj_lab = factor(SJ_LABELS[field_sj], levels = SJ_LABELS),
    modal_sj_lab = factor(SJ_LABELS[modal_sj], levels = SJ_LABELS),
    agree = field_sj == modal_sj
  )

n_compared <- nrow(compare)
if (n_compared == 0) {
  stop(
    "No (site, year) pairs had both a field survey AND a matching annual CSV ",
    "- is ", ANNUAL_INPUT_DIR, " populated for the survey years (",
    paste(survey_years, collapse = ", "), ")?"
  )
}

readr::write_csv(
  compare %>% select(fid, survey_year, field_sj, modal_sj, n_pixels = n_total, agree),
  file.path(OUTPUT_DIR, "01_sj_score_validation_site_year_detail.csv")
)

confusion <- compare %>%
  count(field_sj_lab, modal_sj_lab, name = "n") %>%
  tidyr::complete(field_sj_lab, modal_sj_lab, fill = list(n = 0))

readr::write_csv(confusion, file.path(OUTPUT_DIR, "02_sj_score_validation_confusion_matrix.csv"))

overall_agreement <- mean(compare$agree)

# Cohen's kappa - agreement beyond what row/column marginals alone would
# predict by chance. Computed directly (no extra package) since this
# project's existing modules don't depend on one for it.
tab <- table(compare$field_sj_lab, compare$modal_sj_lab)
po <- sum(diag(tab)) / sum(tab)
pe <- sum((rowSums(tab) / sum(tab)) * (colSums(tab) / sum(tab)))
kappa <- (po - pe) / (1 - pe)

summary_stats <- tibble::tibble(
  n_site_years_compared = n_compared,
  n_distinct_sites = dplyr::n_distinct(compare$fid),
  survey_years = paste(survey_years, collapse = ";"),
  overall_agreement_pct = round(100 * overall_agreement, 1),
  cohens_kappa = round(kappa, 3)
)
readr::write_csv(summary_stats, file.path(OUTPUT_DIR, "03_sj_score_validation_summary.csv"))

cat("\n============================================================\n")
cat("SJ SCORE VALIDATION: field sj vs. modal AE-SJ class, same year\n")
cat("============================================================\n")
cat("Site-years compared: ", n_compared, " (", dplyr::n_distinct(compare$fid), " sites)\n", sep = "")
cat("Overall agreement: ", round(100 * overall_agreement, 1), "%\n", sep = "")
cat("Cohen's kappa: ", round(kappa, 3), "\n\n", sep = "")
print(tab)


# =============================================================================
# 5. PLOT
# =============================================================================

confusion_plot_dat <- confusion %>%
  group_by(field_sj_lab) %>%
  mutate(row_total = sum(n), row_pct = if_else(row_total > 0, n / row_total, NA_real_)) %>%
  ungroup()

p_confusion <- ggplot(confusion_plot_dat, aes(x = modal_sj_lab, y = field_sj_lab, fill = row_pct)) +
  geom_tile(colour = SKOG_THEME$surface, linewidth = 1.4) +
  geom_text(
    aes(label = ifelse(n == 0, "0", skog_num(n)),
        colour = row_pct > 0.5),
    fontface = "bold", size = 4.2
  ) +
  scale_fill_gradientn(
    colours = c(SKOG_THEME$none, SKOG_THEME$cell2, SKOG_THEME$cell3, SKOG_THEME$cell4),
    limits = c(0, 1), na.value = SKOG_THEME$none,
    labels = scales::percent_format(accuracy = 1),
    name = "Share of that\nfield-SJ row"
  ) +
  scale_colour_manual(values = c(`TRUE` = SKOG_THEME$surface, `FALSE` = SKOG_THEME$ink), guide = "none") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(
    title = "Field-surveyed SJ class vs. modal AE-SJ class, same site-year",
    subtitle = paste0(
      skog_num(n_compared), " site-years across ", dplyr::n_distinct(compare$fid),
      " MI sites, ", min(survey_years), "-", max(survey_years), " · ",
      round(100 * overall_agreement, 1), "% overall agreement · Cohen's κ = ", round(kappa, 3)
    ),
    x = "Modal AE-SJ class (this pipeline, that year)",
    y = "Field-surveyed SJ class (MI ground truth)",
    caption = wrap_caption(paste0(
      "Rows: the site's independently field-assessed SJ class in its survey year. Columns: the modal (most ",
      "common) per-pixel AE-SJ class this pipeline assigns that same site in that same year, from n_SJ1_like..",
      "n_SJ4_like in the annual CSVs. A perfectly diagonal matrix would mean the pixel classification exactly ",
      "recovers the field assessment every time it was checked against it directly."
    ))
  ) +
  theme_report() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.title.x = element_text(margin = margin(b = 8))
  )

print(p_confusion)

ggsave(
  file.path(OUTPUT_DIR, "04_sj_score_validation_confusion_matrix.png"),
  plot = p_confusion, width = 7.5, height = 6.5, dpi = 350, bg = "white"
)

cat("\nOutput directory:\n", OUTPUT_DIR, "\n", sep = "")
