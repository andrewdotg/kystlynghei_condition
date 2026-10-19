# =============================================================================
# GRAZING MODULE
# =============================================================================
#
# Effect of grazing (BT) on the AE-SJ trend.
#
# RUN ORDER: source in the same R session after
# R/modules/temporal_trajectories.R. It uses objects built there (mi_sf,
# ae_mi, development_excluded, state_trends, OUTPUT_DIR, find_bt_field(),
# trajectory_step_qc()) and the theme from R/theme.R (run_all.R does this).
# Outputs go to OUTPUT_DIR/grazing/, and those of
# run_grazing_narrow_comparison() to OUTPUT_DIR.
GRAZING_OUTPUT_DIR <- file.path(OUTPUT_DIR, "grazing")
dir.create(GRAZING_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
# =============================================================================


# =============================================================================
# ADD-ON: GRAZING MOSAIC / HETEROGENEITY DIAGNOSTICS
# =============================================================================
#
# Comparing mean trajectories asks whether grazing groups differ in their
# average AE-SJ level. These four diagnostics
# instead ask whether they differ in STRUCTURE - within-site spread,
# composition, and pixel-level movement - which the grazing/heterogeneity
# literature suggests can carry a signal even where group means look similar
# (spatial variability metrics keep differentiating grazing levels after
# mean-based metrics saturate; see e.g. UAV/lidar rangeland-heterogeneity
# studies). Read alongside a caveat worth keeping in view throughout: BT
# (beitetrykk) is not an independent animal-count measure - NiN's own
# documentation notes it is "svært vanskelig å tallfeste på en
# observatøruavhengig måte" and is instead read off observable vegetation
# and ground properties. So BT and the AE-SJ score are both, at root,
# readings of vegetation structure. Anything found here is consistent with a
# grazing effect, not independent proof of one - the burn/tiltak layer
# elsewhere in this project is a cleaner treatment variable if a more
# causal comparison is wanted.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Read BT directly from the original MI attributes
# -----------------------------------------------------------------------------

bt_field_grazing <- find_bt_field(mi_sf)

if (is.na(bt_field_grazing)) {
  stop("Could not locate the BT field in mi_sf.")
}


# Local parser that explicitly supports BT8.
parse_bt_grazing <- function(x) {
  
  z <- trimws(as.character(x))
  
  out <- suppressWarnings(as.integer(z))
  
  missing <- is.na(out)
  
  if (any(missing)) {
    
    extracted <- stringr::str_extract(
      z[missing],
      "(?<![0-9])[1-6](?![0-9])"
    )
    
    out[missing] <- suppressWarnings(
      as.integer(extracted)
    )
  }
  
  out
}
grazing_lookup <- mi_sf %>%
  sf::st_drop_geometry() %>%
  transmute(
    fid = trimws(as.character(fid)),
    bt_grazing = parse_bt_grazing(.data[[bt_field_grazing]])
  ) %>%
  distinct(
    fid,
    .keep_all = TRUE
  ) %>%
  mutate(
    
    grazing_group = case_when(
      
      bt_grazing %in% c(2L, 3L) ~
        "God",
      
      bt_grazing %in% c(1L, 4L) ~
        "Moderat",
      
      bt_grazing == 5L ~
        "Dårlig",
      
      bt_grazing == 6L ~
        "Svært dårlig",
      
      TRUE ~ NA_character_
    ),
    
    grazing_group = factor(
      grazing_group,
      levels = c(
        "God",
        "Moderat",
        "Dårlig",
        "Svært dårlig"
      )
    )
  )

# -----------------------------------------------------------------------------
# 2. Build grazing-analysis annual dataset
# -----------------------------------------------------------------------------
#
# Start again from ae_mi so BT5 and BT8 have not been removed.
#
# Persistent-development sites are excluded as in the main analysis.
# -----------------------------------------------------------------------------

grazing_dat <- ae_mi %>%
  
  select(
    -any_of(
      c(
        "bt"
      )
    )
  ) %>%
  
  left_join(
    grazing_lookup,
    by = "fid"
  ) %>%
  
  filter(
    !is.na(fid),
    fid != "",
    !is.na(grazing_group),
    
    !fid %in% development_excluded,
    
    is.finite(AE_SJ_p50),
    is.finite(AE_SJ_count),
    AE_SJ_count >= MIN_PIXELS_STATE
  )


# -----------------------------------------------------------------------------
# 3. Apply the same annual-step QC used in the main temporal analysis
# -----------------------------------------------------------------------------

grazing_trajectory_qc <- trajectory_step_qc(
  grazing_dat
)


grazing_bad_step_fids <- grazing_trajectory_qc %>%
  filter(
    implausible_step
  ) %>%
  pull(
    fid
  )


grazing_dat <- grazing_dat %>%
  filter(
    !fid %in% grazing_bad_step_fids
  )


# -----------------------------------------------------------------------------
# 4. Summarise annual AE-SJ by grazing group
# -----------------------------------------------------------------------------
#
# Each polygon contributes one annual site-median AE-SJ value.
# We then calculate the mean across polygons in each grazing group.
#
# CI describes variation in the estimated group mean across sites.
# -----------------------------------------------------------------------------

grazing_trajectory <- grazing_dat %>%
  
  group_by(
    grazing_group,
    year
  ) %>%
  
  summarise(
    
    n_fids =
      n_distinct(
        fid
      ),
    
    mean_AE_SJ =
      mean(
        AE_SJ_p50,
        na.rm = TRUE
      ),
    
    sd_AE_SJ =
      sd(
        AE_SJ_p50,
        na.rm = TRUE
      ),
    
    se_AE_SJ =
      sd_AE_SJ /
      sqrt(
        n_fids
      ),
    
    ci_lower =
      mean_AE_SJ -
      1.96 * se_AE_SJ,
    
    ci_upper =
      mean_AE_SJ +
      1.96 * se_AE_SJ,
    
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 5. Group sizes
# -----------------------------------------------------------------------------

grazing_n <- grazing_dat %>%
  distinct(
    fid,
    grazing_group
  ) %>%
  count(
    grazing_group,
    name = "n_fids"
  )


print(
  grazing_n,
  n = Inf
)


readr::write_csv(
  grazing_trajectory,
  file.path(
    OUTPUT_DIR,
    "15_AE_SJ_annual_trajectory_by_grazing_condition.csv"
  )
)


# -----------------------------------------------------------------------------
# 6. Colours
# -----------------------------------------------------------------------------
#
# Same good -> poor visual language used elsewhere in the analysis.
# -----------------------------------------------------------------------------

GRAZING_COLOURS <- c(
  
  "God" =
    SKOG_THEME$g3,
  
  "Moderat" =
    SKOG_THEME$g1,
  
  "Dårlig" =
    SKOG_THEME$p1,
  
  "Svært dårlig" =
    SKOG_THEME$p4
)


# -----------------------------------------------------------------------------
# 7. Labels including sample size
# -----------------------------------------------------------------------------

grazing_labels <- grazing_n %>%
  mutate(
    label = paste0(
      grazing_group,
      " (n = ",
      skog_num(n_fids),
      ")"
    )
  )


grazing_label_values <- setNames(
  grazing_labels$label,
  as.character(
    grazing_labels$grazing_group
  )
)


# =============================================================================
# ADD-ON: NO GRAZING (BT1) VS SOME GRAZING (BT2-3) - NARROW SUCCESSION TEST
# =============================================================================
#
# A single, focused hypothesis test: does the overall AE-SJ succession trend
# (sj_change) differ between sites with no observed grazing sign (BT1)
# and sites with light-to-moderate grazing (BT2-3)? Deliberately restricted
# to the continuous trend, not composition/entropy - composition is the
# metric most exposed to the confound noted earlier: physical management
# carried out specifically to sustain grazing (burning, clearing) shapes
# composition directly, on top of whatever grazing itself does. The trend is
# one step further downstream, though still not a fully independent outcome
# - see the caveat in the plot caption below.
#
# Unlike the wider grazing add-on above, this reuses state_trends directly:
# BT1-4 are never excluded from the main analysis (only BT5/BT6 are, via
# EXCLUDE_BT), so this population is already sitting there. No need to
# rebuild a parallel trajectory pipeline from ae_mi for this narrower test.
# =============================================================================

grazing_narrow_dat <- state_trends %>%
  filter(direction_eligible, bt %in% c(1L, 2L, 3L)) %>%
  mutate(
    grazing_narrow = case_when(
      bt == 1L ~ "No grazing (BT1)",
      bt %in% c(2L, 3L) ~ "Grazing present (BT2-3)"
    ),
    grazing_narrow = factor(grazing_narrow, levels = c("No grazing (BT1)", "Grazing present (BT2-3)"))
  )

grazing_narrow_n <- grazing_narrow_dat %>% count(grazing_narrow, name = "n_fids")
cat("\nNo grazing (BT1) vs grazing present (BT2-3): sample sizes\n")
print(grazing_narrow_n)

# Cross-check against tilstand: BT2/3 map to "God" and BT1 maps to "Moderat"
# in the official primary-variable table, so some of any difference found
# below is expected to double as a tilstand difference - this table doesn't
# correct for that, just shows how much overlap there is so you can judge
# whether a tilstand-stratified follow-up is worth asking for.
grazing_narrow_by_tilstand <- grazing_narrow_dat %>%
  filter(!is.na(tilstand)) %>%
  count(tilstand, grazing_narrow, name = "n_fids")
cat("\nSame comparison, cross-tabulated against tilstand:\n")
print(grazing_narrow_by_tilstand, n = Inf)
readr::write_csv(grazing_narrow_by_tilstand, file.path(GRAZING_OUTPUT_DIR, "20_no_grazing_vs_grazing_by_tilstand_n.csv"))

# -----------------------------------------------------------------------------
# Statistical test
# -----------------------------------------------------------------------------
# Wilcoxon rank-sum (Mann-Whitney) rather than a t-test: sj_change is
# not close to normal (bounded near the SJ1/SJ4 reference endpoints, long
# tails from the small number of strongly-changing sites), so a rank-based
# test is the safer default, the same as for the other grazing comparisons.
# conf.int = TRUE also returns the Hodges-Lehmann
# estimate: a robust, non-parametric estimate of the typical (median-level)
# difference between the two groups, with a confidence interval - reported
# directly on the figure below rather than just the p-value, since a p-value
# alone doesn't say how big the difference actually is.
#
# This treats each site as one independent observation, which it is (one
# summary trend value per site, not repeated annual rows) - but it does not
# account for possible spatial autocorrelation between nearby sites. Worth
# account for possible spatial autocorrelation between nearby sites; that
# would need a mixed model or a spatially-blocked permutation test.

grazing_narrow_test <- wilcox.test(
  sj_change ~ grazing_narrow,
  data = grazing_narrow_dat,
  conf.int = TRUE
)

print(grazing_narrow_test)

grazing_narrow_summary <- grazing_narrow_dat %>%
  group_by(grazing_narrow) %>%
  summarise(
    n_fids = n(),
    median_trend_change = median(sj_change, na.rm = TRUE),
    mean_trend_change = mean(sj_change, na.rm = TRUE),
    q25 = quantile(sj_change, 0.25, na.rm = TRUE),
    q75 = quantile(sj_change, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

grazing_narrow_test_summary <- tibble::tibble(
  comparison = "No grazing (BT1) vs grazing present (BT2-3)",
  n_no_grazing = sum(grazing_narrow_dat$grazing_narrow == "No grazing (BT1)"),
  n_grazing_present = sum(grazing_narrow_dat$grazing_narrow == "Grazing present (BT2-3)"),
  wilcoxon_W = unname(grazing_narrow_test$statistic),
  p_value = grazing_narrow_test$p.value,
  hodges_lehmann_shift = unname(grazing_narrow_test$estimate),
  ci_low = grazing_narrow_test$conf.int[1],
  ci_high = grazing_narrow_test$conf.int[2]
)

readr::write_csv(grazing_narrow_summary, file.path(GRAZING_OUTPUT_DIR, "20_no_grazing_vs_grazing_trend_summary.csv"))
readr::write_csv(grazing_narrow_test_summary, file.path(GRAZING_OUTPUT_DIR, "20_no_grazing_vs_grazing_test_result.csv"))

cat("\nSummary:\n")
print(grazing_narrow_summary)
cat("\nTest result:\n")
print(grazing_narrow_test_summary)


# -----------------------------------------------------------------------------
# Visualisation
# -----------------------------------------------------------------------------
# Red/blue kept consistent with the rest of the deck's direction language
# (succession-leaning = red, restoration/stability-leaning = blue), applied
# here to the two grazing groups rather than to triage tiers. Dashed lines
# are the same moderate/strong triage thresholds used everywhere else in
# this analysis, so the room can see where each group's distribution sits
# relative to a scale they've already seen.

grazing_narrow_labels <- setNames(
  paste0(grazing_narrow_n$grazing_narrow, "\n(n = ", skog_num(grazing_narrow_n$n_fids), ")"),
  as.character(grazing_narrow_n$grazing_narrow)
)

GRAZING_NARROW_COLOURS <- c(
  "No grazing (BT1)" = SKOG_THEME$p3,
  "Grazing present (BT2-3)" = SKOG_THEME$g3
)

p_grazing_narrow <- ggplot(grazing_narrow_dat, aes(x = grazing_narrow, y = sj_change, fill = grazing_narrow)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = SKOG_THEME$ink3) +
  geom_hline(
    yintercept = c(-large_threshold, -moderate_threshold, moderate_threshold, large_threshold),
    linetype = "dashed", linewidth = 0.3, colour = SKOG_THEME$rule2
  ) +
  geom_violin(width = 0.7, alpha = 0.5, colour = NA, trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.15, colour = SKOG_THEME$ink2, linewidth = 0.4, fill = SKOG_THEME$surface) +
  scale_fill_manual(values = GRAZING_NARROW_COLOURS, guide = "none") +
  scale_x_discrete(labels = grazing_narrow_labels) +
  labs(
    title = "Does succession differ between ungrazed and grazed sites?",
    subtitle = "Estimated AE-SJ change over the monitoring period, BT1 (no grazing sign) vs BT2-3 (grazing present)",
    x = NULL,
    y = "Estimated AE-SJ change over the monitoring period",
    caption = wrap_caption(paste0(
      "Wilcoxon rank-sum test: p = ", signif(grazing_narrow_test$p.value, 2),
      "; estimated shift (Hodges-Lehmann) = ", signif(grazing_narrow_test$estimate, 3),
      " AE-SJ units, 95% CI [", signif(grazing_narrow_test$conf.int[1], 3), ", ",
      signif(grazing_narrow_test$conf.int[2], 3), "]. Dashed lines mark the same moderate/strong triage ",
      "thresholds used elsewhere in this analysis. BT2-3 correspond to 'God' and BT1 to 'Moderat' in the ",
      "official tilstand framework, so part of any difference here will overlap with a tilstand difference too."
    ))
  ) +
  theme_report()

# =============================================================================
# ADD-ON: SENSITIVITY OF THE GRAZING COMPARISON TO BASELINE SJ CLASS
# =============================================================================
#
# Restricts the BT1-vs-BT2/3 comparison to sites that STARTED in SJ1/SJ2
# territory - where sheep grazing/burning is the management tool the
# framework diagram actually pairs with those stages, as opposed to SJ3
# (shrub removal) and SJ4 (tree removal), where grazing pressure plausibly
# stops being the operative driver of the AE-SJ trend at all.
#
# Filters on the INDEPENDENTLY FIELD-SURVEYED 7RA-SJ class, not on sj_base
# (the AE-SJ model's own baseline estimate). Using the model's own score to
# decide who enters a test of the model's own trend risks the same kind of
# shared-measurement-basis issue already flagged for BT vs SJ: both are
# reads of vegetation structure, so a systematic AE-SJ offset between
# grazed and ungrazed sites at the same TRUE field stage could admit subtly
# different populations by BT group for reasons unrelated to actual
# successional stage. The field record avoids that - only the field survey
# decides who's in the population; the AE-SJ trend inside it is untouched.
#
# Also filters on BASELINE class, not on whether a site ever entered
# SJ3/SJ4 during monitoring. Filtering on the full trajectory would
# specifically remove the most successional sites - the ones this test is
# trying to detect a difference in - and bias the result toward finding
# nothing regardless of whether a real effect exists. Starting-state
# filtering restricts which sites enter the test, not what they're allowed
# to do once they're in it.
# =============================================================================

find_sj_field <- function(x) {
  nm <- names(x)
  if ("X7RA.SJ" %in% nm) return("X7RA.SJ")
  key <- tolower(gsub("[^a-zA-Z0-9]+", "", nm))
  preferred_keys <- c("x7rasj", "7rasj", "rasj", "sj", "rapidsuccession", "rasuksesjon")
  idx <- which(key %in% preferred_keys)
  if (length(idx) == 0) return(NA_character_)
  nm[idx[1]]
}

parse_sj <- function(x) {
  z <- trimws(as.character(x))
  out <- suppressWarnings(as.integer(z))
  missing <- is.na(out)
  if (any(missing)) {
    extracted <- stringr::str_extract(z[missing], "(?<![0-9])[1-4](?![0-9])")
    out[missing] <- suppressWarnings(as.integer(extracted))
  }
  out
}

sj_field_mi <- find_sj_field(mi_sf)

if (is.na(sj_field_mi)) {
  stop(
    "Could not find a field-surveyed 7RA-SJ column in mi_sf. Run names(mi_sf) ",
    "and check the actual column name, then add it to find_sj_field()'s ",
    "preferred_keys above."
  )
}

field_sj_lookup <- mi_sf %>%
  sf::st_drop_geometry() %>%
  transmute(fid = trimws(as.character(fid)), field_sj = parse_sj(.data[[sj_field_mi]])) %>%
  distinct(fid, .keep_all = TRUE)

cat("\nField-surveyed 7RA-SJ: coverage and distribution\n")
print(field_sj_lookup %>% count(field_sj, name = "n_fids"))

grazing_narrow_dat <- grazing_narrow_dat %>%
  select(-any_of("field_sj")) %>%
  left_join(field_sj_lookup, by = "fid")

cat(
  "\nSites in the BT1/BT2-3 comparison with no matched field_sj value: ",
  sum(is.na(grazing_narrow_dat$field_sj)), " / ", nrow(grazing_narrow_dat), "\n",
  sep = ""
)

grazing_narrow_field_sj12  <- grazing_narrow_dat %>% filter(!is.na(field_sj), field_sj %in% c(1L, 2L))
grazing_narrow_field_sj123 <- grazing_narrow_dat %>% filter(!is.na(field_sj), field_sj %in% c(1L, 2L, 3L))

cat(
  "\nSample sizes as the field-surveyed baseline restriction tightens:\n",
  "  All BT1-3 sites (any baseline SJ):          ", nrow(grazing_narrow_dat), "\n",
  "  Field-surveyed SJ1-SJ3 only (SJ4 excluded): ", nrow(grazing_narrow_field_sj123), "\n",
  "  Field-surveyed SJ1-SJ2 only:                ", nrow(grazing_narrow_field_sj12), "\n",
  sep = ""
)

# Reusable version of the same test/plot built for the unrestricted
# comparison, so this isn't tripling the same 40 lines of code for each
# restricted population.
run_grazing_narrow_comparison <- function(dat, title_suffix, file_suffix, save_plot = TRUE) {
  n_tab <- dat %>% count(grazing_narrow, name = "n_fids")
  
  test <- wilcox.test(sj_change ~ grazing_narrow, data = dat, conf.int = TRUE)
  
  summary_tab <- dat %>%
    group_by(grazing_narrow) %>%
    summarise(
      n_fids = n(),
      median_trend_change = median(sj_change, na.rm = TRUE),
      mean_trend_change = mean(sj_change, na.rm = TRUE),
      q25 = quantile(sj_change, 0.25, na.rm = TRUE),
      q75 = quantile(sj_change, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  test_summary <- tibble::tibble(
    comparison = paste0("No grazing (BT1) vs grazing present (BT2-3), ", title_suffix),
    n_no_grazing = sum(dat$grazing_narrow == "No grazing (BT1)"),
    n_grazing_present = sum(dat$grazing_narrow == "Grazing present (BT2-3)"),
    wilcoxon_W = unname(test$statistic),
    p_value = test$p.value,
    hodges_lehmann_shift = unname(test$estimate),
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2]
  )
  
  labels <- setNames(
    paste0(n_tab$grazing_narrow, "\n(n = ", skog_num(n_tab$n_fids), ")"),
    as.character(n_tab$grazing_narrow)
  )
  
  p <- ggplot(dat, aes(x = grazing_narrow, y = sj_change, fill = grazing_narrow)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = SKOG_THEME$ink3) +
    geom_hline(
      yintercept = c(-large_threshold, -moderate_threshold, moderate_threshold, large_threshold),
      linetype = "dashed", linewidth = 0.3, colour = SKOG_THEME$rule2
    ) +
    geom_violin(width = 0.7, alpha = 0.5, colour = NA, trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.alpha = 0.15, colour = SKOG_THEME$ink2, linewidth = 0.4, fill = SKOG_THEME$surface) +
    scale_fill_manual(values = GRAZING_NARROW_COLOURS, guide = "none") +
    scale_x_discrete(labels = labels) +
    labs(
      title = "Does succession differ between ungrazed and grazed sites?",
      subtitle = paste0("Estimated AE-SJ change over the monitoring period - ", title_suffix),
      x = NULL,
      y = "Estimated AE-SJ change over the monitoring period",
      caption = wrap_caption(paste0(
        "Wilcoxon rank-sum test: p = ", signif(test$p.value, 2),
        "; estimated shift (Hodges-Lehmann) = ", signif(test$estimate, 3),
        " AE-SJ units, 95% CI [", signif(test$conf.int[1], 3), ", ", signif(test$conf.int[2], 3), "]. ",
        "Restricted to sites whose independently field-surveyed 7RA-SJ class falls in this range at the time ",
        "of the MI survey - sites are still free to move beyond it over the monitoring period, as measured by ",
        "AE-SJ. The field survey year may predate the 2018 monitoring baseline. BT2-3 correspond to 'God' and ",
        "BT1 to 'Moderat' in the official tilstand framework, so part of any difference here will overlap ",
        "with a tilstand difference too."
      ))
    ) +
    theme_report()
  
  readr::write_csv(summary_tab, file.path(OUTPUT_DIR, paste0("21_no_grazing_vs_grazing_trend_summary_", file_suffix, ".csv")))
  readr::write_csv(test_summary, file.path(OUTPUT_DIR, paste0("21_no_grazing_vs_grazing_test_result_", file_suffix, ".csv")))
  # save_plot = FALSE skips writing the figure; the plot is returned either
  # way (result$plot).
  if (save_plot) {
    ggsave(
      filename = file.path(OUTPUT_DIR, paste0("21_no_grazing_vs_grazing_trend_", file_suffix, ".png")),
      plot = p, width = 8, height = 7, dpi = 350, bg = "white"
    )
  }

  list(test = test, summary = summary_tab, test_summary = test_summary, plot = p)
}

result_field_sj123 <- run_grazing_narrow_comparison(
  grazing_narrow_field_sj123, "field-surveyed SJ1-SJ3 sites, SJ4 excluded", "field_sj4_excluded",
  save_plot = FALSE
)
result_field_sj12 <- run_grazing_narrow_comparison(
  grazing_narrow_field_sj12, "field-surveyed SJ1-SJ2 sites only", "field_sj12_only",
  save_plot = FALSE
)


# All three versions side by side, so it's visible at a glance whether the
# effect strengthens, weakens, or holds steady as the field-surveyed
# baseline restriction tightens.
grazing_narrow_sensitivity <- bind_rows(
  grazing_narrow_test_summary %>% mutate(population = "All BT1-3 sites (any baseline SJ)"),
  result_field_sj123$test_summary %>% mutate(population = "Field-surveyed SJ1-SJ3 only (SJ4 excluded)"),
  result_field_sj12$test_summary %>% mutate(population = "Field-surveyed SJ1-SJ2 only")
) %>%
  select(population, everything())

readr::write_csv(grazing_narrow_sensitivity, file.path(GRAZING_OUTPUT_DIR, "21_no_grazing_vs_grazing_sensitivity_to_baseline_sj.csv"))

cat("\nEffect size across field-surveyed baseline-SJ restrictions:\n")
print(grazing_narrow_sensitivity)



# =============================================================================
# ADD-ON: TILSTAND-STRATIFIED GRAZING COMPARISON
# =============================================================================
#
# The pooled BT1-vs-BT2/3 result is confounded with tilstand: "God" tilstand
# is 98.3% grazed sites in this dataset, while "Dårlig" and "Svært redusert"
# are ~81% ungrazed - BT is one of the two primary variables that jointly
# define tilstand, so this isn't incidental. Re-running the same comparison
# SEPARATELY within each tilstand class isolates whatever remains once
# starting condition is held fixed - a test of whether grazing explains
# something tilstand doesn't already explain, rather than restating it.
#
# Built on the field-surveyed SJ1/SJ2 population (grazing_narrow_field_sj12)
# established above - the cleanest population from the earlier sensitivity
# check - not a new branch to reconcile against it later.
# =============================================================================

MIN_STRATUM_N <- 30  # same floor MIN_ENDPOINT_POINTS uses in the AE-SJ model script

grazing_tilstand_n <- grazing_narrow_field_sj12 %>%
  filter(!is.na(tilstand)) %>%
  count(tilstand, grazing_narrow, name = "n_fids") %>%
  tidyr::pivot_wider(names_from = grazing_narrow, values_from = n_fids, values_fill = 0)

cat("\nSample sizes by tilstand, field-surveyed SJ1/SJ2 population:\n")
print(grazing_tilstand_n)

grazing_tilstand_eligible <- grazing_tilstand_n %>%
  filter(
    `No grazing (BT1)` >= MIN_STRATUM_N,
    `Grazing present (BT2-3)` >= MIN_STRATUM_N
  ) %>%
  pull(tilstand) %>%
  as.character()

grazing_tilstand_excluded <- setdiff(as.character(grazing_tilstand_n$tilstand), grazing_tilstand_eligible)

cat(
  "\nTilstand classes tested (>= ", MIN_STRATUM_N, " sites on both sides): ",
  paste(grazing_tilstand_eligible, collapse = ", "), "\n",
  "Tilstand classes excluded (too few sites on at least one side): ",
  if (length(grazing_tilstand_excluded) > 0) paste(grazing_tilstand_excluded, collapse = ", ") else "none",
  "\n", sep = ""
)

run_stratum_test <- function(dat, stratum_label) {
  test <- wilcox.test(sj_change ~ grazing_narrow, data = dat, conf.int = TRUE)
  tibble::tibble(
    tilstand = stratum_label,
    n_no_grazing = sum(dat$grazing_narrow == "No grazing (BT1)"),
    n_grazing_present = sum(dat$grazing_narrow == "Grazing present (BT2-3)"),
    wilcoxon_W = unname(test$statistic),
    p_value = test$p.value,
    hodges_lehmann_shift = unname(test$estimate),
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2]
  )
}

grazing_tilstand_results <- purrr::map_dfr(
  grazing_tilstand_eligible,
  function(lvl) {
    run_stratum_test(
      grazing_narrow_field_sj12 %>% filter(as.character(tilstand) == lvl),
      lvl
    )
  }
)

# Pooled estimate shown alongside for reference - result_field_sj12 is the
# same test from the sensitivity check above, not rerun here.
grazing_tilstand_results <- bind_rows(
  grazing_tilstand_results,
  tibble::tibble(
    tilstand = "All tilstand (pooled)",
    n_no_grazing = result_field_sj12$test_summary$n_no_grazing,
    n_grazing_present = result_field_sj12$test_summary$n_grazing_present,
    wilcoxon_W = result_field_sj12$test_summary$wilcoxon_W,
    p_value = result_field_sj12$test_summary$p_value,
    hodges_lehmann_shift = result_field_sj12$test_summary$hodges_lehmann_shift,
    ci_low = result_field_sj12$test_summary$ci_low,
    ci_high = result_field_sj12$test_summary$ci_high
  )
) %>%
  mutate(
    ci_excludes_zero = ci_low > 0 | ci_high < 0,
    tilstand = factor(tilstand, levels = c(rev(TILSTAND_LEVELS), "All tilstand (pooled)"))
  )

readr::write_csv(
  grazing_tilstand_results %>% mutate(tilstand = as.character(tilstand)),
  file.path(GRAZING_OUTPUT_DIR, "23_no_grazing_vs_grazing_by_tilstand.csv")
)

cat("\nGrazing comparison, stratified by tilstand (field-surveyed SJ1/SJ2 population):\n")
print(grazing_tilstand_results)

p_grazing_tilstand <- ggplot(
  grazing_tilstand_results,
  aes(x = hodges_lehmann_shift, y = tilstand, colour = ci_excludes_zero)
) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = SKOG_THEME$ink3) +
  geom_vline(
    xintercept = c(-moderate_threshold, moderate_threshold),
    linetype = "dashed", linewidth = 0.3, colour = SKOG_THEME$rule2
  ) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.6) +
  geom_point(size = 2.6) +
  geom_text(
    aes(label = paste0("n = ", skog_num(n_no_grazing), " vs ", skog_num(n_grazing_present))),
    hjust = 0, nudge_y = 0.28, size = gg_text_size("eyebrow"), colour = SKOG_THEME$ink3, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("TRUE" = SKOG_THEME$p3, "FALSE" = SKOG_THEME$neutral2),
    labels = c("TRUE" = "95% CI excludes zero", "FALSE" = "95% CI includes zero"),
    name = NULL
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Does the grazing difference hold up within each tilstand class?",
    subtitle = "Hodges-Lehmann shift in AE-SJ trend (no grazing minus grazing present), field-surveyed SJ1/SJ2 sites only",
    x = "Estimated shift in AE-SJ trend (AE-SJ units, no grazing vs grazing present)",
    y = NULL,
    caption = wrap_caption(paste0(
      "Each row is the same BT1-vs-BT2/3 Wilcoxon comparison as the pooled result, run separately within that ",
      "tilstand class. A negative shift means BT1 (no grazing) sites trended further toward SJ4 than BT2-3 ",
      "(grazing present) sites, holding starting tilstand fixed. Dashed lines mark the moderate triage threshold. ",
      "Tilstand classes with fewer than ", MIN_STRATUM_N, " sites on either side are not shown - see console ",
      "output or the CSV for which were excluded and why."
    ))
  ) +
  theme_report() +
  theme(legend.position = "bottom", plot.margin = margin(18, 22, 14, 90))



# =============================================================================
# ADD-ON: SINGLE CLEAN TEST OF THE STATED HYPOTHESIS
# =============================================================================
# Grazed sites (BT2-3) move more slowly through SJ1->SJ2->SJ3 succession
# than ungrazed sites (BT1) - tested as ANNUAL RATE (slope), not total
# change over a possibly-varying number of years.
#
# Population: field-surveyed baseline SJ1, SJ2 or SJ3 (SJ4 excluded), and
# tilstand God or Moderat only (Dårlig/Svært redusert excluded, since their
# condition score is itself partly defined by BT).
#
# One test, one population. Requires field_sj joined onto
# grazing_narrow_dat (field-surveyed filtering above).
# =============================================================================

grazing_slope_dat <- grazing_narrow_dat %>%
  filter(
    !is.na(field_sj), field_sj %in% c(1L, 2L, 3L),
    !is.na(tilstand), tilstand %in% c("God", "Moderat")
  )

cat(
  "\nPopulation: field-surveyed SJ1-SJ3, tilstand God/Moderat only\n",
  "  No grazing (BT1):        ", sum(grazing_slope_dat$grazing_narrow == "No grazing (BT1)"), "\n",
  "  Grazing present (BT2-3): ", sum(grazing_slope_dat$grazing_narrow == "Grazing present (BT2-3)"), "\n",
  sep = ""
)

grazing_slope_test <- wilcox.test(sj_slope ~ grazing_narrow, data = grazing_slope_dat, conf.int = TRUE)
print(grazing_slope_test)

grazing_slope_summary <- grazing_slope_dat %>%
  group_by(grazing_narrow) %>%
  summarise(
    n_fids = n(),
    median_slope = median(sj_slope, na.rm = TRUE),
    mean_slope = mean(sj_slope, na.rm = TRUE),
    q25 = quantile(sj_slope, 0.25, na.rm = TRUE),
    q75 = quantile(sj_slope, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

grazing_slope_test_summary <- tibble::tibble(
  comparison = "No grazing (BT1) vs grazing present (BT2-3), field-surveyed SJ1-SJ3, tilstand God/Moderat",
  n_no_grazing = sum(grazing_slope_dat$grazing_narrow == "No grazing (BT1)"),
  n_grazing_present = sum(grazing_slope_dat$grazing_narrow == "Grazing present (BT2-3)"),
  wilcoxon_W = unname(grazing_slope_test$statistic),
  p_value = grazing_slope_test$p.value,
  hodges_lehmann_shift = unname(grazing_slope_test$estimate),
  ci_low = grazing_slope_test$conf.int[1],
  ci_high = grazing_slope_test$conf.int[2]
)

readr::write_csv(grazing_slope_summary, file.path(GRAZING_OUTPUT_DIR, "24_grazing_slope_summary.csv"))
readr::write_csv(grazing_slope_test_summary, file.path(GRAZING_OUTPUT_DIR, "24_grazing_slope_test_result.csv"))

cat("\nSummary:\n"); print(grazing_slope_summary)
cat("\nTest result:\n"); print(grazing_slope_test_summary)

grazing_slope_n <- grazing_slope_dat %>% count(grazing_narrow, name = "n_fids")
grazing_slope_labels <- setNames(
  paste0(grazing_slope_n$grazing_narrow, "\n(n = ", skog_num(grazing_slope_n$n_fids), ")"),
  as.character(grazing_slope_n$grazing_narrow)
)

p_grazing_slope <- ggplot(grazing_slope_dat, aes(x = grazing_narrow, y = sj_slope, fill = grazing_narrow)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = SKOG_THEME$ink3) +
  geom_violin(width = 0.7, alpha = 0.5, colour = NA, trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.15, colour = SKOG_THEME$ink2, linewidth = 0.4, fill = SKOG_THEME$surface) +
  scale_fill_manual(values = GRAZING_NARROW_COLOURS, guide = "none") +
  scale_x_discrete(labels = grazing_slope_labels) +
  labs(
    title = "Do grazed sites move more slowly through SJ1-SJ2-SJ3 succession?",
    subtitle = "Annual AE-SJ slope, field-surveyed baseline SJ1-SJ3, tilstand God/Moderat only",
    x = NULL,
    y = "AE-SJ slope (units per year)\nMore negative = faster movement toward SJ4",
    caption = wrap_caption(paste0(
      "Wilcoxon rank-sum test: p = ", signif(grazing_slope_test$p.value, 2),
      "; estimated shift (Hodges-Lehmann) = ", signif(grazing_slope_test$estimate, 3),
      " AE-SJ units/year, 95% CI [", signif(grazing_slope_test$conf.int[1], 3), ", ",
      signif(grazing_slope_test$conf.int[2], 3), "]. Restricted to sites field-surveyed as SJ1, SJ2 or SJ3 at ",
      "baseline (SJ4 excluded) and mapped as God or Moderat tilstand (Dårlig/Svært redusert excluded)."
    ))
  ) +
  theme_report()


# =============================================================================
# ADD-ON: MEDIAN COMPARISON, ISOLATED FROM THE OVERLAPPING DISTRIBUTIONS
# =============================================================================
#
# The violin plot above is an honest picture of a real result: the within-
# group spread genuinely dwarfs the between-group shift, so the effect will
# never look dramatic there no matter how real it is. This isolates just
# the central-tendency comparison - each group's median with a bootstrap
# confidence interval - so the size of the shift is what's on screen,
# rather than competing visually with the overlap in the tails.
# =============================================================================

set.seed(42)
N_BOOT <- 2000

bootstrap_median_ci <- function(x, n_boot = N_BOOT) {
  boot_medians <- replicate(n_boot, median(sample(x, length(x), replace = TRUE)))
  tibble::tibble(
    median = median(x),
    ci_low = quantile(boot_medians, 0.025, names = FALSE),
    ci_high = quantile(boot_medians, 0.975, names = FALSE)
  )
}

grazing_median_ci <- grazing_slope_dat %>%
  group_by(grazing_narrow) %>%
  summarise(bootstrap_median_ci(sj_slope), n_fids = n(), .groups = "drop")

readr::write_csv(grazing_median_ci, file.path(GRAZING_OUTPUT_DIR, "24_grazing_slope_median_bootstrap_ci.csv"))

cat("\nBootstrapped median AE-SJ slope, with 95% CI:\n")
print(grazing_median_ci)

grazing_median_labels <- setNames(
  paste0(grazing_median_ci$grazing_narrow, "\n(n = ", skog_num(grazing_median_ci$n_fids), ")"),
  as.character(grazing_median_ci$grazing_narrow)
)

# Median AE-SJ change with its 95 % confidence interval by grazing group
# (the figure used in the report).
p_grazing_median_ci <- ggplot(grazing_median_ci, aes(x = grazing_narrow, y = median, colour = grazing_narrow)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = SKOG_THEME$ink3) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), linewidth = 1.1) +
  geom_point(size = 4) +
  scale_colour_manual(values = GRAZING_NARROW_COLOURS, guide = "none") +
  scale_x_discrete(labels = grazing_median_labels) +
  labs(
    title = "Grazed vs ungrazed sites: median succession rate",
    subtitle = "Median annual AE-SJ slope with bootstrap 95% CI",
    x = NULL,
    y = "Median AE-SJ slope (units per year)\nMore negative = faster movement toward SJ4"
  ) +
  theme_report()

print(p_grazing_median_ci)

ggsave(
  filename = file.path(GRAZING_OUTPUT_DIR, "24_grazing_slope_median_ci.png"),
  plot = p_grazing_median_ci, width = 7, height = 6, dpi = 350, bg = "white"
)


