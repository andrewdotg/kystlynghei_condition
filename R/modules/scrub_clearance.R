# =============================================================================
# KYSTLYNGHEI MI - SCRUB/TREE CLEARANCE DETECTION (AE-SJ TRAJECTORY METHOD)
# =============================================================================
#
# The MI analogue of R/naturTiltak/naturTiltak_scrub.R's validation: instead
# of checking the AE-SJ method against independently reported naturTiltak
# clearance actions, this identifies MI kystlynghei sites where the AE-SJ
# trajectory itself looks like a scrub/tree-clearance response - SJ3-like or
# SJ4-like pixels moving to SJ1-like/SJ2-like structure - across the whole
# mapped MI population, not just tiltak-logged sites.
#
# RUN ORDER: source in the same R session right after
# R/modules/temporal_trajectories.R. This module uses objects built there
# (state_dat, ref_dat, recovery_sites) and the plotting helpers from
# R/theme.R that it loads (run_all.R does this).
# =============================================================================

# Outputs go to outputs/scrub_clearance/.
OUTPUT_DIR <- file.path(OUTPUTS_ROOT, "scrub_clearance")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative"), recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# ADD-ON: RESTORATIVE SITES WITH SJ3/SJ4 -> SJ1/SJ2 PIXEL TRANSITIONS
# =============================================================================
#
# PURPOSE
# -------
# Identify restorative MI sites where at least some matched pixels have moved:
#
#   SJ3-like -> SJ1-like or SJ2-like
#   SJ4-like -> SJ1-like or SJ2-like
#
# This is intended as the MI analogue of the independently validated
# scrub/tree-clearance response:
#
#   SJ3-like source structure = consistent with scrub / woody-regrowth removal
#   SJ4-like source structure = consistent with tree / late-successional removal
#
# IMPORTANT
# ---------
# Site selection has TWO requirements:
#
#   1. The polygon is already independently classified as a Moderate or Strong
#      restorative trajectory in the main temporal analysis.
#
#   2. At least one matched 10-m pixel crosses from SJ3/SJ4-like baseline
#      structure into SJ1/SJ2-like final structure.
#
# Thus a single anomalous pixel cannot make a non-restorative site enter this
# population.
#
# Baseline transition period = 2018-2019
# Final transition period    = 2024-2025
#
# Pixel area is approximately:
#   10 m x 10 m = 100 m2 = 0.01 ha
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------

PIXEL_AREA_HA <- 0.01


# Minimum number of SJ3/SJ4 -> SJ1/SJ2 pixels required.
#
# Leave at 1 initially so the analysis reports every qualifying site.
# The sensitivity table below also reports >=5 and >=10 pixels.

MIN_SJ34_TO_SJ12_PIXELS <- 1


# -----------------------------------------------------------------------------
# 2. DERIVE THE LATE -> EARLY TRANSITION METRICS
# -----------------------------------------------------------------------------

# Scrub/tree removal and the burn-coincident split come from python/
# mi_aesj_annual_products.py: scrub = SJ3-source, non-spruce pixels; trees =
# SJ4-source or spruce pixels. Scrub pixels inside any mapped burn and tree
# pixels inside a probable-wildfire burn are burn-coincident.

recovery_late_to_early_all <-
  recovery_sites %>%

  dplyr::filter(
    transition_qc_ok
  ) %>%

  dplyr::mutate(
    
    # ---------------------------------------------------------------
    # SJ3 source pixels
    # ---------------------------------------------------------------
    
    n_SJ3_to_SJ12 =
      SJ3_to_SJ1 +
      SJ3_to_SJ2,
    
    
    # ---------------------------------------------------------------
    # SJ4 source pixels
    # ---------------------------------------------------------------
    
    n_SJ4_to_SJ12 =
      SJ4_to_SJ1 +
      SJ4_to_SJ2,
    
    
    # ---------------------------------------------------------------
    # Combined transition:
    #
    # SJ3/SJ4 -> SJ1/SJ2
    # ---------------------------------------------------------------
    
    n_SJ34_to_SJ12 =
      n_SJ3_to_SJ12 +
      n_SJ4_to_SJ12,
    
    
    # ---------------------------------------------------------------
    # Baseline number of matched pixels in SJ3
    # ---------------------------------------------------------------
    
    n_SJ3_baseline =
      SJ3_to_SJ1 +
      SJ3_to_SJ2 +
      SJ3_to_SJ3 +
      SJ3_to_SJ4,
    
    
    # ---------------------------------------------------------------
    # Baseline number of matched pixels in SJ4
    # ---------------------------------------------------------------
    
    n_SJ4_baseline =
      SJ4_to_SJ1 +
      SJ4_to_SJ2 +
      SJ4_to_SJ3 +
      SJ4_to_SJ4,
    
    
    n_SJ34_baseline =
      n_SJ3_baseline +
      n_SJ4_baseline,
    
    
    # ---------------------------------------------------------------
    # Approximate areas
    # ---------------------------------------------------------------
    
    area_SJ3_to_SJ12_ha =
      n_SJ3_to_SJ12 *
      PIXEL_AREA_HA,
    
    area_SJ4_to_SJ12_ha =
      n_SJ4_to_SJ12 *
      PIXEL_AREA_HA,
    
    area_SJ34_to_SJ12_ha =
      n_SJ34_to_SJ12 *
      PIXEL_AREA_HA,
    
    
    # ---------------------------------------------------------------
    # Fraction of ALL matched site pixels making this transition
    # ---------------------------------------------------------------
    
    prop_site_SJ34_to_SJ12 =
      dplyr::if_else(
        
        is.finite(
          n_transition_pixels
        ) &
          n_transition_pixels > 0,
        
        n_SJ34_to_SJ12 /
          n_transition_pixels,
        
        NA_real_
        
      ),
    
    
    # ---------------------------------------------------------------
    # Of baseline SJ3 pixels, how many ended in SJ1/SJ2?
    # ---------------------------------------------------------------
    
    prop_SJ3_source_to_SJ12 =
      dplyr::if_else(
        
        n_SJ3_baseline > 0,
        
        n_SJ3_to_SJ12 /
          n_SJ3_baseline,
        
        NA_real_
        
      ),
    
    
    # ---------------------------------------------------------------
    # Of baseline SJ4 pixels, how many ended in SJ1/SJ2?
    # ---------------------------------------------------------------
    
    prop_SJ4_source_to_SJ12 =
      dplyr::if_else(
        
        n_SJ4_baseline > 0,
        
        n_SJ4_to_SJ12 /
          n_SJ4_baseline,
        
        NA_real_
        
      ),
    
    
    # ---------------------------------------------------------------
    # Of all baseline SJ3 + SJ4 pixels, what fraction reached SJ1/SJ2?
    # ---------------------------------------------------------------
    
    prop_SJ34_source_to_SJ12 =
      dplyr::if_else(
        
        n_SJ34_baseline > 0,
        
        n_SJ34_to_SJ12 /
          n_SJ34_baseline,
        
        NA_real_
        
      ),
    
    
    # ---------------------------------------------------------------
    # Which late-successional source contributed?
    # ---------------------------------------------------------------
    
    late_source_type =
      dplyr::case_when(
        
        n_SJ3_to_SJ12 > 0 &
          n_SJ4_to_SJ12 > 0 ~
          "SJ3 + SJ4",
        
        n_SJ3_to_SJ12 > 0 ~
          "SJ3 only",
        
        n_SJ4_to_SJ12 > 0 ~
          "SJ4 only",
        
        TRUE ~
          "None"
        
      ),
    
    
    late_source_type =
      factor(
        
        late_source_type,
        
        levels = c(
          "SJ3 only",
          "SJ4 only",
          "SJ3 + SJ4",
          "None"
        )
        
      ),
    
    
    # ---------------------------------------------------------------
    # A much stricter site-level interpretation:
    #
    # Dominant class SJ3/SJ4 at baseline and SJ1/SJ2 in the final period?
    #
    # This is NOT used for primary site selection below.
    # ---------------------------------------------------------------
    
    dominant_SJ34_to_SJ12 =
      dominant_base %in%
      c(
        "SJ3-like",
        "SJ4-like"
      ) &
      dominant_final %in%
      c(
        "SJ1-like",
        "SJ2-like"
      )
    
  )


# -----------------------------------------------------------------------------
# 3. SELECT THE TARGET RESTORATIVE POPULATION
# -----------------------------------------------------------------------------

recovery_SJ34_to_SJ12_sites <-
  recovery_late_to_early_all %>%
  
  dplyr::filter(
    
    is.finite(
      n_SJ34_to_SJ12
    ),
    
    n_SJ34_to_SJ12 >=
      MIN_SJ34_TO_SJ12_PIXELS
    
  ) %>%
  
  dplyr::arrange(
    
    dplyr::desc(
      n_SJ34_to_SJ12
    )
    
  )


# -----------------------------------------------------------------------------
# 4. NUMBER OF SITES
# -----------------------------------------------------------------------------

n_recovery_total <-
  dplyr::n_distinct(
    recovery_sites$fid
  )


n_recovery_transition_qc <-
  dplyr::n_distinct(
    recovery_late_to_early_all$fid
  )


n_recovery_SJ34_to_SJ12 <-
  dplyr::n_distinct(
    recovery_SJ34_to_SJ12_sites$fid
  )


cat(
  "\n============================================================\n",
  "RESTORATIVE SITES WITH SJ3/SJ4 -> SJ1/SJ2 TRANSITIONS\n",
  "============================================================\n\n",
  
  "All Moderate/Strong restorative sites: ",
  n_recovery_total,
  "\n",
  
  "Restorative sites passing transition QC: ",
  n_recovery_transition_qc,
  "\n",
  
  "Restorative sites with >= ",
  MIN_SJ34_TO_SJ12_PIXELS,
  " SJ3/SJ4 -> SJ1/SJ2 pixel(s): ",
  n_recovery_SJ34_to_SJ12,
  "\n",
  
  "Share of transition-QC restorative sites: ",
  round(
    100 *
      n_recovery_SJ34_to_SJ12 /
      n_recovery_transition_qc,
    1
  ),
  "%\n\n",
  
  sep = ""
)


# -----------------------------------------------------------------------------
# 5. SENSITIVITY TO THE NUMBER OF TRANSITION PIXELS
# -----------------------------------------------------------------------------
#
# Useful because one 10-m pixel is a permissive definition.
#
# 1 pixel  = approximately 0.01 ha
# 5 pixels = approximately 0.05 ha
# 10 pixels = approximately 0.10 ha
# -----------------------------------------------------------------------------

SJ34_transition_sensitivity <-
  tibble::tibble(
    
    minimum_pixels =
      c(
        1,
        5,
        10
      )
    
  ) %>%
  
  dplyr::mutate(
    
    minimum_area_ha =
      minimum_pixels *
      PIXEL_AREA_HA,
    
    n_sites =
      vapply(
        
        minimum_pixels,
        
        function(x) {
          
          sum(
            recovery_late_to_early_all$n_SJ34_to_SJ12 >= x,
            na.rm = TRUE
          )
          
        },
        
        integer(1)
        
      ),
    
    pct_transition_QC_recovery_sites =
      100 *
      n_sites /
      n_recovery_transition_qc
    
  )


cat(
  "Sensitivity to minimum number of late-to-early pixels:\n"
)


print(
  SJ34_transition_sensitivity,
  n = Inf
)


# -----------------------------------------------------------------------------
# 6. SITE COUNTS BY TILSTAND AND SOURCE TYPE
# -----------------------------------------------------------------------------

SJ34_site_counts_by_tilstand <-
  recovery_SJ34_to_SJ12_sites %>%
  
  dplyr::count(
    
    tilstand,
    late_source_type,
    
    name =
      "n_sites"
    
  ) %>%
  
  tidyr::complete(
    
    tilstand =
      factor(
        TILSTAND_LEVELS,
        levels =
          TILSTAND_LEVELS
      ),
    
    late_source_type =
      factor(
        c(
          "SJ3 only",
          "SJ4 only",
          "SJ3 + SJ4"
        ),
        levels = c(
          "SJ3 only",
          "SJ4 only",
          "SJ3 + SJ4",
          "None"
        )
      ),
    
    fill =
      list(
        n_sites = 0L
      )
    
  )


cat(
  "\nSites by tilstand and source transition:\n"
)


print(
  SJ34_site_counts_by_tilstand,
  n = Inf
)


# -----------------------------------------------------------------------------
# 7. TOTAL APPROXIMATE AREA OF SJ3/SJ4 -> SJ1/SJ2 PIXELS
# -----------------------------------------------------------------------------
#
# These are summed matched-pixel areas.
#
# If MI polygons overlap spatially, this is summed polygon-pixel area rather
# than guaranteed unique landscape area.
# -----------------------------------------------------------------------------

SJ34_area_summary <-
  recovery_SJ34_to_SJ12_sites %>%
  
  dplyr::summarise(
    
    n_sites =
      dplyr::n_distinct(
        fid
      ),
    
    n_SJ3_to_SJ12_pixels =
      sum(
        n_SJ3_to_SJ12,
        na.rm = TRUE
      ),
    
    n_SJ4_to_SJ12_pixels =
      sum(
        n_SJ4_to_SJ12,
        na.rm = TRUE
      ),
    
    n_SJ34_to_SJ12_pixels =
      sum(
        n_SJ34_to_SJ12,
        na.rm = TRUE
      ),
    
    area_SJ3_to_SJ12_ha =
      sum(
        area_SJ3_to_SJ12_ha,
        na.rm = TRUE
      ),
    
    area_SJ4_to_SJ12_ha =
      sum(
        area_SJ4_to_SJ12_ha,
        na.rm = TRUE
      ),
    
    area_SJ34_to_SJ12_ha =
      sum(
        area_SJ34_to_SJ12_ha,
        na.rm = TRUE
      ),

    # Burn-coincident vs clearance, split per pixel in
    # python/mi_aesj_annual_products.py with the dNBR burn mask: scrub inside
    # any burn and trees inside a probable-wildfire burn are burn-coincident.
    # The same split is given for the spruce subset below.
    area_SJ34_to_SJ12_burnt_ha =
      sum(
        area_SJ34_to_SJ12_burnt_ha,
        na.rm = TRUE
      ),

    area_SJ34_to_SJ12_clearance_ha =
      sum(
        area_SJ34_to_SJ12_clearance_ha,
        na.rm = TRUE
      ),

    n_SJ34_to_SJ12_picea_pixels =
      sum(
        n_SJ34_to_SJ12_picea,
        na.rm = TRUE
      ),

    area_SJ34_to_SJ12_picea_ha =
      sum(
        area_SJ34_to_SJ12_picea_ha,
        na.rm = TRUE
      ),

    area_SJ34_to_SJ12_picea_burnt_ha =
      sum(
        area_SJ34_to_SJ12_picea_burnt_ha,
        na.rm = TRUE
      ),

    area_SJ34_to_SJ12_picea_clearance_ha =
      sum(
        area_SJ34_to_SJ12_picea_clearance_ha,
        na.rm = TRUE
      ),

    median_area_per_site_ha =
      stats::median(
        area_SJ34_to_SJ12_ha,
        na.rm = TRUE
      ),
    
    mean_area_per_site_ha =
      mean(
        area_SJ34_to_SJ12_ha,
        na.rm = TRUE
      ),
    
    median_prop_site_pixels =
      stats::median(
        prop_site_SJ34_to_SJ12,
        na.rm = TRUE
      ),
    
    mean_prop_site_pixels =
      mean(
        prop_site_SJ34_to_SJ12,
        na.rm = TRUE
      ),
    
    median_prop_baseline_SJ34_reaching_SJ12 =
      stats::median(
        prop_SJ34_source_to_SJ12,
        na.rm = TRUE
      ),
    
    mean_prop_baseline_SJ34_reaching_SJ12 =
      mean(
        prop_SJ34_source_to_SJ12,
        na.rm = TRUE
      )
    
  )


cat(
  "\nSJ3/SJ4 -> SJ1/SJ2 area summary:\n"
)


print(
  SJ34_area_summary,
  n = Inf
)


# -----------------------------------------------------------------------------
# 8. AREA BY TILSTAND
# -----------------------------------------------------------------------------

SJ34_area_by_tilstand <-
  recovery_SJ34_to_SJ12_sites %>%
  
  dplyr::group_by(
    tilstand
  ) %>%
  
  dplyr::summarise(
    
    n_sites =
      dplyr::n_distinct(
        fid
      ),
    
    n_transition_pixels =
      sum(
        n_SJ34_to_SJ12,
        na.rm = TRUE
      ),
    
    transition_area_ha =
      sum(
        area_SJ34_to_SJ12_ha,
        na.rm = TRUE
      ),
    
    median_area_per_site_ha =
      stats::median(
        area_SJ34_to_SJ12_ha,
        na.rm = TRUE
      ),
    
    median_prop_site_pixels =
      stats::median(
        prop_site_SJ34_to_SJ12,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  )


# -----------------------------------------------------------------------------
# 9. EXACT TRANSITION PATHWAYS
# -----------------------------------------------------------------------------

SJ34_pathway_summary <-
  tibble::tibble(
    
    transition =
      c(
        "SJ3 -> SJ1",
        "SJ3 -> SJ2",
        "SJ4 -> SJ1",
        "SJ4 -> SJ2"
      ),
    
    n_pixels =
      c(
        
        sum(
          recovery_SJ34_to_SJ12_sites$SJ3_to_SJ1,
          na.rm = TRUE
        ),
        
        sum(
          recovery_SJ34_to_SJ12_sites$SJ3_to_SJ2,
          na.rm = TRUE
        ),
        
        sum(
          recovery_SJ34_to_SJ12_sites$SJ4_to_SJ1,
          na.rm = TRUE
        ),
        
        sum(
          recovery_SJ34_to_SJ12_sites$SJ4_to_SJ2,
          na.rm = TRUE
        )
        
      )
    
  ) %>%
  
  dplyr::mutate(
    
    area_ha =
      n_pixels *
      PIXEL_AREA_HA,
    
    proportion_of_late_to_early_area =
      n_pixels /
      sum(
        n_pixels
      )
    
  )


cat(
  "\nExact SJ3/SJ4 -> SJ1/SJ2 pathways:\n"
)


print(
  SJ34_pathway_summary,
  n = Inf
)


# -----------------------------------------------------------------------------
# 10. STRICT DOMINANT-CLASS SHIFT COUNT
# -----------------------------------------------------------------------------
#
# This is substantially stricter than the main definition.
#
# A site qualifies only if its dominant baseline class was SJ3/SJ4 AND its
# dominant final class was SJ1/SJ2.
#
# Report this as a secondary indicator, not the main pixel-transition result.
# -----------------------------------------------------------------------------

n_dominant_SJ34_to_SJ12 <-
  recovery_SJ34_to_SJ12_sites %>%
  
  dplyr::filter(
    dominant_SJ34_to_SJ12
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n_distinct(
        fid
      )
    
  ) %>%
  
  dplyr::pull(
    n
  )


cat(
  "\nRestorative sites whose DOMINANT class moved SJ3/SJ4 -> SJ1/SJ2: ",
  n_dominant_SJ34_to_SJ12,
  "\n",
  sep = ""
)


# =============================================================================
# 11. REPLICATE RECOVERY COMPOSITION FIGURE FOR THIS SUBSET
# =============================================================================

SJ34_recovery_composition_summary <-
  build_period_composition_summary(
    
    recovery_SJ34_to_SJ12_sites,
    
    group_vars =
      "tilstand"
    
  )


p_SJ34_recovery_composition <-
  make_composition_change_plot(
    
    SJ34_recovery_composition_summary,
    
    title =
      "Restorative sites with late-to-early structural transitions",
    
    subtitle =
      paste0(
        "Sites containing matched pixels moving from SJ3/SJ4-like in 2018-2019 ",
        "to SJ1/SJ2-like in 2024-2025"
      )
    
  )




# =============================================================================
# 12. REPLICATE RECOVERY TRAJECTORY FIGURE FOR THIS SUBSET
# =============================================================================

p_SJ34_recovery_trajectories <-
  make_trajectory_plot(
    
    state_dat,
    
    recovery_SJ34_to_SJ12_sites,
    
    ref_dat,
    
    line_colour =
      SKOG_THEME$g3,
    
    title =
      "Restorative sites with SJ3/SJ4 to SJ1/SJ2 transitions",
    
    subtitle =
      paste0(
        "Annual AE-SJ trajectories; top ",
        min(
          SITE_GRAPH_TOP_N_PER_DIRECTION,
          n_recovery_SJ34_to_SJ12
        ),
        " sites by estimated restorative change"
      ),
    
    top_n =
      min(
        SITE_GRAPH_TOP_N_PER_DIRECTION,
        n_recovery_SJ34_to_SJ12
      )
    
  )




# =============================================================================
# 13. WRITE RESULTS
# =============================================================================

readr::write_csv(
  
  recovery_SJ34_to_SJ12_sites,
  
  file.path(
    file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "15_restorative_sites_SJ34_to_SJ12.csv")
  )
  
)


readr::write_csv(
  
  SJ34_transition_sensitivity,
  
  file.path(
    file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "16_SJ34_to_SJ12_site_count_sensitivity.csv")
  )
  
)


readr::write_csv(
  
  SJ34_site_counts_by_tilstand %>%
    dplyr::mutate(
      tilstand =
        as.character(
          tilstand
        ),
      late_source_type =
        as.character(
          late_source_type
        )
    ),
  
  file.path(
    file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "17_SJ34_to_SJ12_counts_by_tilstand.csv")
  )
  
)


readr::write_csv(
  
  SJ34_area_summary,
  
  file.path(
    file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "18_SJ34_to_SJ12_area_summary.csv")
  )
  
)


readr::write_csv(
  
  SJ34_area_by_tilstand,
  
  file.path(
    file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "19_SJ34_to_SJ12_area_by_tilstand.csv")
  )
  
)


readr::write_csv(
  
  SJ34_pathway_summary,
  
  file.path(
    file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "20_SJ34_to_SJ12_transition_pathways.csv")
  )
  
)


readr::write_csv(
  
  SJ34_recovery_composition_summary,
  
  file.path(
    file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "21_SJ34_to_SJ12_composition_summary.csv")
  )
  
)


# =============================================================================
# 14. SAVE FIGURES
# =============================================================================
#


# =============================================================================
# 15. FINAL CONSOLE SUMMARY
# =============================================================================

cat(
  "\n============================================================\n",
  "SJ3/SJ4 -> SJ1/SJ2 RESTORATION SUMMARY\n",
  "============================================================\n",
  
  "Qualifying restorative sites: ",
  n_recovery_SJ34_to_SJ12,
  "\n",
  
  "Strict dominant-class SJ3/SJ4 -> SJ1/SJ2 sites: ",
  n_dominant_SJ34_to_SJ12,
  "\n",
  
  "Matched pixels making SJ3/SJ4 -> SJ1/SJ2 transition: ",
  round(
    SJ34_area_summary$n_SJ34_to_SJ12_pixels,
    1
  ),
  "\n",
  
  "Approximate summed transition area: ",
  round(
    SJ34_area_summary$area_SJ34_to_SJ12_ha,
    2
  ),
  " ha\n",
  
  "  From SJ3: ",
  round(
    SJ34_area_summary$area_SJ3_to_SJ12_ha,
    2
  ),
  " ha\n",
  
  "  From SJ4: ",
  round(
    SJ34_area_summary$area_SJ4_to_SJ12_ha,
    2
  ),
  " ha\n",
  
  sep = ""
)


p_SJ34_recovery_trajectories <-
  make_trajectory_plot(
    
    state_dat,
    
    recovery_SJ34_to_SJ12_sites,
    
    ref_dat,
    
    line_colour =
      SKOG_THEME$g3,
    
    title =
      "Sites moving from SJ3/SJ4-like to SJ1/SJ2-like conditions",
    
    subtitle =
      paste0(
        "Annual site-median AE-SJ trajectories; ",
        "baseline = 2018-2019, final = 2024-2025"
      ),

    top_n =
      min(
        SITE_GRAPH_TOP_N_PER_DIRECTION,
        nrow(
          recovery_SJ34_to_SJ12_sites
        )
      )

  )


# =============================================================================
# COMPOSITION FIGURE:
# SJ3/SJ4 -> SJ1/SJ2 RESTORATIVE SITES, ALL TILSTAND COMBINED
# =============================================================================
#
# Each qualifying FID contributes equally.
#
# Baseline = mean pixel composition across 2018-2019
# Final    = mean pixel composition across 2024-2025
#
# No stratification by tilstand.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Build pre/post composition table
# -----------------------------------------------------------------------------

SJ34_recovery_composition_overall <-
  dplyr::bind_rows(
    
    recovery_SJ34_to_SJ12_sites %>%
      
      dplyr::transmute(
        
        fid,
        
        period =
          "2018-2019",
        
        `SJ1-like` =
          p1_base,
        
        `SJ2-like` =
          p2_base,
        
        `SJ3-like` =
          p3_base,
        
        `SJ4-like` =
          p4_base
        
      ),
    
    
    recovery_SJ34_to_SJ12_sites %>%
      
      dplyr::transmute(
        
        fid,
        
        period =
          "2024-2025",
        
        `SJ1-like` =
          p1_final,
        
        `SJ2-like` =
          p2_final,
        
        `SJ3-like` =
          p3_final,
        
        `SJ4-like` =
          p4_final
        
      )
    
  ) %>%
  
  tidyr::pivot_longer(
    
    cols =
      dplyr::all_of(
        SJ_ZONE_LEVELS
      ),
    
    names_to =
      "sj_zone",
    
    values_to =
      "proportion"
    
  ) %>%
  
  dplyr::mutate(
    
    period =
      factor(
        period,
        levels = c(
          "2018-2019",
          "2024-2025"
        )
      ),
    
    sj_zone =
      factor(
        sj_zone,
        levels =
          SJ_ZONE_LEVELS
      )
    
  )


# -----------------------------------------------------------------------------
# 2. Equal-site-weighted mean composition
# -----------------------------------------------------------------------------

SJ34_recovery_composition_summary <-
  SJ34_recovery_composition_overall %>%
  
  dplyr::group_by(
    period,
    sj_zone
  ) %>%
  
  dplyr::summarise(
    
    n_fids =
      dplyr::n_distinct(
        fid
      ),
    
    mean_proportion =
      mean(
        proportion,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  )


# Labels only for segments >= 5%.

SJ34_recovery_composition_labels <-
  SJ34_recovery_composition_summary %>%
  
  dplyr::filter(
    mean_proportion >= 0.05
  )


# -----------------------------------------------------------------------------
# 3. Number of qualifying sites
# -----------------------------------------------------------------------------

n_SJ34_recovery_sites <-
  dplyr::n_distinct(
    recovery_SJ34_to_SJ12_sites$fid
  )


# -----------------------------------------------------------------------------
# 4. Plot
# -----------------------------------------------------------------------------

p_SJ34_recovery_composition <-
  ggplot2::ggplot(
    
    SJ34_recovery_composition_summary,
    
    ggplot2::aes(
      x = period,
      y = mean_proportion,
      fill = sj_zone
    )
    
  ) +
  
  ggplot2::geom_col(
    width = 0.62
  ) +
  
  ggplot2::geom_text(
    
    data =
      SJ34_recovery_composition_labels,
    
    ggplot2::aes(
      label =
        percent_axis(
          mean_proportion
        )
    ),
    
    position =
      ggplot2::position_stack(
        vjust = 0.5
      ),
    
    colour =
      SKOG_THEME$surface,
    
    fontface =
      "bold",
    
    size =
      gg_text_size(
        "on_bar"
      )
    
  ) +
  
  ggplot2::scale_fill_manual(
    
    values =
      SJ_ZONE_COLOURS,
    
    drop =
      FALSE,
    
    name =
      "SJ Class"
    
  ) +
  
  ggplot2::scale_y_continuous(
    
    limits = c(
      0,
      1
    ),
    
    breaks =
      seq(
        0,
        1,
        0.2
      ),
    
    labels =
      percent_axis,
    
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0
        )
      )
    
  ) +
  
  ggplot2::labs(
    
    title =
      "Sites moving from SJ3/SJ4-like to SJ1/SJ2-like conditions",
    
    subtitle =
      paste0(
        "Mean pixel composition before and after the monitoring period; ",
        "n = ",
        skog_num(
          n_SJ34_recovery_sites
        ),
        " restorative sites"
      ),
    
    x =
      NULL,
    
    y =
      "Mean proportion of site pixels",
    
    caption =
      wrap_caption(
        paste0(
          "Each site contributes equally. Baseline is the mean composition ",
          "across 2018-2019; final is the mean composition across 2024-2025. ",
          "SJ Classes are AE-SJ score zones rather than field-observed classes."
        )
      )
    
  ) +
  
  theme_report() +
  
  ggplot2::theme(
    legend.position =
      "top"
  )


print(
  p_SJ34_recovery_composition
)


# -----------------------------------------------------------------------------
# 5. Save
# -----------------------------------------------------------------------------

ggplot2::ggsave(
  
  filename =
    file.path(
      file.path(OUTPUT_DIR, "SJ34_to_SJ12_restorative", "14_SJ34_to_SJ12_composition_change.png")
    ),
  
  plot =
    p_SJ34_recovery_composition,
  
  width =
    7.5,
  
  height =
    6,
  
  dpi =
    350,
  
  bg =
    "white"
  
)



