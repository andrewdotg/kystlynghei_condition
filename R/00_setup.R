# =============================================================================
# 00_setup.R - project root resolution + shared package loading
# =============================================================================
#
# Sourced by every script in this repo before anything else. Resolves TWO
# roots, because this repo (kystlynghei-condition/) sits inside, but is
# distinct from, the larger project folder the raw input data live in:
#
#   REPO_ROOT - this git repo. Captured as getwd() before anything else
#               changes it - safe, because the very fact that
#               source("R/00_setup.R") found this file means the working
#               directory was already the repo root when it was called.
#               Used for source()-ing other files in R/ and for OUTPUT_DIR
#               below, so results always land inside the repo regardless of
#               where the data root turns out to be.
#
#   PROJECT_ROOT - the larger Andrew_kystlynghei project folder, where the
#               raw input CSVs and rasters live (e.g. the field calibration
#               points). Input paths in the scripts (e.g.
#               "equalArea2/GEEpointStats/", "deploy/BurningAnalysis/...") are
#               relative to PROJECT_ROOT, which is why the working directory
#               is changed to it below.
#
# PROJECT_ROOT is found via, in order: KYSTLYNGHEI_PROJECT_ROOT (set this if
# you're running from a plain clone that isn't inside the original P: drive
# project folder), the P:\ mapped-drive path, the RStudio Server Linux mount.
#
# CALL IT GUARDED, NOT BARE: every script that needs this file sources it as
#   if (!exists("REPO_ROOT")) source("R/00_setup.R")
# not a bare source("R/00_setup.R"). Reason: this file's own setwd() below
# moves the working directory away from REPO_ROOT, so a second bare
# source("R/00_setup.R") later in the same session (e.g. run_all.R sourcing
# several scripts back to back, or running scripts one after another in the
# same RStudio console without restarting R) would fail outright - R can't
# find "R/00_setup.R" relative to a working directory that is no longer
# REPO_ROOT. The guard makes this file run exactly once per session: the
# first script sources it for real (cwd is still REPO_ROOT at that point,
# e.g. because you opened this folder as an RStudio Project), every script
# after that just reuses the REPO_ROOT/PROJECT_ROOT/OUTPUTS_ROOT it already
# set. Everything else that sources another file in R/ (theme.R,
# ae_sj_core.R, or a module sourcing this file itself) uses a
# REPO_ROOT-anchored file.path(), not a bare relative path, for the same
# reason.
# =============================================================================

# Headless Linux: use the ragg device, so ggplot text renders with the
# report fonts when no display is attached (the base png()/X11 device fails
# with "invalid font type").
options(bitmapType = "cairo")
if (requireNamespace("ragg", quietly = TRUE)) {
  options(device = function(...) ragg::agg_png(..., width = 800, height = 600, units = "px", res = 100))
}

REPO_ROOT <- getwd()

resolve_project_root <- function() {
  candidates <- c(
    Sys.getenv("KYSTLYNGHEI_PROJECT_ROOT", unset = NA),
    "P:/156004_08_effektovervaking_2026_2030_opsjon_2_fjernmaling/Andrew_kystlynghei",
    "/data/P-Prosjekter2/156004_08_effektovervaking_2026_2030_opsjon_2_fjernmaling/Andrew_kystlynghei"
  )
  candidates <- candidates[!is.na(candidates)]
  hit <- candidates[dir.exists(candidates)]

  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))

  # Nothing matched - fall back to whatever the working directory already
  # was (e.g. this repo was opened as its own .Rproj right inside the
  # original project folder, so REPO_ROOT already doubles as PROJECT_ROOT).
  REPO_ROOT
}

PROJECT_ROOT <- resolve_project_root()
setwd(PROJECT_ROOT)

# Every module writes its results under here (each in its own subfolder),
# instead of back out into the original scattered deploy/ Development/
# ScrubClearence/ GeethenComp/ folders. Gitignored - see .gitignore - but
# having everything one script produces land in one predictable place inside
# the repo is what a future Quarto document will read from.
OUTPUTS_ROOT <- file.path(REPO_ROOT, "outputs")

# Loads packages, stopping with a clear message if any are missing.
required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) {
    stop(
      "Install these packages before running this script: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}
