# ================================================================
# DEVELOPMENT / URBANISATION LOSS MODULE
# VISUALISE DYNAMIC WORLD NATURE LOSS TO DEVELOPMENT
# ================================================================
#
# Three parts: development-loss summary figures; joining the results to the
# MI point shapefile; cumulative loss by fylke.

if (!exists("REPO_ROOT")) source("R/00_setup.R")
# ggrepel labels the fylke figure.
required_packages(c("tidyverse", "scales", "patchwork", "sf", "ggrepel"))

# All outputs go to OUTPUT_DIR.
OUTPUT_DIR <- file.path(OUTPUTS_ROOT, "development_loss")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# LOAD DATA
# ------------------------------------------------

# Change to your downloaded filename
dat <- read_csv("Development/MI_DW_persistent_development.csv")

glimpse(dat)


# ------------------------------------------------
# CLEAN DATA
# ------------------------------------------------

dat <- dat %>%
  mutate(
    across(
      c(
        dev_px_16_20,
        dev_ha_16_20,
        dev_px_21_25,
        dev_ha_21_25,
        dev_px_total,
        dev_ha_total
      ),
      ~replace_na(as.numeric(.x), 0)
    )
  )


# ================================================================
# BASIC SUMMARY
# ================================================================

summary_stats <- dat %>%
  summarise(
    n_features = n(),
    
    features_with_loss =
      sum(dev_ha_total > 0),
    
    pct_features_with_loss =
      mean(dev_ha_total > 0) * 100,
    
    loss_2016_2020_ha =
      sum(dev_ha_16_20),
    
    loss_2021_2025_ha =
      sum(dev_ha_21_25),
    
    total_loss_ha =
      sum(dev_ha_total),
    
    total_pixels =
      sum(dev_px_total)
  )

print(summary_stats)


# ================================================================
# 1. TOTAL DEVELOPMENT LOSS BY PERIOD
# ================================================================

period_summary <- dat %>%
  summarise(
    `2016–2020` = sum(dev_ha_16_20),
    `2021–2025` = sum(dev_ha_21_25)
  ) %>%
  pivot_longer(
    everything(),
    names_to = "period",
    values_to = "area_ha"
  )


p_period <- ggplot(
  period_summary,
  aes(
    x = period,
    y = area_ha
  )
) +
  geom_col(
    width = 0.65
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f ha", area_ha)
    ),
    vjust = -0.5
  ) +
  labs(
    title = "Detected nature loss to development",
    subtitle = "Persistent Dynamic World built transitions",
    x = NULL,
    y = "Detected area (ha)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank()
  )


# ================================================================
# 2. DISTRIBUTION AMONG FEATURES
# ================================================================

p_distribution <- dat %>%
  filter(
    dev_ha_total > 0
  ) %>%
  ggplot(
    aes(
      x = dev_ha_total
    )
  ) +
  geom_histogram(
    bins = 30
  ) +
  labs(
    title = "Distribution of detected development loss",
    subtitle = "Only features with detected loss",
    x = "Nature lost to development (ha)",
    y = "Number of features"
  ) +
  theme_minimal(base_size = 12)



# ================================================================
# 3. TOP FEATURES
# ================================================================

# Determine the most useful feature ID column automatically.
#
# Change this manually if your identifying column has another name.
possible_id_fields <- c(
  "fid",
  "id",
  "system:index",
  "lokalId",
  "lokalid"
)

id_field <- possible_id_fields[
  possible_id_fields %in% names(dat)
][1]


if (is.na(id_field)) {
  
  dat <- dat %>%
    mutate(
      plot_id = row_number()
    )
  
} else {
  
  dat <- dat %>%
    mutate(
      plot_id = .data[[id_field]]
    )
}


top20 <- dat %>%
  arrange(
    desc(dev_ha_total)
  ) %>%
  slice_head(
    n = 20
  ) %>%
  mutate(
    plot_id = factor(
      plot_id,
      levels = rev(plot_id)
    )
  )


p_top <- ggplot(
  top20,
  aes(
    x = plot_id,
    y = dev_ha_total
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Features with greatest detected development",
    x = "Feature",
    y = "Detected nature loss (ha)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank()
  )



# ================================================================
# 4. DEVELOPMENT BY PERIOD FOR TOP FEATURES
# ================================================================

top_period <- top20 %>%
  select(
    plot_id,
    dev_ha_16_20,
    dev_ha_21_25
  ) %>%
  pivot_longer(
    c(
      dev_ha_16_20,
      dev_ha_21_25
    ),
    names_to = "period",
    values_to = "area_ha"
  ) %>%
  mutate(
    period = recode(
      period,
      dev_ha_16_20 = "2016–2020",
      dev_ha_21_25 = "2021–2025"
    )
  )


p_top_period <- ggplot(
  top_period,
  aes(
    x = plot_id,
    y = area_ha,
    fill = period
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Timing of development in highest-loss features",
    x = "Feature",
    y = "Detected nature loss (ha)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank()
  )



# ================================================================
# 5. PIXEL COUNT VS AREA
# ================================================================
#
# Mostly a sanity check that the exported pixel-count and actual
# pixel-area measures behave as expected.
# ================================================================

p_pixels <- dat %>%
  filter(
    dev_px_total > 0
  ) %>%
  ggplot(
    aes(
      x = dev_px_total,
      y = dev_ha_total
    )
  ) +
  geom_point(
    alpha = 0.5
  ) +
  labs(
    title = "Detected pixels versus calculated area",
    x = "Detected Dynamic World pixels",
    y = "Detected area (ha)"
  ) +
  theme_minimal(base_size = 12)



# ================================================================
# 6. FEATURES WITH DETECTED DEVELOPMENT
# ================================================================

development_features <- dat %>%
  filter(
    dev_ha_total > 0
  ) %>%
  arrange(
    desc(dev_ha_total)
  ) %>%
  select(
    plot_id,
    dev_px_16_20,
    dev_ha_16_20,
    dev_px_21_25,
    dev_ha_21_25,
    dev_px_total,
    dev_ha_total,
    everything()
  )


print(
  development_features,
  n = 30
)


# Export a convenient ranked table.
write_csv(
  development_features,
  file.path(OUTPUT_DIR, "MI_DW_development_ranked.csv")
)


# ================================================================
# 7. SIMPLE SUMMARY FIGURE
# ================================================================

summary_figure <-
  p_period +
  p_distribution +
  p_top +
  plot_layout(
    widths = c(1, 1, 1.5)
  )



# (library(tidyverse); library(scales); library(patchwork); library(sf) -
# now loaded once at the top of the file, see the header note.)


# ================================================================
# JOIN DEVELOPMENT RESULTS TO POINT SHAPEFILE
# ================================================================

# ------------------------------------------------
# 1. CHOOSE FOLDER CONTAINING POINT SHAPEFILES
# ------------------------------------------------
#
# Folder (relative to PROJECT_ROOT) containing the candidate point layers.
# ------------------------------------------------

point_folder <- "equalArea2/MI_SJtrend"


# ------------------------------------------------
# 2. FIND AVAILABLE SHAPEFILES
# ------------------------------------------------

point_files <- list.files(
  point_folder,
  pattern = "\\.shp$",
  full.names = TRUE,
  ignore.case = TRUE
)


if (length(point_files) == 0) {
  stop("No shapefiles found in point_folder.")
}


# Print available files.
cat("\nAvailable shapefiles:\n\n")

for (i in seq_along(point_files)) {
  cat(
    i,
    ": ",
    basename(point_files[i]),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------
# 3. PICK THE POINT LAYER
# ------------------------------------------------
#
# OPTION A — manually specify its number after looking
# at the list printed above.
#
# For example:
# point_layer_number <- 3
#
# ------------------------------------------------

point_layer_number <- 3   # <-- CHANGE THIS


point_file <- point_files[
  point_layer_number
]


cat(
  "\nUsing point layer:\n",
  point_file,
  "\n"
)


# ------------------------------------------------
# 4. READ POINT LAYER
# ------------------------------------------------

points <- st_read(
  point_file,
  quiet = FALSE
)


# ------------------------------------------------
# 5. CHECK THAT IT REALLY IS A POINT LAYER
# ------------------------------------------------

print(
  unique(
    st_geometry_type(points)
  )
)


if (!all(
  st_geometry_type(points) %in%
  c("POINT", "MULTIPOINT")
)) {
  
  warning(
    "Selected layer is not entirely POINT/MULTIPOINT geometry."
  )
}


# ------------------------------------------------
# 6. CHECK FID EXISTS
# ------------------------------------------------

if (!"fid" %in% names(points)) {
  stop(
    "Selected point layer does not contain a column called 'fid'."
  )
}


if (!"fid" %in% names(dat)) {
  stop(
    "Development CSV does not contain a column called 'fid'."
  )
}


# ------------------------------------------------
# 7. MAKE FID TYPES MATCH
# ------------------------------------------------
#
# Character is safest for joins because one source may
# have imported fid as numeric and another as text.
# ------------------------------------------------

points <- points %>%
  mutate(
    fid = as.character(fid)
  )


dat <- dat %>%
  mutate(
    fid = as.character(fid)
  )


# ------------------------------------------------
# 8. CHECK FOR DUPLICATE FIDS
# ------------------------------------------------

duplicate_csv_fids <- dat %>%
  count(fid) %>%
  filter(n > 1)


if (nrow(duplicate_csv_fids) > 0) {
  
  print(
    duplicate_csv_fids
  )
  
  warning(
    "Duplicate fid values found in the development CSV."
  )
}


duplicate_point_fids <- points %>%
  st_drop_geometry() %>%
  count(fid) %>%
  filter(n > 1)


if (nrow(duplicate_point_fids) > 0) {
  
  print(
    duplicate_point_fids
  )
  
  warning(
    "Duplicate fid values found in the point layer."
  )
}


# ------------------------------------------------
# 9. KEEP DEVELOPMENT VARIABLES FOR JOIN
# ------------------------------------------------
#
# This avoids accidentally duplicating lots of columns that
# already exist in the point shapefile.
# ------------------------------------------------

development_join <- dat %>%
  select(
    fid,
    
    dev_px_16_20,
    dev_ha_16_20,
    
    dev_px_21_25,
    dev_ha_21_25,
    
    dev_px_total,
    dev_ha_total
  )


# ------------------------------------------------
# 10. JOIN TO POINTS
# ------------------------------------------------
#
# left_join means ALL points are retained.
# Points without a development result get NA initially.
# ------------------------------------------------

points_development <- points %>%
  left_join(
    development_join,
    by = "fid"
  )


# ------------------------------------------------
# 11. REPLACE MISSING DEVELOPMENT VALUES WITH ZERO
# ------------------------------------------------

points_development <- points_development %>%
  mutate(
    across(
      c(
        dev_px_16_20,
        dev_ha_16_20,
        dev_px_21_25,
        dev_ha_21_25,
        dev_px_total,
        dev_ha_total
      ),
      ~replace_na(.x, 0)
    )
  )


# ------------------------------------------------
# 12. CHECK JOIN SUCCESS
# ------------------------------------------------

join_check <- points_development %>%
  st_drop_geometry() %>%
  summarise(
    
    n_points = n(),
    
    n_with_development =
      sum(dev_px_total > 0),
    
    total_dev_pixels =
      sum(dev_px_total),
    
    total_dev_ha =
      sum(dev_ha_total)
    
  )


print(
  join_check
)


# FIDs in the CSV that were not found in the point layer.
missing_from_points <- anti_join(
  dat %>% select(fid),
  points %>%
    st_drop_geometry() %>%
    select(fid),
  by = "fid"
)


cat(
  "\nNumber of CSV fids not found in point layer:",
  nrow(missing_from_points),
  "\n"
)


if (nrow(missing_from_points) > 0) {
  
  print(
    missing_from_points
  )
  
}


# ================================================================
# WRITE POINT SHAPEFILE
# ================================================================

# The joined point shapefile goes to OUTPUT_DIR.
output_shp <- file.path(
  OUTPUT_DIR,
  "MI_DW_development_points.shp"
)


st_write(
  points_development,
  output_shp,
  delete_layer = TRUE
)


cat(
  "\nSaved point shapefile to:\n",
  output_shp,
  "\n"
)

# ================================================================
# 8. CUMULATIVE AREA LOSS OVER TIME
# ================================================================
#
# Only two coarse periods are available (2016-2020, 2021-2025), so
# this is three checkpoints, not a smooth curve. The straight lines
# between checkpoints imply steady loss within each period, which
# isn't actually known - real timing within a period could be more
# front- or back-loaded. Noted in the subtitle so it isn't read as
# more precise than it is.
# ================================================================

cumulative_loss <- tibble(
  year =
    c(2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025),
  
  cumulative_ha =
    c(
      0,
      sum(dat$dev_ha_16_20),
      sum(dat$dev_ha_16_20) + sum(dat$dev_ha_21_25)
    )
)


p_cumulative <- ggplot(
  cumulative_loss,
  aes(
    x = year,
    y = cumulative_ha
  )
) +
  geom_area(
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f ha", cumulative_ha)
    ),
    vjust = -1
  ) +
  scale_x_continuous(
    breaks = cumulative_loss$year
  ) +
  labs(
    title = "Cumulative nature loss to development",
    subtitle = "Running total of persistent Dynamic World built transitions; straight lines between the two available periods, not a measured trajectory",
    x = NULL,
    y = "Cumulative detected area (ha)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank()
  )


# ================================================================
# ADD-ON: CUMULATIVE DEVELOPMENT AREA LOSS BY FYLKE
# ================================================================
#
# Reads the annual Dynamic World development CSV and joins each FID
# to fylke using the shapefile field `filkesna_4`.
#
# Output:
#   - annual development loss by fylke
#   - cumulative development loss by fylke
#   - cumulative area-loss line graph, coloured by fylke
#
# The CSV contains annual area fields:
#   dev_ha_2018 ... dev_ha_2025
# ================================================================


# (library(tidyverse); library(sf); library(scales) - now loaded once at the
# top of the file, see the header note.)


# ------------------------------------------------
# 1. INPUTS
# ------------------------------------------------

ANNUAL_DEV_CSV <-
  "Development/MI_DW_persistent_development_annual.csv"
# MI polygons with fylke (fields fid, filkesna_4), in the shared Data folder
# next to the project folder.
FYLKE_SHP <- file.path(PROJECT_ROOT, "..", "Data", "Kystlynghei", "MI_fixed_fylke.shp")


# ------------------------------------------------
# 2. READ ANNUAL DEVELOPMENT DATA
# ------------------------------------------------

dev_annual <- readr::read_csv(
  ANNUAL_DEV_CSV,
  show_col_types = FALSE
) %>%
  mutate(
    fid = as.character(fid)
  )


# ------------------------------------------------
# 3. CHECK REQUIRED ANNUAL FIELDS
# ------------------------------------------------

annual_ha_fields <- paste0(
  "dev_ha_",
  2018:2025
)

missing_annual_fields <- setdiff(
  annual_ha_fields,
  names(dev_annual)
)

if (length(missing_annual_fields) > 0) {
  stop(
    "Annual development CSV is missing: ",
    paste(
      missing_annual_fields,
      collapse = ", "
    )
  )
}


# Convert annual loss fields to numeric and replace NA with 0.
dev_annual <- dev_annual %>%
  mutate(
    across(
      all_of(annual_ha_fields),
      ~ replace_na(
        suppressWarnings(as.numeric(.x)),
        0
      )
    )
  )


# ------------------------------------------------
# 4. READ FYLKE LOOKUP SHAPEFILE
# ------------------------------------------------

fylke_sf <- sf::st_read(
  FYLKE_SHP,
  quiet = TRUE
)


if (!"fid" %in% names(fylke_sf)) {
  stop(
    "Fylke shapefile does not contain a field called 'fid'."
  )
}

if (!"fylkesna_4" %in% names(fylke_sf)) {
  stop(
    "Fylke shapefile does not contain a field called 'filkesna_4'."
  )
}


# ------------------------------------------------
# 5. CREATE ONE FID -> FYLKE LOOKUP
# ------------------------------------------------

fylke_lookup <- fylke_sf %>%
  sf::st_drop_geometry() %>%
  transmute(
    fid = as.character(fid),
    fylke = trimws(as.character(fylkesna_4))
  ) %>%
  filter(
    !is.na(fid),
    fid != "",
    !is.na(fylke),
    fylke != ""
  ) %>%
  distinct(
    fid,
    .keep_all = TRUE
  )


# ------------------------------------------------
# 6. JOIN FYLKE TO DEVELOPMENT DATA
# ------------------------------------------------

dev_with_fylke <- dev_annual %>%
  left_join(
    fylke_lookup,
    by = "fid"
  )


# Quick join check.
cat(
  "\nAnnual development FIDs: ",
  n_distinct(dev_annual$fid),
  "\nFIDs with fylke: ",
  n_distinct(dev_with_fylke$fid[!is.na(dev_with_fylke$fylke)]),
  "\nFIDs without fylke: ",
  n_distinct(dev_with_fylke$fid[is.na(dev_with_fylke$fylke)]),
  "\n",
  sep = ""
)


# Show unmatched FIDs if there are any.
missing_fylke <- dev_with_fylke %>%
  filter(
    is.na(fylke)
  ) %>%
  distinct(
    fid
  )

if (nrow(missing_fylke) > 0) {
  
  cat(
    "\nFirst FIDs without fylke match:\n"
  )
  
  print(
    head(
      missing_fylke,
      20
    )
  )
}


# ------------------------------------------------
# 7. CONVERT ANNUAL FIELDS TO LONG FORMAT
# ------------------------------------------------

dev_fylke_long <- dev_with_fylke %>%
  filter(
    !is.na(fylke)
  ) %>%
  select(
    fid,
    fylke,
    all_of(annual_ha_fields)
  ) %>%
  pivot_longer(
    cols = all_of(annual_ha_fields),
    names_to = "year",
    values_to = "annual_loss_ha"
  ) %>%
  mutate(
    year = as.integer(
      stringr::str_extract(
        year,
        "[0-9]{4}"
      )
    )
  )


# ------------------------------------------------
# 8. SUM ANNUAL LOSS BY FYLKE
# ------------------------------------------------

fylke_annual_loss <- dev_fylke_long %>%
  group_by(
    fylke,
    year
  ) %>%
  summarise(
    annual_loss_ha = sum(
      annual_loss_ha,
      na.rm = TRUE
    ),
    n_fids_with_loss = n_distinct(
      fid[annual_loss_ha > 0]
    ),
    .groups = "drop"
  ) %>%
  arrange(
    fylke,
    year
  )


# ------------------------------------------------
# 9. CALCULATE CUMULATIVE AREA LOSS
# ------------------------------------------------

fylke_cumulative_loss <- fylke_annual_loss %>%
  group_by(
    fylke
  ) %>%
  arrange(
    year,
    .by_group = TRUE
  ) %>%
  mutate(
    cumulative_loss_ha = cumsum(
      annual_loss_ha
    )
  ) %>%
  ungroup()


# Optional zero baseline at 2017 so all lines clearly start at zero.
fylke_baseline <- fylke_cumulative_loss %>%
  distinct(
    fylke
  ) %>%
  mutate(
    year = 2017L,
    annual_loss_ha = 0,
    n_fids_with_loss = 0L,
    cumulative_loss_ha = 0
  )


fylke_cumulative_plot_dat <- bind_rows(
  fylke_baseline,
  fylke_cumulative_loss
) %>%
  arrange(
    fylke,
    year
  )


# ------------------------------------------------
# 10. SAVE SUMMARY TABLES
# ------------------------------------------------

readr::write_csv(
  fylke_annual_loss,
  file.path(OUTPUT_DIR, "MI_DW_annual_development_loss_by_fylke.csv")
)

readr::write_csv(
  fylke_cumulative_loss,
  file.path(OUTPUT_DIR, "MI_DW_cumulative_development_loss_by_fylke.csv")
)




# ------------------------------------------------
# 12. CUMULATIVE AREA LOSS GRAPH
# ------------------------------------------------
# Put the labels in one tidy right-hand column
fylke_end_labels <- fylke_cumulative_plot_dat %>%
  group_by(fylke) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  mutate(
    label_x = 2025.4,
    label = paste0(fylke, "  ", sprintf("%.1f ha", cumulative_loss_ha))
  )

p_cumulative_fylke <- ggplot(
  fylke_cumulative_plot_dat,
  aes(
    x = year,
    y = cumulative_loss_ha,
    colour = fylke,
    group = fylke
  )
) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 1.8
  ) +
  ggrepel::geom_text_repel(
    data = fylke_end_labels,
    aes(
      x = label_x,
      y = cumulative_loss_ha,
      label = label,
      colour = fylke
    ),
    inherit.aes = FALSE,
    hjust = 0,
    direction = "y",
    box.padding = 0.45,
    point.padding = 0.25,
    force = 2.5,
    force_pull = 0.15,
    min.segment.length = 0,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 42,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = 2017:2025,
    expand = expansion(
      mult = c(0.01, 0.24)
    )
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      accuracy = 1,
      big.mark = " "
    ),
    expand = expansion(
      mult = c(0, 0.06)
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    title = "Cumulative nature loss to development by fylke",
    x = NULL,
    y = "Cumulative detected area loss (ha)",
    colour = "Fylke",
    caption = "2017 is a zero baseline; detected annual losses begin in 2018."
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 170, 10, 10)
  )


print(
  p_cumulative_fylke
)


# ------------------------------------------------
# 13. SAVE FIGURE
# ------------------------------------------------

ggsave(
  file.path(OUTPUT_DIR, "MI_DW_cumulative_area_loss_by_fylke.png"),
  p_cumulative_fylke,
  width = 10,
  height = 10,
  dpi = 350,
  bg = "white"
)
