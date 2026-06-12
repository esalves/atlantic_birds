# explore_biotime.R
# Download BioTIME v2.0 from Zenodo and assess coverage of Brazilian
# terrestrial arthropods for 1990–2018, as a potential replacement or
# supplement for the PREDICTS arthropod node in the atlantic birds SEM.
#
# Runtime: first run ~5 min (166 MB RDS download + load); subsequent runs use cache.
# Output:  biotime_brazil_arthropods.rds  — filtered site × year table
#          biotime_summary.txt            — console-style summary
#
# Usage: source("explore_biotime.R")  or  Rscript explore_biotime.R

library(dplyr)
library(ggplot2)

# ── paths ─────────────────────────────────────────────────────────────────────
CACHE_DIR  <- "biotime_cache"
QUERY_FILE <- file.path(CACHE_DIR, "BioTIMEQuery.rds")    # v2.0 ships as .rds
META_FILE  <- file.path(CACHE_DIR, "BioTIMEMetadata.csv")
OUT_FILE   <- "biotime_brazil_arthropods.rds"
dir.create(CACHE_DIR, showWarnings = FALSE)

# ── Brazil bounding box ───────────────────────────────────────────────────────
BR_LAT <- c(-33.8, 5.3)
BR_LON <- c(-73.9, -28.8)

# ── 1. Discover files on Zenodo (v2.0, record 15222193) ──────────────────────
message("Querying Zenodo record 15222193 for BioTIME v2.0 files…")
zenodo_api <- "https://zenodo.org/api/records/15222193"
rec <- jsonlite::fromJSON(url(zenodo_api))
files_df <- as.data.frame(rec$files)
message("Files available in record:")
print(files_df[, c("key", "size")])

# ── 2. Download helper (skips if already cached) ──────────────────────────────
download_if_missing <- function(url, dest) {
  if (file.exists(dest)) {
    message("  already cached: ", basename(dest))
    return(invisible(dest))
  }
  message("  downloading: ", basename(dest), " (", round(file.size(dest) / 1e6, 1), " MB)…")
  download.file(url, dest, mode = "wb", quiet = FALSE)
  invisible(dest)
}

# Match files to local destinations.
# v2.0 query is an .rds; metadata is a .csv.
pick_url <- function(files_df, pattern) {
  row <- files_df[grep(pattern, files_df$key, ignore.case = TRUE), ]
  if (nrow(row) == 0) stop("No file matching '", pattern, "' in Zenodo record.")
  if (nrow(row) > 1) { row <- row[1, ]; message("Multiple matches; using: ", row$key[1]) }
  row$links$self   # Zenodo self-link for direct download
}

query_url <- pick_url(files_df, "query.*\\.rds")
meta_url  <- pick_url(files_df, "metadata")

download_if_missing(query_url, QUERY_FILE)
download_if_missing(meta_url,  META_FILE)

# ── 3. Load metadata (small) ──────────────────────────────────────────────────
message("\nLoading metadata…")
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
message("  ", nrow(meta), " studies in metadata")

# ── 4. Load the query RDS and filter to Brazil ────────────────────────────────
# v2.0 ships as a compressed RDS (~166 MB on disk). Load it fully then filter;
# it expands to a data.table or data.frame in memory.
message("\nLoading BioTIME query RDS (166 MB compressed – may take ~1 min)…")
q_all <- readRDS(QUERY_FILE)
message("  ", nrow(q_all), " total records loaded")
message("  columns: ", paste(names(q_all), collapse = ", "))

message("\nFiltering to Brazil bounding box…")
brazil_rows <- q_all %>%
  filter(
    LATITUDE  >= BR_LAT[1], LATITUDE  <= BR_LAT[2],
    LONGITUDE >= BR_LON[1], LONGITUDE <= BR_LON[2]
  )
rm(q_all); gc()

message("  ", nrow(brazil_rows), " records in Brazil bounding box")
message("  ", n_distinct(brazil_rows$STUDY_ID), " studies")

# ── 5. Join metadata for realm + taxa information ─────────────────────────────
brazil_rows <- brazil_rows %>%
  left_join(
    meta %>% select(STUDY_ID, REALM, TAXA, ORGANISMS, TITLE, BIOME_MAP,
                    starts_with("CENT_"), starts_with("NUMBER_")),
    by = "STUDY_ID"
  )

# ── 6. Filter: terrestrial realm, invertebrate/arthropod taxa ─────────────────
message("\nRealm × Taxa breakdown for Brazil:")
print(brazil_rows %>% count(REALM, TAXA, sort = TRUE))

arth_brazil <- brazil_rows %>%
  filter(
    REALM == "Terrestrial",
    grepl("invertebrate|arthropod|insect|beetle|butterfly|moth|bug|fly|bee|wasp|ant|spider",
          TAXA, ignore.case = TRUE) |
    grepl("invertebrate|arthropod|insect|beetle|butterfly|moth|bug|fly|bee|wasp|ant|spider",
          ORGANISMS, ignore.case = TRUE)
  )

message("\nAfter terrestrial + invertebrate filter: ", nrow(arth_brazil), " records")
message("Studies: ", n_distinct(arth_brazil$STUDY_ID))

# ── 7. Temporal coverage ──────────────────────────────────────────────────────
message("\nYear range in arthropod subset:")
print(range(arth_brazil$YEAR, na.rm = TRUE))

target_range <- arth_brazil %>%
  filter(YEAR >= 1990, YEAR <= 2018)
message("Records within 1990–2018: ", nrow(target_range))
message("Studies within 1990–2018: ", n_distinct(target_range$STUDY_ID))

# Per-study temporal span
study_span <- arth_brazil %>%
  group_by(STUDY_ID) %>%
  summarise(
    yr_min   = min(YEAR),
    yr_max   = max(YEAR),
    n_years  = n_distinct(YEAR),
    n_sites  = n_distinct(paste(LATITUDE, LONGITUDE)),
    n_records = n(),
    has_abundance = any(!is.na(ABUNDANCE)),
    has_biomass   = any(!is.na(BIOMASS)),
    .groups = "drop"
  ) %>%
  left_join(meta %>% select(STUDY_ID, TAXA, ORGANISMS, TITLE, BIOME_MAP), by = "STUDY_ID") %>%
  arrange(desc(n_years))

message("\nStudy temporal spans (sorted by number of years sampled):")
print(study_span, n = 30)

# Studies that overlap with 1990–2018 AND have multiple years
useful <- study_span %>%
  filter(yr_min <= 2018, yr_max >= 1990, n_years >= 3)
message("\nStudies overlapping 1990–2018 with ≥3 sampled years: ", nrow(useful))
print(useful)

# ── 8. Abundance vs biomass availability ──────────────────────────────────────
message("\nAbundance/biomass data:")
cat("  Has ABUNDANCE values:", sum(!is.na(arth_brazil$ABUNDANCE)), "\n")
cat("  Has BIOMASS values:  ", sum(!is.na(arth_brazil$BIOMASS)), "\n")

# ── 9. Spatial distribution ───────────────────────────────────────────────────
message("\nSpatial extent of arthropod records:")
cat("  Latitude range: ", paste(range(arth_brazil$LATITUDE, na.rm=TRUE), collapse=" – "), "\n")
cat("  Longitude range:", paste(range(arth_brazil$LONGITUDE, na.rm=TRUE), collapse=" – "), "\n")

# ── 10. Plots ─────────────────────────────────────────────────────────────────
dir.create("biotime_plots", showWarnings = FALSE)

# Records per year
p_years <- arth_brazil %>%
  filter(YEAR >= 1985, YEAR <= 2025) %>%
  count(YEAR) %>%
  ggplot(aes(x = YEAR, y = n)) +
  geom_col(fill = "steelblue") +
  geom_vline(xintercept = c(1990, 2018), linetype = "dashed", colour = "red") +
  labs(title = "BioTIME: Brazil terrestrial arthropod records per year",
       subtitle = "Red dashed lines = 1990–2018 bird dataset window",
       x = "Year", y = "N records") +
  theme_minimal()
ggsave("biotime_plots/records_per_year.png", p_years, width = 8, height = 4)
message("Saved: biotime_plots/records_per_year.png")

# Per-study year × site heatmap (for studies with ≥3 years)
if (nrow(useful) > 0) {
  heat_dat <- arth_brazil %>%
    filter(STUDY_ID %in% useful$STUDY_ID) %>%
    mutate(site_id = paste0(round(LATITUDE, 2), "_", round(LONGITUDE, 2))) %>%
    distinct(STUDY_ID, YEAR, site_id) %>%
    count(STUDY_ID, YEAR)

  p_heat <- ggplot(heat_dat, aes(x = YEAR, y = as.factor(STUDY_ID), fill = log1p(n))) +
    geom_tile() +
    geom_vline(xintercept = c(1990, 2018), linetype = "dashed", colour = "red") +
    scale_fill_viridis_c(name = "log(n sites)") +
    labs(title = "BioTIME: Brazil arthropod studies – temporal coverage",
         subtitle = "Each row = one study; red dashes = 1990–2018 window",
         x = "Year", y = "Study ID") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7))
  ggsave("biotime_plots/study_coverage_heatmap.png", p_heat, width = 10, height = max(4, nrow(useful) * 0.4))
  message("Saved: biotime_plots/study_coverage_heatmap.png")
}

# Site map
p_map <- ggplot(arth_brazil %>% distinct(LATITUDE, LONGITUDE, STUDY_ID),
                aes(x = LONGITUDE, y = LATITUDE, colour = as.factor(STUDY_ID))) +
  geom_point(alpha = 0.6, size = 1.5) +
  coord_fixed(xlim = BR_LON, ylim = BR_LAT) +
  labs(title = "BioTIME: Brazil terrestrial arthropod sampling sites",
       x = "Longitude", y = "Latitude", colour = "Study") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("biotime_plots/site_map.png", p_map, width = 7, height = 7)
message("Saved: biotime_plots/site_map.png")

# ── 11. Save filtered dataset ──────────────────────────────────────────────────
saveRDS(
  list(
    records      = arth_brazil,
    study_spans  = study_span,
    useful_studies = useful
  ),
  OUT_FILE
)
message("\nSaved filtered data to: ", OUT_FILE)

# ── 12. Summary verdict ───────────────────────────────────────────────────────
message("\n=== SUMMARY ===")
message("Total Brazil terrestrial arthropod records:  ", nrow(arth_brazil))
message("Total studies:                               ", n_distinct(arth_brazil$STUDY_ID))
message("Records within 1990–2018:                    ", nrow(target_range))
message("Studies with ≥3 years in 1990–2018 window:  ", nrow(useful))
if (nrow(useful) > 0) {
  message("\nUseful studies for temporal analysis:")
  useful %>%
    select(STUDY_ID, yr_min, yr_max, n_years, n_sites, TAXA, TITLE) %>%
    print(n = Inf)
  message("\nVerdict: BioTIME HAS usable Brazil arthropod time-series.")
  message("Next step: compare site coordinates with bird localities; check abundance column units.")
} else {
  message("\nVerdict: BioTIME has insufficient Brazil arthropod coverage for this analysis.")
  message("Recommendation: stick with PREDICTS + explore ATLANTIC ANTS / PELD data.")
}
