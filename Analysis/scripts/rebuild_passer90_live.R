# rebuild_passer90_live.R
# ---------------------------------------------------------------------------
# Regenerates Analysis/data/derived/passer90.rda as the LIVE-ONLY analytical
# sample, identically to the `newdataframes` chunk of atlantic_birds_ms.Rmd
# (same filters, mutations, column selection) but restricted to live-measured
# birds (Status == "live").
#
# WHY (2026-06): museum specimens shrink/shift on preservation (museum tail
# reads ~+5%, mass ~-1.8%) and are concentrated in the early years, so they can
# confound the temporal trend (see body_mass_descriptive.qmd, "museum vs live").
# Status == "live" is applied BEFORE the n>=30 / range>=5 species filter, so the
# sample stays internally consistent (every species has >=30 live records).
# Effect: 89 -> 73 species, 15,332 -> 12,571 records; scaled_yr / scaled_lat are
# re-standardised on the live-only sample.
#
# This is the data fed to atlantic_parallel.R (wing), atlantic_parallel_bill.R,
# atlantic_diet_interaction.R, climate_extraction.R (-> passer90_climate.rds,
# used by atlantic_drmsem.R) and update_descriptive_stats.R, all of which simply
# load(passer90.rda). Run this ONCE before those scripts. The mass and
# sensitivity scripts re-read the raw csv and carry the same Status filter
# inline, so they do not depend on this file.
#
# Run:  Rscript Analysis/scripts/rebuild_passer90_live.R
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(readr) })

# --- robust path resolution (repo layout: Analysis/{scripts,data,output}) ---
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts")))
      return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/. Run from inside the atlantic_birds repo.")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)

LIVE_ONLY <- TRUE   # the entire purpose of this script; FALSE reproduces the old 89-spp sample

birds <- read_csv(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"),
                  guess_max = 70000, show_col_types = FALSE)

passer90 <- birds %>%
  filter(Year >= 1990, Age == "Adult", Sex != "Unknown",
         Order == "Passeriformes",
         AtlanticForests_20km_Buffer == "inside the 20 km polygon")
if (LIVE_ONLY) passer90 <- passer90 %>% filter(Status == "live")

passer90 <- passer90 %>%
  mutate(
    conc.wing.length    = coalesce(Wing_length_right.mm., Wing_length_left.mm., Wing_length.mm.),
    ln_conc_wing_length = log(conc.wing.length),
    ln_wing_length      = log(Wing_length_right.mm.),
    scaled_lat          = scale(Latitude_decimal_degrees * -1),
    scaled_yr           = scale(Year),
    Binomial            = gsub(" ", "_", Binomial),
    spp                 = Binomial
  ) %>%
  group_by(Binomial) %>%
  filter(n() >= 30) %>%
  mutate(min = min(Year), max = max(Year), range = max - min) %>%
  filter(range >= 5) %>%
  ungroup() %>%
  dplyr::select(spp, Binomial, Sex, Status, Year, scaled_yr, Annual_mean_temperature,
                Latitude_decimal_degrees, scaled_lat, Longitude_decimal_degrees,
                Wing_length_right.mm., ln_wing_length, conc.wing.length,
                ln_conc_wing_length, Bill_width.mm.) %>%
  as.data.frame()

stopifnot(!LIVE_ONLY || all(passer90$Status == "live"))
cat(sprintf("passer90 (%s): %d records, %d species; wing non-NA %d, bill non-NA %d\n",
            if (LIVE_ONLY) "LIVE-only" else "all-status",
            nrow(passer90), length(unique(passer90$Binomial)),
            sum(!is.na(passer90$conc.wing.length)), sum(!is.na(passer90$Bill_width.mm.))))

save(passer90, file = derived_path("passer90.rda"))
cat("Wrote", derived_path("passer90.rda"), "\n")
