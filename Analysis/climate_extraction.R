# climate_extraction.R
# ---------------------------------------------------------------------------
# Attaches TIME-RESOLVED climate to every bird record, so body size can be
# modelled against the temperature and precipitation each individual actually
# experienced, rather than against "year" as a proxy.
#
# WHY (see CODE_REVIEW.md §B): the `Annual_mean_temperature` and
# `Annual_rainfall` fields in ATLANTIC BIRD TRAITS are STATIC per-locality
# climatologies (97% of localities carry a single value across all sampling
# years), so they cannot support a temporal test. We instead use TerraClimate
# (monthly, 1958-present, ~4 km; Abatzoglou et al. 2018).
#
# REPRODUCIBLE DOWNLOAD: climate is pulled straight from R with the climateR
# package (Johnson, mikejohnson51/climateR), which subsets TerraClimate on the
# server via OPeNDAP — only the needed cells/months are transferred, so there is
# NO multi-GB raster download and nothing to place by hand.
# install: pak::pak("mikejohnson51/climateR")   (or remotes::install_github)
#
# Output: passer90_climate.rds = passer90 + per-record rec_tmean / rec_ppt
# (and scaled versions), plus the reconciled `species_name` used by the model.
# ---------------------------------------------------------------------------

library(climateR)
library(sf)
library(dplyr)

# --- robust path resolution (find repo files regardless of getwd()) --------
.find_file <- function(fname, subdir = "Analysis") {
  dirs <- getwd(); d <- getwd()
  for (i in 1:8) { d <- dirname(d); dirs <- c(dirs, d) }
  cand <- unique(c(file.path(dirs, fname), file.path(dirs, subdir, fname)))
  hit  <- cand[file.exists(cand)]
  if (!length(hit)) stop("Could not locate '", fname,
      "'. Run from inside the atlantic_birds repo, or setwd() to the Analysis/ folder.")
  normalizePath(hit[1], winslash = "/")
}
ANALYSIS_DIR <- dirname(.find_file("passer90.rda"))
apath <- function(...) file.path(ANALYSIS_DIR, ...)

load(apath("passer90.rda"))

# Match the model's species labels (reconciled eBird names, underscored) so the
# saved frame's `species_name` lines up with the phylogeny used in the models.
passer90$species_name <- gsub("_", " ", passer90$Binomial)
ebird_synonyms <- c(
  "Antilophia galeata"       = "Chiroxiphia galeata",
  "Tachyphonus cristatus"    = "Loriotus cristatus",
  "Pyrrhocoma ruficeps"      = "Thlypopsis pyrrhocoma",
  "Pyriglena pernambucensis" = "Pyriglena leuconota",
  "Tangara sayaca"           = "Thraupis sayaca",
  "Tangara cayana"           = "Stilpnia cayana",
  "Tiaris fuliginosus"       = "Asemospiza fuliginosa"
)
hit <- passer90$species_name %in% names(ebird_synonyms)
passer90$species_name[hit] <- ebird_synonyms[passer90$species_name[hit]]
passer90$species_name <- gsub(" ", "_", passer90$species_name)

# --- unique sampling localities -> sf points -------------------------------
locs <- passer90 %>%
  filter(!is.na(Longitude_decimal_degrees), !is.na(Latitude_decimal_degrees)) %>%
  distinct(Longitude_decimal_degrees, Latitude_decimal_degrees) %>%
  mutate(loc_id = dplyr::row_number())
locs_sf <- sf::st_as_sf(locs,
                        coords = c("Longitude_decimal_degrees", "Latitude_decimal_degrees"),
                        crs = 4326)

# --- pull monthly TerraClimate for those points, 1990-2018 -----------------
# For point AOIs climateR returns the time series of the intersecting cell;
# passing ID keeps the per-locality identifier in the output.
clim_raw <- getTerraClim(
  AOI       = locs_sf,
  varname   = c("tmax", "tmin", "ppt"),
  startDate = "1990-01-01",
  endDate   = "2018-12-31",
  ID        = "loc_id"
)
# Defensive: inspect names(clim_raw) once if columns differ from the below.
# Expected columns: loc_id, date, tmax, tmin, ppt (monthly rows per locality).

clim_year <- clim_raw %>%
  mutate(year  = as.integer(format(as.Date(date), "%Y")),
         tmean = (tmax + tmin) / 2) %>%
  group_by(loc_id, year) %>%
  summarise(rec_tmean = mean(tmean, na.rm = TRUE),
            rec_ppt   = mean(ppt,   na.rm = TRUE),  # mean monthly precip; use sum() for annual total
            .groups = "drop")

# --- join climate back onto each record (by locality + capture year) -------
passer90 <- passer90 %>%
  left_join(locs, by = c("Longitude_decimal_degrees", "Latitude_decimal_degrees")) %>%
  left_join(clim_year, by = c("loc_id", "Year" = "year"))

# standardise climate predictors on a common scale (as in the lab manuscript)
passer90$scaled_tmean <- as.numeric(scale(passer90$rec_tmean))
passer90$scaled_ppt   <- as.numeric(scale(passer90$rec_ppt))
passer90$spp          <- passer90$species_name   # match the model's grouping terms

saveRDS(passer90, apath("passer90_climate.rds"))
message("Saved passer90_climate.rds with time-resolved TerraClimate per record.")

# ---------------------------------------------------------------------------
# Then the `model_climate` chunk (or a script) fits body size against ACTUAL
# climate, keeping year to separate climate tracking from other temporal change:
#
#   brm(conc.wing.length ~ 1 + Sex + scaled_tmean + scaled_ppt + scaled_lat +
#         scaled_yr + (1 + scaled_tmean || spp) + (1 | gr(species_name, cov = A)),
#       data = passer90_climate, data2 = list(A = A), ...)
#
# (A comes from the phylogeny chunk / atlantic_parallel.R, built on the same
#  reconciled species_name labels.)
# ---------------------------------------------------------------------------
