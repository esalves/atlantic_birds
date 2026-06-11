# climate_extraction.R
# ---------------------------------------------------------------------------
# Attaches TIME-RESOLVED climate to every bird record, so body size can be
# modelled against the temperature and precipitation each individual actually
# experienced, rather than against "year" as a proxy.
#
# WHY (see CODE_REVIEW.md §B): the `Annual_mean_temperature` and
# `Annual_rainfall` fields in ATLANTIC BIRD TRAITS are STATIC per-locality
# climatologies (97% of localities carry a single value across all sampling
# years), so they cannot support a temporal test. The lab's Bergmann/migration
# manuscript builds climate from WorldClim v2.1 rasters; WorldClim is a
# 1970-2000 long-term average and is fine for a SPATIAL/exposure question, but
# for our TEMPORAL question we need a monthly time series. We therefore use
# TerraClimate (monthly, 1958-present, ~4 km; Abatzoglou et al. 2018) — CHELSA
# monthly (Karger et al. 2017) or CRU TS / ERA5-Land are drop-in alternatives.
#
# Output: passer90 with columns rec_tmax, rec_tmin, rec_tmean (deg C) and
# rec_ppt (mm/month) for the record's coordinates and YEAR (and month, if the
# capture Date is available), saved as passer90_climate.rds.
# ---------------------------------------------------------------------------

library(terra)
library(dplyr)
library(lubridate)

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

# Need a capture month where possible; fall back to annual summary otherwise.
# (The ABT "Date" column is dropped during wrangling — re-join it from the raw
#  CSV if you want month-level matching. Here we use annual summaries.)
pts <- vect(passer90,
            geom = c("Longitude_decimal_degrees", "Latitude_decimal_degrees"),
            crs  = "EPSG:4326")

years <- sort(unique(na.omit(passer90$Year)))
years <- years[years >= 1990 & years <= 2018]

# ---------------------------------------------------------------------------
# TerraClimate yearly files (one NetCDF per variable per year):
#   tmax, tmin (deg C, monthly), ppt (mm, monthly).
# Download from https://www.climatologylab.org/terraclimate.html or via
# the THREDDS server. Set the directory that holds the *.nc files:
TC_DIR <- apath("terraclimate")   # <-- put the downloaded NetCDFs in Analysis/terraclimate/

read_var_year <- function(var, yr) {
  f <- file.path(TC_DIR, sprintf("TerraClimate_%s_%d.nc", var, yr))
  if (!file.exists(f)) stop("missing: ", f)
  rast(f)                       # 12 monthly layers
}

extract_year <- function(yr) {
  idx  <- which(passer90$Year == yr)
  p    <- pts[idx, ]
  tmax <- read_var_year("tmax", yr); tmin <- read_var_year("tmin", yr)
  ppt  <- read_var_year("ppt",  yr)
  tmean_m <- (tmax + tmin) / 2
  data.frame(
    row_id   = idx,
    rec_tmean = rowMeans(terra::extract(tmean_m, p)[ , -1], na.rm = TRUE),
    rec_tmax  = apply(terra::extract(tmax,  p)[ , -1], 1, max,  na.rm = TRUE),
    rec_tmin  = apply(terra::extract(tmin,  p)[ , -1], 1, min,  na.rm = TRUE),
    rec_ppt   = rowMeans(terra::extract(ppt, p)[ , -1], na.rm = TRUE)
  )
}

clim <- bind_rows(lapply(years, extract_year)) |> arrange(row_id)
passer90$rec_tmean <- NA_real_; passer90$rec_tmax <- NA_real_
passer90$rec_tmin  <- NA_real_; passer90$rec_ppt  <- NA_real_
passer90[clim$row_id, c("rec_tmean","rec_tmax","rec_tmin","rec_ppt")] <-
  clim[, c("rec_tmean","rec_tmax","rec_tmin","rec_ppt")]

# standardise climate predictors on a common scale (as in the lab manuscript)
passer90$scaled_tmean <- as.numeric(scale(passer90$rec_tmean))
passer90$scaled_ppt   <- as.numeric(scale(passer90$rec_ppt))

saveRDS(passer90, apath("passer90_climate.rds"))
message("Saved passer90_climate.rds with time-resolved climate per record.")

# ---------------------------------------------------------------------------
# Then refit the primary body-size model against ACTUAL climate (and keep year
# to separate climate tracking from other temporal change), e.g.:
#
#   brm(conc.wing.length ~ 1 + Sex + scaled_tmean + scaled_ppt + scaled_lat +
#         scaled_yr + (1 + scaled_tmean || spp) + (1 | gr(Binomial, cov = A)),
#       data = passer90, data2 = list(A = A), ...)
#
# This directly tests Bergmann-style temperature tracking and adds precipitation
# as the second climatic axis used in the lab manuscript.
# ---------------------------------------------------------------------------
