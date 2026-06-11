# climate_extraction.R
# ---------------------------------------------------------------------------
# Attaches TIME-RESOLVED climate to every bird record using WorldClim's
# HISTORICAL MONTHLY WEATHER data, so body size can be modelled against the
# temperature each individual actually experienced, rather than against "year"
# as a proxy.
#
# WHY WorldClim (not TerraClimate): the previous OPeNDAP-based TerraClimate
# extraction streamed from the U Idaho THREDDS host, which hangs or goes down.
# WorldClim historical monthly weather (CRU-TS 4.09 downscaled) are plain
# static files on a fast CDN (geodata.ucdavis.edu), so the flakiness disappears.
#
# IMPORTANT — WHICH WorldClim PRODUCT: this uses the *historical monthly weather*
# series (monthly tmin/tmax, 1950-2024, CRU-TS 4.09 downscaled with WorldClim
# 2.1; https://www.worldclim.org/data/monthlywth.html), which is TIME-RESOLVED.
# It is NOT the WorldClim v2.1 "historical climate" normals (1970-2000 averages)
# — those are static in time.
#
# RESOLUTION: 2.5 arc-minutes (~4.6 km).
#   Set RES <- "5m" or "10m" for a quarter / sixteenth of the download size.
#
# DOWNLOAD: WorldClim serves GLOBAL rasters (no server-side country subset). We
# download the decadal zips covering 1990-2018 once (tmin, tmax x three decades),
# then crop to a bounding box around the bird localities and extract per record.
# At 2.5m this is a few GB total; cached, so re-running is cheap.
#
# CAVEAT (note in Methods): WorldClim monthly's interannual signal derives from
# CRU-TS 4.09 (station-interpolated, sparse in the interior Atlantic Forest), so
# the temporal signal is somewhat smoothed — same caliber of limitation as
# TerraClimate, and fine for a broad temporal Bergmann test.
#
# OUTPUT: passer90_climate.rds = passer90 + per-record rec_tmean (and scaled
# version scaled_tmean), the reconciled `species_name`, and `spp` set to it.
#
# Cite: CRU-TS 4.09 (Harris et al. 2020) downscaled with WorldClim 2.1
#       (Fick & Hijmans 2017).
# ---------------------------------------------------------------------------

library(terra)     # install.packages("terra")
library(dplyr)

# --- config -----------------------------------------------------------------
RES      <- "2.5m"                     # "2.5m" (default), "5m", or "10m"
DECADES  <- c("1990-1999", "2000-2009", "2010-2019")  # cover 1990-2018
YEAR_MIN <- 1990L
YEAR_MAX <- 2018L
BBOX_PAD <- 0.5                        # degrees of padding around localities
WC_BASE  <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.09"

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

# --- match the model's species labels (reconciled eBird names, underscored) -
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

# --- unique sampling localities --------------------------------------------
locs <- passer90 %>%
  filter(!is.na(Longitude_decimal_degrees), !is.na(Latitude_decimal_degrees)) %>%
  distinct(Longitude_decimal_degrees, Latitude_decimal_degrees) %>%
  mutate(loc_id = dplyr::row_number())

pts  <- terra::vect(as.data.frame(locs),
                    geom = c("Longitude_decimal_degrees", "Latitude_decimal_degrees"),
                    crs  = "EPSG:4326")
bbox <- terra::ext(
  min(locs$Longitude_decimal_degrees) - BBOX_PAD,
  max(locs$Longitude_decimal_degrees) + BBOX_PAD,
  min(locs$Latitude_decimal_degrees)  - BBOX_PAD,
  max(locs$Latitude_decimal_degrees)  + BBOX_PAD
)

# --- download + cache helpers ----------------------------------------------
WC_DIR <- apath("worldclim")
if (!dir.exists(WC_DIR)) dir.create(WC_DIR, recursive = TRUE)
options(timeout = 3600)

# Resume-capable download (system curl if available, else download.file).
fetch <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 1e6) {
    message("  cached: ", basename(dest)); return(invisible())
  }
  message("  downloading ", basename(dest), " ...")
  if (nzchar(Sys.which("curl"))) {
    ret <- system(sprintf(
      'curl -L -C - --retry 5 --retry-delay 20 --connect-timeout 60 --max-time 3600 -o "%s" "%s"',
      dest, url))
  } else {
    ret <- tryCatch(
      utils::download.file(url, dest, mode = "wb", quiet = FALSE),
      error = function(e) 1L)
  }
  if (ret != 0 || !file.exists(dest) || file.size(dest) < 1e6) {
    if (file.exists(dest)) unlink(dest)
    stop("Download failed for ", basename(dest))
  }
}

# Build a cropped monthly SpatRaster for one variable across 1990-2018, plus the
# year/month each layer represents (parsed robustly from the .tif filenames).
get_var_monthly <- function(var) {
  tifs <- character(0)
  for (dec in DECADES) {
    stem    <- sprintf("wc2.1_cruts4.09_%s_%s_%s", RES, var, dec)
    zip_dst <- file.path(WC_DIR, paste0(stem, ".zip"))
    out_dir <- file.path(WC_DIR, stem)
    fetch(sprintf("%s/%s.zip", WC_BASE, stem), zip_dst)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
      utils::unzip(zip_dst, exdir = out_dir)
    }
    tifs <- c(tifs, list.files(out_dir, pattern = "\\.tif$",
                               recursive = TRUE, full.names = TRUE))
  }
  # Parse YYYY-MM from each file; keep only the modelling window.
  ym  <- regmatches(basename(tifs), regexpr("[0-9]{4}-[0-9]{2}", basename(tifs)))
  yr  <- as.integer(substr(ym, 1, 4))
  mo  <- as.integer(substr(ym, 6, 7))
  keep <- !is.na(yr) & yr >= YEAR_MIN & yr <= YEAR_MAX
  tifs <- tifs[keep]; yr <- yr[keep]; mo <- mo[keep]
  ord  <- order(yr, mo)                       # tidy chronological order
  list(r = terra::crop(terra::rast(tifs[ord]), bbox), yr = yr[ord], mo = mo[ord])
}

# Per locality x year aggregate: MEAN of the 12 monthly values.
extract_var_annual <- function(var) {
  message("WorldClim ", var, " (", RES, ") ...")
  s    <- get_var_monthly(var)
  vals <- terra::extract(s$r, pts, ID = FALSE)        # rows = localities, cols = months
  mat  <- as.matrix(vals)
  do.call(rbind, lapply(sort(unique(s$yr)), function(y) {
    cols <- which(s$yr == y)
    data.frame(loc_id = locs$loc_id, year = y,
               m = rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE))
  }))
}

tmax_y <- extract_var_annual("tmax")
tmin_y <- extract_var_annual("tmin")

# --- assemble per locality-year climate ------------------------------------
clim_year <- tmax_y %>% rename(tmax = m) %>%
  left_join(rename(tmin_y, tmin = m), by = c("loc_id", "year")) %>%
  mutate(rec_tmean = (tmax + tmin) / 2) %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
  select(loc_id, year, rec_tmean)

# --- join climate back onto each record (by locality + capture year) -------
passer90 <- passer90 %>%
  left_join(locs, by = c("Longitude_decimal_degrees", "Latitude_decimal_degrees")) %>%
  left_join(clim_year, by = c("loc_id", "Year" = "year"))

passer90$scaled_tmean <- as.numeric(scale(passer90$rec_tmean))
passer90$spp          <- passer90$species_name   # match the model's grouping terms

saveRDS(passer90, apath("passer90_climate.rds"))
message("Saved passer90_climate.rds with time-resolved WorldClim (",
        RES, ", CRU-TS 4.09 downscaled) per record.")

# Quick sanity check: time-resolved temperature MUST vary within a locality
# across years (unlike the static Annual_mean_temperature). Expect > 0.
.chk <- passer90 %>%
  filter(!is.na(rec_tmean)) %>%
  group_by(loc_id) %>%
  summarise(n_yr = n_distinct(Year), rng = diff(range(rec_tmean)), .groups = "drop") %>%
  filter(n_yr > 1)
message(sprintf(
  "Within-locality temperature range across years: median %.3f C, max %.3f C (n=%d multi-year localities).",
  stats::median(.chk$rng), max(.chk$rng), nrow(.chk)))
