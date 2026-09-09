# climate_extraction.R
# ---------------------------------------------------------------------------
# WHAT: Two things.
#   (1) PER RECORD (unchanged output schema): attaches TIME-RESOLVED climate to
#       every bird record from WorldClim's HISTORICAL MONTHLY WEATHER data, so
#       body size can be modelled against the temperature each individual
#       actually experienced rather than against "year" as a proxy.
#       -> data/derived/passer90_climate.rds
#   (2) PER LOCALITY (added 2026-09, revision Phase 5, REVISION_PLAN.md §1
#       "§3.4" and §3 "Phase 5"): for every sampled coordinate locality, the
#       annual mean temperature, the warm-quarter maximum temperature and a
#       12-month SPEI for 1990-2018, a linear trend 1995-2018 per locality and
#       pooled across localities, and a figure (Fig. S-climate).
#       -> output/climate_trends.rds, figures/climate_trends.png
#
# WHY WorldClim (not TerraClimate): the previous OPeNDAP-based TerraClimate
# extraction streamed from the U Idaho THREDDS host, which hangs or goes down.
# WorldClim historical monthly weather (CRU-TS 4.09 downscaled) are plain
# static files on a fast CDN (geodata.ucdavis.edu), so the flakiness disappears.
#
# IMPORTANT — WHICH WorldClim PRODUCT: this uses the *historical monthly weather*
# series (monthly tmin/tmax/prec, 1950-2024, CRU-TS 4.09 downscaled with
# WorldClim 2.1; https://www.worldclim.org/data/monthlywth.html), which is
# TIME-RESOLVED. It is NOT the WorldClim v2.1 "historical climate" normals
# (1970-2000 averages) — those are static in time.
#
# RESOLUTION: 2.5 arc-minutes (~4.6 km).
#   Set RES <- "5m" or "10m" for a quarter / sixteenth of the download size.
#
# DOWNLOAD: WorldClim serves GLOBAL rasters (no server-side country subset). We
# download the decadal zips covering 1990-2018 once (tmin, tmax, prec x three
# decades), then crop to a bounding box around the bird localities and extract
# per locality. At 2.5m this is a few GB total; cached under data/raw/worldclim,
# so re-running is cheap. The script NEVER downloads when the unzipped decade
# directory is already present.
#
# SPEI: the plan asks for SPEI from the TerraClimate files in
# data/raw/terraclimate/. SPEI needs monthly precipitation (ppt) and potential
# evapotranspiration (pet). At the time of writing that directory holds ONLY
# TerraClimate_tmax_1990.nc and a truncated TerraClimate_tmax_1991.nc, so the
# TerraClimate route is coded (spei_terraclimate(), expects
# TerraClimate_{ppt,pet}_YYYY.nc for every year) but cannot run. When those
# files are absent the script computes SPEI-12 from the WorldClim monthly series
# itself: PET by Hargreaves-Samani (FAO-56 eq. 52; tmin, tmax, extraterrestrial
# radiation from latitude), water balance D = P - PET, 12-month sums, and per
# calendar month and locality a GENERALIZED LOGISTIC distribution fitted from
# sample L-moments (unbiased PWMs; Hosking 1990; Hosking & Wallis 1997 A.7),
# standardised to N(0,1). This is the parameterisation SPEI::spei() uses
# (lmom::pelglo / cdfglo); the three-parameter log-logistic of Vicente-Serrano
# et al. (2010) is its kappa < 0 branch. The strict log-logistic (first run of
# this script) cannot be fitted when the sample L-skewness is negative
# (beta = 1/tau3 <= 0) and failed in 23 % of locality x month cells, dropping
# 71 of 357 localities from the SPEI trend; the generalized logistic covers both
# signs of skew by reflection, and the number of reflected (kappa > 0) cells is
# recorded. The reference period is the whole 1990-2018 window (29 years; a
# 30-year normal is not available inside the cached files). The output records
# which source was used (`spei_source`). If the SPEI package is installed the
# in-house calculation is cross-checked against SPEI::spei() on one locality.
#
# TREND INFERENCE: localities share the same year-to-year weather, so a
# random-slope model that treats 357 localities as independent replicates
# (`y ~ decade + (1 + decade || loc_id)`) gives a trend SE that is far too small
# (t ~ 47 for Tmean). The headline pooled trend therefore adds a year random
# intercept, `+ (1 | year)`, which absorbs the shared anomaly and leaves ~24
# effective time points; the naive rows are kept for comparison, and an OLS on
# the across-locality regional mean series (with the lag-1 residual
# autocorrelation) is reported as the simplest honest check.
#
# CAVEATS (note in Methods): WorldClim monthly's interannual signal derives from
# CRU-TS 4.09 (station-interpolated, sparse in the interior Atlantic Forest), so
# the temporal signal is smoothed — same caliber of limitation as TerraClimate,
# fine for describing regional trends. Exposure-window mismatch: capture-year
# temperature cannot cause a fixed adult wing length; the locality trends are
# CONTEXT for the study region and period, not a per-individual exposure.
# COASTAL MASK: 98 of 455 coordinate localities (772 records; 91 in Espirito
# Santo around Guarapari, 4 in SC, 2 in PB, 1 in RJ) fall in cells the WorldClim
# land mask leaves empty, ~10 km from the nearest valid cell. They are NA in
# both outputs (as in every previous build of passer90_climate.rds; no
# nearest-cell fill is applied so the per-record file is unchanged in content)
# and are listed with their distance to the nearest valid cell in
# climate_trends.rds$no_raster_localities.
#
# OUTPUT (1): passer90_climate.rds = passer90 (all columns of the current
#   build) + loc_id, rec_tmean (annual mean of (tmax+tmin)/2 in the capture
#   year at the record's coordinates), scaled_tmean, the reconciled
#   `species_name`, and `spp` set to it. Schema unchanged from the previous
#   version except that passer90 now carries the 42 columns of the 2026-09 build.
# OUTPUT (2): climate_trends.rds — list: annual (locality x year table of
#   tmean, tmax_warmq, prec_annual, spei12_dec, spei12_mean), locality_slopes
#   (OLS per locality, 1995-2018, per decade), slope_summary, pooled (lme4
#   random-slope model per variable: unweighted / record-weighted, each with
#   and without the `(1 | year)` term; use the "+ year RE" rows), regional_mean
#   (year series of the across-locality mean), regional_trend (OLS on that
#   series), record_level_cor (cor(Year, rec_tmean) at the records, total and
#   split within / between locality), no_raster_localities, n_localities,
#   n_records_by_locality, spei_source, spei_* diagnostics, settings, session.
#   climate_trends.png — four panels (Tmean, warm-quarter Tmax, SPEI-12,
#   distribution of locality slopes).
#
# INPUTS:  data/derived/passer90.rda; data/raw/worldclim/ (cached rasters);
#          data/raw/terraclimate/ (optional ppt/pet for SPEI)
# RUN:     Rscript Analysis/scripts/climate_extraction.R
#          (from the repo root, Analysis/, or Analysis/scripts/; ~5-10 min with
#          cached rasters, dominated by cropping 3 x 348 global layers)
#
# Cite: CRU-TS 4.09 (Harris et al. 2020) downscaled with WorldClim 2.1
#       (Fick & Hijmans 2017); SPEI: Vicente-Serrano et al. 2010; Beguería &
#       Vicente-Serrano (SPEI R package); L-moments: Hosking 1990, Hosking &
#       Wallis 1997; Hargreaves & Samani 1985; Allen et al. 1998 (FAO-56).
# Session (2026-09-09, run locally in full on the cached rasters): R 4.6.0,
#   terra 1.9.27, dplyr 1.2.1, lme4 2.0.1, ggplot2 4.0.3, patchwork 1.3.2,
#   tidyr 1.3.2, SPEI 1.8.1 (cross-check only).
# ---------------------------------------------------------------------------

suppressMessages({
  library(terra)     # install.packages("terra")
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(ggplot2)
  library(patchwork)
})
set.seed(20260909)   # nothing stochastic below; recorded for the convention

# --- config -----------------------------------------------------------------
RES      <- "2.5m"                     # "2.5m" (default), "5m", or "10m"
DECADES  <- c("1990-1999", "2000-2009", "2010-2019")  # cover 1990-2018
YEAR_MIN <- 1990L
YEAR_MAX <- 2018L
TREND_MIN <- 1995L                     # trend window (the bird record is 1995-2018)
TREND_MAX <- 2018L
BBOX_PAD <- 0.5                        # degrees of padding around localities
WC_BASE  <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.09"
SPEI_SCALE <- 12L                      # months of accumulation for the SPEI

# --- robust path resolution (repo layout: Analysis/{scripts,data,output,figures}) ---
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts")))
      return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/ (need Analysis/data/derived and Analysis/scripts). ",
       "Run from inside the atlantic_birds repo.")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
if (!dir.exists(fig_path())) dir.create(fig_path(), recursive = TRUE)

.pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
session_info <- list(date = as.character(Sys.Date()), R = R.version.string,
                     packages = sapply(c("terra", "dplyr", "tidyr", "lme4", "ggplot2",
                                         "patchwork", "SPEI"), .pkg_ver))

load(derived_path("passer90.rda"))
message(sprintf("passer90: %d records, %d species (live-only build)",
                nrow(passer90), dplyr::n_distinct(passer90$Binomial)))

# --- match the model's species labels (reconciled eBird names, underscored) -
passer90$species_name <- gsub("_", " ", passer90$Binomial)
ebird_synonyms <- c(
  "Antilophia galeata"       = "Chiroxiphia galeata",
  "Tachyphonus cristatus"    = "Loriotus cristatus",
  "Pyrrhocoma ruficeps"      = "Thlypopsis pyrrhocoma",
  "Pyriglena pernambucensis" = "Pyriglena leuconota",
  "Tangara sayaca"           = "Thraupis sayaca",
  "Tangara cayana"           = "Stilpnia cayana",
  "Tiaris fuliginosus"       = "Asemospiza fuliginosa",
  "Tangara palmarum"         = "Thraupis palmarum",    # Palm Tanager stays in Thraupis
  "Tangara peruviana"        = "Stilpnia peruviana",   # Black-backed Tanager (eBird Stilpnia split)
  "Dixiphia pipra"           = "Pseudopipra pipra"     # White-crowned Manakin (SACC 876)
  # Herpsilochmus sellowi (Caatinga Antwren) left unmapped on purpose: absent
  # from the clootl eBird/Clements taxonomy used by the phylo step; relabelling
  # it here would only desync the grouping label from that decision.
)
hit <- passer90$species_name %in% names(ebird_synonyms)
passer90$species_name[hit] <- ebird_synonyms[passer90$species_name[hit]]
passer90$species_name <- gsub(" ", "_", passer90$species_name)

# --- unique sampling localities --------------------------------------------
locs <- passer90 %>%
  filter(!is.na(Longitude_decimal_degrees), !is.na(Latitude_decimal_degrees)) %>%
  distinct(Longitude_decimal_degrees, Latitude_decimal_degrees) %>%
  mutate(loc_id = dplyr::row_number())
message(sprintf("%d unique coordinate localities", nrow(locs)))

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
WC_DIR <- raw_path("worldclim")
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
# An already-unzipped decade directory is used as is (no zip needed).
get_var_monthly <- function(var) {
  tifs <- character(0)
  for (dec in DECADES) {
    stem    <- sprintf("wc2.1_cruts4.09_%s_%s_%s", RES, var, dec)
    zip_dst <- file.path(WC_DIR, paste0(stem, ".zip"))
    out_dir <- file.path(WC_DIR, stem)
    if (!dir.exists(out_dir) || length(list.files(out_dir, pattern = "\\.tif$", recursive = TRUE)) < 120) {
      fetch(sprintf("%s/%s.zip", WC_BASE, stem), zip_dst)
      if (!dir.exists(out_dir)) dir.create(out_dir)
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

# Locality x month matrix for one variable (rows = localities in `locs` order).
extract_var_monthly <- function(var) {
  message("WorldClim ", var, " (", RES, ") ...")
  s    <- get_var_monthly(var)
  vals <- terra::extract(s$r, pts, ID = FALSE)        # rows = localities, cols = months
  list(mat = as.matrix(vals), yr = s$yr, mo = s$mo)
}

# Per locality x year aggregate: MEAN of the 12 monthly values (identical to the
# previous extract_var_annual(); now derived from the retained monthly matrix).
annual_mean_from_monthly <- function(m) {
  do.call(rbind, lapply(sort(unique(m$yr)), function(y) {
    cols <- which(m$yr == y)
    data.frame(loc_id = locs$loc_id, year = y,
               m = rowMeans(m$mat[, cols, drop = FALSE], na.rm = TRUE))
  }))
}

tmax_m <- extract_var_monthly("tmax")
tmin_m <- extract_var_monthly("tmin")
stopifnot(identical(tmax_m$yr, tmin_m$yr), identical(tmax_m$mo, tmin_m$mo))
tmax_y <- annual_mean_from_monthly(tmax_m)
tmin_y <- annual_mean_from_monthly(tmin_m)

# ============================================================================
# (1) PER-RECORD OUTPUT — unchanged schema
# ============================================================================
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

saveRDS(passer90, derived_path("passer90_climate.rds"))
message(sprintf("Saved passer90_climate.rds (%d records, %d species) with time-resolved WorldClim (%s, CRU-TS 4.09 downscaled) per record; rec_tmean NA for %d records.",
                nrow(passer90), dplyr::n_distinct(passer90$Binomial), RES, sum(is.na(passer90$rec_tmean))))

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

# ============================================================================
# (2) LOCALITY-LEVEL CLIMATE TRENDS 1995-2018 (added 2026-09)
# ============================================================================
message("\n=== Locality-level climate series and trends ===")
n_by_loc <- passer90 %>% count(loc_id, name = "n_records")

# --- annual mean temperature (same quantity as rec_tmean) -------------------
tmean_annual <- clim_year %>% rename(tmean = rec_tmean)

# --- warm-quarter maximum temperature ---------------------------------------
# Mean monthly Tmax over the warmest consecutive three-month window whose
# CENTRAL month falls in the year (so a DJF warm quarter straddling New Year is
# assigned to the January year). Southern-hemisphere summers are DJF, which is
# why a calendar-year quarter split would be wrong here.
running3 <- function(x) {                      # centred 3-month running mean
  n <- length(x); out <- rep(NA_real_, n)
  if (n >= 3) out[2:(n - 1)] <- (x[1:(n - 2)] + x[2:(n - 1)] + x[3:n]) / 3
  out
}
tmax_run <- t(apply(tmax_m$mat, 1, running3))  # localities x months
tmax_warmq <- do.call(rbind, lapply(sort(unique(tmax_m$yr)), function(y) {
  cols <- which(tmax_m$yr == y)
  data.frame(loc_id = locs$loc_id, year = y,
             tmax_warmq = apply(tmax_run[, cols, drop = FALSE], 1, function(v)
               if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)))
}))

# --- precipitation (WorldClim prec, mm/month) --------------------------------
prec_m <- extract_var_monthly("prec")
stopifnot(identical(prec_m$yr, tmax_m$yr))
prec_annual <- do.call(rbind, lapply(sort(unique(prec_m$yr)), function(y) {
  cols <- which(prec_m$yr == y)
  data.frame(loc_id = locs$loc_id, year = y,
             prec_annual = rowSums(prec_m$mat[, cols, drop = FALSE], na.rm = FALSE))
}))

# --- SPEI --------------------------------------------------------------------
# Hargreaves-Samani PET (mm/month) from monthly tmin/tmax and latitude.
# Ra (extraterrestrial radiation, MJ m-2 d-1) from FAO-56 eq. 21-25 at the
# mid-month day of year; PET (mm/d) = 0.0023 * 0.408 * Ra * (Tmean + 17.8) * sqrt(Tmax - Tmin).
hargreaves_pet <- function(tmin, tmax, lat_deg, yr, mo) {
  # tmin/tmax: localities x months; lat_deg: vector per locality
  stopifnot(ncol(tmin) == length(yr), length(lat_deg) == nrow(tmin))
  days_in <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  leap    <- (yr %% 4 == 0 & yr %% 100 != 0) | yr %% 400 == 0
  ndays   <- days_in[mo] + ifelse(leap & mo == 2, 1, 0)
  doy_mid <- sapply(seq_along(mo), function(i) {   # day of year at mid-month
    d <- days_in; if (leap[i]) d[2] <- 29
    sum(d[seq_len(mo[i] - 1)]) + d[mo[i]] / 2
  })
  lat <- lat_deg * pi / 180
  Ra  <- matrix(NA_real_, nrow(tmin), ncol(tmin))
  for (j in seq_along(mo)) {
    dr    <- 1 + 0.033 * cos(2 * pi * doy_mid[j] / 365)
    delta <- 0.409 * sin(2 * pi * doy_mid[j] / 365 - 1.39)
    ws    <- acos(pmin(1, pmax(-1, -tan(lat) * tan(delta))))
    Ra[, j] <- (24 * 60 / pi) * 0.0820 * dr *
      (ws * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(ws))
  }
  tmean <- (tmin + tmax) / 2
  trange <- pmax(tmax - tmin, 0)
  pet_d <- 0.0023 * 0.408 * Ra * (tmean + 17.8) * sqrt(trange)
  pet_d <- pmax(pet_d, 0)
  sweep(pet_d, 2, ndays, `*`)                   # mm/month
}

# SPEI from a monthly water-balance matrix (localities x months): k-month sums,
# then per calendar month a generalized logistic distribution fitted from the
# sample L-moments (unbiased PWMs, Hosking 1990; parameters as in Hosking &
# Wallis 1997, Appendix A.7 — the same estimator as lmom::pelglo, which
# SPEI::spei(distribution = "log-Logistic", fit = "ub-pwm") calls), then the
# standard normal quantile of the fitted CDF. The generalized logistic with
# shape kappa is the three-parameter log-logistic of Vicente-Serrano et al.
# (2010) when kappa < 0 (positive L-skewness, beta = -1/kappa) and its
# reflection when kappa > 0; a strict log-logistic cannot be fitted to a
# negatively skewed sample, which is why the previous version of this
# function lost 23 % of the locality x calendar-month cells.
spei_from_balance <- function(D, yr, mo, scale = 12L) {
  n <- ncol(D)
  acc <- matrix(NA_real_, nrow(D), n)
  for (j in scale:n) acc[, j] <- rowSums(D[, (j - scale + 1):j, drop = FALSE])
  out <- matrix(NA_real_, nrow(D), n)
  n_fail <- 0L; n_reflected <- 0L
  # Generalized-logistic parameters (xi, alpha, kappa) from sorted sample xs.
  glo_pars <- function(xs) {
    nn <- length(xs); ii <- seq_len(nn)
    b0 <- mean(xs)                                             # unbiased PWMs b_r
    b1 <- sum((ii - 1) / (nn - 1) * xs) / nn
    b2 <- sum((ii - 1) * (ii - 2) / ((nn - 1) * (nn - 2)) * xs) / nn
    l1 <- b0; l2 <- 2 * b1 - b0; l3 <- 6 * b2 - 6 * b1 + b0    # L-moments
    if (!is.finite(l2) || l2 <= 0) return(NULL)
    t3 <- l3 / l2                                              # L-skewness, |t3| < 1
    if (!is.finite(t3) || abs(t3) >= 1) return(NULL)
    k <- -t3
    if (abs(k) < 1e-8) {                                       # logistic limit
      a <- l2; xi <- l1
    } else {
      a  <- l2 * sin(k * pi) / (k * pi)
      xi <- l1 - a * (1 / k - pi / sin(k * pi))
    }
    if (!is.finite(a) || a <= 0) return(NULL)
    list(xi = xi, a = a, k = k)
  }
  # CDF: F(x) = 1 / (1 + exp(-y)), y = -log(1 - k (x - xi) / a) / k (k != 0).
  # Outside the finite bound of the support F is 1 (k > 0, upper bound) or 0
  # (k < 0, lower bound).
  glo_cdf <- function(x, p) {
    if (abs(p$k) < 1e-8) {
      y <- (x - p$xi) / p$a
    } else {
      arg <- 1 - p$k * (x - p$xi) / p$a
      y <- ifelse(arg > 0, -log(pmax(arg, 1e-300)) / p$k, if (p$k > 0) Inf else -Inf)
    }
    1 / (1 + exp(-y))
  }
  for (m in 1:12) {
    cols <- which(mo == m & seq_len(n) >= scale)
    for (i in seq_len(nrow(D))) {
      x <- acc[i, cols]
      ok <- !is.na(x); if (sum(ok) < 8) { n_fail <- n_fail + 1L; next }
      p <- glo_pars(sort(x[ok]))
      if (is.null(p)) { n_fail <- n_fail + 1L; next }
      if (p$k > 0) n_reflected <- n_reflected + 1L
      Fx <- pmin(pmax(glo_cdf(x, p), 1e-6), 1 - 1e-6)     # keep qnorm finite
      out[i, cols] <- qnorm(Fx)
    }
  }
  attr(out, "n_failed_cells")    <- n_fail        # locality x calendar-month fits that failed (all-NA cells or < 8 years)
  attr(out, "n_reflected_cells") <- n_reflected   # cells with kappa > 0 (negative L-skewness; not fittable by a strict log-logistic)
  out
}

# TerraClimate route (preferred by the plan): expects TerraClimate_ppt_YYYY.nc
# and TerraClimate_pet_YYYY.nc (global 1/24 deg, 12 monthly layers each) for
# every year in YEAR_MIN:YEAR_MAX under data/raw/terraclimate/.
TC_DIR <- raw_path("terraclimate")
tc_files <- function(var) file.path(TC_DIR, sprintf("TerraClimate_%s_%d.nc", var, YEAR_MIN:YEAR_MAX))
tc_ok <- function(f) file.exists(f) & file.size(f) > 1e6
tc_available <- all(tc_ok(tc_files("ppt"))) && all(tc_ok(tc_files("pet")))
spei_terraclimate <- function() {
  read_var <- function(var) {
    r <- terra::rast(tc_files(var))              # 12 layers per file, chronological
    vals <- terra::extract(terra::crop(r, bbox), pts, ID = FALSE)
    as.matrix(vals)
  }
  P <- read_var("ppt"); PET <- read_var("pet")
  yr <- rep(YEAR_MIN:YEAR_MAX, each = 12); mo <- rep(1:12, times = length(YEAR_MIN:YEAR_MAX))
  list(spei = spei_from_balance(P - PET, yr, mo, SPEI_SCALE), yr = yr, mo = mo,
       source = "TerraClimate ppt - pet (Penman-Monteith PET, Abatzoglou et al. 2018)")
}
spei_worldclim <- function() {
  PET <- hargreaves_pet(tmin_m$mat, tmax_m$mat, locs$Latitude_decimal_degrees, tmax_m$yr, tmax_m$mo)
  list(spei = spei_from_balance(prec_m$mat - PET, prec_m$yr, prec_m$mo, SPEI_SCALE),
       yr = prec_m$yr, mo = prec_m$mo, pet = PET,
       source = "WorldClim 2.1 / CRU-TS 4.09 monthly prec - Hargreaves-Samani PET (tmin, tmax, latitude)")
}
if (tc_available) {
  message("SPEI: TerraClimate ppt and pet files present -> using TerraClimate.")
  sp <- spei_terraclimate()
} else {
  present <- list.files(TC_DIR, pattern = "\\.nc$")
  message("SPEI: TerraClimate ppt/pet files NOT present in data/raw/terraclimate/ (found: ",
          if (length(present)) paste(present, collapse = ", ") else "nothing",
          "). Falling back to WorldClim prec + Hargreaves PET.")
  sp <- spei_worldclim()
}
spei_source <- sp$source
spei_failed_cells    <- attr(sp$spei, "n_failed_cells")
spei_reflected_cells <- attr(sp$spei, "n_reflected_cells")
message(sprintf("SPEI-%d computed for %d localities; of %d locality x calendar-month generalized-logistic fits, %d had kappa > 0 (negative L-skewness) and %d failed (left NA; expected = 12 x localities without raster data).",
                SPEI_SCALE, nrow(sp$spei), 12L * nrow(sp$spei), spei_reflected_cells, spei_failed_cells))
# Optional cross-check of the in-house SPEI against the SPEI package (one locality).
spei_check <- NULL
if (requireNamespace("SPEI", quietly = TRUE) && !tc_available) {
  i <- which.max(n_by_loc$n_records[match(locs$loc_id, n_by_loc$loc_id)])
  D <- prec_m$mat[i, ] - sp$pet[i, ]
  ref <- tryCatch(as.numeric(SPEI::spei(ts(D, start = c(YEAR_MIN, 1), frequency = 12),
                                        scale = SPEI_SCALE, verbose = FALSE)$fitted),
                  error = function(e) NULL)
  if (!is.null(ref)) {
    ok <- is.finite(ref) & is.finite(sp$spei[i, ])
    spei_check <- list(loc_id = locs$loc_id[i], r = cor(ref[ok], sp$spei[i, ok]),
                       max_abs_diff = max(abs(ref[ok] - sp$spei[i, ok])))
    message(sprintf("SPEI cross-check vs SPEI::spei (loc %d): r = %.4f, max |diff| = %.3f",
                    spei_check$loc_id, spei_check$r, spei_check$max_abs_diff))
  }
}
spei_annual <- do.call(rbind, lapply(sort(unique(sp$yr)), function(y) {
  cols <- which(sp$yr == y)
  dec  <- which(sp$yr == y & sp$mo == 12)
  data.frame(loc_id = locs$loc_id, year = y,
             spei12_dec  = if (length(dec)) sp$spei[, dec] else NA_real_,   # calendar-year balance
             spei12_mean = rowMeans(sp$spei[, cols, drop = FALSE], na.rm = TRUE))
}))
spei_annual$spei12_mean[!is.finite(spei_annual$spei12_mean)] <- NA_real_

# --- assemble the locality x year table -------------------------------------
annual <- tmean_annual %>%
  left_join(tmax_warmq, by = c("loc_id", "year")) %>%
  left_join(prec_annual, by = c("loc_id", "year")) %>%
  left_join(spei_annual, by = c("loc_id", "year")) %>%
  left_join(locs, by = "loc_id") %>%
  left_join(n_by_loc, by = "loc_id") %>%
  arrange(loc_id, year)

# --- trends: per locality OLS (1995-2018) and pooled random-slope models -----
trend_vars <- c("tmean", "tmax_warmq", "spei12_dec", "prec_annual")
trend_dat <- annual %>% filter(year >= TREND_MIN, year <= TREND_MAX) %>%
  mutate(year_c = year - mean(c(TREND_MIN, TREND_MAX)), decade = year_c / 10)

ols_slope <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 10) return(c(slope = NA_real_, se = NA_real_, n = sum(ok)))
  f <- lm(y[ok] ~ x[ok]); cf <- summary(f)$coefficients
  c(slope = cf[2, 1], se = cf[2, 2], n = sum(ok))
}
locality_slopes <- bind_rows(lapply(trend_vars, function(v) {
  trend_dat %>% group_by(loc_id, n_records, Longitude_decimal_degrees, Latitude_decimal_degrees) %>%
    summarise(res = list(ols_slope(.data[[v]], decade)), .groups = "drop") %>%
    mutate(variable = v, slope_per_decade = sapply(res, `[[`, "slope"),
           se = sapply(res, `[[`, "se"), n_years = sapply(res, `[[`, "n")) %>%
    select(-res)
}))
slope_summary <- locality_slopes %>% group_by(variable) %>%
  summarise(n_localities = sum(is.finite(slope_per_decade)),
            mean = mean(slope_per_decade, na.rm = TRUE),
            median = median(slope_per_decade, na.rm = TRUE),
            q025 = quantile(slope_per_decade, 0.025, na.rm = TRUE),
            q975 = quantile(slope_per_decade, 0.975, na.rm = TRUE),
            prop_positive = mean(slope_per_decade > 0, na.rm = TRUE),
            record_weighted_mean = weighted.mean(slope_per_decade, n_records, na.rm = TRUE),
            .groups = "drop")
message("Per-locality OLS slopes per decade, 1995-2018:")
print(as.data.frame(slope_summary), digits = 3)

# Pooled: y ~ decade + (1 + decade || loc_id) [+ (1 | year)]; the fixed slope
# is the mean locality trend. WITHOUT the year term the SE only respects
# locality-level heterogeneity and treats the localities as independent
# replicates of the same 24 years, which they are not (they share the regional
# weather), so t-values of ~47 result. WITH `(1 | year)` the shared annual
# anomaly is absorbed and the SE reflects ~24 effective time points; those rows
# (weighting ending in "+ year RE") are the ones to quote. The record-weighted
# variant weights localities by the number of bird records (climate at the
# places where most birds were measured).
pooled_fit <- function(v, weighted = FALSE, year_re = FALSE) {
  d <- trend_dat %>% filter(is.finite(.data[[v]])) %>% mutate(year_f = factor(year))
  rhs <- "decade + (1 + decade || loc_id)"
  if (year_re) rhs <- paste(rhs, "+ (1 | year_f)")
  fml <- as.formula(paste(v, "~", rhs))
  f <- if (weighted) {
    d$w <- d$n_records / mean(d$n_records)
    lmer(fml, data = d, REML = TRUE, weights = w, control = lmerControl(optimizer = "bobyqa"))
  } else {
    lmer(fml, data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  }
  cf <- summary(f)$coefficients
  vc <- as.data.frame(VarCorr(f))
  data.frame(variable = v,
             weighting = paste0(if (weighted) "record-weighted" else "unweighted",
                                if (year_re) " + year RE" else ""),
             year_random_effect = year_re,
             intercept = cf["(Intercept)", "Estimate"],
             slope_per_decade = cf["decade", "Estimate"], se = cf["decade", "Std. Error"],
             t = cf["decade", "t value"],
             lower = cf["decade", "Estimate"] - 1.96 * cf["decade", "Std. Error"],
             upper = cf["decade", "Estimate"] + 1.96 * cf["decade", "Std. Error"],
             sd_locality_slope = vc$sdcor[vc$var1 %in% "decade"][1],
             sd_year = if (year_re) vc$sdcor[vc$grp == "year_f"][1] else NA_real_,
             sd_residual = vc$sdcor[vc$grp == "Residual"][1],
             n_localities = dplyr::n_distinct(d$loc_id), n_years = dplyr::n_distinct(d$year),
             n_obs = nrow(d), singular = isSingular(f), row.names = NULL)
}
pooled <- bind_rows(lapply(trend_vars, function(v)
  bind_rows(pooled_fit(v, FALSE, FALSE), pooled_fit(v, TRUE, FALSE),
            pooled_fit(v, FALSE, TRUE),  pooled_fit(v, TRUE, TRUE))))
message("Pooled random-slope trends per decade, 1995-2018 (quote the '+ year RE' rows):")
print(pooled %>% select(variable, weighting, slope_per_decade, se, t, lower, upper,
                        sd_locality_slope, sd_year, n_localities, singular), digits = 3)

# Total change over the trend window implied by the pooled slope
pooled <- pooled %>% mutate(change_over_window = slope_per_decade * (TREND_MAX - TREND_MIN) / 10,
                            change_lower = lower * (TREND_MAX - TREND_MIN) / 10,
                            change_upper = upper * (TREND_MAX - TREND_MIN) / 10)

# Regional mean series (across localities with data, each locality once)
regional_mean <- annual %>% filter(is.finite(tmean)) %>% group_by(year) %>%
  summarise(across(all_of(c("tmean", "tmax_warmq", "prec_annual", "spei12_dec", "spei12_mean")),
                   ~ mean(.x, na.rm = TRUE)), n_localities = dplyr::n_distinct(loc_id), .groups = "drop")

# OLS trend on the regional-mean series (one value per year; the simplest
# inference that does not pretend localities are independent). t-based CI;
# lag-1 residual autocorrelation reported so the reader can judge the SE.
regional_trend <- bind_rows(lapply(trend_vars, function(v) {
  d <- regional_mean %>% filter(year >= TREND_MIN, year <= TREND_MAX, is.finite(.data[[v]]))
  d$decade <- (d$year - mean(c(TREND_MIN, TREND_MAX))) / 10
  f  <- lm(reformulate("decade", v), data = d)
  cf <- summary(f)$coefficients; r <- residuals(f)
  q  <- qt(0.975, f$df.residual)
  data.frame(variable = v, slope_per_decade = cf["decade", 1], se = cf["decade", 2],
             t = cf["decade", 3], p = cf["decade", 4],
             lower = cf["decade", 1] - q * cf["decade", 2], upper = cf["decade", 1] + q * cf["decade", 2],
             change_over_window = cf["decade", 1] * (TREND_MAX - TREND_MIN) / 10,
             n_years = nrow(d), resid_lag1_autocor = cor(r[-1], r[-length(r)]), row.names = NULL)
}))
message("OLS trend on the regional-mean series, 1995-2018 (per decade):")
print(regional_trend, digits = 3)

# Correlation of year with tmean AT THE RECORDS (is "later" also "warmer" in the
# bird sample?), split into the within-locality part (the warming a fixed site
# experienced) and the between-locality part (later sampling happening at
# warmer places: REVISION_NOTES_P0 reports mean latitude drifting from -24.8 to
# -21.9 and r(year, lat) = 0.17).
rec_cor <- passer90 %>% filter(!is.na(rec_tmean)) %>%
  group_by(loc_id) %>%
  mutate(yr_w = Year - mean(Year), t_w = rec_tmean - mean(rec_tmean),
         yr_b = mean(Year), t_b = mean(rec_tmean)) %>% ungroup() %>%
  summarise(r_year_tmean = cor(Year, rec_tmean),
            r_within_locality = cor(yr_w, t_w),
            r_between_locality_record_weighted = cor(yr_b, t_b),
            n = n(), n_localities = dplyr::n_distinct(loc_id))
.loc_means <- passer90 %>% filter(!is.na(rec_tmean)) %>% group_by(loc_id) %>%
  summarise(yr_b = mean(Year), t_b = mean(rec_tmean), .groups = "drop")
rec_cor$r_between_locality_unweighted <- cor(.loc_means$yr_b, .loc_means$t_b)
message(sprintf("Record-level cor(Year, rec_tmean) = %.3f (n = %d); within-locality %.3f; between-locality %.3f (record-weighted) / %.3f (each locality once)",
                rec_cor$r_year_tmean, rec_cor$n, rec_cor$r_within_locality,
                rec_cor$r_between_locality_record_weighted, rec_cor$r_between_locality_unweighted))

# Localities that fall outside the WorldClim land mask (all months NA): where
# they are, how many records they carry, distance to the nearest valid cell.
.first_tif <- list.files(file.path(WC_DIR, sprintf("wc2.1_cruts4.09_%s_tmax_%s", RES, DECADES[1])),
                         pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)[1]
.r1 <- terra::crop(terra::rast(.first_tif), bbox)
no_raster <- locs %>%
  filter(!(loc_id %in% annual$loc_id[is.finite(annual$tmean)])) %>%
  left_join(n_by_loc, by = "loc_id") %>%
  left_join(passer90 %>% distinct(loc_id, State, Municipality) %>%
              group_by(loc_id) %>% slice(1) %>% ungroup(), by = "loc_id")
no_raster$km_to_nearest_valid_cell <- vapply(seq_len(nrow(no_raster)), function(i) {
  x <- no_raster$Longitude_decimal_degrees[i]; y <- no_raster$Latitude_decimal_degrees[i]
  w  <- terra::crop(.r1, terra::ext(x - 0.5, x + 0.5, y - 0.5, y + 0.5))
  lp <- terra::as.points(w, na.rm = TRUE)
  if (nrow(lp) == 0) return(NA_real_)
  p1 <- terra::vect(data.frame(x = x, y = y), geom = c("x", "y"), crs = "EPSG:4326")
  min(terra::distance(p1, lp)) / 1000
}, numeric(1))
message(sprintf("Localities outside the WorldClim land mask: %d (%d records); median %.1f km to the nearest valid cell. By state: %s",
                nrow(no_raster), sum(no_raster$n_records), stats::median(no_raster$km_to_nearest_valid_cell, na.rm = TRUE),
                paste(names(table(no_raster$State)), table(no_raster$State), sep = "=", collapse = ", ")))

# --- figure -----------------------------------------------------------------
# Anomalies relative to each locality's 1995-2018 mean, so localities at
# different absolute temperatures share a panel.
anom <- trend_dat %>% group_by(loc_id) %>%
  mutate(tmean_anom = tmean - mean(tmean, na.rm = TRUE),
         tmax_warmq_anom = tmax_warmq - mean(tmax_warmq, na.rm = TRUE)) %>% ungroup()
reg_anom <- anom %>% group_by(year) %>%
  summarise(tmean_anom = mean(tmean_anom, na.rm = TRUE),
            tmax_warmq_anom = mean(tmax_warmq_anom, na.rm = TRUE),
            spei12_dec = mean(spei12_dec, na.rm = TRUE), .groups = "drop")
# Headline pooled trend for the figure: unweighted WITH the year random effect
# (the SE that respects the shared regional weather).
.pp <- function(v) pooled %>% filter(variable == v, weighting == "unweighted + year RE")
# `intercept0` = value of the pooled line at the window midpoint: 0 for the
# anomaly panels (anomalies are centred per locality on the same window), the
# fitted intercept for SPEI, which is plotted on its own scale.
panel <- function(var_anom, var_pooled, ylab, title, intercept0 = 0) {
  p <- .pp(var_pooled)
  ggplot(anom, aes(x = year, y = .data[[var_anom]], group = loc_id)) +
    geom_line(alpha = 0.06, colour = "steelblue") +
    geom_line(data = reg_anom, aes(x = year, y = .data[[var_anom]], group = 1), colour = "black", linewidth = 0.9) +
    geom_point(data = reg_anom, aes(x = year, y = .data[[var_anom]], group = 1), colour = "black", size = 1.4) +
    geom_abline(intercept = intercept0 - p$slope_per_decade * mean(c(TREND_MIN, TREND_MAX)) / 10,
                slope = p$slope_per_decade / 10, colour = "firebrick", linewidth = 0.8, linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    labs(x = "Year", y = ylab, title = title,
         subtitle = sprintf("Pooled trend %+.2f [%+.2f, %+.2f] per decade; %d localities",
                            p$slope_per_decade, p$lower, p$upper, p$n_localities)) +
    theme_classic(base_size = 10)
}
p1 <- panel("tmean_anom", "tmean", "Annual mean temperature anomaly (C)",
            "Annual mean temperature at the sampled localities")
p2 <- panel("tmax_warmq_anom", "tmax_warmq", "Warm-quarter Tmax anomaly (C)",
            "Warmest-quarter maximum temperature")
p3 <- panel("spei12_dec", "spei12_dec", "SPEI-12 (December)",
            "12-month SPEI (calendar-year water balance)",
            intercept0 = .pp("spei12_dec")$intercept) +
  labs(caption = paste0("SPEI source: ", spei_source))
p4 <- ggplot(locality_slopes %>% filter(variable %in% c("tmean", "tmax_warmq")),
             aes(x = slope_per_decade, fill = variable)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity", colour = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c(tmean = "steelblue", tmax_warmq = "firebrick"),
                    labels = c(tmean = "annual mean T", tmax_warmq = "warm-quarter Tmax"), name = NULL) +
  labs(x = "Per-locality OLS slope 1995-2018 (C per decade)", y = "Localities",
       title = "Distribution of locality trends") +
  theme_classic(base_size = 10) +
  theme(legend.position = "inside", legend.position.inside = c(0.25, 0.85))
p_clim <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Climate at the bird sampling localities, 1995-2018",
                  subtitle = sprintf("WorldClim 2.1 historical monthly weather (CRU-TS 4.09 downscaled, %s); thin lines = %d coordinate localities with data (%d more fall outside the land mask)\nblack = across-locality mean; dashed = pooled trend from the locality random-slope model with a year random effect (95%% CI in panel subtitles)",
                                     RES, dplyr::n_distinct(anom$loc_id[is.finite(anom$tmean)]), nrow(no_raster)),
                  theme = theme(plot.subtitle = element_text(size = 8)))
ggsave(fig_path("climate_trends.png"), p_clim, width = 11, height = 8, dpi = 150)
message("Saved: ", fig_path("climate_trends.png"))

# --- save --------------------------------------------------------------------
climate_trends <- list(
  annual            = annual,
  locality_slopes   = locality_slopes,
  slope_summary     = slope_summary,
  pooled            = pooled,
  regional_mean     = regional_mean,
  regional_trend    = regional_trend,
  record_level_cor  = rec_cor,
  n_localities      = nrow(locs),
  n_localities_with_data = nrow(locs) - nrow(no_raster),
  no_raster_localities = no_raster,
  n_records_by_locality = n_by_loc,
  spei_source       = spei_source,
  spei_terraclimate_available = tc_available,
  spei_failed_cells = spei_failed_cells,
  spei_reflected_cells = spei_reflected_cells,
  spei_check        = spei_check,
  settings = list(resolution = RES, year_window = c(YEAR_MIN, YEAR_MAX),
                  trend_window = c(TREND_MIN, TREND_MAX), spei_scale = SPEI_SCALE,
                  warm_quarter = "max over centred 3-month running mean of monthly Tmax, window centre in year",
                  spei_reference_period = paste(YEAR_MIN, YEAR_MAX, sep = "-"),
                  spei_distribution = "generalized logistic (Hosking & Wallis 1997 A.7) from unbiased-PWM L-moments, per calendar month and locality; = SPEI::spei(distribution = 'log-Logistic', fit = 'ub-pwm')",
                  pooled_headline = "weighting == 'unweighted + year RE'"),
  session = session_info
)
saveRDS(climate_trends, out_path("climate_trends.rds"))
message("Saved: ", out_path("climate_trends.rds"))
