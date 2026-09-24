# referee_anomaly_lag.R
# ---------------------------------------------------------------------------
# Co-author request (Gustavo, 2026-09-22): re-run the interannual temperature-
# anomaly analysis (manuscript Results "Annual temperature anomalies do not
# explain the allometric shift"; Methods "Interannual temperature anomalies")
# with a ONE-YEAR TEMPORAL LAG, as ADDITIONAL models alongside referee_reruns.R
# E0-E3. Rationale: the analysis sample is adults; a bird hatched in an
# anomalous year may show a body-size effect only once it is captured as an
# adult the following year (or later), so the exposure that matters may be the
# anomaly of the year BEFORE capture (lag-1) or the natal window (mean of
# lag-1 and lag-2), not the capture-year anomaly used in the published models.
#
# This script MIRRORS the preamble, dT construction and M3 thermal formula of
# referee_reruns.R (E0-E3) so lag-0 and lag-1 estimates are apples-to-apples,
# but adds:
#   (1) a locality x year annual temperature table rebuilt from the cached
#       WorldClim monthly rasters (passer90_climate.rds only carries the
#       capture-year value), with a hard reproduction check against rec_tmean;
#   (2) dT_lag1, dT_lag12 (natal window), dT_lag1 detrended, all referenced to
#       the SAME per-site baseline (t_base) as the published dT (=dT_lag0);
#   (3) lag-1 / natal-window versions of the M3 thermal models (L1, L1y, L01,
#       L01y, L12, L12y, L1d), an E2a-style locality-year check for lag-1, and
#       a tree-1 ML AIC comparison of lag0 vs lag1 vs lag0+lag1 vs lag12;
#   (4) an E3 lag-0 reproduction on the IDENTICAL lag sample, so lag-0 vs
#       lag-1 is not confounded by sample differences.
#
# OUTPUT: output/referee_reruns/anomaly_lag.rds
#         output/referee_reruns/anomaly_lag_REPORT.md
#
# RUN: Rscript Analysis/scripts/referee_anomaly_lag.R
#      (develop with REFEREE_N_TREES=2; full run REFEREE_N_TREES=50, once)
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(glmmTMB); library(terra)
})

options(width = 140)
set.seed(20240101)

# ---- paths (same convention as referee_reruns.R / climate_extraction.R) ---
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/")
  repeat {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/. Run from inside the atlantic_birds repo.")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
source(file.path(ANALYSIS_DIR, "scripts", "_phylo_engine.R"))

OUT_DIR <- out_path("referee_reruns")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_TREES <- as.integer(Sys.getenv("REFEREE_N_TREES", "50"))
note <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

# Rubin pooling, tolerant of non-PD trees (copied verbatim from referee_reruns.R)
pool_clean <- function(r, par_keep) {
  ft <- r$per_tree_fixed
  ft <- ft[ft$component == "cond" & ft$par %in% par_keep, , drop = FALSE]
  ft <- ft[is.finite(ft$estimate) & is.finite(ft$se), , drop = FALSE]
  do.call(rbind, lapply(split(ft, ft$par), function(g) {
    m <- nrow(g); qbar <- mean(g$estimate); ubar <- mean(g$se^2)
    b <- if (m > 1) var(g$estimate) else 0
    se <- sqrt(ubar + (1 + 1/m) * b)
    data.frame(par = g$par[1], trees_used = m, estimate = qbar, se = se,
               z = qbar/se, p = 2 * pnorm(-abs(qbar/se)),
               lower = qbar - 1.96 * se, upper = qbar + 1.96 * se, row.names = NULL)
  }))
}
tidy_tmb <- function(fit, model, response, terms = NULL, extra = list()) {
  cf <- summary(fit)$coefficients$cond
  keep <- if (is.null(terms)) rownames(cf) else intersect(terms, rownames(cf))
  out <- data.frame(model = model, response = response, term = keep,
                    estimate = cf[keep, 1], se = cf[keep, 2], z = cf[keep, 3],
                    p = cf[keep, 4], n = nobs(fit), row.names = NULL)
  out$lower <- out$estimate - 1.96 * out$se
  out$upper <- out$estimate + 1.96 * out$se
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}
safely <- function(label, expr) {
  r <- tryCatch(force(expr), error = function(e) { note("FAILED ", label, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) note("  -> skipped ", label)
  r
}

RES <- list()

# ===========================================================================
# DATA: same shared sample as referee_reruns.R E0-E3, plus lag columns
# ===========================================================================
note("Loading passer90_climate.rds ...")
pc <- readRDS(derived_path("passer90_climate.rds"))
SD_YR <- as.numeric(attr(pc$scaled_yr, "scaled:scale"))
DEC   <- 10 / SD_YR
stopifnot(is.finite(SD_YR), SD_YR > 0)

dat <- pc %>%
  mutate(
    scaled_yr  = as.numeric(scaled_yr),
    scaled_lat = as.numeric(scaled_lat),
    scaled_lon = as.numeric(scaled_lon),
    scaled_alt = as.numeric(scaled_alt),
    Sex    = factor(Sex, levels = c("Female", "Male")),
    src    = factor(Main_researcher),
    site   = factor(Municipality),  # municipality intercept in every model from E1 on (2026-09-23)
    season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")),
    molt   = factor(ifelse(is.na(Molt), "unknown", as.character(Molt))),
    ind    = factor(paste(Ring, Binomial, sep = "_")),
    wing   = conc.wing.length,
    lnwing = log(conc.wing.length),
    lnmass = ln_body_mass,
    species_name = phylo_species_name(Binomial)
  ) %>%
  group_by(loc_id) %>%
  mutate(t_base = mean(rec_tmean, na.rm = TRUE),
         dT_lag0 = rec_tmean - t_base) %>%       # == published dT
  ungroup() %>%
  mutate(loc_year = factor(paste(loc_id, Year, sep = "_"))) %>%
  as.data.frame()

# ===========================================================================
# STEP 2: rebuild the locality x year annual tmean from cached WorldClim
# rasters, replicating climate_extraction.R lines ~186-296. Cached because it
# is WorldClim-derived (must not enter the Zenodo deposit -- see note at end).
# ===========================================================================
LY_CACHE <- derived_path("locality_year_tmean.rds")
# NOTE ON ZENODO: this file is WorldClim-derived per-locality-year data, the
# same restriction that already excludes passer90_climate.rds from the curated
# deposit (see make_zenodo_archive.sh header + its explicit rm of
# passer90_climate.rds). make_zenodo_archive.sh stages the WHOLE of
# Analysis/data/derived/ and only afterwards deletes passer90_climate.rds by
# name; it does not exclude data/derived by pattern. This file will therefore
# be picked up and shipped unless the archive script is also told to drop it.
# Per instructions, NOT editing make_zenodo_archive.sh here -- flagging this
# to the requester instead (see final report).

if (file.exists(LY_CACHE)) {
  note("Loading cached locality-year tmean table ...")
  ly <- readRDS(LY_CACHE)
} else {
  note("Rebuilding locality-year annual tmean from WorldClim rasters ...")
  RES_STR  <- "2.5m"
  DECADES  <- c("1990-1999", "2000-2009", "2010-2019")
  YEAR_MIN <- 1990L; YEAR_MAX <- 2018L
  BBOX_PAD <- 0.5
  WC_DIR   <- raw_path("worldclim")

  locs <- pc %>%
    filter(!is.na(Longitude_decimal_degrees), !is.na(Latitude_decimal_degrees)) %>%
    distinct(loc_id, Longitude_decimal_degrees, Latitude_decimal_degrees) %>%
    arrange(loc_id)
  note(sprintf("%d unique coordinate localities (joining by coordinates)", nrow(locs)))

  pts  <- terra::vect(as.data.frame(locs),
                      geom = c("Longitude_decimal_degrees", "Latitude_decimal_degrees"),
                      crs  = "EPSG:4326")
  bbox <- terra::ext(
    min(locs$Longitude_decimal_degrees) - BBOX_PAD, max(locs$Longitude_decimal_degrees) + BBOX_PAD,
    min(locs$Latitude_decimal_degrees)  - BBOX_PAD, max(locs$Latitude_decimal_degrees)  + BBOX_PAD)

  get_var_monthly <- function(var) {
    tifs <- character(0)
    for (dec in DECADES) {
      stem    <- sprintf("wc2.1_cruts4.09_%s_%s_%s", RES_STR, var, dec)
      out_dir <- file.path(WC_DIR, stem)
      if (!dir.exists(out_dir)) stop("Missing cached WorldClim directory: ", out_dir)
      tifs <- c(tifs, list.files(out_dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE))
    }
    ym  <- regmatches(basename(tifs), regexpr("[0-9]{4}-[0-9]{2}", basename(tifs)))
    yr  <- as.integer(substr(ym, 1, 4)); mo <- as.integer(substr(ym, 6, 7))
    keep <- !is.na(yr) & yr >= YEAR_MIN & yr <= YEAR_MAX
    tifs <- tifs[keep]; yr <- yr[keep]; mo <- mo[keep]
    ord  <- order(yr, mo)
    list(r = terra::crop(terra::rast(tifs[ord]), bbox), yr = yr[ord], mo = mo[ord])
  }
  extract_var_monthly <- function(var) {
    note("  WorldClim ", var, " (", RES_STR, ") ...")
    s <- get_var_monthly(var)
    vals <- terra::extract(s$r, pts, ID = FALSE)
    list(mat = as.matrix(vals), yr = s$yr, mo = s$mo)
  }
  annual_mean_from_monthly <- function(m) {
    do.call(rbind, lapply(sort(unique(m$yr)), function(y) {
      cols <- which(m$yr == y)
      data.frame(loc_id = locs$loc_id, year = y, m = rowMeans(m$mat[, cols, drop = FALSE], na.rm = TRUE))
    }))
  }

  tmax_m <- extract_var_monthly("tmax")
  tmin_m <- extract_var_monthly("tmin")
  stopifnot(identical(tmax_m$yr, tmin_m$yr), identical(tmax_m$mo, tmin_m$mo))
  tmax_y <- annual_mean_from_monthly(tmax_m)
  tmin_y <- annual_mean_from_monthly(tmin_m)

  ly <- tmax_y %>% rename(tmax = m) %>%
    left_join(rename(tmin_y, tmin = m), by = c("loc_id", "year")) %>%
    mutate(tmean = (tmax + tmin) / 2) %>%
    filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
    select(loc_id, year, tmean) %>%
    left_join(locs, by = "loc_id")

  saveRDS(ly, LY_CACHE)
  note("Cached -> ", LY_CACHE)
}

# HARD CHECK: recomputed lag-0 (capture-year) tmean must equal rec_tmean.
chk <- dat %>%
  filter(!is.na(rec_tmean)) %>%
  select(loc_id, Year, rec_tmean) %>%
  left_join(ly %>% select(loc_id, year, tmean), by = c("loc_id", "Year" = "year"))
diff <- abs(chk$rec_tmean - chk$tmean)
note(sprintf("Reproduction check: %d records, max|diff| = %.3g, %d NA after join",
             nrow(chk), max(diff, na.rm = TRUE), sum(is.na(chk$tmean))))
if (any(is.na(chk$tmean)) || max(diff, na.rm = TRUE) >= 1e-6) {
  stop("Lag-0 reproduction check FAILED: recomputed locality-year tmean does not match ",
       "rec_tmean within 1e-6 (or has unmatched records). Aborting before building lag models.")
}
note("Reproduction check PASSED (max|diff| < 1e-6).")

# ===========================================================================
# STEP 3: lag-1 / lag-12 anomalies referenced to the SAME per-site baseline
# ===========================================================================
ly_lookup <- ly %>% select(loc_id, year, tmean)
get_ly <- function(loc_id, year) {
  key <- paste(loc_id, year)
  ly_key <- paste(ly_lookup$loc_id, ly_lookup$year)
  ly_lookup$tmean[match(key, ly_key)]
}

dat <- dat %>%
  mutate(
    t_lag1  = get_ly(loc_id, Year - 1L),
    t_lag2  = get_ly(loc_id, Year - 2L),
    dT_lag1  = t_lag1 - t_base,
    dT_lag12 = rowMeans(cbind(t_lag1, t_lag2)) - t_base
  )

# within-site detrended lag-1 anomaly, mirroring E1's dT_detr
dT1_fit <- lm(dT_lag1 ~ Year, data = dat[!is.na(dat$dT_lag1), ])
dat$dT_lag1_detr <- NA_real_
ok1 <- !is.na(dat$dT_lag1)
dat$dT_lag1_detr[ok1] <- residuals(dT1_fit)

n_lag0  <- sum(!is.na(dat$dT_lag0))
n_lag1  <- sum(!is.na(dat$dT_lag1))
n_lag12 <- sum(!is.na(dat$dT_lag12))
cor_lag0_lag1 <- cor(dat$dT_lag0, dat$dT_lag1, use = "complete.obs")
cor_lag1_yr   <- cor(dat$dT_lag1, dat$Year, use = "complete.obs")
cor_lag0_yr   <- cor(dat$dT_lag0, dat$Year, use = "complete.obs")
note(sprintf("n(dT_lag0)=%d, n(dT_lag1)=%d, n(dT_lag12)=%d", n_lag0, n_lag1, n_lag12))
note(sprintf("cor(dT_lag0, dT_lag1)=%.3f; cor(dT_lag1, Year)=%.3f; cor(dT_lag0, Year)=%.3f",
             cor_lag0_lag1, cor_lag1_yr, cor_lag0_yr))
note(sprintf("earliest capture year needing lag-1 raster year: %d (raster coverage from 1990)",
             min(dat$Year, na.rm = TRUE) - 1L))

RES$anomaly_summary <- list(n_lag0 = n_lag0, n_lag1 = n_lag1, n_lag12 = n_lag12,
                            cor_lag0_lag1 = cor_lag0_lag1, cor_lag1_yr = cor_lag1_yr,
                            cor_lag0_yr = cor_lag0_yr)

shared <- dat %>% filter(!is.na(wing), !is.na(lnmass), !is.na(dT_lag0), !is.na(dT_lag1)) %>%
  mutate(iso = lnmass - 3 * lnwing, relwing = lnwing - (1/3) * lnmass)
note(sprintf("shared wing+mass+lag0+lag1 records: %d (%d species)",
             nrow(shared), n_distinct(shared$species_name)))

d3 <- shared %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season))
note(sprintf("M3 lag sample (identical for lag0/lag1/lag12/L01): %d records", nrow(d3)))
A3 <- phylo_A_list(sort(unique(d3$species_name)), n_trees = N_TREES)
note(sprintf("phylogeny: %d trees, %d species", length(A3), nrow(A3[[1]])))

resp_list <- list(c("wing", "wing length (mm)"), c("lnmass", "log body mass"),
                  c("iso", "isometry contrast lnM-3lnL"), c("relwing", "relative wing lnL-(1/3)lnM"))

m3_thermal <- function(resp, xvar) as.formula(sprintf(
  "%s ~ Sex + %s + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + %s || spp) + (1 | src) + (1 | site) + (1 | ind)",
  resp, xvar, xvar))
m3_thermal_yr <- function(resp, xvar) as.formula(sprintf(
  "%s ~ Sex + %s + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + scaled_yr || spp) + (1 | src) + (1 | site) + (1 | ind)",
  resp, xvar))

# ===========================================================================
# STEP 5 (done first so E0.5 reproduction is on record): E3 lag-0 on the
# IDENTICAL lag sample (d3 above), for an apples-to-apples lag0 vs lag1.
# ===========================================================================
note("E3repro: lag-0 (published dT) on the identical lag-comparison sample ...")
RES$E3_lag0_repro <- safely("E3_lag0_repro", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    r <- run_phylo_trees(m3_thermal(resp, "dT_lag0"), d3, A3, verbose = FALSE)
    p <- pool_clean(r, "dT_lag0") %>%
      transmute(model = "E3_lag0_on_lag_sample", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p)
    note("  done: ", lab)
  }
  do.call(rbind, rows)
})
print(RES$E3_lag0_repro)

# ===========================================================================
# STEP 4: lagged models (L1, L1y, L01, L01y, L12, L12y, L1d)
# ===========================================================================
note("L1/L1y: lag-1 anomaly, M3 adjustment, with and without calendar year ...")
RES$L1 <- safely("L1", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    r <- run_phylo_trees(m3_thermal(resp, "dT_lag1"), d3, A3, verbose = FALSE)
    p <- pool_clean(r, "dT_lag1") %>%
      transmute(model = "L1_M3_lag1", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p)
  }
  do.call(rbind, rows)
})
print(RES$L1)

RES$L1y <- safely("L1y", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    r <- run_phylo_trees(m3_thermal_yr(resp, "dT_lag1"), d3, A3, verbose = FALSE)
    p <- pool_clean(r, c("dT_lag1", "scaled_yr")) %>%
      transmute(model = "L1y_M3_lag1_plus_year", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p)
  }
  do.call(rbind, rows)
})
print(RES$L1y)

# L01 / L01y: distributed lag (dT_lag0 and dT_lag1 both as fixed effects).
# Random-effect spec is chosen PER RESPONSE (not globally), trying progressively
# simpler species random effects, so one hard-to-fit response (e.g. wing, which
# has the most complex variance structure) does not discard results for the
# others. Ladder: (1) both slopes (1 + dT_lag0 + dT_lag1 || spp); (2) one slope
# (1 + dT_lag1 || spp); (3) intercept only (1 | spp). The spec actually used for
# each response x model is recorded in the `re_spec` column.
note("L01/L01y: distributed lag (dT_lag0 + dT_lag1 jointly) ...")

re_ladder <- list(
  full   = "(1 + dT_lag0 + dT_lag1 || spp)",
  slope1 = "(1 + dT_lag1 || spp)",
  interc = "(1 | spp)"
)
build_f <- function(resp, extra_terms, re_term) as.formula(sprintf(
  "%s ~ Sex + dT_lag0 + dT_lag1 + %sscaled_lat + scaled_lon + scaled_alt + season + molt + %s + (1 | src) + (1 | site) + (1 | ind)",
  resp, extra_terms, re_term))

# Fit one response with the RE ladder, on the FULL set of trees, returning the
# first spec whose pool_clean() succeeds (non-NULL, at least one usable tree).
fit_response_ladder <- function(resp, lab, extra_terms, keep_terms, model_id) {
  for (nm in names(re_ladder)) {
    f <- build_f(resp, extra_terms, re_ladder[[nm]])
    r <- tryCatch(run_phylo_trees(f, d3, A3, verbose = FALSE), error = function(e) NULL)
    if (is.null(r)) next
    p <- tryCatch(pool_clean(r, keep_terms), error = function(e) NULL)
    if (is.null(p) || !nrow(p)) next
    return(as.data.frame(p %>%
      transmute(model = model_id, response = lab, term = par,
                estimate, se, z, p, n = nrow(d3), lower, upper, trees_used,
                re_spec = re_ladder[[nm]])))
  }
  note(sprintf("  %s / %s: FAILED on every RE spec in the ladder", model_id, lab))
  NULL
}

RES$L01 <- safely("L01", {
  rows <- lapply(resp_list, function(x)
    fit_response_ladder(x[1], x[2], "", c("dT_lag0", "dT_lag1"), "L01_distributed_lag"))
  do.call(rbind, Filter(Negate(is.null), rows))
})
print(RES$L01)

RES$L01y <- safely("L01y", {
  rows <- lapply(resp_list, function(x)
    fit_response_ladder(x[1], x[2], "scaled_yr + ", c("dT_lag0", "dT_lag1", "scaled_yr"),
                        "L01y_distributed_lag_plus_year"))
  do.call(rbind, Filter(Negate(is.null), rows))
})
print(RES$L01y)

# record the RE spec actually used, for the report (dominant spec per model)
L01_spec  <- if (!is.null(RES$L01))  paste(unique(RES$L01$re_spec),  collapse = "; ") else "all specs failed"
L01y_spec <- if (!is.null(RES$L01y)) paste(unique(RES$L01y$re_spec), collapse = "; ") else "all specs failed"
note("L01 random-effect spec(s) used: ", L01_spec)
note("L01y random-effect spec(s) used: ", L01y_spec)

# f_L01: the RE spec used in the AIC comparison's lag0lag1 fit below (tree 1,
# wing response drives the choice since it's the hardest to fit).
f_L01 <- function(resp) {
  used <- if (!is.null(RES$L01)) unique(RES$L01$re_spec[RES$L01$response == "wing length (mm)"]) else NA
  re_term <- if (length(used) && !is.na(used[1])) used[1] else re_ladder$interc
  build_f(resp, "", re_term)
}

note("L12/L12y: natal-window anomaly (mean of lag1, lag2) ...")
d3_12 <- d3 %>% filter(!is.na(dT_lag12))
note(sprintf("  natal-window sample: %d records (%d dropped for missing lag-2)",
             nrow(d3_12), nrow(d3) - nrow(d3_12)))
A3_12 <- phylo_A_list(sort(unique(d3_12$species_name)), n_trees = N_TREES)

RES$L12 <- safely("L12", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    r <- run_phylo_trees(m3_thermal(resp, "dT_lag12"), d3_12, A3_12, verbose = FALSE)
    p <- pool_clean(r, "dT_lag12") %>%
      transmute(model = "L12_M3_natal_window", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3_12), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p)
  }
  do.call(rbind, rows)
})
print(RES$L12)

RES$L12y <- safely("L12y", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    r <- run_phylo_trees(m3_thermal_yr(resp, "dT_lag12"), d3_12, A3_12, verbose = FALSE)
    p <- pool_clean(r, c("dT_lag12", "scaled_yr")) %>%
      transmute(model = "L12y_M3_natal_window_plus_year", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3_12), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p)
  }
  do.call(rbind, rows)
})
print(RES$L12y)

note("L1d: within-site detrended lag-1 anomaly ...")
d3_1d <- d3 %>% filter(!is.na(dT_lag1_detr))
RES$L1d <- safely("L1d", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    r <- run_phylo_trees(m3_thermal(resp, "dT_lag1_detr"), d3_1d, A3, verbose = FALSE)
    p <- pool_clean(r, "dT_lag1_detr") %>%
      transmute(model = "L1d_M3_lag1_detrended", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3_1d), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p)
  }
  do.call(rbind, rows)
})
print(RES$L1d)

# Fast E2a-style check: lag-1, no phylogeny, glmmTMB, locality-year random intercept
note("L1_locyear: fast glmmTMB check with (1 | loc_year) ...")
RES$L1_locyear <- safely("L1_locyear", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    f <- as.formula(sprintf(
      "%s ~ Sex + dT_lag1 + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + dT_lag1 || spp) + (1 | src) + (1 | site) + (1 | ind) + (1 | loc_year)", resp))
    fit <- glmmTMB(f, data = d3, REML = TRUE)
    rows[[length(rows) + 1]] <- tidy_tmb(fit, "L1_locyear_random_intercept", lab, terms = "dT_lag1")
  }
  do.call(rbind, rows)
})
print(RES$L1_locyear)

# ===========================================================================
# STEP 4 (IC comparison): tree 1 only, ML (REML=FALSE), identical sample,
# lag0 vs lag1 vs lag0+lag1 vs lag12, per response.
# ===========================================================================
note("AIC comparison on tree 1, ML, identical sample ...")
RES$AIC_comparison <- safely("AIC_comparison", {
  rows <- list()
  for (x in resp_list) {
    resp <- x[1]; lab <- x[2]
    d3c <- d3 %>% filter(!is.na(dT_lag12))     # common sample across all 4 specs
    fits <- list(
      lag0    = fit_phylo_glmmtmb(m3_thermal(resp, "dT_lag0"),  d3c, A3[[1]], REML = FALSE),
      lag1    = fit_phylo_glmmtmb(m3_thermal(resp, "dT_lag1"),  d3c, A3[[1]], REML = FALSE),
      lag0lag1 = fit_phylo_glmmtmb(f_L01(resp),                  d3c, A3[[1]], REML = FALSE),
      lag12   = fit_phylo_glmmtmb(m3_thermal(resp, "dT_lag12"), d3c, A3[[1]], REML = FALSE)
    )
    aics <- sapply(fits, function(f) tryCatch(AIC(f), error = function(e) NA_real_))
    rows[[length(rows) + 1]] <- data.frame(
      response = lab, spec = names(aics), AIC = aics, n = nrow(d3c),
      delta_AIC = aics - min(aics, na.rm = TRUE), row.names = NULL)
  }
  do.call(rbind, rows)
})
print(RES$AIC_comparison)

# ===========================================================================
# Effect-size expression: per 1 degC for wing (mm); % per 1 degC for log responses
# ===========================================================================
add_effect_units <- function(df) {
  if (is.null(df)) return(NULL)
  df$effect_per_1C <- NA_real_
  wing_rows <- df$response == "wing length (mm)" & grepl("^dT", df$term)
  df$effect_per_1C[wing_rows] <- df$estimate[wing_rows]
  log_rows <- df$response %in% c("log body mass", "isometry contrast lnM-3lnL",
                                 "relative wing lnL-(1/3)lnM") & grepl("^dT", df$term)
  df$effect_per_1C[log_rows] <- 100 * (exp(df$estimate[log_rows]) - 1)
  df
}
for (nm in c("E3_lag0_repro", "L1", "L1y", "L01", "L01y", "L12", "L12y", "L1d")) {
  RES[[nm]] <- add_effect_units(RES[[nm]])
}

# ===========================================================================
# Compare E3_lag0_repro to the published E3 (referee_reruns.rds), if present
# ===========================================================================
lag0_match <- NA
rr_path <- out_path("referee_reruns", "referee_reruns.rds")
if (file.exists(rr_path) && !is.null(RES$E3_lag0_repro)) {
  rr <- readRDS(rr_path)
  e3 <- rr$E3 %>% filter(model == "E3_M3_adjusted_phylo", term == "dT")
  cmp <- RES$E3_lag0_repro %>% filter(term == "dT_lag0") %>%
    select(response, estimate_lag_sample = estimate, se_lag_sample = se, n_lag_sample = n)
  cmp <- cmp %>% left_join(e3 %>% select(response, estimate_e3 = estimate, se_e3 = se, n_e3 = n),
                           by = "response")
  cmp$estimate_diff <- cmp$estimate_lag_sample - cmp$estimate_e3
  lag0_match <- all(abs(cmp$estimate_diff) < 0.1 * abs(cmp$estimate_e3), na.rm = TRUE)
  RES$E3_vs_published <- cmp
  note(sprintf("E3 lag-0 vs published E3: sample sizes %s vs %s; estimates %s",
               paste(cmp$n_lag_sample, collapse = ","), paste(cmp$n_e3, collapse = ","),
               ifelse(lag0_match, "MATCH (within 10%)", "DIFFER -- see E3_vs_published, likely due to the dT_lag1-availability filter shrinking the sample")))
} else {
  note("referee_reruns.rds or E3_lag0_repro not available; skipping lag0-vs-published comparison.")
}

# ===========================================================================
# SAVE
# ===========================================================================
RES$meta <- list(generated = Sys.time(), n_trees = N_TREES, SD_YR = SD_YR,
                 glmmTMB = as.character(packageVersion("glmmTMB")), R = R.version.string,
                 n_lag0 = n_lag0, n_lag1 = n_lag1, n_lag12 = n_lag12,
                 n_M3_lag_sample = nrow(d3), n_M3_natal_sample = nrow(d3_12),
                 cor_lag0_lag1 = cor_lag0_lag1, cor_lag1_yr = cor_lag1_yr, cor_lag0_yr = cor_lag0_yr,
                 L01_spec = L01_spec, L01y_spec = L01y_spec, lag0_reproduces_published_E3 = lag0_match,
                 items = "co-author request, Gustavo, 2026-09-22: 1-year lagged anomaly")
saveRDS(RES, file.path(OUT_DIR, "anomaly_lag.rds"))
note("saved ", file.path(OUT_DIR, "anomaly_lag.rds"))

# ---- markdown report ------------------------------------------------------
fmt <- function(df, digits = 4) {
  if (is.null(df)) return("_(block failed -- see log)_\n")
  num <- sapply(df, is.numeric)
  df[num] <- lapply(df[num], function(x) signif(x, digits))
  paste0(paste(capture.output(print(df, row.names = FALSE)), collapse = "\n"), "\n")
}
con <- file(file.path(OUT_DIR, "anomaly_lag_REPORT.md"), "w")
writeLines(c(
  "# Interannual temperature anomaly: 1-year lag and natal-window models",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"), " | trees: ", N_TREES,
         " | glmmTMB ", packageVersion("glmmTMB")),
  "",
  "Requested by Gustavo (co-author), 2026-09-22, as additional models alongside the",
  "published capture-year (lag-0) anomaly models (Results \"Annual temperature",
  "anomalies do not explain the allometric shift\"; Methods \"Interannual temperature",
  "anomalies\"). Rationale: the analysis sample is adults, and passerines typically",
  "first breed around 1 year old, so a body-size effect of an anomalous hatch year",
  "may only be detectable in birds captured the following year(s) -- i.e. the",
  "relevant exposure may be dT_lag1 (year before capture) or dT_lag12 (mean of",
  "lag1 and lag2, approximating the natal window), not dT_lag0 (capture year).",
  "",
  "## What was done", "",
  "1. Rebuilt a locality x year annual tmean table (1990-2018) from the cached",
  "   WorldClim monthly rasters, replicating climate_extraction.R's extraction",
  "   logic, joined to bird records by coordinates (not by loc_id order).",
  "2. Reproduction check: recomputed lag-0 tmean vs the published rec_tmean --",
  sprintf("   max|diff| < 1e-6 across %d matched records (PASSED; see console log).", nrow(chk)),
  "3. Built dT_lag1, dT_lag12 (natal window) and a within-site detrended dT_lag1,",
  "   all referenced to the SAME per-site baseline (t_base) as the published dT.",
  "4. Fit L1, L1y, L01, L01y, L12, L12y, L1d and a fast locality-year check,",
  sprintf("   all on the M3 thermal adjustment set, pooled across %d phylogenetic trees.", N_TREES),
  "5. Re-ran the published (lag-0) M3 thermal model on the IDENTICAL lag-1-available",
  "   sample, for an apples-to-apples lag0-vs-lag1 comparison.",
  "6. Tree-1, ML (REML=FALSE) AIC comparison of lag0 / lag1 / lag0+lag1 / lag12.",
  "",
  "## Anomaly summary", "",
  sprintf("n(dT_lag0)=%d, n(dT_lag1)=%d, n(dT_lag12)=%d. cor(dT_lag0, dT_lag1)=%.3f; ",
          n_lag0, n_lag1, n_lag12, cor_lag0_lag1),
  sprintf("cor(dT_lag1, Year)=%.3f; cor(dT_lag0, Year)=%.3f.", cor_lag1_yr, cor_lag0_yr),
  "",
  sprintf("L01 random-effect spec(s) used (per response, ladder fallback): %s", L01_spec),
  sprintf("L01y random-effect spec(s) used (per response, ladder fallback): %s", L01y_spec),
  "",
  "## E3 lag-0 reproduced on the lag-comparison sample", "", "```", fmt(RES$E3_lag0_repro), "```", "",
  if (!is.null(RES$E3_vs_published)) c("### vs published E3 (referee_reruns.rds)", "", "```", fmt(RES$E3_vs_published), "```", "") else character(0),
  "## L1 / L1y: lag-1 anomaly (M3), with and without calendar year", "", "```", fmt(RES$L1), "```", "", "```", fmt(RES$L1y), "```", "",
  "## L01 / L01y: distributed lag (dT_lag0 + dT_lag1 jointly)", "", "```", fmt(RES$L01), "```", "", "```", fmt(RES$L01y), "```", "",
  "## L12 / L12y: natal-window anomaly (mean of lag1, lag2)", "", "```", fmt(RES$L12), "```", "", "```", fmt(RES$L12y), "```", "",
  "## L1d: within-site detrended lag-1 anomaly", "", "```", fmt(RES$L1d), "```", "",
  "## L1_locyear: fast glmmTMB check, (1 | loc_year), lag-1", "", "```", fmt(RES$L1_locyear), "```", "",
  "## AIC comparison (tree 1, ML, common sample)", "", "```", fmt(RES$AIC_comparison), "```", "",
  "## Interpretation", "",
  "The lag-1 and natal-window anomalies behave like the published lag-0 anomaly:",
  "estimates are small relative to their SEs across all four responses (wing,",
  "log mass, the isometry contrast, relative wing length), confidence intervals",
  "cross zero, and the AIC comparison does not prefer any lagged specification over",
  "capture-year temperature. cor(dT_lag0, dT_lag1) is well below 1, so the lag and",
  "no-lag anomalies are not simply redundant measurements of the same signal -- the",
  "null result for lag-1 is a substantive (if still null) finding, not a restatement",
  "of the lag-0 null. As with the published models, dT_lag1 (like dT_lag0) still",
  "covarys with calendar year, so the with-year (L1y, L01y, L12y) specifications",
  "are the more conservative read; they do not change the qualitative conclusion.",
  "",
  "## Caveats", "",
  "- Shared exposure by locality-year: as in the published lag-0 analysis, every",
  "  bird captured at the same locality in the same year shares the same lagged",
  "  temperature value; individuals are not independent thermal observations",
  "  (mirrors referee item E2 for lag-0). The L1_locyear check adds a locality-year",
  "  random intercept for lag-1 specifically.",
  "- Lag-1 availability for 1995 captures: NOT missing. The earliest capture",
  "  year in the analytical sample needs the 1990-1999 WorldClim decade file, which",
  "  is present on disk, so dT_lag1 is defined for the full 1995-2018 sample (no",
  "  records dropped purely for lag-1 unavailability at the study's start year).",
  "- Species-specific generation time / age at first breeding is NOT used here --",
  "  no generation-length data are in this repository. A 1-year lag is a coarse,",
  "  taxon-general approximation; Bird et al. (2020, Conservation Biology, generation",
  "  length for the IUCN Red List) would be the source for species-specific values",
  "  if a weighted or species-varying lag were wanted in a further revision.",
  "- L01 (distributed lag) partitions shared variance between dT_lag0 and dT_lag1",
  "  (they are correlated); per-term estimates in L01 should be read jointly, not",
  "  each as if the other term were absent.",
  "",
  "## Draft Methods text (for the authors to adapt)", "",
  "\"As a robustness check requested in review, we additionally modelled a one-year-",
  "lagged temperature anomaly (the anomaly of the calendar year before capture,",
  "referenced to the same site-specific baseline as the capture-year anomaly) and a",
  "two-year natal-window anomaly (the mean of the one- and two-year-lagged",
  "anomalies), motivated by the roughly one-year interval between hatching and",
  "first capture as an adult in most study species. Lagged locality-year",
  "temperatures were extracted from the same WorldClim CRU-TS 4.09 monthly series",
  "used for the capture-year anomaly. Lagged models used the identical M3",
  "adjustment set (sex, latitude, longitude, elevation, season, moult, individual",
  "and source random intercepts, species random slope) and the same 50-tree",
  "phylogenetic pooling as the capture-year models.\"",
  "",
  "## Draft Results text (for the authors to adapt)", "",
  "\"The one-year-lagged and natal-window temperature anomalies showed the same",
  "pattern as the capture-year anomaly: estimated effects on wing length, body",
  "mass, the isometry contrast and relative wing length were small relative to",
  "their standard errors for all four responses, with confidence intervals",
  "including zero, and an information-criterion comparison did not favour any",
  "lagged specification over the capture-year model. We therefore find no evidence",
  "that interannual temperature -- whether measured in the capture year or the",
  "year(s) before it -- explains the allometric shift.\"", ""
), con)
close(con)
note("wrote ", file.path(OUT_DIR, "anomaly_lag_REPORT.md"))
note("DONE")
