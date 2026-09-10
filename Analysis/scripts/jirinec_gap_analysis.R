# jirinec_gap_analysis.R
# ===========================================================================
# EXPLORATORY SIDE-REPORT — NOT part of the manuscript pipeline.
#
# WHAT: tests the seven analyses that Jirinec et al. (2021, Sci Adv 7:eabk1743)
#       performed on the BDFFP Amazonian data and that the Atlantic Forest
#       manuscript does NOT currently do. Nothing here overwrites any file the
#       manuscript reads: every output goes to output/jirinec_gaps/ and
#       figures/jirinec_gaps/.
#
# GAPS TESTED
#   G1  Lagged climate covariates in the MEAN models (their central mechanistic
#       result: the year trend in mass vanishes conditional on lagged
#       temperature; precipitation models fit better than temperature models).
#       We currently use scaled_tmean only in drmSEM and the dispersion models.
#   G2  Species-level mass slopes (they report per-species tallies; we only have
#       a species forest plot for wing).
#   G3  mass:wing ratio as a response (their headline metric, 69% of species).
#   G4  Phylogenetic signal in the SLOPES (Moran's I / Pagel's lambda on the
#       per-species year coefficients). We only have phylo signal on intercepts.
#   G5  Vertical foraging stratum as a moderator (the ONE trait that gave them a
#       signal). We only test diet-invertebrate proportion.
#   G6  Design-matched restricted tier (their single-site / single-protocol
#       design approximated by long-span contributor blocks).
#   G7  Their symmetric coverage filter (>= 5 records BOTH before and after the
#       median year), which our >= 30 records / >= 5-year span filter does not
#       guarantee.
#
# KEY DESIGN NOTE FOR G1 --------------------------------------------------
#   Jirinec had ONE site, so their climate coefficients are, by construction,
#   pure WITHIN-site temporal anomalies. Our records span 357 localities across
#   ~15 degrees of latitude, so raw temperature is dominated by the spatial
#   gradient (lowland hot vs montane cool) -- a Bergmann-in-space effect, not a
#   temporal response. To make our estimate comparable to theirs we decompose
#   every climate variable into
#       between = locality mean over 1990-2018   (spatial gradient)
#       within  = value - locality mean          (the temporal anomaly)
#   and read the temporal climate response off the WITHIN terms.
#   Their seasonal lags (dry season of capture / previous wet / previous dry)
#   are not transferable: our captures run year-round across a biome with no
#   common dry season. We use ANNUAL lags 0 / 1 / 2 instead and say so.
#
# RUN
#   Rscript Analysis/scripts/jirinec_gap_analysis.R              # lme4 tier, fast
#   Rscript Analysis/scripts/jirinec_gap_analysis.R --trees 50   # + phylo tier
# ===========================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(lme4); library(glmmTMB); library(ape)
})

.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/.")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
gap_path     <- function(...) file.path(ANALYSIS_DIR, "output", "jirinec_gaps", ...)
gfig_path    <- function(...) file.path(ANALYSIS_DIR, "figures", "jirinec_gaps", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
dir.create(gap_path(),  showWarnings = FALSE, recursive = TRUE)
dir.create(gfig_path(), showWarnings = FALSE, recursive = TRUE)

args    <- commandArgs(trailingOnly = TRUE)
get_flag <- function(f) {
  i <- grep(paste0("^", f, "(=|$)"), args)[1]
  if (is.na(i)) return(NULL)
  v <- sub(paste0("^", f, "=?"), "", args[i]); if (!nzchar(v)) v <- args[i + 1]; v
}
N_TREES <- if (!is.null(get_flag("--trees"))) as.integer(get_flag("--trees")) else 0L
DO_PHYLO <- N_TREES > 0L
if (DO_PHYLO) source(script_path("_phylo_engine.R"))

res <- list(generated = Sys.time(), n_trees = N_TREES,
            session = c(R = R.version.string,
                        sapply(c("lme4", "glmmTMB", "ape", "phytools"),
                               function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA))))

# ---------------------------------------------------------------------------
# 1. Data: the manuscript's analytical sample + record-level climate lags
# ---------------------------------------------------------------------------
d0 <- readRDS(derived_path("passer90_climate.rds"))
stopifnot(all(d0$Status == "live"), all(d0$known_sex))

SD_YR  <- 5.019925                       # attr(passer90$scaled_yr, "scaled:scale")
CTR_YR <- 2009.487
DEC    <- 10 / SD_YR                     # SD-years per decade

d <- d0 %>%
  mutate(scaled_yr  = as.numeric(scaled_yr),
         scaled_lat = as.numeric(scaled_lat),
         scaled_lon = as.numeric(scaled_lon),
         scaled_alt = as.numeric(scaled_alt),
         Sex        = factor(Sex, levels = c("Female", "Male")),
         src        = factor(Main_researcher),
         site       = factor(Municipality),
         ring       = factor(ifelse(is.na(Ring), paste0("unringed_", ID_ABT),
                                    paste(Ring, Binomial, sep = "|"))),
         season     = factor(season, levels = c("DJF", "MAM", "JJA", "SON")),
         spp        = factor(Binomial),
         wing       = conc.wing.length)

# --- climate: locality x year series, annual lags 0/1/2, within/between split --
ann <- readRDS(out_path("climate_trends.rds"))$annual %>%
  select(loc_id, year, tmean, prec_annual, spei12_mean) %>%
  filter(!is.na(tmean))

# between-locality means over the full 1990-2018 window
loc_mean <- ann %>% group_by(loc_id) %>%
  summarise(tmean_loc = mean(tmean), prec_loc = mean(prec_annual),
            spei_loc  = mean(spei12_mean), .groups = "drop")

# within-locality anomalies, then lag them by calendar year
anom <- ann %>% left_join(loc_mean, by = "loc_id") %>%
  transmute(loc_id, year,
            t_an = tmean - tmean_loc,
            p_an = prec_annual - prec_loc,
            s_an = spei12_mean - spei_loc)

lag_tbl <- function(k) anom %>% mutate(year = year + k) %>%
  rename(!!paste0("t_an_l", k) := t_an, !!paste0("p_an_l", k) := p_an,
         !!paste0("s_an_l", k) := s_an)

clim <- lag_tbl(0) %>%
  left_join(lag_tbl(1), by = c("loc_id", "year")) %>%
  left_join(lag_tbl(2), by = c("loc_id", "year")) %>%
  left_join(loc_mean,   by = "loc_id")

d <- d %>% left_join(clim, by = c("loc_id" = "loc_id", "Year" = "year"))

# standardise anomalies and between-locality means (SD units, so coefficients compare)
an_cols  <- c("t_an_l0","t_an_l1","t_an_l2","p_an_l0","p_an_l1","p_an_l2","s_an_l0","s_an_l1","s_an_l2")
loc_cols <- c("tmean_loc","prec_loc","spei_loc")
an_sd <- sapply(d[an_cols],  sd, na.rm = TRUE)
lc_c  <- sapply(d[loc_cols], mean, na.rm = TRUE); lc_s <- sapply(d[loc_cols], sd, na.rm = TRUE)
for (v in an_cols)  d[[paste0("z_", v)]] <- d[[v]] / an_sd[[v]]
for (v in loc_cols) d[[paste0("z_", v)]] <- (d[[v]] - lc_c[[v]]) / lc_s[[v]]

res$climate_scaling <- list(anomaly_sd = an_sd, loc_center = lc_c, loc_scale = lc_s,
  note = paste("Anomalies are within-locality deviations from the 1990-2018 locality mean,",
               "divided by their SD across records. z_tmean_loc etc. are the between-locality",
               "(spatial) components. Lags are ANNUAL (l0 = capture year, l1 = previous year,",
               "l2 = two years before), NOT the seasonal lags of Jirinec et al."))

clim_ok <- complete.cases(d[paste0("z_", an_cols)]) & complete.cases(d[paste0("z_", loc_cols)])
res$climate_coverage <- data.frame(
  n_records_total = nrow(d), n_records_climate = sum(clim_ok),
  pct = round(100 * sum(clim_ok) / nrow(d), 1),
  n_mass_climate = sum(clim_ok & !is.na(d$ln_body_mass)),
  n_wing_climate = sum(clim_ok & !is.na(d$wing)),
  n_spp_climate  = n_distinct(d$spp[clim_ok]),
  n_loc_climate  = n_distinct(d$loc_id[clim_ok]))
message("[data] climate-complete records: ", sum(clim_ok), " / ", nrow(d),
        " (mass ", res$climate_coverage$n_mass_climate, ", wing ", res$climate_coverage$n_wing_climate, ")")

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------
CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_lmer <- function(formula, data, REML = TRUE) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    tryCatch(lmer(formula, data = data, REML = REML, control = CTRL),
             error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL }),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  attr(fit, "msgs") <- paste(unique(msgs), collapse = " | "); fit
}

# tidy fixed effects of interest from an lmer fit
tidy_lmer <- function(fit, terms_keep, trait, model, note = "") {
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients
  do.call(rbind, lapply(intersect(terms_keep, rownames(cf)), function(tm) {
    est <- cf[tm, "Estimate"]; se <- cf[tm, "Std. Error"]
    data.frame(trait = trait, model = model, term = tm, estimate = est, se = se,
               t = cf[tm, "t value"], lower = est - 1.96 * se, upper = est + 1.96 * se,
               n = nobs(fit), AIC = AIC(fit), singular = isSingular(fit),
               note = note, stringsAsFactors = FALSE)
  }))
}

# per-decade / percent-per-decade units for year terms
add_units <- function(tab, trait_mean, log_scale) {
  yr <- tab$term %in% c("scaled_yr", "yr_within_src", "yr_mean_src")
  pct <- function(x) if (log_scale) 100 * (exp(x) - 1) else 100 * x / trait_mean
  tab$per_decade       <- ifelse(yr, tab$estimate * DEC, NA_real_)
  tab$per_decade_lower <- ifelse(yr, tab$lower * DEC, NA_real_)
  tab$per_decade_upper <- ifelse(yr, tab$upper * DEC, NA_real_)
  tab$pct_per_decade       <- ifelse(yr, pct(tab$per_decade), NA_real_)
  tab$pct_per_decade_lower <- ifelse(yr, pct(tab$per_decade_lower), NA_real_)
  tab$pct_per_decade_upper <- ifelse(yr, pct(tab$per_decade_upper), NA_real_)
  tab
}

# ===========================================================================
# G1. Lagged climate covariates in the MEAN models
# ---------------------------------------------------------------------------
# Their result: for MASS the year trend disappears once lagged temperature is
# conditioned on; it persists for wing and mass:wing. Precipitation models fit
# better than temperature models in all nine comparisons.
# Here: does conditioning on within-locality climate anomalies change our
# (already null) year slope, and do Atlantic Forest birds respond to
# interannual climate anomalies at all?
# ===========================================================================
message("\n=== G1: lagged climate covariates in the mean models ===")

dc <- d[clim_ok, ]

# collinearity check, as they did (their max VIF < 2.01)
vif_check <- function(vars, data) {
  X <- data[complete.cases(data[vars]), vars, drop = FALSE]
  sapply(vars, function(v) {
    f <- as.formula(paste(v, "~", paste(setdiff(vars, v), collapse = " + ")))
    1 / (1 - summary(lm(f, data = X))$r.squared)
  })
}
res$G1_vif <- list(
  temperature   = vif_check(c("scaled_yr", "z_t_an_l0", "z_t_an_l1", "z_t_an_l2", "z_tmean_loc", "scaled_lat"), dc),
  precipitation = vif_check(c("scaled_yr", "z_p_an_l0", "z_p_an_l1", "z_p_an_l2", "z_prec_loc",  "scaled_lat"), dc),
  spei          = vif_check(c("scaled_yr", "z_s_an_l0", "z_s_an_l1", "z_s_an_l2", "z_spei_loc",  "scaled_lat"), dc))
res$G1_year_climate_cor <- c(
  r_yr_tanom_l0 = cor(dc$scaled_yr, dc$z_t_an_l0, use = "complete.obs"),
  r_yr_panom_l0 = cor(dc$scaled_yr, dc$z_p_an_l0, use = "complete.obs"),
  r_yr_sanom_l0 = cor(dc$scaled_yr, dc$z_s_an_l0, use = "complete.obs"))

TERMS_T <- c("z_t_an_l0", "z_t_an_l1", "z_t_an_l2", "z_tmean_loc")
TERMS_P <- c("z_p_an_l0", "z_p_an_l1", "z_p_an_l2", "z_prec_loc")
TERMS_S <- c("z_s_an_l0", "z_s_an_l1", "z_s_an_l2", "z_spei_loc")

# base RHS at two control tiers, mirroring the manuscript's ladder
rhs_M1 <- "Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) + (1 | src) + (1 | site)"
rhs_M3 <- paste("Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season +",
                "(1 + scaled_yr || spp) + (1 | src) + (1 | site) + (1 | ring)")

g1_rows <- list(); g1_aic <- list()
for (tr in c("ln_body_mass", "wing")) {
  log_sc <- tr == "ln_body_mass"
  for (tier in c("M1", "M3")) {
    rhs <- if (tier == "M1") rhs_M1 else rhs_M3
    dt  <- dc[!is.na(dc[[tr]]), ]
    if (tier == "M3") dt <- dt[!is.na(dt$season) & !is.na(dt$scaled_alt), ]
    dt$y <- dt[[tr]]
    tmean_raw <- if (log_sc) mean(exp(dt$y)) else mean(dt$y)

    specs <- list(
      year_only   = paste("y ~", rhs),
      temperature = paste("y ~", rhs, "+", paste(TERMS_T, collapse = " + ")),
      precipitation = paste("y ~", rhs, "+", paste(TERMS_P, collapse = " + ")),
      spei        = paste("y ~", rhs, "+", paste(TERMS_S, collapse = " + ")))

    for (nm in names(specs)) {
      f   <- as.formula(specs[[nm]])
      fit <- fit_lmer(f, dt, REML = TRUE)
      keep <- c("scaled_yr", TERMS_T, TERMS_P, TERMS_S)
      tb  <- tidy_lmer(fit, keep, tr, paste0(tier, "_", nm), note = attr(fit, "msgs"))
      if (!is.null(tb)) g1_rows[[length(g1_rows) + 1]] <- add_units(tb, tmean_raw, log_sc)
      # ML refit for the AIC comparison between climate variants (fixed effects differ)
      fitml <- fit_lmer(f, dt, REML = FALSE)
      if (!is.null(fitml)) g1_aic[[length(g1_aic) + 1]] <- data.frame(
        trait = tr, tier = tier, spec = nm, AIC_ML = AIC(fitml), n = nobs(fitml))
    }
    message("  [G1] ", tr, " ", tier, " done (n = ", nrow(dt), ")")
  }
}
res$G1_table <- bind_rows(g1_rows)
res$G1_aic   <- bind_rows(g1_aic) %>% group_by(trait, tier) %>%
  mutate(dAIC = AIC_ML - min(AIC_ML)) %>% ungroup() %>% arrange(trait, tier, AIC_ML)

# headline: how much does the year slope move when climate is conditioned on?
res$G1_year_shift <- res$G1_table %>% filter(term == "scaled_yr") %>%
  select(trait, model, estimate, lower, upper, pct_per_decade,
         pct_per_decade_lower, pct_per_decade_upper, per_decade,
         per_decade_lower, per_decade_upper, n)
print(as.data.frame(res$G1_year_shift %>% select(trait, model, pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n)))
cat("\n-- AIC (ML) across climate variants --\n"); print(as.data.frame(res$G1_aic))
cat("\n-- climate anomaly coefficients (M3 tier) --\n")
print(as.data.frame(res$G1_table %>% filter(term != "scaled_yr", grepl("^M3", model)) %>%
                      select(trait, model, term, estimate, lower, upper, t)))

# ===========================================================================
# G2. Species-level MASS slopes (they report per-species tallies; we only have
#     a species forest plot for wing)
# ===========================================================================
message("\n=== G2: species-level mass slopes ===")

species_slopes_lmer <- function(fit, term = "scaled_yr", grp = "spp") {
  re <- ranef(fit, condVar = TRUE)
  comp <- Filter(function(x) term %in% names(x), re)
  stopifnot(length(comp) >= 1)
  comp <- comp[[1]]
  # lme4 stores postVar either as one array (correlated terms) or, for `||`
  # models, as a NAMED LIST of 1x1xN arrays keyed by term.
  cv <- attr(comp, "postVar")
  se_blup <- if (is.list(cv)) {
    stopifnot(term %in% names(cv)); sqrt(as.numeric(cv[[term]]))
  } else if (!is.null(cv) && length(dim(cv)) == 3) {
    j <- match(term, names(comp)); sqrt(cv[j, j, ])
  } else NA_real_
  fe <- fixef(fit)[[term]]; se_fe <- summary(fit)$coefficients[term, "Std. Error"]
  data.frame(spp = rownames(comp), ranef_yr = comp[[term]], se_ranef = se_blup,
             fixed = fe, se_fixed = se_fe, row.names = NULL) %>%
    mutate(slope = fixed + ranef_yr, se_total = sqrt(se_ranef^2 + se_fixed^2),
           lo = slope - 1.96 * se_total, hi = slope + 1.96 * se_total)
}

g2 <- list()
for (tr in c("ln_body_mass", "wing")) {
  dt <- d[!is.na(d[[tr]]), ]; dt$y <- dt[[tr]]
  log_sc <- tr == "ln_body_mass"
  tm <- if (log_sc) mean(exp(dt$y)) else mean(dt$y)
  for (tier in c("M0", "M3")) {
    f <- if (tier == "M0") as.formula("y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)")
         else as.formula(paste("y ~", rhs_M3))
    dd <- if (tier == "M3") dt[!is.na(dt$season) & !is.na(dt$scaled_alt), ] else dt
    fit <- fit_lmer(f, dd)
    if (is.null(fit)) next
    sl <- species_slopes_lmer(fit) %>%
      mutate(trait = tr, model = tier,
             per_decade = slope * DEC, lo_dec = lo * DEC, hi_dec = hi * DEC,
             pct_per_decade = if (log_sc) 100 * (exp(per_decade) - 1) else 100 * per_decade / tm,
             pct_lo = if (log_sc) 100 * (exp(lo_dec) - 1) else 100 * lo_dec / tm,
             pct_hi = if (log_sc) 100 * (exp(hi_dec) - 1) else 100 * hi_dec / tm)
    g2[[paste0(tr, "_", tier)]] <- sl
  }
}
res$G2_species_slopes <- bind_rows(g2)

# the tally in Jirinec's terms: sign of point estimate, and CI excluding zero
res$G2_tally <- res$G2_species_slopes %>% group_by(trait, model) %>%
  summarise(n_spp = n(),
            n_negative_mean = sum(slope < 0),
            n_positive_mean = sum(slope > 0),
            n_CI_negative   = sum(hi < 0),
            n_CI_positive   = sum(lo > 0),
            n_CI_spans_zero = sum(lo <= 0 & hi >= 0),
            pct_CI_negative = round(100 * sum(hi < 0) / n(), 1),
            median_pct_per_decade = round(median(pct_per_decade), 3), .groups = "drop")
print(as.data.frame(res$G2_tally))

# ===========================================================================
# G3. mass:wing ratio as a response (their headline metric)
# ===========================================================================
message("\n=== G3: mass:wing ratio ===")
dsh <- d %>% filter(!is.na(ln_body_mass), !is.na(wing)) %>%
  mutate(ln_ratio  = ln_body_mass - log(wing),   # log(mass/wing): trend is %/decade
         raw_ratio = Body_mass.g. / wing)

g3 <- list()
for (resp in c("ln_ratio", "raw_ratio")) {
  dsh$y <- dsh[[resp]]
  tm <- mean(dsh$y)
  for (tier in c("M0", "M1", "M3")) {
    f <- switch(tier,
      M0 = as.formula("y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)"),
      M1 = as.formula(paste("y ~", rhs_M1)),
      M3 = as.formula(paste("y ~", rhs_M3)))
    dd <- if (tier == "M3") dsh[!is.na(dsh$season) & !is.na(dsh$scaled_alt), ] else dsh
    fit <- fit_lmer(f, dd); if (is.null(fit)) next
    tb <- tidy_lmer(fit, "scaled_yr", resp, tier, note = attr(fit, "msgs"))
    g3[[length(g3) + 1]] <- add_units(tb, tm, log_scale = (resp == "ln_ratio"))
    if (tier == "M3") {
      sl <- species_slopes_lmer(fit) %>%
        mutate(trait = resp, model = tier, per_decade = slope * DEC,
               lo_dec = lo * DEC, hi_dec = hi * DEC,
               pct_per_decade = if (resp == "ln_ratio") 100 * (exp(per_decade) - 1) else 100 * per_decade / tm,
               pct_lo = if (resp == "ln_ratio") 100 * (exp(lo_dec) - 1) else 100 * lo_dec / tm,
               pct_hi = if (resp == "ln_ratio") 100 * (exp(hi_dec) - 1) else 100 * hi_dec / tm)
      res$G3_species_slopes <- bind_rows(res$G3_species_slopes, sl)
    }
  }
}
res$G3_table <- bind_rows(g3)
res$G3_tally <- res$G3_species_slopes %>% group_by(trait, model) %>%
  summarise(n_spp = n(), n_negative_mean = sum(slope < 0), n_CI_negative = sum(hi < 0),
            n_CI_positive = sum(lo > 0), n_CI_spans_zero = sum(lo <= 0 & hi >= 0),
            pct_CI_negative = round(100 * sum(hi < 0) / n(), 1), .groups = "drop")
res$G3_n_shared <- nrow(dsh)
print(as.data.frame(res$G3_table %>% select(trait, model, estimate, lower, upper,
                                            pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n)))
print(as.data.frame(res$G3_tally))

# ===========================================================================
# G4. Phylogenetic signal in the SLOPES (they: Moran's I on per-species year
#     coefficients via phylosignal; we currently only have phylo signal on the
#     species INTERCEPTS via the propto term)
# ---------------------------------------------------------------------------
# CAVEAT (applies to them as much as to us): these are shrunk BLUPs carrying
# sampling error, so signal estimated on them is biased downward by shrinkage
# toward a common mean and upward by any phylogeny already inside the fitting
# model. We therefore report BOTH the M0 (weakly shrunk) and M3 slopes, and
# neither model contains a phylogenetic term, so the test is not circular here.
# ===========================================================================
message("\n=== G4: phylogenetic signal in the slopes ===")

Acache <- readRDS(derived_path("phylo_A_50trees.rds"))
ebird_syn <- c("Antilophia galeata"="Chiroxiphia galeata","Tachyphonus cristatus"="Loriotus cristatus",
  "Pyrrhocoma ruficeps"="Thlypopsis pyrrhocoma","Pyriglena pernambucensis"="Pyriglena leuconota",
  "Tangara sayaca"="Thraupis sayaca","Tangara cayana"="Stilpnia cayana",
  "Tangara palmarum"="Thraupis palmarum","Tangara peruviana"="Stilpnia peruviana",
  "Dixiphia pipra"="Pseudopipra pipra","Tiaris fuliginosus"="Asemospiza fuliginosa")
tip_of <- function(binomial) {
  nm <- gsub("_", " ", binomial); hit <- nm %in% names(ebird_syn); nm[hit] <- ebird_syn[nm[hit]]
  gsub(" ", "_", nm)
}

# Pagel's lambda by ML on a correlation matrix (profile out mu and sigma^2)
pagel_lambda <- function(y, A) {
  n <- length(y); one <- rep(1, n); I <- diag(n)
  nll <- function(lam) {
    V <- lam * A + (1 - lam) * I
    ch <- tryCatch(chol(V), error = function(e) NULL); if (is.null(ch)) return(1e10)
    Vi <- chol2inv(ch); ldet <- 2 * sum(log(diag(ch)))
    mu <- as.numeric((t(one) %*% Vi %*% y) / (t(one) %*% Vi %*% one))
    r  <- y - mu; s2 <- as.numeric(t(r) %*% Vi %*% r) / n
    0.5 * (n * log(2 * pi * s2) + ldet + n)
  }
  op <- optimize(nll, c(0, 1), tol = 1e-6)
  c(lambda = op$minimum, logLik = -op$objective, logLik0 = -nll(0),
    LRT = 2 * (-op$objective + nll(0)), p = pchisq(2 * (-op$objective + nll(0)), 1, lower.tail = FALSE))
}

signal_across_trees <- function(slopes_df, label) {
  sl <- slopes_df %>% mutate(tip = tip_of(as.character(spp))) %>% filter(tip %in% Acache$species)
  y  <- setNames(sl$per_decade, sl$tip)
  out <- lapply(seq_len(min(50, Acache$n_trees)), function(i) {
    A  <- Acache$A[[i]][names(y), names(y), drop = FALSE]
    dm <- 2 * (1 - A); diag(dm) <- 0                 # ultrametric, height-1: d = 2(1-A)
    w  <- 1 / dm; diag(w) <- 0; w[!is.finite(w)] <- 0
    mi <- ape::Moran.I(y, w, scaled = TRUE)
    pl <- pagel_lambda(as.numeric(y), A)
    data.frame(tree = i, moran_I = mi$observed, moran_expected = mi$expected,
               moran_sd = mi$sd, moran_p = mi$p.value,
               lambda = pl[["lambda"]], lambda_LRT = pl[["LRT"]], lambda_p = pl[["p"]])
  })
  bind_rows(out) %>% summarise(
    metric = label, n_spp = length(y), n_trees = n(),
    moran_I = mean(moran_I), moran_I_lo = quantile(moran_I, .025), moran_I_hi = quantile(moran_I, .975),
    moran_expected = mean(moran_expected),
    moran_p_median = median(moran_p), moran_p_max = max(moran_p),
    lambda = mean(lambda), lambda_lo = quantile(lambda, .025), lambda_hi = quantile(lambda, .975),
    lambda_p_median = median(lambda_p), lambda_p_max = max(lambda_p))
}

g4 <- list()
for (tr in unique(res$G2_species_slopes$trait)) for (md in unique(res$G2_species_slopes$model)) {
  s <- res$G2_species_slopes %>% filter(trait == tr, model == md)
  if (nrow(s) < 10) next
  g4[[paste(tr, md)]] <- signal_across_trees(s, paste0(tr, " (", md, ")"))
}
if (!is.null(res$G3_species_slopes))
  g4[["ratio"]] <- signal_across_trees(res$G3_species_slopes %>% filter(trait == "ln_ratio"), "ln(mass:wing) (M3)")
res$G4_slope_signal <- bind_rows(g4)
print(as.data.frame(res$G4_slope_signal %>%
  select(metric, n_spp, moran_I, moran_expected, moran_p_median, lambda, lambda_lo, lambda_hi, lambda_p_median)))

# ===========================================================================
# G5. Vertical foraging stratum as a moderator (the one trait that gave them a
#     signal: mass:wing declined for midstory species)
# ===========================================================================
message("\n=== G5: vertical foraging stratum ===")

elton <- suppressMessages(readr::read_tsv(raw_path("BirdFuncDat.txt"), show_col_types = FALSE))
elton$Scientific <- gsub(" ", "_", elton$Scientific)
strat <- elton %>%
  select(Scientific, ground = `ForStrat-ground`, understory = `ForStrat-understory`,
         midhigh = `ForStrat-midhigh`, canopy = `ForStrat-canopy`, aerial = `ForStrat-aerial`) %>%
  filter(Scientific %in% levels(d$spp)) %>%
  mutate(above_understory = (midhigh + canopy + aerial) / 100,
         dom = c("ground","understory","midhigh","canopy","aerial")[
           max.col(across(c(ground, understory, midhigh, canopy, aerial)), ties.method = "first")],
         stratum = factor(ifelse(dom %in% c("ground","understory"), "low", "high"), levels = c("low","high")))
res$G5_species_strata <- strat
res$G5_strata_counts <- strat %>% count(dom) %>% arrange(desc(n))
message("  [G5] matched ", nrow(strat), " / ", nlevels(d$spp), " species to EltonTraits strata")

ds <- d %>% left_join(strat %>% select(Scientific, above_understory, stratum),
                      by = c("Binomial" = "Scientific")) %>%
  filter(!is.na(above_understory)) %>%
  mutate(au_std = as.numeric(scale(above_understory)),
         ln_ratio = ln_body_mass - log(wing))

g5 <- list()
for (tr in c("ln_body_mass", "wing", "ln_ratio")) {
  dt <- ds[!is.na(ds[[tr]]), ]; dt$y <- dt[[tr]]
  log_sc <- tr %in% c("ln_body_mass", "ln_ratio")
  tm <- if (log_sc) mean(exp(dt$y)) else mean(dt$y)
  for (tier in c("M1", "M3")) {
    rhs <- if (tier == "M1") rhs_M1 else rhs_M3
    dd <- if (tier == "M3") dt[!is.na(dt$season) & !is.na(dt$scaled_alt), ] else dt
    for (mod in c("continuous", "categorical")) {
      extra <- if (mod == "continuous") "+ au_std + scaled_yr:au_std" else "+ stratum + scaled_yr:stratum"
      fit <- fit_lmer(as.formula(paste("y ~", rhs, extra)), dd)
      if (is.null(fit)) next
      keep <- c("scaled_yr", "au_std", "scaled_yr:au_std", "stratumhigh", "scaled_yr:stratumhigh")
      tb <- tidy_lmer(fit, keep, tr, paste0(tier, "_", mod), note = attr(fit, "msgs"))
      if (!is.null(tb)) g5[[length(g5) + 1]] <- add_units(tb, tm, log_sc)
    }
  }
}
res$G5_table <- bind_rows(g5)
print(as.data.frame(res$G5_table %>% filter(grepl(":", term)) %>%
                      select(trait, model, term, estimate, lower, upper, t, n)))

# ===========================================================================
# G6. Design-matched restricted tier -- approximating their single-site,
#     single-protocol, single-team design with our longest-running contributors
# ===========================================================================
message("\n=== G6: design-matched restricted tier ===")

src_cov <- d %>% group_by(src) %>%
  summarise(n = n(), n_mass = sum(!is.na(ln_body_mass)), n_wing = sum(!is.na(wing)),
            yr_min = min(Year), yr_max = max(Year), span = max(Year) - min(Year),
            n_years = n_distinct(Year), n_site = n_distinct(site), n_spp = n_distinct(spp),
            .groups = "drop") %>% arrange(desc(span), desc(n))
res$G6_contributor_coverage <- src_cov

SPAN_MIN <- 8L; N_MIN <- 150L
keep_src <- src_cov %>% filter(span >= SPAN_MIN, n >= N_MIN) %>% pull(src)
res$G6_criteria <- list(span_min = SPAN_MIN, n_min = N_MIN,
                        n_contributors_kept = length(keep_src),
                        contributors = as.character(keep_src))
d6 <- d %>% filter(src %in% keep_src) %>% mutate(ln_ratio = ln_body_mass - log(wing))
res$G6_n <- c(n_records = nrow(d6), n_spp = n_distinct(d6$spp),
              n_mass = sum(!is.na(d6$ln_body_mass)), n_wing = sum(!is.na(d6$wing)),
              n_src = length(keep_src), yr_min = min(d6$Year), yr_max = max(d6$Year))
message("  [G6] ", length(keep_src), " contributors, ", nrow(d6), " records")

fit_ladder <- function(data, tag) {
  out <- list()
  for (tr in c("ln_body_mass", "wing", "ln_ratio")) {
    dt <- data[!is.na(data[[tr]]), ]; if (nrow(dt) < 200) next
    dt$y <- dt[[tr]]; log_sc <- tr %in% c("ln_body_mass", "ln_ratio")
    tm <- if (log_sc) mean(exp(dt$y)) else mean(dt$y)
    for (tier in c("M0", "M1", "M3")) {
      f <- switch(tier,
        M0 = as.formula("y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)"),
        M1 = as.formula(paste("y ~", rhs_M1)),
        M3 = as.formula(paste("y ~", rhs_M3)))
      dd <- if (tier == "M3") dt[!is.na(dt$season) & !is.na(dt$scaled_alt), ] else dt
      fit <- fit_lmer(f, dd); if (is.null(fit)) next
      tb <- tidy_lmer(fit, "scaled_yr", tr, paste0(tag, "_", tier), note = attr(fit, "msgs"))
      if (!is.null(tb)) out[[length(out) + 1]] <- add_units(tb, tm, log_sc)
    }
  }
  bind_rows(out)
}
res$G6_table <- fit_ladder(d6, "restricted")
print(as.data.frame(res$G6_table %>% select(trait, model, estimate, lower, upper,
                                            pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n)))

# ===========================================================================
# G7. Their symmetric coverage filter: >= 5 records BOTH before and after the
#     median year. Our filter (>= 30 records, >= 5-year span) does not
#     guarantee a species is present in both halves of the record.
# ===========================================================================
message("\n=== G7: symmetric coverage filter ===")
MED_YR <- median(d$Year)
cov_spp <- d %>% group_by(spp) %>%
  summarise(n = n(), n_before = sum(Year <= MED_YR), n_after = sum(Year > MED_YR),
            n_mass_before = sum(Year <= MED_YR & !is.na(ln_body_mass)),
            n_mass_after  = sum(Year >  MED_YR & !is.na(ln_body_mass)),
            n_wing_before = sum(Year <= MED_YR & !is.na(wing)),
            n_wing_after  = sum(Year >  MED_YR & !is.na(wing)), .groups = "drop") %>%
  mutate(passes_jirinec_mass = n_mass_before >= 5 & n_mass_after >= 5,
         passes_jirinec_wing = n_wing_before >= 5 & n_wing_after >= 5)
res$G7_median_year <- MED_YR
res$G7_species_coverage <- cov_spp
res$G7_summary <- data.frame(
  median_year = MED_YR, n_spp_total = nrow(cov_spp),
  n_spp_pass_mass = sum(cov_spp$passes_jirinec_mass),
  n_spp_fail_mass = sum(!cov_spp$passes_jirinec_mass),
  n_spp_pass_wing = sum(cov_spp$passes_jirinec_wing),
  n_spp_fail_wing = sum(!cov_spp$passes_jirinec_wing))
print(res$G7_summary)

d7m <- d %>% filter(spp %in% cov_spp$spp[cov_spp$passes_jirinec_mass]) %>%
  mutate(ln_ratio = ln_body_mass - log(wing))
res$G7_table <- fit_ladder(d7m, "symmetric")
print(as.data.frame(res$G7_table %>% select(trait, model, estimate, lower, upper,
                                            pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n)))

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
res$elapsed_min <- round(as.numeric(difftime(Sys.time(), res$generated, units = "mins")), 2)
saveRDS(res, gap_path("jirinec_gaps_results.rds"))
message("\n[write] ", gap_path("jirinec_gaps_results.rds"), "  (", res$elapsed_min, " min)")
