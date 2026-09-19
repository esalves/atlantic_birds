# referee_reruns.R
# ---------------------------------------------------------------------------
# Co-author review re-runs commissioned 2026-09-18 (see REVIEW_RESPONSE_PLAN.md).
# Addresses referee items E1-E4, E6-E7:
#
#   E1  Thermal models: annual temperature anomaly still covaries with calendar
#       year over a warming study period. Refit with year included simultaneously,
#       and with a within-site DETRENDED anomaly, to separate the long-term trend
#       from year-to-year thermal deviation.              (Mizuno c84)
#   E2  The anomaly is shared by every bird caught at the same locality in the
#       same year, so individuals are not independent thermal observations.
#       Add a locality-year random intercept and cluster-robust SEs. (Mizuno c85)
#   E3  Thermal models used a more parsimonious adjustment set than the main
#       temporal models. Refit at the full M3 level across 50 trees. (Mizuno c86)
#   E4  Same question for the bivariate isometry model, which omitted longitude,
#       elevation, season and moult. Refit at M3. Also fits the isometry contrast
#       DIRECTLY as a response (ln M - 3 ln L), which gives an exact SE and
#       sidesteps the delta-method slope covariance.       (Mizuno c81, c82)
#   E6  The secular series bridges museum skins and live captures, although the
#       primary analysis excludes skins because preparation shrinks them. Report
#       the source-adjusted vs unadjusted slope, test year x Status, and run a
#       museum-only sensitivity on the spanning cohort.     (Mizuno c41)
#   E7  The secular series pools all age classes although the primary analysis
#       excludes juveniles. Add Age as a covariate and run an adult-only
#       sensitivity.                                        (Mizuno c42)
#
# NOTE ON REPRODUCTION: the script that produced the reported thermal-coupling
# estimates (beta = +0.0212, t = 2.17 on passer90) does NOT exist in the repository --
# only its figure and its numbers in output/secular_expansion/REPORT.md. Block E0
# below refits that model from passer90_climate.rds so the reported number has a
# reproducible source; every E1-E3 variant is a modification of E0.
#
# OUTPUT: output/referee_reruns/referee_reruns.rds  (all tables)
#         output/referee_reruns/REPORT.md           (human-readable summary)
#
# RUN: Rscript Analysis/scripts/referee_reruns.R
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(glmmTMB); library(lme4); library(readr)
  library(sandwich); library(clubSandwich)
})

options(width = 140)
set.seed(20240101)

# ---- paths (same convention as the other scripts) -------------------------
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

# Year scaling is read from the data below, once passer90_climate.rds is loaded,
# rather than hardcoded (see SCRIPT_CONSTANTS_AUDIT.md).
SD_YR <- NA_real_
DEC   <- NA_real_
N_TREES <- as.integer(Sys.getenv("REFEREE_N_TREES", "50"))

RES <- list()   # collected result tables
note <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

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

# Rubin pooling that tolerates trees whose Hessian was not positive-definite.
# The shared engine's pool_rubin_df() takes mean(se^2) without na.rm, so a single
# failed tree turns the pooled SE into NA. Here we drop non-finite per-tree fits
# and record how many were used.
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

safely <- function(label, expr) {
  r <- tryCatch(force(expr), error = function(e) { note("FAILED ", label, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) note("  -> skipped ", label)
  r
}

# ===========================================================================
# DATA: modern analytical sample with climate
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
    season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")),
    molt   = factor(ifelse(is.na(Molt), "unknown", as.character(Molt))),
    ind    = factor(paste(Ring, Binomial, sep = "_")),
    wing   = conc.wing.length,
    lnwing = log(conc.wing.length),
    lnmass = ln_body_mass,
    species_name = phylo_species_name(Binomial)
  ) %>%
  # local annual temperature anomaly: record-year site temperature minus the
  # site's own long-run mean over the sampled years (site-specific baseline)
  group_by(loc_id) %>%
  mutate(t_base = mean(rec_tmean, na.rm = TRUE),
         dT     = rec_tmean - t_base) %>%
  ungroup() %>%
  mutate(loc_year = factor(paste(loc_id, Year, sep = "_"))) %>%
  as.data.frame()

# E1: within-site DETRENDED anomaly -- remove the part of dT explained by year,
# so the remaining variation is genuine interannual deviation.
dT_fit  <- lm(dT ~ Year, data = dat[!is.na(dat$dT), ])
dat$dT_detr <- NA_real_
ok <- !is.na(dat$dT)
dat$dT_detr[ok] <- residuals(dT_fit)

shared <- dat %>% filter(!is.na(wing), !is.na(lnmass), !is.na(dT)) %>%
  mutate(iso = lnmass - 3 * lnwing,           # exact log-scale isometry contrast
         relwing = lnwing - (1/3) * lnmass)   # relative wing length
note(sprintf("shared wing+mass+climate records: %d (%d species, %d locality-years)",
             nrow(shared), n_distinct(shared$species_name), nlevels(droplevels(shared$loc_year))))

cat(sprintf("\ndT: sd=%.3f, range [%.2f, %.2f]; cor(dT, Year)=%.3f; cor(dT_detr, Year)=%.3f\n",
            sd(shared$dT), min(shared$dT), max(shared$dT),
            cor(shared$dT, shared$Year), cor(shared$dT_detr, shared$Year)))

A_list <- phylo_A_list(sort(unique(shared$species_name)), n_trees = N_TREES)
note(sprintf("phylogeny: %d trees, %d species", length(A_list), nrow(A_list[[1]])))

# ===========================================================================
# E0. Reconstruct the parsimonious thermal models reported in REPORT.md 5.5
# ===========================================================================
note("E0: reconstructing the parsimonious thermal models reported in REPORT.md 5.5 ...")

par_rhs <- function(resp, xvar) {
  as.formula(sprintf("%s ~ Sex + %s + scaled_lat + (1 + %s || spp) + (1 | src)", resp, xvar, xvar))
}

RES$E0 <- safely("E0", {
  do.call(rbind, lapply(
    list(c("wing", "wing length (mm)"), c("lnmass", "log body mass"),
         c("iso", "isometry contrast lnM-3lnL"), c("relwing", "relative wing lnL-(1/3)lnM")),
    function(x) {
      fit <- glmmTMB(par_rhs(x[1], "dT"), data = shared, REML = TRUE)
      tidy_tmb(fit, "E0_parsimonious", x[2], terms = "dT")
    }))
})
print(RES$E0)

# ===========================================================================
# E1. Anomaly + calendar year simultaneously, and detrended anomaly
# ===========================================================================
note("E1: separating the secular trend from interannual thermal deviation ...")

RES$E1 <- safely("E1", {
  rows <- list()
  for (x in list(c("wing", "wing length (mm)"), c("lnmass", "log body mass"),
                 c("iso", "isometry contrast lnM-3lnL"), c("relwing", "relative wing lnL-(1/3)lnM"))) {
    resp <- x[1]; lab <- x[2]
    # (a) anomaly AND calendar year in the same model
    f_a <- as.formula(sprintf(
      "%s ~ Sex + dT + scaled_yr + scaled_lat + (1 + scaled_yr || spp) + (1 | src)", resp))
    fa <- glmmTMB(f_a, data = shared, REML = TRUE)
    rows[[length(rows) + 1]] <- tidy_tmb(fa, "E1a_anomaly_plus_year", lab, terms = c("dT", "scaled_yr"))
    # (b) within-site detrended anomaly
    fb <- glmmTMB(par_rhs(resp, "dT_detr"), data = shared, REML = TRUE)
    rows[[length(rows) + 1]] <- tidy_tmb(fb, "E1b_detrended_anomaly", lab, terms = "dT_detr")
  }
  do.call(rbind, rows)
})
print(RES$E1)

# ===========================================================================
# E2. Locality-year clustering of the temperature effect
# ===========================================================================
note("E2: locality-year clustering ...")

RES$E2 <- safely("E2", {
  rows <- list()
  for (x in list(c("wing", "wing length (mm)"), c("lnmass", "log body mass"),
                 c("iso", "isometry contrast lnM-3lnL"), c("relwing", "relative wing lnL-(1/3)lnM"))) {
    resp <- x[1]; lab <- x[2]
    # (a) locality-year random intercept -- the shared-exposure unit
    f <- as.formula(sprintf(
      "%s ~ Sex + dT + scaled_lat + (1 + dT || spp) + (1 | src) + (1 | loc_year)", resp))
    fit <- glmmTMB(f, data = shared, REML = TRUE)
    rows[[length(rows) + 1]] <- tidy_tmb(fit, "E2a_locyear_random_intercept", lab, terms = "dT")
    # (b) cluster-robust SE by locality-year on the parsimonious fit (lme4 + clubSandwich)
    cr <- tryCatch({
      fl <- lmer(as.formula(sprintf("%s ~ Sex + dT + scaled_lat + (1 + dT || spp) + (1 | src)", resp)),
                 data = shared, REML = TRUE,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      ct <- clubSandwich::coef_test(fl, vcov = "CR2", cluster = shared$loc_year, test = "Satterthwaite")
      ct <- ct[ct$Coef == "dT", ]
      data.frame(model = "E2b_clusterrobust_CR2", response = lab, term = "dT",
                 estimate = ct$beta, se = ct$SE, z = ct$tstat, p = ct$p_Satt,
                 n = nobs(fl), lower = ct$beta - 1.96 * ct$SE, upper = ct$beta + 1.96 * ct$SE,
                 row.names = NULL)
    }, error = function(e) { note("  CR2 failed for ", resp, ": ", conditionMessage(e)); NULL })
    if (!is.null(cr)) rows[[length(rows) + 1]] <- cr
  }
  do.call(rbind, rows)
})
print(RES$E2)

# ===========================================================================
# E3. Thermal models at the full M3 adjustment level, across 50 trees
# ===========================================================================
note("E3: thermal models at M3 adjustment across ", N_TREES, " trees ...")

m3_thermal <- function(resp, xvar) as.formula(sprintf(
  "%s ~ Sex + %s + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + %s || spp) + (1 | src) + (1 | ind)",
  resp, xvar, xvar))

RES$E3 <- safely("E3", {
  d3 <- shared %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season))
  note(sprintf("  M3 thermal sample: %d records", nrow(d3)))
  A3 <- phylo_A_list(sort(unique(d3$species_name)), n_trees = N_TREES)
  rows <- list()
  for (x in list(c("wing", "wing length (mm)"), c("lnmass", "log body mass"),
                 c("iso", "isometry contrast lnM-3lnL"), c("relwing", "relative wing lnL-(1/3)lnM"))) {
    resp <- x[1]; lab <- x[2]
    r <- run_phylo_trees(m3_thermal(resp, "dT"), d3, A3, verbose = FALSE)
    p <- pool_clean(r, "dT") %>%
      transmute(model = "E3_M3_adjusted_phylo", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p)
    # same model, plus calendar year (E1 + E3 combined: the strictest test)
    r2 <- run_phylo_trees(
      as.formula(sprintf("%s ~ Sex + dT + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + scaled_yr || spp) + (1 | src) + (1 | ind)", resp)),
      d3, A3, verbose = FALSE)
    p2 <- pool_clean(r2, c("dT", "scaled_yr")) %>%
      transmute(model = "E3b_M3_adjusted_phylo_plus_year", response = lab, term = par,
                estimate, se, z, p, n = nrow(d3), lower, upper, trees_used)
    rows[[length(rows) + 1]] <- as.data.frame(p2)
    note("  done: ", lab)
  }
  do.call(rbind, rows)
})
print(RES$E3)

# ===========================================================================
# E4. Isometry contrast under the full M3 adjustment set
# ===========================================================================
note("E4: isometry contrast at M3 adjustment ...")

RES$E4 <- safely("E4", {
  d4 <- shared %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season))
  A4 <- phylo_A_list(sort(unique(d4$species_name)), n_trees = N_TREES)
  rows <- list()

  # (a) contrast fitted DIRECTLY as a response -- exact SE, no delta method
  for (spec in list(
    list(id = "E4a_iso_parsimonious",
         f = iso ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) + (1 | src)),
    list(id = "E4b_iso_M3",
         f = iso ~ Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season + molt +
               (1 + scaled_yr || spp) + (1 | src) + (1 | ind)))) {
    r <- run_phylo_trees(spec$f, d4, A4, verbose = FALSE)
    p <- pool_clean(r, "scaled_yr") %>%
      transmute(model = spec$id, response = "isometry contrast lnM-3lnL", term = par,
                estimate, se, z, p, n = nrow(d4), lower, upper, trees_used,
                pct_decade = 100 * estimate * DEC,
                pct_decade_lo = 100 * lower * DEC, pct_decade_hi = 100 * upper * DEC)
    rows[[length(rows) + 1]] <- as.data.frame(p)
    note("  done: ", spec$id)
  }

  # (b) separate wing / mass slopes under M3, for comparison with the reported
  #     delta-method contrast (Mizuno c82: how was the covariance obtained?)
  for (resp in c("lnwing", "lnmass")) {
    r <- run_phylo_trees(
      as.formula(sprintf("%s ~ Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + scaled_yr || spp) + (1 | src) + (1 | ind)", resp)),
      d4, A4, verbose = FALSE)
    p <- pool_clean(r, "scaled_yr") %>%
      transmute(model = "E4c_separate_M3_slopes", response = resp, term = par,
                estimate, se, z, p, n = nrow(d4), lower, upper, trees_used,
                pct_decade = 100 * (exp(estimate * DEC) - 1))
    rows[[length(rows) + 1]] <- as.data.frame(p)
  }
  dplyr::bind_rows(rows)
})
print(RES$E4)

# ===========================================================================
# E6 / E7. Secular museum series: source calibration and age composition
# ===========================================================================
note("E6/E7: secular museum series ...")

RES$secular <- safely("E6/E7", {
  birds_raw <- read_csv(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"),
                        guess_max = 70000, show_col_types = FALSE)
  CTR_YR <- as.numeric(attr(pc$scaled_yr, "scaled:center"))   # read from the data, not hardcoded
  d_base <- birds_raw %>%
    filter(AtlanticForests_20km_Buffer == "inside the 20 km polygon", !is.na(Year)) %>%
    mutate(conc_wing = coalesce(Wing_length_right.mm., Wing_length_left.mm., Wing_length.mm.),
           Binomial  = factor(gsub(" ", "_", Binomial)),
           Status    = factor(ifelse(Status %in% c("live", "museum"), Status, NA_character_)),
           site      = factor(Municipality),
           Age3      = factor(ifelse(is.na(Age) | Age == "Unknown", "Unknown",
                                     ifelse(Age == "Adult", "Adult", "Immature")),
                              levels = c("Adult", "Immature", "Unknown")),
           scaled_yr = (Year - CTR_YR) / SD_YR)

  # spanning cohort, as in secular_museum_expansion.R Tier 6
  spp_pre  <- d_base %>% filter(Status == "museum", Year < 1990, !is.na(conc_wing)) %>% pull(Binomial) %>% unique()
  spp_post <- d_base %>% filter(Status == "live",  Year >= 1995, !is.na(conc_wing)) %>% pull(Binomial) %>% unique()
  d6 <- d_base %>% filter(Binomial %in% intersect(spp_pre, spp_post),
                          Status %in% c("live", "museum"), Year >= 1880, !is.na(conc_wing))
  MW <- mean(d6$conc_wing)
  note(sprintf("  spanning cohort: %d records, %d species (%d museum / %d live); Age: %s",
               nrow(d6), n_distinct(d6$Binomial), sum(d6$Status == "museum"), sum(d6$Status == "live"),
               paste(names(table(d6$Age3)), table(d6$Age3), sep = "=", collapse = ", ")))

  CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  grab <- function(fit, id, desc, data) {
    if (is.null(fit)) return(NULL)
    cf <- summary(fit)$coefficients
    if (!"scaled_yr" %in% rownames(cf)) return(NULL)
    e <- cf["scaled_yr", 1]; s <- cf["scaled_yr", 2]
    data.frame(model = id, sample = desc, term = "scaled_yr", estimate = e, se = s,
               t = cf["scaled_yr", 3], n = nobs(fit),
               mm_decade = e * DEC, pct_decade = 100 * e * DEC / MW,
               pct_lo = 100 * (e - 1.96 * s) * DEC / MW,
               pct_hi = 100 * (e + 1.96 * s) * DEC / MW, row.names = NULL)
  }
  ff <- function(f, d) tryCatch(lmer(f, data = d, control = CTRL), error = function(e) NULL)

  rows <- list()
  # spec as run in secular_museum_expansion.R (Status intercept offset only)
  rows[[1]] <- grab(ff(conc_wing ~ scaled_yr + Status + (1 + scaled_yr || Binomial) + (1 | site), d6),
                    "S0_asrun_status_offset", "spanning cohort, Status intercept", d6)
  # E6a: NO source adjustment at all -- how much does the offset matter?
  rows[[2]] <- grab(ff(conc_wing ~ scaled_yr + (1 + scaled_yr || Binomial) + (1 | site), d6),
                    "E6a_no_source_adjustment", "spanning cohort, no Status term", d6)
  # E6b: year x Status interaction -- does the slope itself differ by source?
  f6b <- ff(conc_wing ~ scaled_yr * Status + (1 + scaled_yr || Binomial) + (1 | site), d6)
  if (!is.null(f6b)) {
    cf <- summary(f6b)$coefficients
    rn <- grep("scaled_yr:Status", rownames(cf), value = TRUE)
    rows[[3]] <- data.frame(model = "E6b_year_x_status_interaction", sample = "spanning cohort",
                            term = c("scaled_yr", rn), estimate = cf[c("scaled_yr", rn), 1],
                            se = cf[c("scaled_yr", rn), 2], t = cf[c("scaled_yr", rn), 3],
                            n = nobs(f6b), mm_decade = cf[c("scaled_yr", rn), 1] * DEC,
                            pct_decade = 100 * cf[c("scaled_yr", rn), 1] * DEC / MW,
                            pct_lo = NA_real_, pct_hi = NA_real_, row.names = NULL)
  }
  # E6c: museum specimens only (no bridging across sources at all)
  dm <- d6 %>% filter(Status == "museum")
  rows[[4]] <- grab(ff(conc_wing ~ scaled_yr + (1 + scaled_yr || Binomial) + (1 | site), dm),
                    "E6c_museum_only", sprintf("spanning cohort, museum only (n=%d)", nrow(dm)), dm)
  # E6d: live captures only
  dl <- d6 %>% filter(Status == "live")
  rows[[5]] <- grab(ff(conc_wing ~ scaled_yr + (1 + scaled_yr || Binomial) + (1 | site), dl),
                    "E6d_live_only", sprintf("spanning cohort, live only (n=%d)", nrow(dl)), dl)
  # E7a: Age as a covariate alongside the source offset
  rows[[6]] <- grab(ff(conc_wing ~ scaled_yr + Status + Age3 + (1 + scaled_yr || Binomial) + (1 | site), d6),
                    "E7a_age_covariate", "spanning cohort, Status + Age3", d6)
  # E7b: adults only, matching the primary analysis filter
  da <- d6 %>% filter(Age3 == "Adult")
  rows[[7]] <- grab(ff(conc_wing ~ scaled_yr + Status + (1 + scaled_yr || Binomial) + (1 | site), da),
                    "E7b_adults_only", sprintf("spanning cohort, adults only (n=%d)", nrow(da)), da)
  # E7c: adults only AND age-composition check -- is Age composition trending?
  age_trend <- glm(I(Age3 == "Adult") ~ Year, family = binomial, data = d6)
  note(sprintf("  P(Adult) ~ Year: beta=%.4f, p=%.4g",
               coef(age_trend)[2], summary(age_trend)$coefficients[2, 4]))

  attr_out <- dplyr::bind_rows(Filter(Negate(is.null), rows))
  attr(attr_out, "age_trend") <- summary(age_trend)$coefficients
  attr_out
})
print(RES$secular)

# ===========================================================================
# SAVE
# ===========================================================================
RES$meta <- list(generated = Sys.time(), n_trees = N_TREES, SD_YR = SD_YR,
                 glmmTMB = as.character(packageVersion("glmmTMB")),
                 R = R.version.string,
                 items = "E1-E4, E6-E7 (co-author review, 2026-09-18)")
saveRDS(RES, file.path(OUT_DIR, "referee_reruns.rds"))
note("saved ", file.path(OUT_DIR, "referee_reruns.rds"))

# ---- markdown report ------------------------------------------------------
fmt <- function(df, digits = 4) {
  if (is.null(df)) return("_(block failed -- see log)_\n")
  num <- sapply(df, is.numeric)
  df[num] <- lapply(df[num], function(x) signif(x, digits))
  paste0(paste(capture.output(print(df, row.names = FALSE)), collapse = "\n"), "\n")
}
con <- file(file.path(OUT_DIR, "REPORT.md"), "w")
writeLines(c(
  "# Referee re-runs (E1-E4, E6-E7)",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"), " | trees: ", N_TREES,
         " | glmmTMB ", packageVersion("glmmTMB")),
  "",
  "Commissioned from the co-author review of 2026-09-18. See `REVIEW_RESPONSE_PLAN.md`.",
  "",
  "**Reproduction note.** The script that produced the reported thermal-coupling",
  "estimates was not in the repository; block E0 refits that model from",
  "`passer90_climate.rds` so the reported numbers have a reproducible source.",
  "",
  "## E0. Published parsimonious thermal models (reproduction)", "", "```", fmt(RES$E0), "```", "",
  "## E1. Anomaly with calendar year, and detrended anomaly (Mizuno c84)", "", "```", fmt(RES$E1), "```", "",
  "## E2. Locality-year clustering (Mizuno c85)", "", "```", fmt(RES$E2), "```", "",
  "## E3. Thermal models at M3 adjustment across trees (Mizuno c86)", "", "```", fmt(RES$E3), "```", "",
  "## E4. Isometry contrast at M3 adjustment (Mizuno c81, c82)", "", "```", fmt(RES$E4), "```", "",
  "## E6/E7. Secular series: source calibration and age (Mizuno c41, c42)", "", "```", fmt(RES$secular), "```", ""
), con)
close(con)
note("wrote ", file.path(OUT_DIR, "REPORT.md"))
note("DONE")
