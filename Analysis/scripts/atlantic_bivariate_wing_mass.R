# atlantic_bivariate_wing_mass.R
# ---------------------------------------------------------------------------
# WHAT: Isometry test on the records that carry BOTH wing length and body mass
# (7,577 live, known-sex, adult records; 73 species). Two brms models fitted
# per phylogenetic tree and pooled across trees with Rubin's rules:
#
#   (1) BIVARIATE  mvbind(wing, lnmass) ~ 1 + Sex + scaled_yr + scaled_lat +
#                    (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)) +
#                    (1 | src) + (1 | site),   set_rescor(TRUE)
#       -> each response has its own species intercept + year slope, its own
#          phylogenetic term, contributor (Main_researcher) and site
#          (Municipality) intercepts; residuals are correlated across the two
#          responses within a record. The DERIVED quantity is the difference of
#          the two year slopes expressed in % per decade, with its posterior
#          interval, plus the departure from isometry (observed mass change
#          minus the (1 + w)^3 - 1 expectation).
#   (2) RELATIVE WING  wing ~ 1 + Sex + scaled_yr + scaled_lat + ln_body_mass +
#                    (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)) +
#                    (1 | src) + (1 | site)
#       -> the year slope of wing length conditional on body mass, and the
#          within-model allometric coefficient of log mass.
#
# WHY (REVISION_PLAN.md §1 "§2.5 Mass and allometry", Phase 2 item 1): the
# manuscript claims wings shorten while body mass does not. The referee asked
# for that contrast to be tested on the same individuals with correlated
# errors, so that "the difference in year slopes gets its own posterior", with
# the provenance and site terms of Phase 1, and rephrased as departure from
# isometry rather than wing loading.
#
# ENGINES (2026-09-09): the phylogenetic tier has two engines.
#   --engine glmmTMB   DEFAULT for --trees (and for a bare run). glmmTMB's `propto`
#                 covariance structure (Williams, McGillycuddy, Drobniak, Bolker,
#                 Warton & Nakagawa 2025, bioRxiv 10.64898/2025.12.20.695312) via
#                 the shared engine _phylo_engine.R, on the IDENTICAL 50 trees of
#                 the published brms analysis (phylo_A_list()), Rubin-pooled.
#                 glmmtmb_validation_wing.R shows it reproduces the published
#                 brms wing model to 2-3 decimals in 24 s. Fits, on the shared
#                 subset, M0 / M1 / M1cc / M2 x {wing, log mass, wing | log mass}
#                 univariate phylogenetic models, the slope difference in % per
#                 decade with a delta-method interval that ASSUMES the two
#                 estimates are independent, and a JOINT long-format model
#                 (record x trait rows) with trait-specific phylogenetic species
#                 effects (two propto terms), whose slope difference carries the
#                 estimated covariance. Writes output/bivariate_phylo_results.rds.
#                 Minutes on a laptop. No brms, no Stan.
#   --engine brms      The original Bayesian path below (mvbind + set_rescor(TRUE)
#                 bivariate model and the relative-wing model), 50 trees, SAMPLING
#                 settings: the cross-check for Totoro. Hours. Writes
#                 output/bivariate_results.rds and models/bivariate_wing_mass_multiphylo.rda.
#   --smoke is always brms (it is a brms wiring test); --fast-lme4 is always lme4.
#
# MODES (command-line flags; default = --engine glmmTMB, 50 trees):
#   --fast-lme4   lme4 (REML) analogues on the SAME shared subset, no phylogeny:
#                 wing and log-mass univariate models (baseline; + src + site;
#                 + lon + alt + season), the relative-wing model, and a stacked
#                 glmmTMB joint model (shared record intercept, trait-specific
#                 residual variance) on the 100*ln scale, which gives the slope
#                 difference WITH its covariance. Runs in ~1-2 minutes.
#                 Writes output/bivariate_fast_lme4.rds.
#   --smoke       brms SMOKE TEST ONLY: 1 tree, chains = 2, iter = 400,
#                 warmup = 200, rstan. Checks that the models compile and run
#                 and that the parameter names / derived quantities are wired
#                 correctly. Its numbers are NOT results. Writes
#                 output/bivariate_smoke_results.rds and
#                 output/models/bivariate_smoke.rda (never the full-run files).
#   --subsample N (with --smoke only) fit on a random subset of N shared records
#                 so the compile / parameter-name / derived-quantity wiring can be
#                 checked in minutes on a laptop. Never a result.
#   --trees N     number of trees for the phylogenetic tier (default 50, as in
#                 atlantic_parallel.R; both engines). N = 1 is the "single-tree
#                 validation" tier of the plan.
#
# UNITS: scaled_yr = (Year - 2009.487) / 5.0199; per decade = beta * 10 / 5.0199.
# Wing % per decade = 100 * per_decade / mean(wing on the shared subset); mass
# % per decade = 100 * (exp(per_decade) - 1).
#
# INPUTS:  data/derived/passer90.rda.
#          glmmTMB engine: scripts/_phylo_engine.R and its cache
#          data/derived/phylo_A_50trees.rds (the 50 species correlation matrices
#          of the published brms fits, models/brm0_multiphylo.rda; clootl
#          fallback inside the engine); output/bivariate_fast_lme4.rds (optional,
#          for the lme4-tier comparison table).
#          brms engine: clootl / AvesData tree cloud (McTavish et al. 2025) via
#          prepR4pcm, as in atlantic_parallel.R. If AVESDATA_PATH is unset and
#          data/raw/AvesDataLite-main exists it is used (no download); otherwise
#          clootl::get_avesdata_repo(raw_path()) as before.
# OUTPUTS: output/bivariate_phylo_results.rds       (--engine glmmTMB, default)
#          output/bivariate_results.rds             (--engine brms full run)
#          output/models/bivariate_wing_mass_multiphylo.rda (brms full run: fits,
#                                                    Rubin summaries, clootl version)
#          output/bivariate_fast_lme4.rds            (--fast-lme4)
#          output/bivariate_smoke_results.rds, output/models/bivariate_smoke.rda (--smoke)
#
# RUN:  Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --fast-lme4
#       Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --trees 50                 # glmmTMB (default engine), minutes
#       Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --smoke [--subsample 800]  # brms wiring test
#       Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --engine brms --trees 50   # Totoro, hours
# Session used locally (2026-09-09): R 4.6.0, brms 2.23.0, rstan 2.32.7, lme4
# 2.0.1, glmmTMB 1.1.14, prepR4pcm 0.5.0.9000, clootl 0.1.4, ape 5.8.1,
# MCMCglmm 2.36, phytools 2.5.2, posterior 1.7.0, future.apply 1.20.2.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
MODE_FAST  <- "--fast-lme4" %in% args
MODE_SMOKE <- "--smoke" %in% args
N_TREES    <- if (any(grepl("^--trees(=|$)", args))) {
  i <- grep("^--trees", args)[1]
  v <- sub("^--trees=?", "", args[i]); if (!nzchar(v)) v <- args[i + 1]
  as.integer(v)
} else 50L
if (MODE_SMOKE) N_TREES <- 1L
# --subsample N : fit the brms models on a random subset of N shared records
# (wiring / compile check only; NEVER a result). Only honoured together with --smoke.
SUBSAMPLE_N <- if (any(grepl("^--subsample(=|$)", args))) {
  i <- grep("^--subsample", args)[1]
  v <- sub("^--subsample=?", "", args[i]); if (!nzchar(v)) v <- args[i + 1]
  as.integer(v)
} else NA_integer_
if (!is.na(SUBSAMPLE_N) && !MODE_SMOKE) stop("--subsample is only allowed with --smoke")
# --engine glmmTMB|brms : engine of the phylogenetic tier. Default glmmTMB (propto,
# _phylo_engine.R). --smoke is always brms (it tests the brms wiring); --fast-lme4 is lme4.
ENGINE <- if (any(grepl("^--engine(=|$)", args))) {
  i <- grep("^--engine", args)[1]
  v <- sub("^--engine=?", "", args[i]); if (!nzchar(v)) v <- args[i + 1]
  match.arg(tolower(v), c("glmmtmb", "brms"))
} else "glmmtmb"
ENGINE <- c(glmmtmb = "glmmTMB", brms = "brms")[[ENGINE]]
if (MODE_SMOKE) ENGINE <- "brms"
if (MODE_FAST)  ENGINE <- "lme4"

suppressMessages({ library(dplyr) })
set.seed(20240101)   # same seed as atlantic_parallel.R (tree subsample)

# --- robust path resolution (repo layout: Analysis/{scripts,data,output,figures}) ---
# Locate the Analysis/ directory regardless of getwd() (script may be run from the
# repo root, Analysis/, or Analysis/scripts/), then build paths into its subfolders.
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
ANALYSIS_DIR <- .find_analysis_dir()                 # = .../atlantic_birds/Analysis
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
source(file.path(ANALYSIS_DIR, "scripts", "_sampling_config.R"))  # SAMPLING settings
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
# Shared phylogenetic engine (glmmTMB propto; needs the path helpers above):
# phylo_species_name(), phylo_A_list(), run_phylo_trees(), pool_rubin_df(), ...
if (ENGINE == "glmmTMB") source(file.path(ANALYSIS_DIR, "scripts", "_phylo_engine.R"))

# ---------------------------------------------------------------------------
# 1. Shared-record subset
# ---------------------------------------------------------------------------
load(derived_path("passer90.rda"))         # cleaned data (rebuild_passer90_live.R)
stopifnot(all(passer90$Status == "live"), all(passer90$known_sex))
SD_YR  <- as.numeric(attr(passer90$scaled_yr, "scaled:scale"))    # 5.0199
DEC    <- 10 / SD_YR                                              # SD-years per decade
YEARS_RECORD <- 23                                                # 1995 -> 2018, as effect_scale.rds

dat <- passer90 %>%
  filter(!is.na(conc.wing.length), !is.na(ln_body_mass)) %>%
  mutate(scaled_yr  = as.numeric(scaled_yr), scaled_lat = as.numeric(scaled_lat),
         Sex  = factor(Sex, levels = c("Female", "Male")),
         src  = factor(Main_researcher), site = factor(Municipality),
         season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")),
         wing   = conc.wing.length,        # short names -> clean brms parameter names
         lnmass = ln_body_mass) %>%
  as.data.frame()
stopifnot(!anyNA(dat$src), !anyNA(dat$site), !anyNA(dat$scaled_lat))
MEAN_WING <- mean(dat$wing)
cat(sprintf("Shared wing + mass subset: %d records, %d species, %d contributors, %d municipalities; mean wing %.2f mm, mean mass %.2f g\n",
            nrow(dat), n_distinct(dat$spp), nlevels(dat$src), nlevels(dat$site), MEAN_WING, mean(exp(dat$lnmass))))
stopifnot(nrow(dat) == 7577L)               # REVISION_PLAN.md / REVISION_NOTES_P0.md

# Conversions shared by all modes ------------------------------------------
pct_wing_decade <- function(b) 100 * b * DEC / MEAN_WING            # linear model on mm
pct_mass_decade <- function(b) 100 * (exp(b * DEC) - 1)             # model on ln g
# departure from isometry: observed mass % change minus the (1+w)^3-1 expectation
iso_gap_decade  <- function(bw, bm) pct_mass_decade(bm) - 100 * ((1 + bw * DEC / MEAN_WING)^3 - 1)
# log-scale isometry contrast per SD-year: d ln(mass) - 3 d ln(wing) (0 under isometry)
iso_contrast    <- function(bw, bm) bm - 3 * bw / MEAN_WING

# ===========================================================================
# MODE A: --fast-lme4
# ===========================================================================
if (MODE_FAST) {
  suppressMessages({ library(lme4); library(glmmTMB) })
  t0 <- Sys.time()
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  tidy_fixef <- function(fit, model, response) {
    cf <- summary(fit)$coefficients
    data.frame(model = model, response = response, term = rownames(cf), estimate = cf[, 1], se = cf[, 2],
               t = cf[, 3], lower = cf[, 1] - 1.96 * cf[, 2], upper = cf[, 1] + 1.96 * cf[, 2],
               n = nobs(fit), singular = isSingular(fit), row.names = NULL, stringsAsFactors = FALSE)
  }
  dcc <- dat %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season))
  # within-species-centred log mass for the relative-wing sensitivity
  dat <- dat %>% group_by(spp) %>% mutate(lnmass_c = lnmass - mean(lnmass)) %>% ungroup() %>% as.data.frame()
  dcc <- dcc %>% group_by(spp) %>% mutate(lnmass_c = lnmass - mean(lnmass)) %>% ungroup() %>% as.data.frame()

  rhs0 <- "Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)"
  rhs1 <- paste(rhs0, "+ (1 | src) + (1 | site)")
  rhs2 <- paste(rhs1, "+ scaled_lon + scaled_alt + season")
  specs <- list(
    list(id = "M0_baseline",   rhs = rhs0, data = dat),
    list(id = "M1_src_site",   rhs = rhs1, data = dat),
    list(id = "M1cc_src_site", rhs = rhs1, data = dcc),
    list(id = "M2_geo_season", rhs = rhs2, data = dcc))
  uni <- list()
  for (s in specs) {
    fw <- lmer(as.formula(paste("wing ~", s$rhs)),   data = s$data, REML = TRUE, control = ctrl)
    fm <- lmer(as.formula(paste("lnmass ~", s$rhs)), data = s$data, REML = TRUE, control = ctrl)
    fr <- lmer(as.formula(paste("wing ~ lnmass +", s$rhs)),   data = s$data, REML = TRUE, control = ctrl)
    frc <- lmer(as.formula(paste("wing ~ lnmass_c +", s$rhs)), data = s$data, REML = TRUE, control = ctrl)
    uni[[s$id]] <- rbind(tidy_fixef(fw, s$id, "wing"), tidy_fixef(fm, s$id, "lnmass"),
                         tidy_fixef(fr, s$id, "wing | lnmass"), tidy_fixef(frc, s$id, "wing | lnmass_c (species-centred)"))
    cat(sprintf("[%s] wing b_yr = %.3f (t %.2f); lnmass b_yr = %.4f (t %.2f); wing|lnmass b_yr = %.3f (t %.2f), b_lnmass = %.2f\n",
                s$id, fixef(fw)["scaled_yr"], coef(summary(fw))["scaled_yr", 3],
                fixef(fm)["scaled_yr"], coef(summary(fm))["scaled_yr", 3],
                fixef(fr)["scaled_yr"], coef(summary(fr))["scaled_yr", 3], fixef(fr)["lnmass"]))
  }
  uni <- do.call(rbind, uni); rownames(uni) <- NULL

  # Derived contrast from the two univariate fits, ASSUMING independence of the
  # two estimates (conservative if the residual correlation is positive, because
  # Var(a - b) = Va + Vb - 2 Cov and Cov > 0 is then ignored). The stacked
  # glmmTMB model below provides the covariance-correct version.
  contrast_indep <- function(model) {
    w <- uni[uni$model == model & uni$response == "wing" & uni$term == "scaled_yr", ]
    m <- uni[uni$model == model & uni$response == "lnmass" & uni$term == "scaled_yr", ]
    pw <- pct_wing_decade(w$estimate); pm <- pct_mass_decade(m$estimate)
    se_pw <- 100 * w$se * DEC / MEAN_WING
    se_pm <- 100 * exp(m$estimate * DEC) * m$se * DEC          # delta method
    dif <- pw - pm; se_dif <- sqrt(se_pw^2 + se_pm^2)
    data.frame(model = model, n = w$n,
               wing_pct_decade = pw, wing_pct_decade_lower = pw - 1.96 * se_pw, wing_pct_decade_upper = pw + 1.96 * se_pw,
               mass_pct_decade = pm, mass_pct_decade_lower = pm - 1.96 * se_pm, mass_pct_decade_upper = pm + 1.96 * se_pm,
               diff_pct_decade = dif, diff_lower = dif - 1.96 * se_dif, diff_upper = dif + 1.96 * se_dif,
               isometric_mass_pct_decade_expected = 100 * ((1 + pw / 100)^3 - 1),
               iso_gap_pct_decade = iso_gap_decade(w$estimate, m$estimate),
               stringsAsFactors = FALSE)
  }
  contrasts_indep <- do.call(rbind, lapply(unique(uni$model), contrast_indep))

  # --- stacked joint model (glmmTMB): both responses on the 100*ln scale ------
  # Rows = record x trait; trait-specific fixed effects; independent species
  # intercepts and year slopes per trait; independent contributor and site
  # intercepts per trait; a SHARED record-level intercept (1 | rec) plus a
  # trait-specific residual variance (dispformula = ~ trait). The shared
  # intercept induces a positive within-record correlation between the two
  # responses (the analogue of set_rescor(TRUE), constrained to one common
  # component), so the two year slopes are estimated jointly and their
  # difference has a covariance-correct SE. Tried and rejected (2026-09-09; see
  # REVISION_NOTES_P2.md): dispformula = ~0 with us(0 + trait | rec) — the
  # textbook multivariate trick — did not converge here (non-PD Hessian) under
  # REML or ML, on raw or unit-SD responses. Slopes are in % per SD-year.
  stack_fit <- function(data, rhs_extra, id) {
    keep <- c("Sex", "scaled_yr", "scaled_lat", "spp", "src", "site", "scaled_lon", "scaled_alt", "season")
    long <- rbind(
      data.frame(rec = data$ID_ABT, trait = "wing",   y = 100 * log(data$wing), data[, keep]),
      data.frame(rec = data$ID_ABT, trait = "lnmass", y = 100 * data$lnmass,    data[, keep]))
    long$trait <- factor(long$trait, levels = c("wing", "lnmass")); long$rec <- factor(long$rec)
    f <- as.formula(paste0("y ~ 0 + trait + trait:(Sex + scaled_yr + scaled_lat", rhs_extra, ") + ",
                           "diag(0 + trait + trait:scaled_yr | spp) + diag(0 + trait | src) + diag(0 + trait | site) + (1 | rec)"))
    t1 <- Sys.time()
    fit <- tryCatch(glmmTMB(f, data = long, dispformula = ~ trait, REML = TRUE), error = function(e) e)
    if (inherits(fit, "error")) return(list(id = id, error = conditionMessage(fit)))
    cf <- summary(fit)$coefficients$cond
    V  <- vcov(fit)$cond
    iw <- grep("traitwing:scaled_yr", rownames(cf)); im <- grep("traitlnmass:scaled_yr", rownames(cf))
    bw <- cf[iw, 1]; bm <- cf[im, 1]
    L  <- rep(0, nrow(cf)); L[iw] <- 1; L[im] <- -1                    # difference wing - mass (% per SD-yr)
    d_est <- sum(L * cf[, 1]); d_se <- sqrt(as.numeric(t(L) %*% V %*% L))
    L3 <- rep(0, nrow(cf)); L3[im] <- 1; L3[iw] <- -3                  # isometry contrast: mass - 3*wing (0 under isometry)
    i_est <- sum(L3 * cf[, 1]); i_se <- sqrt(as.numeric(t(L3) %*% V %*% L3))
    # implied within-record correlation: s2_rec / sqrt((s2_rec + s2_w)(s2_rec + s2_m))
    s2_rec <- as.numeric(VarCorr(fit)$cond$rec[1, 1])
    disp   <- exp(fixef(fit)$disp)                                    # dispersion model: log(sigma^2) ~ trait (0 + trait -> per-trait)
    sig2   <- if (length(disp) == 2 && all(grepl("trait", names(disp)))) disp else c(disp[1], disp[1] * disp[2])
    names(sig2) <- c("wing", "lnmass")
    r_rec  <- s2_rec / sqrt((s2_rec + sig2[["wing"]]) * (s2_rec + sig2[["lnmass"]]))
    list(id = id, n_records = nlevels(long$rec), converged = fit$fit$convergence == 0, pdHess = fit$sdr$pdHess,
         elapsed_min = as.numeric(difftime(Sys.time(), t1, units = "mins")),
         wing_pct_per_sdyr = bw, wing_se = cf[iw, 2], mass_pct_per_sdyr = bm, mass_se = cf[im, 2],
         cov_wing_mass_slopes = V[iw, im], cor_wing_mass_slopes = V[iw, im] / sqrt(V[iw, iw] * V[im, im]),
         wing_pct_decade = bw * DEC, wing_pct_decade_lower = (bw - 1.96 * cf[iw, 2]) * DEC, wing_pct_decade_upper = (bw + 1.96 * cf[iw, 2]) * DEC,
         mass_pct_decade = bm * DEC, mass_pct_decade_lower = (bm - 1.96 * cf[im, 2]) * DEC, mass_pct_decade_upper = (bm + 1.96 * cf[im, 2]) * DEC,
         diff_pct_decade = d_est * DEC, diff_lower = (d_est - 1.96 * d_se) * DEC, diff_upper = (d_est + 1.96 * d_se) * DEC, diff_z = d_est / d_se,
         iso_contrast_pct_decade = i_est * DEC, iso_lower = (i_est - 1.96 * i_se) * DEC, iso_upper = (i_est + 1.96 * i_se) * DEC, iso_z = i_est / i_se,
         record_level_var = s2_rec, residual_var = sig2, implied_within_record_cor = r_rec,
         fixef = data.frame(term = rownames(cf), estimate = cf[, 1], se = cf[, 2], z = cf[, 3], row.names = NULL))
  }
  cat("Fitting stacked glmmTMB joint models (this is the slow part of --fast-lme4) ...\n")
  stacked <- list(
    M1_src_site   = stack_fit(dat, "", "M1_src_site"),
    M2_geo_season = stack_fit(dcc, " + scaled_lon + scaled_alt + season", "M2_geo_season"))
  for (s in stacked) {
    if (!is.null(s$error)) { cat(sprintf("[stacked %s] ERROR: %s\n", s$id, s$error)); next }
    cat(sprintf("[stacked %s] wing %.2f %%/decade [%.2f, %.2f]; mass %.2f %%/decade [%.2f, %.2f]; diff (wing - mass) %.2f [%.2f, %.2f] z = %.2f; isometry contrast (mass - 3*wing) %.2f [%.2f, %.2f] z = %.2f; implied within-record r = %.3f; converged+pdHess %s (%.1f min)\n",
                s$id, s$wing_pct_decade, s$wing_pct_decade_lower, s$wing_pct_decade_upper,
                s$mass_pct_decade, s$mass_pct_decade_lower, s$mass_pct_decade_upper,
                s$diff_pct_decade, s$diff_lower, s$diff_upper, s$diff_z,
                s$iso_contrast_pct_decade, s$iso_lower, s$iso_upper, s$iso_z,
                s$implied_within_record_cor, s$converged && s$pdHess, s$elapsed_min))
  }

  fast <- list(
    generated = as.character(Sys.time()), mode = "fast-lme4",
    n_records = nrow(dat), n_species = n_distinct(dat$spp), n_src = nlevels(dat$src), n_site = nlevels(dat$site),
    n_complete_case = nrow(dcc), mean_wing_mm = MEAN_WING, mean_mass_g = mean(exp(dat$lnmass)),
    sd_year = SD_YR, sd_years_per_decade = DEC,
    univariate = uni, contrasts_independent = contrasts_indep, stacked_glmmTMB = stacked,
    estimation = "lme4::lmer REML (bobyqa), Wald 95% CI; glmmTMB REML with dispformula = ~ trait and a shared (1 | record) intercept for the joint model; no phylogenetic term",
    r_version = R.version.string, lme4_version = as.character(packageVersion("lme4")),
    glmmTMB_version = as.character(packageVersion("glmmTMB")),
    elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(fast, out_path("bivariate_fast_lme4.rds"))
  cat("\n=== univariate year slopes on the shared subset ===\n")
  print(uni %>% filter(term %in% c("scaled_yr", "lnmass", "lnmass_c")) %>%
          transmute(model, response, term, estimate = signif(estimate, 4), se = signif(se, 3), t = round(t, 2),
                    lower = signif(lower, 4), upper = signif(upper, 4), n, singular) %>% as.data.frame(), row.names = FALSE)
  cat("\n=== slope contrasts, independence assumption ===\n"); print(contrasts_indep, digits = 3)
  cat("\nWrote", out_path("bivariate_fast_lme4.rds"), sprintf("(%.1f min)\n", fast$elapsed_min))
  quit(save = "no", status = 0)
}

# ===========================================================================
# MODE C: --engine glmmTMB (default): phylogenetic tier across N_TREES trees
# ===========================================================================
# Univariate phylogenetic models on the shared subset, the delta-method slope
# contrast, and a joint long-format model with trait-specific phylogenetic
# effects. Every fit goes through _phylo_engine.R (fit_phylo_glmmtmb appends
# propto(0 + species_name | g, A) to the formula; run_phylo_trees loops over the
# trees and pools with Rubin's rules). The joint model is the one exception: it
# needs TWO propto terms (one per trait), which the engine does not build, so it
# is fitted here with glmmTMB directly on the engine's matrices and pooled with
# the engine's pool_rubin_df().
if (ENGINE == "glmmTMB") {
  suppressMessages({ library(glmmTMB) })
  t0 <- Sys.time()
  lme4_tier <- if (file.exists(out_path("bivariate_fast_lme4.rds"))) readRDS(out_path("bivariate_fast_lme4.rds")) else NULL
  if (is.null(lme4_tier)) message("bivariate_fast_lme4.rds not found: the lme4 comparison table will be empty (run --fast-lme4 first).")

  # Species names -> tree tips (eBird synonyms). The non-phylogenetic species
  # term (spp) and the phylogenetic term (species_name) use the same tip name,
  # exactly as in the published brms model. Herpsilochmus_sellowi is absent from
  # the tree (engine note); it has no wing + mass record, the filter is a guard.
  dat$species_name <- phylo_species_name(dat$Binomial)
  dat <- dat[dat$species_name != "Herpsilochmus_sellowi", , drop = FALSE]
  dat$spp <- dat$species_name
  dat <- dat %>% group_by(spp) %>% mutate(lnmass_c = lnmass - mean(lnmass)) %>% ungroup() %>% as.data.frame()
  dcc <- dat %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season)) %>% as.data.frame()
  species_all <- sort(unique(dat$species_name))
  A_list <- phylo_A_list(species_all, n_trees = N_TREES)   # the 50 published trees (cached), subset to these species
  N_TREES <- length(A_list)
  # each fit uses only the species present in its data frame (complete-case subsets may lose species)
  A_for <- function(data) { s <- sort(unique(as.character(data$species_name))); lapply(A_list, function(A) A[s, s, drop = FALSE]) }
  cat(sprintf("glmmTMB phylogenetic tier: %d trees, %d species, %d records (complete-case %d)\n",
              N_TREES, length(species_all), nrow(dat), nrow(dcc)))

  # ---- convergence repair around the engine loop ------------------------------
  # With BOTH a non-phylogenetic species intercept (1 | spp) and the phylogenetic
  # term (the published structure), the REML surface has two modes: the correct
  # one (phylogenetic SD ~ 23 mm, iid species SD ~ 0) and a worse local optimum
  # (phylogenetic SD ~ 0.1, iid SD ~ 16, REML objective +24 units, Hessian not
  # positive-definite) that glmmTMB's default all-zero start reaches on some trees
  # (4 of the first 10 for the M0 wing model on this subset, 2026-09-09 probe;
  # the year slope differs by 0.002 between modes). run_phylo_trees() is run
  # unchanged; afterwards every tree without a positive-definite Hessian is
  # refitted from the free variance parameters of a converged tree of the SAME
  # model (or a generic seed: small iid species SD, large phylogenetic SD),
  # placed into that tree's own full theta vector (the propto entries that encode
  # A are tree-specific and must not be shared). The refit replaces the original
  # only if it converged AND its objective is not worse. Both counts are stored.
  # (Candidate addition to _phylo_engine.R; kept here because that file is shared.)
  # A usable tree = positive-definite Hessian AND finite fixed-effect SEs. Trees still
  # unusable after two seeded restarts (seed from a converged tree, then a generic
  # seed) are EXCLUDED from the Rubin pooling — their SEs are undefined, and pooling
  # them would return NA — and listed in $excluded_trees; n_converged counts the
  # trees actually pooled. The point estimates of excluded trees stay in $per_tree_fixed.
  run_phylo_robust <- function(formula, data, A_list, ...) {
    r <- run_phylo_trees(formula, data, A_list, keep_fits = TRUE, verbose = FALSE, ...)
    usable_tree <- function(i) isTRUE(r$varcomp$converged[r$varcomp$tree == i]) &&
      all(is.finite(r$per_tree_fixed$se[r$per_tree_fixed$tree == i]))
    r$varcomp$usable <- vapply(r$varcomp$tree, usable_tree, logical(1))
    r$n_converged_first_pass <- sum(r$varcomp$usable)
    bad <- which(!r$varcomp$usable); r$retried_trees <- bad; r$repaired_trees <- integer(0)
    sp_sd <- sd(tapply(data[[all.vars(formula)[1]]], as.character(data$species_name), mean), na.rm = TRUE)
    for (i in bad) {
      good <- which(r$varcomp$usable)
      s <- fit_phylo_glmmtmb(formula, data, A_list[[i]], doFit = FALSE, ...)
      th0 <- s$parameters$theta; free <- which(!is.na(s$mapArg$theta))
      seeds <- list()
      if (length(good)) { seed <- r$fits[[good[1]]]$fit$par; th <- th0; th[free] <- seed[names(seed) == "theta"]; seeds$converged_tree <- th }
      th <- th0; th[free] <- 0; th[free[1]] <- log(0.01 * sp_sd); th[free[length(free)]] <- log(sp_sd); seeds$generic <- th
      old_obj <- r$fits[[i]]$fit$objective
      for (th in seeds) {
        fit2 <- tryCatch(fit_phylo_glmmtmb(formula, data, A_list[[i]], start = list(theta = th), ...), error = function(e) NULL)
        if (is.null(fit2)) next
        td <- tidy_phylo_fit(fit2, tree = i)
        ok2 <- isTRUE(fit2$sdr$pdHess) && all(is.finite(td$fixed$se))
        # tolerance 0.01: a same-mode refit differs from the original by optimizer round-off
        # (~1e-4), the two modes by ~24 objective units
        if (ok2 && (is.na(old_obj) || fit2$fit$objective <= old_obj + 0.01)) {
          r$fits[[i]] <- fit2
          r$per_tree_fixed <- rbind(r$per_tree_fixed[r$per_tree_fixed$tree != i, ], td$fixed)
          r$varcomp[r$varcomp$tree == i, names(td$varcomp)] <- td$varcomp
          r$varcomp$usable[r$varcomp$tree == i] <- TRUE
          r$repaired_trees <- c(r$repaired_trees, i)
          break
        }
      }
    }
    r$per_tree_fixed <- r$per_tree_fixed[order(r$per_tree_fixed$tree), ]
    r$varcomp$objective <- sapply(r$fits, function(f) f$fit$objective)
    r$varcomp$retried <- r$varcomp$tree %in% bad
    r$excluded_trees <- r$varcomp$tree[!r$varcomp$usable]
    keep <- r$per_tree_fixed$tree %in% r$varcomp$tree[r$varcomp$usable]
    r$pooled <- pool_rubin_df(r$per_tree_fixed[keep, ])
    r$varcomp_summary <- r$varcomp %>% filter(usable) %>%
      summarise(across(where(is.numeric) & !tree & !objective, list(mean = ~mean(.x, na.rm = TRUE), lo = ~quantile(.x, .025, na.rm = TRUE), hi = ~quantile(.x, .975, na.rm = TRUE))))
    r$n_converged <- sum(r$varcomp$usable)
    r$pooling <- "Rubin's rules over usable trees only (positive-definite Hessian, finite SEs); see excluded_trees"
    r$fits <- NULL
    r
  }

  # ---- univariate phylogenetic models -----------------------------------------
  rhs0 <- "Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)"
  rhs1 <- paste(rhs0, "+ (1 | src) + (1 | site)")
  rhs2 <- paste(rhs1, "+ scaled_lon + scaled_alt + season")
  specs <- list(
    list(id = "M0_baseline",   rhs = rhs0, data = dat),
    list(id = "M1_src_site",   rhs = rhs1, data = dat),
    list(id = "M1cc_src_site", rhs = rhs1, data = dcc),
    list(id = "M2_geo_season", rhs = rhs2, data = dcc))
  responses <- c(wing = "wing ~", lnmass = "lnmass ~", `wing | lnmass` = "wing ~ lnmass +")
  uni_pooled <- list(); uni_varcomp <- list(); uni_pertree <- list(); timing <- list()
  for (s in specs) for (rn in names(responses)) {
    f <- as.formula(paste(responses[[rn]], s$rhs))
    t_m <- Sys.time()
    r <- run_phylo_robust(f, s$data, A_for(s$data))
    secs_m <- as.numeric(Sys.time() - t_m, units = "secs")          # includes the repair refits
    n_conv <- r$n_converged
    key <- paste(s$id, rn, sep = " :: ")
    uni_pooled[[key]]  <- cbind(model = s$id, response = rn, as.data.frame(r$pooled), n = nrow(s$data),
                                n_species = length(unique(s$data$species_name)), n_trees = r$n_trees, n_converged = n_conv,
                                n_converged_first_pass = r$n_converged_first_pass, stringsAsFactors = FALSE)
    uni_varcomp[[key]] <- cbind(model = s$id, response = rn, as.data.frame(r$varcomp_summary), stringsAsFactors = FALSE)
    uni_pertree[[key]] <- list(fixed = r$per_tree_fixed, varcomp = r$varcomp)
    timing[[key]] <- data.frame(model = s$id, response = rn, n = nrow(s$data), n_trees = r$n_trees, n_converged = n_conv,
                                n_converged_first_pass = r$n_converged_first_pass, n_retried = length(r$retried_trees),
                                n_repaired = length(r$repaired_trees), secs = secs_m, secs_per_tree = secs_m / r$n_trees, stringsAsFactors = FALSE)
    yr <- r$pooled %>% filter(component == "cond", par == "scaled_yr")
    cat(sprintf("[%s | %s] b_yr = %.4f [%.4f, %.4f] z = %.2f | phylo SD %.3f, phylo prop %.3f | %d/%d trees pdHess (first pass %d, repaired %d) | %.0f s\n",
                s$id, rn, yr$estimate, yr$lower, yr$upper, yr$z, r$varcomp_summary$sd_phylo_mean,
                r$varcomp_summary$phylo_prop_mean, n_conv, r$n_trees, r$n_converged_first_pass, length(r$repaired_trees), secs_m))
    if (rn == "wing | lnmass") {
      bm <- r$pooled %>% filter(component == "cond", par == "lnmass")
      cat(sprintf("           allometric b(ln mass) = %.3f [%.3f, %.3f]\n", bm$estimate, bm$lower, bm$upper))
    }
  }
  uni_pooled  <- as.data.frame(bind_rows(uni_pooled))    # bind_rows: M0 has no src / site variance columns
  uni_varcomp <- as.data.frame(bind_rows(uni_varcomp))
  timing      <- as.data.frame(bind_rows(timing))

  # ---- slope difference from the univariate fits, INDEPENDENCE ASSUMED --------
  # The wing and mass models are fitted separately, so Cov(b_wing, b_mass) is not
  # estimated; the delta-method SE below sets it to zero. That is conservative
  # for the difference if the two slope estimates are positively correlated (the
  # joint model below estimates that correlation; in the lme4 tier it was 0.04).
  # Wing %/decade is linear in b (100 * b * DEC / mean wing); mass %/decade is
  # 100 * (exp(b * DEC) - 1), whose delta-method SE is 100 * exp(b * DEC) * DEC * se.
  contrast_indep_phylo <- function(model) {
    w <- uni_pooled %>% filter(model == !!model, response == "wing",   component == "cond", par == "scaled_yr")
    m <- uni_pooled %>% filter(model == !!model, response == "lnmass", component == "cond", par == "scaled_yr")
    pw <- pct_wing_decade(w$estimate); se_pw <- 100 * w$se * DEC / MEAN_WING
    pm <- pct_mass_decade(m$estimate); se_pm <- 100 * exp(m$estimate * DEC) * m$se * DEC
    dif <- pw - pm; se_dif <- sqrt(se_pw^2 + se_pm^2)
    iso <- iso_contrast(w$estimate, m$estimate)                     # d ln mass - 3 d ln wing, per SD-yr (0 under isometry)
    se_iso <- sqrt(m$se^2 + 9 * w$se^2 / MEAN_WING^2)
    data.frame(model = model, n = w$n, n_trees = w$n_trees, assumption = "independent wing and mass estimates (Cov = 0)",
               wing_pct_decade = pw, wing_pct_decade_lower = pw - 1.96 * se_pw, wing_pct_decade_upper = pw + 1.96 * se_pw,
               mass_pct_decade = pm, mass_pct_decade_lower = pm - 1.96 * se_pm, mass_pct_decade_upper = pm + 1.96 * se_pm,
               diff_pct_decade = dif, diff_se = se_dif, diff_lower = dif - 1.96 * se_dif, diff_upper = dif + 1.96 * se_dif, diff_z = dif / se_dif,
               iso_contrast_pct_decade = 100 * iso * DEC, iso_lower = 100 * (iso - 1.96 * se_iso) * DEC,
               iso_upper = 100 * (iso + 1.96 * se_iso) * DEC, iso_z = iso / se_iso,
               isometric_mass_pct_decade_expected = 100 * ((1 + pw / 100)^3 - 1),
               iso_gap_pct_decade = iso_gap_decade(w$estimate, m$estimate),
               stringsAsFactors = FALSE)
  }
  contrasts_indep <- do.call(rbind, lapply(unique(uni_pooled$model), contrast_indep_phylo))

  # ---- joint long-format model with trait-specific phylogenetic effects -------
  # Rows = record x trait, both responses on the 100 * ln scale (slopes in % per
  # SD-year). Trait-specific fixed effects; trait-specific species year slopes
  # diag(0 + trait:scaled_yr | spp); trait-specific contributor and site
  # intercepts; a shared record intercept (1 | rec) inducing the within-record
  # correlation; trait-specific residual variance (dispformula = ~ trait); and
  # TWO phylogenetic terms, propto(0 + species_name:is_w | g, A_w) for wing rows
  # and propto(0 + species_name:is_m | g, A_m) for mass rows (is_w / is_m are 0/1
  # indicators, so each term acts on one trait only and has its own variance;
  # propto matches the matrix by the design-column names, hence the renamed
  # dimnames). The non-phylogenetic trait-specific species INTERCEPTS are
  # omitted: with the phylogenetic term present they are estimated at ~0 in
  # every univariate fit (validation and this script), and keeping them gave a
  # REML fit whose Hessian was positive-definite but whose vcov was not
  # invertible (NA standard errors); dropping them gives finite SEs (probe of
  # 2026-09-09, 1 tree: identical slopes to 4 decimals). The species year slopes
  # are kept. The slope difference (wing - mass) and the isometry contrast
  # (mass - 3 * wing) use the estimated covariance of the two year slopes.
  make_long <- function(data) {
    keep <- c("Sex", "scaled_yr", "scaled_lat", "spp", "src", "site", "scaled_lon", "scaled_alt", "season", "species_name")
    long <- rbind(
      data.frame(rec = data$ID_ABT, trait = "wing",   y = 100 * log(data$wing), data[, keep]),
      data.frame(rec = data$ID_ABT, trait = "lnmass", y = 100 * data$lnmass,    data[, keep]))
    long$trait <- factor(long$trait, levels = c("wing", "lnmass")); long$rec <- factor(long$rec)
    long$is_w <- as.numeric(long$trait == "wing"); long$is_m <- 1 - long$is_w
    long$g <- factor(1)
    long
  }
  # seed_theta: free variance parameters of a converged tree of the same model, used to
  # restart a tree whose default-start fit has no positive-definite Hessian (same
  # repair logic as run_phylo_robust; the propto entries of theta encode A and are
  # taken from THIS tree's own doFit = FALSE structure).
  joint_fit_one <- function(long, A, rhs_extra, tree, seed_theta = NULL) {
    A <- A[levels(long$species_name), levels(long$species_name)]
    A_w <- A; dimnames(A_w) <- rep(list(paste0("species_name", rownames(A), ":is_w")), 2)
    A_m <- A; dimnames(A_m) <- rep(list(paste0("species_name", rownames(A), ":is_m")), 2)
    f <- as.formula(paste0("y ~ 0 + trait + trait:(Sex + scaled_yr + scaled_lat", rhs_extra, ") + ",
                           "diag(0 + trait:scaled_yr | spp) + diag(0 + trait | src) + diag(0 + trait | site) + (1 | rec) + ",
                           "propto(0 + species_name:is_w | g, A_w) + propto(0 + species_name:is_m | g, A_m)"))
    environment(f) <- list2env(list(A_w = A_w, A_m = A_m), parent = environment())
    t1 <- Sys.time()
    start <- NULL
    if (!is.null(seed_theta)) {
      s <- glmmTMB(f, data = long, dispformula = ~ trait, REML = TRUE, doFit = FALSE)
      th <- s$parameters$theta; free <- which(!is.na(s$mapArg$theta))
      if (length(free) == length(seed_theta)) { th[free] <- seed_theta; start <- list(theta = th) }
    }
    fit <- tryCatch(glmmTMB(f, data = long, dispformula = ~ trait, REML = TRUE, start = start), error = function(e) e)
    if (inherits(fit, "error")) return(list(error = conditionMessage(fit), tree = tree))
    cf <- summary(fit)$coefficients$cond; V <- vcov(fit)$cond
    ok <- isTRUE(fit$sdr$pdHess) && all(is.finite(cf[, 2])) && all(is.finite(V))
    theta_free <- fit$fit$par[names(fit$fit$par) == "theta"]
    iw <- grep("^traitwing:scaled_yr$", rownames(cf)); im <- grep("^traitlnmass:scaled_yr$", rownames(cf))
    bw <- cf[iw, 1]; bm <- cf[im, 1]
    d_se <- sqrt(V[iw, iw] + V[im, im] - 2 * V[iw, im])                 # wing - mass
    i_se <- sqrt(V[im, im] + 9 * V[iw, iw] - 6 * V[iw, im])             # mass - 3 * wing
    vc <- VarCorr(fit)$cond
    sd_of <- function(nm) { x <- vc[[nm]]; if (is.null(x)) NA_real_ else unname(attr(x, "stddev")[1]) }
    s2_rec <- sd_of("rec")^2
    disp <- exp(fixef(fit)$disp)                                        # log(sigma^2) ~ trait: (Intercept), traitlnmass
    sig2 <- c(wing = unname(disp[1]), lnmass = unname(disp[1] * disp[2]))
    list(tree = tree, pdHess = isTRUE(fit$sdr$pdHess), ok = ok, objective = fit$fit$objective, theta_free = theta_free,
         convergence = fit$fit$convergence, secs = as.numeric(Sys.time() - t1, units = "secs"), seeded = !is.null(start),
         fixed = data.frame(tree = tree, component = "cond", par = rownames(cf), estimate = cf[, 1], se = cf[, 2], row.names = NULL),
         derived = data.frame(tree = tree, component = "joint",
                              par = c("wing_pct_per_sdyr", "mass_pct_per_sdyr", "diff_pct_per_sdyr", "iso_contrast_pct_per_sdyr"),
                              estimate = c(bw, bm, bw - bm, bm - 3 * bw),
                              se = c(cf[iw, 2], cf[im, 2], d_se, i_se), row.names = NULL),
         varcomp = data.frame(tree = tree, cov_slopes = V[iw, im], cor_slopes = V[iw, im] / sqrt(V[iw, iw] * V[im, im]),
                              sd_phylo_wing = sd_of("g"), sd_phylo_lnmass = sd_of("g.1"),
                              sd_spp_yr_wing = unname(attr(vc$spp, "stddev")["traitwing:scaled_yr"]),
                              sd_spp_yr_lnmass = unname(attr(vc$spp, "stddev")["traitlnmass:scaled_yr"]),
                              sd_src_wing = unname(attr(vc$src, "stddev")[1]), sd_src_lnmass = unname(attr(vc$src, "stddev")[2]),
                              sd_site_wing = unname(attr(vc$site, "stddev")[1]), sd_site_lnmass = unname(attr(vc$site, "stddev")[2]),
                              sd_rec = sqrt(s2_rec), sigma_wing = sqrt(sig2[["wing"]]), sigma_lnmass = sqrt(sig2[["lnmass"]]),
                              implied_within_record_cor = s2_rec / sqrt((s2_rec + sig2[["wing"]]) * (s2_rec + sig2[["lnmass"]])),
                              pdHess = isTRUE(fit$sdr$pdHess), ok = ok, objective = fit$fit$objective, seeded = !is.null(start)))
  }
  # Loop over trees; trees whose fit is not usable (no positive-definite Hessian or
  # non-finite SEs) are refitted from the free variance parameters of the first
  # converged tree and kept only if that refit is usable and its objective is not
  # worse (tolerance 0.01). Trees still unusable after the retry are EXCLUDED from
  # the Rubin pooling (their fixed-effect SEs are undefined) and counted.
  joint_run <- function(data, rhs_extra, id) {
    long <- make_long(data)
    sp <- sort(unique(as.character(long$species_name))); long$species_name <- factor(long$species_name, levels = sp)
    t1 <- Sys.time(); out <- vector("list", N_TREES)
    for (i in seq_len(N_TREES)) {
      out[[i]] <- joint_fit_one(long, A_list[[i]], rhs_extra, i)
      if (i %% 10 == 0 || i == N_TREES) message("[joint ", id, "] tree ", i, "/", N_TREES, " ", round(as.numeric(Sys.time() - t1, units = "secs")), " s")
    }
    usable <- function(x) is.null(x$error) && isTRUE(x$ok)
    n_first_pass <- sum(vapply(out, usable, logical(1)))
    bad <- which(!vapply(out, usable, logical(1))); repaired <- integer(0)
    good <- which(vapply(out, usable, logical(1)))
    if (length(bad) && length(good)) {
      seed <- out[[good[1]]]$theta_free
      for (i in bad) {
        r2 <- joint_fit_one(long, A_list[[i]], rhs_extra, i, seed_theta = seed)
        old_obj <- if (is.null(out[[i]]$error)) out[[i]]$objective else NA_real_
        if (usable(r2) && (is.na(old_obj) || r2$objective <= old_obj + 0.01)) { out[[i]] <- r2; repaired <- c(repaired, i) }
      }
      message("[joint ", id, "] retried ", length(bad), " tree(s) with a seeded start; repaired ", length(repaired))
    }
    keep <- vapply(out, usable, logical(1))
    errs <- unlist(lapply(out, function(x) x$error))
    if (!any(keep)) return(list(id = id, n_records = nlevels(long$rec), n_trees = N_TREES, n_fitted = sum(!vapply(out, function(x) !is.null(x$error), logical(1))),
                                n_converged = 0L, n_converged_first_pass = n_first_pass, errors = errs, secs = as.numeric(Sys.time() - t1, units = "secs")))
    ok <- out[keep]
    fixed <- do.call(rbind, lapply(ok, `[[`, "fixed")); derived <- do.call(rbind, lapply(ok, `[[`, "derived"))
    varcomp <- do.call(rbind, lapply(out[vapply(out, function(x) is.null(x$error), logical(1))], `[[`, "varcomp"))   # all fitted trees, incl. excluded
    pooled_derived <- pool_rubin_df(derived) %>%
      mutate(pct_decade = estimate * DEC, pct_decade_lower = lower * DEC, pct_decade_upper = upper * DEC)
    list(id = id, n_records = nlevels(long$rec), n_species = length(sp), n_trees = N_TREES,
         n_fitted = sum(vapply(out, function(x) is.null(x$error), logical(1))),
         n_converged = sum(keep), n_converged_first_pass = n_first_pass, n_retried = length(bad), n_repaired = length(repaired),
         excluded_trees = which(!keep), errors = errs, secs = as.numeric(Sys.time() - t1, units = "secs"),
         formula = paste0("y ~ 0 + trait + trait:(Sex + scaled_yr + scaled_lat", rhs_extra, ") + diag(0 + trait:scaled_yr | spp) + diag(0 + trait | src) + diag(0 + trait | site) + (1 | rec) + propto(0 + species_name:is_w | g, A_w) + propto(0 + species_name:is_m | g, A_m); dispformula = ~ trait; REML"),
         pooling = "Rubin's rules over the trees with a usable fit (positive-definite Hessian, finite SEs); excluded_trees lists the rest",
         pooled_fixed = pool_rubin_df(fixed), pooled_derived = pooled_derived,
         varcomp_summary = varcomp %>% filter(ok) %>%
           summarise(across(where(is.numeric) & !tree & !objective, list(mean = ~mean(.x, na.rm = TRUE), lo = ~quantile(.x, .025, na.rm = TRUE), hi = ~quantile(.x, .975, na.rm = TRUE)))),
         per_tree = list(fixed = fixed, derived = derived, varcomp = varcomp))
  }
  cat("Fitting the joint long-format phylogenetic model (two propto terms) ...\n")
  joint <- list(M1_src_site   = joint_run(dat, "", "M1_src_site"),
                M2_geo_season = joint_run(dcc, " + scaled_lon + scaled_alt + season", "M2_geo_season"))
  for (j in joint) {
    if (j$n_converged == 0) { cat(sprintf("[joint %s] no usable tree: %s\n", j$id, paste(unique(j$errors), collapse = " | "))); next }
    pd <- j$pooled_derived
    g <- function(p) pd[pd$par == p, ]
    cat(sprintf("[joint %s] n = %d records; %d/%d trees usable and pooled (first pass %d, repaired %d, excluded %d); %.0f s\n  wing %.2f %%/decade [%.2f, %.2f]; mass %.2f [%.2f, %.2f]; diff (wing - mass) %.2f [%.2f, %.2f] z = %.2f; isometry contrast (mass - 3 wing) %.2f [%.2f, %.2f] z = %.2f\n  slope correlation %.3f; phylo SD wing %.1f / mass %.1f (100 ln scale); within-record r %.3f\n",
                j$id, j$n_records, j$n_converged, j$n_trees, j$n_converged_first_pass, j$n_repaired, length(j$excluded_trees), j$secs,
                g("wing_pct_per_sdyr")$pct_decade, g("wing_pct_per_sdyr")$pct_decade_lower, g("wing_pct_per_sdyr")$pct_decade_upper,
                g("mass_pct_per_sdyr")$pct_decade, g("mass_pct_per_sdyr")$pct_decade_lower, g("mass_pct_per_sdyr")$pct_decade_upper,
                g("diff_pct_per_sdyr")$pct_decade, g("diff_pct_per_sdyr")$pct_decade_lower, g("diff_pct_per_sdyr")$pct_decade_upper, g("diff_pct_per_sdyr")$z,
                g("iso_contrast_pct_per_sdyr")$pct_decade, g("iso_contrast_pct_per_sdyr")$pct_decade_lower, g("iso_contrast_pct_per_sdyr")$pct_decade_upper, g("iso_contrast_pct_per_sdyr")$z,
                j$varcomp_summary$cor_slopes_mean, j$varcomp_summary$sd_phylo_wing_mean, j$varcomp_summary$sd_phylo_lnmass_mean,
                j$varcomp_summary$implied_within_record_cor_mean))
  }

  # ---- comparison with the lme4 tier on disk (same subset, same specifications) ----
  comparison <- NULL
  if (!is.null(lme4_tier)) {
    l <- lme4_tier$univariate %>% filter(term %in% c("scaled_yr", "lnmass")) %>%
      transmute(model, response, par = term, lme4_estimate = estimate, lme4_se = se, lme4_lower = lower, lme4_upper = upper, lme4_n = n)
    p <- uni_pooled %>% filter(component == "cond", par %in% c("scaled_yr", "lnmass")) %>%
      transmute(model, response, par, phylo_estimate = estimate, phylo_se = se, phylo_lower = lower, phylo_upper = upper, phylo_n = n, n_trees, n_converged)
    comparison <- inner_join(l, p, by = c("model", "response", "par")) %>%
      mutate(delta_estimate = phylo_estimate - lme4_estimate, se_ratio_phylo_over_lme4 = phylo_se / lme4_se,
             lme4_excludes_zero = lme4_lower > 0 | lme4_upper < 0, phylo_excludes_zero = phylo_lower > 0 | phylo_upper < 0,
             verdict_changes = lme4_excludes_zero != phylo_excludes_zero)
    lc <- lme4_tier$contrasts_independent %>% transmute(model, lme4_diff = diff_pct_decade, lme4_diff_lower = diff_lower, lme4_diff_upper = diff_upper)
    pc <- contrasts_indep %>% transmute(model, phylo_diff = diff_pct_decade, phylo_diff_lower = diff_lower, phylo_diff_upper = diff_upper)
    comparison_contrasts <- inner_join(lc, pc, by = "model")
    js <- lme4_tier$stacked_glmmTMB
    comparison_joint <- do.call(rbind, lapply(names(joint), function(id) {
      j <- joint[[id]]; s <- js[[id]]
      if (is.null(s) || !is.null(s$error) || j$n_converged == 0) return(NULL)
      pd <- j$pooled_derived; g <- function(p) pd[pd$par == p, ]
      data.frame(model = id,
                 lme4_joint_diff = s$diff_pct_decade, lme4_joint_lower = s$diff_lower, lme4_joint_upper = s$diff_upper,
                 phylo_joint_diff = g("diff_pct_per_sdyr")$pct_decade, phylo_joint_lower = g("diff_pct_per_sdyr")$pct_decade_lower, phylo_joint_upper = g("diff_pct_per_sdyr")$pct_decade_upper,
                 lme4_joint_iso = s$iso_contrast_pct_decade, lme4_joint_iso_lower = s$iso_lower, lme4_joint_iso_upper = s$iso_upper,
                 phylo_joint_iso = g("iso_contrast_pct_per_sdyr")$pct_decade, phylo_joint_iso_lower = g("iso_contrast_pct_per_sdyr")$pct_decade_lower, phylo_joint_iso_upper = g("iso_contrast_pct_per_sdyr")$pct_decade_upper,
                 lme4_cor_slopes = s$cor_wing_mass_slopes, phylo_cor_slopes = j$varcomp_summary$cor_slopes_mean, stringsAsFactors = FALSE)
    }))
    comparison <- list(univariate = comparison, contrasts_independent = comparison_contrasts, joint = comparison_joint,
                       lme4_generated = lme4_tier$generated)
  }

  results <- list(
    generated = as.character(Sys.time()), engine = "glmmTMB propto (Williams et al. 2025) via _phylo_engine.R; REML; Rubin's rules across trees",
    mode = sprintf("phylogenetic tier, %d trees", N_TREES), n_trees = N_TREES,
    trees = "the 50 clootl trees of the published brms analysis (phylo_A_list(); cache data/derived/phylo_A_50trees.rds built from models/brm0_multiphylo.rda)",
    n_records = nrow(dat), n_species = length(species_all), n_src = nlevels(droplevels(dat$src)), n_site = nlevels(droplevels(dat$site)),
    n_complete_case = nrow(dcc), mean_wing_mm = MEAN_WING, mean_mass_g = mean(exp(dat$lnmass)), sd_year = SD_YR, sd_years_per_decade = DEC,
    specifications = c(M0_baseline = rhs0, M1_src_site = rhs1, M1cc_src_site = paste(rhs1, "[complete cases]"), M2_geo_season = rhs2,
                       phylogenetic_term = "propto(0 + species_name | g, A) appended by fit_phylo_glmmtmb()"),
    univariate = list(pooled = uni_pooled, varcomp = uni_varcomp, per_tree = uni_pertree),
    contrasts_independent = contrasts_indep,
    joint = joint,
    comparison_lme4 = comparison,
    timing = list(univariate = timing, joint = data.frame(model = names(joint), secs = sapply(joint, `[[`, "secs"),
                                                          n_fitted = sapply(joint, `[[`, "n_fitted"),
                                                          n_converged = sapply(joint, function(j) if (is.null(j$n_converged)) 0L else j$n_converged),
                                                          n_converged_first_pass = sapply(joint, function(j) if (is.null(j$n_converged_first_pass)) 0L else j$n_converged_first_pass),
                                                          n_excluded = sapply(joint, function(j) length(j$excluded_trees)), row.names = NULL)),
    definitions = c(
      wing_pct_decade = "100 * b_wing * (10/SD_yr) / mean wing (mm, shared subset)",
      mass_pct_decade = "100 * (exp(b_lnmass * 10/SD_yr) - 1); delta-method SE 100 * exp(b*DEC) * DEC * se",
      diff_pct_decade = "wing_pct_decade - mass_pct_decade; negative = wing shortens relative to mass; contrasts_independent assumes Cov(b_wing, b_mass) = 0, joint uses the estimated covariance",
      iso_contrast = "b_lnmass - 3 * b_wing / mean wing (per SD-yr; 0 under isometry), reported x 100 x DEC as % per decade; joint: (mass - 3 * wing) on the 100 ln scale",
      intervals = "estimate +/- 1.96 * Rubin SE (within-tree Wald variance + (1 + 1/m) between-tree variance)",
      converged = "n_converged = trees with a positive-definite Hessian (glmmTMB $sdr$pdHess)"),
    glmmTMB_version = as.character(packageVersion("glmmTMB")), r_version = R.version.string,
    total_elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(results, out_path("bivariate_phylo_results.rds"))

  cat("\n=== Univariate phylogenetic year slopes (Rubin-pooled) ===\n")
  print(uni_pooled %>% filter(component == "cond", par %in% c("scaled_yr", "lnmass")) %>%
          transmute(model, response, par, estimate = signif(estimate, 4), se = signif(se, 3), lower = signif(lower, 4), upper = signif(upper, 4),
                    z = round(z, 2), n, n_converged) %>% as.data.frame(), row.names = FALSE)
  cat("\n=== Slope contrasts, independence assumption (% per decade) ===\n")
  print(contrasts_indep %>% select(model, n, wing_pct_decade, mass_pct_decade, diff_pct_decade, diff_lower, diff_upper, diff_z,
                                   iso_contrast_pct_decade, iso_lower, iso_upper, iso_z), digits = 3, row.names = FALSE)
  if (!is.null(comparison)) {
    cat("\n=== lme4 tier vs phylogenetic tier (same subset and specification) ===\n")
    print(comparison$univariate %>% transmute(model, response, par, lme4 = signif(lme4_estimate, 4), lme4_se = signif(lme4_se, 3),
                                              phylo = signif(phylo_estimate, 4), phylo_se = signif(phylo_se, 3), delta = signif(delta_estimate, 3),
                                              se_ratio = round(se_ratio_phylo_over_lme4, 3), verdict_changes) %>% as.data.frame(), row.names = FALSE)
  }
  cat(sprintf("\nWrote %s (%d trees; %.1f min total)\n", out_path("bivariate_phylo_results.rds"), N_TREES, results$total_elapsed_min))
  quit(save = "no", status = 0)
}

# ===========================================================================
# MODE B: --engine brms (full run, or --smoke) — the Bayesian cross-check (Totoro)
# ===========================================================================
suppressMessages({
  library(brms); library(ape); library(MCMCglmm); library(prepR4pcm); library(phytools)
  library(future.apply); library(posterior)
})
t_start <- Sys.time()

if (MODE_SMOKE) {
  SAMPLING$chains <- 2L; SAMPLING$iter <- 400L; SAMPLING$warmup <- 200L
  SAMPLING$workers <- 1L; SAMPLING$backend <- "rstan"; SAMPLING$cores <- 1L
  message("[smoke] 1 tree, chains = 2, iter = 400, warmup = 200, rstan — NOT A RESULT")
  if (!is.na(SUBSAMPLE_N) && SUBSAMPLE_N < nrow(dat)) {
    set.seed(20260909)
    dat <- dat[sort(sample.int(nrow(dat), SUBSAMPLE_N)), ]
    dat$src <- droplevels(dat$src); dat$site <- droplevels(dat$site)
    message(sprintf("[smoke] random SUBSAMPLE of %d records (%d species) — wiring check only", nrow(dat), n_distinct(dat$spp)))
  }
}

# clootl/eBird uses scientific names WITHOUT underscores; provide a clean column.
dat$species_name <- gsub("_", " ", dat$Binomial)

# --- ABT (2018) -> current eBird/Clements names (identical to atlantic_parallel.R) ---
ebird_synonyms <- c(
  "Antilophia galeata"       = "Chiroxiphia galeata",   # Antilophia merged into Chiroxiphia (2024)
  "Tachyphonus cristatus"    = "Loriotus cristatus",    # moved to Loriotus
  "Pyrrhocoma ruficeps"      = "Thlypopsis pyrrhocoma", # genus AND epithet changed
  "Pyriglena pernambucensis" = "Pyriglena leuconota",   # now a subspecies of leuconota
  "Tangara sayaca"           = "Thraupis sayaca",
  "Tangara cayana"           = "Stilpnia cayana",        # moved to Stilpnia
  "Tangara palmarum"         = "Thraupis palmarum",      # Palm Tanager stays in Thraupis
  "Tangara peruviana"        = "Stilpnia peruviana",     # moved to Stilpnia (eBird 2020+)
  "Dixiphia pipra"           = "Pseudopipra pipra",      # SACC 876 (2023): Dixiphia -> Pseudopipra
  "Tiaris fuliginosus"       = "Asemospiza fuliginosa"   # Tiaris split -> Asemospiza
)
hit <- dat$species_name %in% names(ebird_synonyms)
dat$species_name[hit] <- ebird_synonyms[dat$species_name[hit]]
spp_data <- unique(dat$species_name)

# --- AvesData repo for the clootl tree cloud --------------------------------
# atlantic_parallel.R downloads it when AVESDATA_PATH is unset. Here we first
# look for the copy already in data/raw (present in this repo), so a local run
# needs no network; otherwise fall back to the download exactly as before.
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
  if (dir.exists(raw_path("AvesDataLite-main"))) {
    Sys.setenv(AVESDATA_PATH = raw_path("AvesDataLite-main"))
  } else {
    clootl::get_avesdata_repo(path = raw_path())
  }
}
.check <- pr_get_tree(spp_data, source = "clootl", n_tree = 1)
if (length(.check$unmatched) > 0) {
  message("Removing ", length(.check$unmatched), " species absent from eBird taxonomy: ",
          paste(.check$unmatched, collapse = ", "))
  spp_data <- setdiff(spp_data, .check$unmatched)
}
stopifnot(length(.check$unmatched) <= 1)
got    <- pr_get_tree(spp_data, source = "clootl", n_tree = 100, cache = TRUE)
trees  <- got$tree                              # multiPhylo (100 dated trees)
clootl_version <- pr_cite_tree(got, format = "text")

# --- Reconcile data species names to the tree tips (as atlantic_parallel.R) --
rec <- reconcile_tree(x = dat, tree = trees[[1]], x_species = "species_name",
                      fuzzy = TRUE, resolve = "flag")
print(reconcile_summary(rec))
aligned <- reconcile_apply(rec, data = dat, tree = trees[[1]],
                           species_col = "species_name", drop_unresolved = TRUE)
dat <- aligned$data
dat$species_name <- gsub(" ", "_", dat$species_name)  # underscore form for brms grouping
dat$spp <- dat$species_name                            # both grouping terms use matched names
cat(sprintf("After reconciliation: %d records, %d species\n", nrow(dat), n_distinct(dat$spp)))

norm_us   <- function(x) gsub(" ", "_", x)
keep_tips <- norm_us(aligned$tree$tip.label)
trees_pruned <- lapply(trees, function(t) {
  t$tip.label <- norm_us(t$tip.label)
  ape::keep.tip(t, intersect(keep_tips, t$tip.label))
})
class(trees_pruned) <- "multiPhylo"
tree_samp <- sample(trees_pruned, N_TREES)      # seed set above; 50 in the full run

# --- helper: per-tree phylogenetic covariance matrix (verbatim, atlantic_parallel.R) ---
make_A <- function(tree) {
  tree$tip.label <- gsub(" ", "_", tree$tip.label)  # match dat$species_name format
  # clootl trees are dated but pruning + floating-point leaves them very slightly
  # non-ultrametric, which inverseA(scale = TRUE) refuses. Coerce to ultrametric
  # (minimal NNLS adjustment) first.
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A   <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

# --- priors: wing as atlantic_parallel.R, log mass as atlantic_parallel_mass.R ---
priors_bi <- c(
  prior(normal(71, 15), class = Intercept, resp = wing),
  prior(normal(3, 1.5),  class = Intercept, resp = lnmass),   # log grams ~ 3
  prior(normal(0, 10),   class = b,  resp = wing),
  prior(normal(0, 10),   class = b,  resp = lnmass),
  prior(cauchy(0, 1),    class = sd, resp = wing),
  prior(cauchy(0, 1),    class = sd, resp = lnmass),
  prior(cauchy(0, 1),    class = sigma, resp = wing),
  prior(cauchy(0, 1),    class = sigma, resp = lnmass),
  prior(lkj(2),          class = rescor))
priors_rel <- c(prior(normal(71, 15), class = Intercept),
                prior(normal(0, 10),  class = b),
                prior(cauchy(0, 1),   class = sd),
                prior(cauchy(0, 1),   class = sigma))

# Group-level terms of the two responses are independent (no |ID| syntax): each
# response has its own species intercept + slope, phylogenetic effect,
# contributor and site intercepts, matching the univariate models of Phases 1
# and 2; only the residuals are correlated (set_rescor(TRUE)).
bf_wing <- bf(wing   ~ 1 + Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) +
                (1 | gr(species_name, cov = A)) + (1 | src) + (1 | site))
bf_mass <- bf(lnmass ~ 1 + Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) +
                (1 | gr(species_name, cov = A)) + (1 | src) + (1 | site))

fit_bivariate <- function(tree) {
  A <- make_A(tree)
  brm(bf_wing + bf_mass + set_rescor(TRUE),
      data = dat, data2 = list(A = A),
      family = gaussian(), prior = priors_bi,
      iter = SAMPLING$iter, warmup = SAMPLING$warmup, chains = SAMPLING$chains, cores = SAMPLING$cores, backend = SAMPLING$backend,
      control = SAMPLING_CONTROL,
      seed = 20240101)
}
fit_relative <- function(tree) {
  A <- make_A(tree)
  brm(wing ~ 1 + Sex + scaled_yr + scaled_lat + lnmass +
        (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)) + (1 | src) + (1 | site),
      data = dat, data2 = list(A = A),
      family = gaussian(), prior = priors_rel,
      iter = SAMPLING$iter, warmup = SAMPLING$warmup, chains = SAMPLING$chains, cores = SAMPLING$cores, backend = SAMPLING$backend,
      control = SAMPLING_CONTROL,
      seed = 20240101)
}

# --- Rubin's rules pooling across trees (verbatim from atlantic_parallel.R) ---
# For each fixed effect: pooled estimate = mean of per-tree posterior means;
# total variance = within-tree var (mean of per-tree posterior variances)
# + between-tree var (variance of per-tree means) inflated by (1 + 1/m).
# (Nakagawa & de Villemereuil 2019, Syst. Biol. 68:632-641.)
pool_rubin <- function(fits, pars = c("b_Intercept","b_SexMale","b_scaled_yr","b_scaled_lat")) {
  m <- length(fits)
  draws <- lapply(fits, function(f) as_draws_df(f)[, pars, drop = FALSE])
  means <- t(sapply(draws, function(d) sapply(d, mean)))
  vars  <- t(sapply(draws, function(d) sapply(d, var)))
  qbar  <- colMeans(means)                 # pooled point estimate
  ubar  <- colMeans(vars)                  # within-imputation variance
  b     <- apply(means, 2, var)            # between-tree variance
  tot   <- ubar + (1 + 1/m) * b            # total variance
  se    <- sqrt(tot)
  data.frame(par = pars, estimate = qbar, se = se,
             lower = qbar - 1.96*se, upper = qbar + 1.96*se)
}
# Same rules applied to a list of draws data frames that already contain the
# derived quantities (so the slope difference is pooled like any parameter).
pool_rubin_draws <- function(draws_list, pars) {
  m <- length(draws_list)
  means <- t(sapply(draws_list, function(d) sapply(d[, pars, drop = FALSE], mean)))
  vars  <- t(sapply(draws_list, function(d) sapply(d[, pars, drop = FALSE], var)))
  if (m == 1) { means <- matrix(means, nrow = 1, dimnames = list(NULL, pars)); vars <- matrix(vars, nrow = 1, dimnames = list(NULL, pars)) }
  qbar <- colMeans(means); ubar <- colMeans(vars)
  b    <- if (m > 1) apply(means, 2, var) else rep(0, length(pars))
  se   <- sqrt(ubar + (1 + 1/m) * b)
  # per-tree posterior quantiles pooled by averaging (reported alongside the Wald-Rubin interval)
  q025 <- colMeans(t(sapply(draws_list, function(d) sapply(d[, pars, drop = FALSE], quantile, 0.025))))
  q975 <- colMeans(t(sapply(draws_list, function(d) sapply(d[, pars, drop = FALSE], quantile, 0.975))))
  data.frame(par = pars, estimate = qbar, se = se, lower = qbar - 1.96 * se, upper = qbar + 1.96 * se,
             mean_q025 = q025, mean_q975 = q975, row.names = NULL)
}
# With a single tree the verbatim pool_rubin() has no between-tree variance to
# estimate (var of one value = NA); route m = 1 through pool_rubin_draws(), which
# sets the between-tree term to zero, so the --trees 1 validation tier still
# returns intervals. For m > 1 the two functions give identical numbers.
pool_fixed <- function(fits, pars) {
  if (length(fits) > 1) return(pool_rubin(fits, pars))
  pool_rubin_draws(lapply(fits, function(f) as.data.frame(as_draws_df(f))[, pars, drop = FALSE]), pars)[, c("par", "estimate", "se", "lower", "upper")]
}
# Derived quantities from one bivariate fit
derive_bivariate <- function(f) {
  d <- as.data.frame(as_draws_df(f))
  bw <- d$b_wing_scaled_yr; bm <- d$b_lnmass_scaled_yr
  data.frame(
    b_wing_scaled_yr = bw, b_lnmass_scaled_yr = bm,
    b_wing_scaled_lat = d$b_wing_scaled_lat, b_lnmass_scaled_lat = d$b_lnmass_scaled_lat,
    b_wing_SexMale = d$b_wing_SexMale, b_lnmass_SexMale = d$b_lnmass_SexMale,
    rescor_wing_lnmass = d$rescor__wing__lnmass,
    wing_mm_decade = bw * DEC,
    wing_pct_decade = pct_wing_decade(bw),
    mass_pct_decade = pct_mass_decade(bm),
    diff_pct_decade = pct_wing_decade(bw) - pct_mass_decade(bm),        # wing % - mass %; < 0 = wing shrinks relative to mass
    iso_gap_pct_decade = iso_gap_decade(bw, bm),                        # observed mass % - isometric expectation
    iso_contrast_per_sdyr = iso_contrast(bw, bm),                       # d ln mass - 3 d ln wing
    wing_pct_record = 100 * bw * (YEARS_RECORD / SD_YR) / MEAN_WING,
    mass_pct_record = 100 * (exp(bm * YEARS_RECORD / SD_YR) - 1),
    isometric_mass_pct_record_expected = 100 * ((1 + bw * (YEARS_RECORD / SD_YR) / MEAN_WING)^3 - 1),
    sd_spp_wing_yr = d$sd_spp__wing_scaled_yr, sd_spp_lnmass_yr = d$sd_spp__lnmass_scaled_yr,
    sd_src_wing = d$sd_src__wing_Intercept, sd_src_lnmass = d$sd_src__lnmass_Intercept,
    sd_site_wing = d$sd_site__wing_Intercept, sd_site_lnmass = d$sd_site__lnmass_Intercept,
    sd_phylo_wing = d$sd_species_name__wing_Intercept, sd_phylo_lnmass = d$sd_species_name__lnmass_Intercept,
    sigma_wing = d$sigma_wing, sigma_lnmass = d$sigma_lnmass)
}
derive_relative <- function(f) {
  d <- as.data.frame(as_draws_df(f))
  data.frame(b_scaled_yr = d$b_scaled_yr, b_lnmass = d$b_lnmass, b_scaled_lat = d$b_scaled_lat, b_SexMale = d$b_SexMale,
             wing_mm_decade_given_mass = d$b_scaled_yr * DEC,
             wing_pct_decade_given_mass = pct_wing_decade(d$b_scaled_yr),
             sd_src = d$sd_src__Intercept, sd_site = d$sd_site__Intercept, sigma = d$sigma)
}
diag_of <- function(f) {
  s <- summary(f)
  fx <- s$fixed
  list(max_rhat = max(fx$Rhat, na.rm = TRUE), min_bulk_ess = min(fx$Bulk_ESS, na.rm = TRUE),
       min_tail_ess = min(fx$Tail_ESS, na.rm = TRUE),
       n_divergent = tryCatch(sum(nuts_params(f, pars = "divergent__")$Value), error = function(e) NA))
}

# --- fits ---------------------------------------------------------------------
plan(multisession, workers = SAMPLING$workers)
cat(sprintf("Fitting the bivariate model on %d tree(s) ...\n", N_TREES))
t1 <- Sys.time()
fits_bi <- future_lapply(tree_samp, fit_bivariate, future.seed = TRUE)
el_bi <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
cat(sprintf("Bivariate fits done in %.1f min\n", el_bi))

# Smoke-test time guard: the whole smoke run should stay within ~20 minutes.
SMOKE_BUDGET_MIN <- 20
skip_relative <- FALSE
if (MODE_SMOKE && as.numeric(difftime(Sys.time(), t_start, units = "mins")) + el_bi > SMOKE_BUDGET_MIN) {
  skip_relative <- TRUE
  message(sprintf("[smoke] elapsed %.1f min; the relative-wing model would exceed the %d-min budget and is SKIPPED",
                  as.numeric(difftime(Sys.time(), t_start, units = "mins")), SMOKE_BUDGET_MIN))
}
fits_rel <- NULL; el_rel <- NA_real_
if (!skip_relative) {
  cat(sprintf("Fitting the relative-wing (wing | log mass) model on %d tree(s) ...\n", N_TREES))
  t2 <- Sys.time()
  fits_rel <- future_lapply(tree_samp, fit_relative, future.seed = TRUE)
  el_rel <- as.numeric(difftime(Sys.time(), t2, units = "mins"))
  cat(sprintf("Relative-wing fits done in %.1f min\n", el_rel))
}
plan(sequential)

# --- pooling -------------------------------------------------------------------
draws_bi <- lapply(fits_bi, derive_bivariate)
pars_bi  <- setdiff(names(draws_bi[[1]]), character(0))
rubin_bi <- pool_rubin_draws(draws_bi, pars_bi)
rubin_bi_fixed <- pool_fixed(fits_bi, pars = c("b_wing_Intercept", "b_lnmass_Intercept", "b_wing_SexMale", "b_lnmass_SexMale",
                                               "b_wing_scaled_yr", "b_lnmass_scaled_yr", "b_wing_scaled_lat", "b_lnmass_scaled_lat",
                                               "rescor__wing__lnmass"))
# posterior probability that the difference is negative (wing shrinks relative to mass), averaged over trees
p_diff_neg   <- mean(sapply(draws_bi, function(d) mean(d$diff_pct_decade < 0)))
p_isogap_pos <- mean(sapply(draws_bi, function(d) mean(d$iso_gap_pct_decade > 0)))   # mass above isometric expectation
diag_bi <- lapply(fits_bi, diag_of)

rubin_rel <- NULL; rubin_rel_fixed <- NULL; diag_rel <- NULL
if (!is.null(fits_rel)) {
  draws_rel <- lapply(fits_rel, derive_relative)
  rubin_rel <- pool_rubin_draws(draws_rel, names(draws_rel[[1]]))
  rubin_rel_fixed <- pool_fixed(fits_rel, pars = c("b_Intercept", "b_SexMale", "b_scaled_yr", "b_scaled_lat", "b_lnmass"))
  diag_rel <- lapply(fits_rel, diag_of)
}

cat("\n=== Bivariate model: Rubin-pooled fixed effects ===\n"); print(rubin_bi_fixed, digits = 4)
cat("\n=== Bivariate model: derived quantities (pooled) ===\n")
print(rubin_bi[rubin_bi$par %in% c("wing_mm_decade", "wing_pct_decade", "mass_pct_decade", "diff_pct_decade",
                                   "iso_gap_pct_decade", "iso_contrast_per_sdyr", "rescor_wing_lnmass",
                                   "wing_pct_record", "mass_pct_record", "isometric_mass_pct_record_expected"), ], digits = 3)
cat(sprintf("P(diff < 0) = %.3f (wing shrinks relative to mass); P(iso gap > 0) = %.3f (mass above isometric expectation)\n",
            p_diff_neg, p_isogap_pos))
if (!is.null(rubin_rel_fixed)) { cat("\n=== Relative-wing model (wing | log mass): pooled fixed effects ===\n"); print(rubin_rel_fixed, digits = 4) }
cat("\nDiagnostics (per tree): max Rhat", sapply(diag_bi, `[[`, "max_rhat"), "; min bulk ESS", sapply(diag_bi, `[[`, "min_bulk_ess"),
    "; divergent", sapply(diag_bi, `[[`, "n_divergent"), "\n")

results <- list(
  generated = as.character(Sys.time()),
  mode = if (MODE_SMOKE) sprintf("SMOKE TEST — NOT A RESULT (1 tree, chains = 2, iter = 400%s)",
                                  if (!is.na(SUBSAMPLE_N)) paste0(", random subsample of ", nrow(dat), " records") else "")
         else sprintf("full run, %d trees", N_TREES),
  is_smoke = MODE_SMOKE, subsample_n = if (is.na(SUBSAMPLE_N)) NA_integer_ else nrow(dat), n_trees = N_TREES,
  sampling = SAMPLING[c("chains", "iter", "warmup", "adapt_delta", "max_treedepth", "backend")],
  n_records = nrow(dat), n_species = n_distinct(dat$spp), n_src = nlevels(droplevels(dat$src)), n_site = nlevels(droplevels(dat$site)),
  mean_wing_mm = MEAN_WING, mean_mass_g = mean(exp(dat$lnmass)), sd_year = SD_YR, sd_years_per_decade = DEC, years_record = YEARS_RECORD,
  bivariate = list(fixed = rubin_bi_fixed, derived = rubin_bi, p_diff_negative = p_diff_neg, p_isogap_positive = p_isogap_pos,
                   diagnostics = diag_bi, elapsed_min = el_bi,
                   formula = "mvbind(wing, lnmass) ~ 1 + Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)) + (1 | src) + (1 | site), set_rescor(TRUE)"),
  relative_wing = list(fixed = rubin_rel_fixed, derived = rubin_rel, diagnostics = diag_rel, elapsed_min = el_rel, skipped = skip_relative,
                       formula = "wing ~ 1 + Sex + scaled_yr + scaled_lat + lnmass + (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)) + (1 | src) + (1 | site)"),
  definitions = c(
    wing_pct_decade = "100 * b_wing_scaled_yr * (10/SD_yr) / mean wing (mm, shared subset)",
    mass_pct_decade = "100 * (exp(b_lnmass_scaled_yr * 10/SD_yr) - 1)",
    diff_pct_decade = "wing_pct_decade - mass_pct_decade; negative = wing shortens relative to mass",
    iso_gap_pct_decade = "mass_pct_decade - 100*((1 + wing fraction per decade)^3 - 1); positive = mass declines LESS than isometry predicts",
    iso_contrast_per_sdyr = "b_lnmass_scaled_yr - 3 * b_wing_scaled_yr / mean wing; 0 under isometry",
    intervals = "estimate +/- 1.96 * Rubin SE (as atlantic_parallel.R); mean_q025/mean_q975 = per-tree posterior 2.5/97.5% quantiles averaged over trees"),
  clootl_version = clootl_version, r_version = R.version.string, brms_version = as.character(packageVersion("brms")),
  total_elapsed_min = as.numeric(difftime(Sys.time(), t_start, units = "mins")))

if (MODE_SMOKE) {
  saveRDS(results, out_path("bivariate_smoke_results.rds"))
  save(fits_bi, fits_rel, rubin_bi, rubin_bi_fixed, rubin_rel, rubin_rel_fixed, clootl_version,
       file = out_path("models", "bivariate_smoke.rda"))
  cat(sprintf("\n[smoke] wrote bivariate_smoke_results.rds and models/bivariate_smoke.rda (total %.1f min). These are NOT results.\n",
              results$total_elapsed_min))
} else {
  saveRDS(results, out_path("bivariate_results.rds"))
  save(fits_bi, fits_rel, rubin_bi, rubin_bi_fixed, rubin_rel, rubin_rel_fixed, clootl_version,
       file = out_path("models", "bivariate_wing_mass_multiphylo.rda"))
  cat(sprintf("\nWrote bivariate_results.rds and models/bivariate_wing_mass_multiphylo.rda (%d trees, total %.1f min)\n",
              N_TREES, results$total_elapsed_min))
}
