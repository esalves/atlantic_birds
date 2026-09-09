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
# MODES (command-line flags; default = full brms run, 50 trees, SAMPLING settings):
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
#   --trees N     number of trees for the full brms run (default 50, as in
#                 atlantic_parallel.R). N = 1 is the "single-tree validation"
#                 tier of the plan.
#
# UNITS: scaled_yr = (Year - 2009.487) / 5.0199; per decade = beta * 10 / 5.0199.
# Wing % per decade = 100 * per_decade / mean(wing on the shared subset); mass
# % per decade = 100 * (exp(per_decade) - 1).
#
# INPUTS:  data/derived/passer90.rda; clootl / AvesData tree cloud (McTavish et
#          al. 2025) via prepR4pcm, as in atlantic_parallel.R. If AVESDATA_PATH
#          is unset and data/raw/AvesDataLite-main exists it is used (no
#          download); otherwise clootl::get_avesdata_repo(raw_path()) as before.
# OUTPUTS: output/bivariate_results.rds             (full run)
#          output/models/bivariate_wing_mass_multiphylo.rda (full run: fits,
#                                                    Rubin summaries, clootl version)
#          output/bivariate_fast_lme4.rds            (--fast-lme4)
#          output/bivariate_smoke_results.rds, output/models/bivariate_smoke.rda (--smoke)
#
# RUN:  Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --fast-lme4
#       Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --smoke [--subsample 800]
#       Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --trees 50   # Totoro
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
# MODE B: brms (full run, or --smoke)
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
