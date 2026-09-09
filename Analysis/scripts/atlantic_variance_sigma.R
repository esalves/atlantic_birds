# atlantic_variance_sigma.R
# ---------------------------------------------------------------------------
# WHAT: Phase 3 of REVISION_PLAN.md — re-analysis of the "rising wing-length
# variance" result with provenance (contributor), site, sex and season controls.
# Four parts, each written to output/variance_results.rds:
#
#   (1) PHYLOGENETIC DISTRIBUTIONAL model (the primary variance evidence, replacing
#       drmSEM): the M3 mean structure with a phylogenetic species intercept, and the
#       plan's sigma sub-model:
#         conc.wing.length ~ Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt +
#                            season + (1 + scaled_yr || spp) + <phylogenetic species
#                            intercept> + (1 | Main_researcher) + (1 | Municipality)
#         sigma ~ scaled_yr [+ scaled_tmean] + Sex + (1 | Main_researcher) + (1 | Municipality)
#       Two engines, chosen with --engine (glmmTMB is the default since 2026-09-09):
#         * glmmTMB: _phylo_engine.R fits the phylogenetic intercept with glmmTMB's
#           propto() covariance structure (Williams, McGillycuddy, Drobniak, Bolker,
#           Warton & Nakagawa 2025, bioRxiv 10.64898/2025.12.20.695312) and the sigma
#           sub-model with dispformula, REML, on the IDENTICAL 50 trees as the
#           published brms wing model (data/derived/phylo_A_50trees.rds, cached from
#           brm0_multiphylo.rda), pooled with Rubin's rules. glmmTMB reproduces the
#           published brms wing model to 2-3 decimals (glmmtmb_validation_wing.R).
#           EVERY tier of the sigma ladder in part 2 (S0-S7) is refitted with the
#           phylogenetic term across the trees (~3 s per fit), so "phylogeny vs no
#           phylogeny" is a like-for-like comparison on the same records. Results:
#           output/variance_phylo_results.rds (+ $phylo pointer/headline slot in
#           variance_results.rds). Headline models: S7 (plan specification, with
#           scaled_tmean, 7,651 records) and S3 (the --no-tmean variant, 8,282 records).
#         * brms: the original Stan path, kept as the Bayesian cross-check for Totoro:
#           bf(<mean with (1 | gr(species_name, cov = A))>, sigma ~ ...) across N
#           clootl trees, priors below, Rubin-pooled exactly as in atlantic_parallel.R;
#           output/models/variance_sigma_multiphylo.rda and variance_results.rds$brms.
#       scaled_tmean (record-year annual mean temperature at
#       the capture coordinates, from P5's climate_extraction.R) is included ONLY
#       if data/derived/passer90_climate.rds carries the ID_ABT key and covers
#       >= 90 % of the complete-case wing sample; the tmean model then runs on the
#       tmean-complete subset (7,651 of 8,282 records; the 631 missing are coastal
#       / island localities outside the WorldClim land mask, see
#       REVISION_NOTES_P3.md s1). --no-tmean forces the 8,282-record model without
#       the term (glmmTMB engine: fits both S3 and S7 whenever tmean is usable, so
#       --no-tmean only drops the tmean tiers). If the climate file is absent or
#       stale (no ID_ABT: the pre-2026-09 89-species build) the term is dropped
#       with a message.
#   (2) FAST glmmTMB proxies of the same question that run in seconds
#       (dispformula = sigma sub-model, no phylogeny): a ladder of models from
#       "sigma ~ year only" (the continuous analogue of the pooled lnCVR) to the
#       full provenance/site/sex/season-controlled model, a Mundlak
#       within-/between-contributor decomposition of the year effect on sigma, a
#       protocol (wing_col) sensitivity, species intercepts in sigma, the
#       tmean-sample pair (S3t = S3 refitted on the tmean-complete records, S7 =
#       S3t + scaled_tmean in sigma, so sample and covariate effects separate),
#       and a two-stage log|residual| cross-check in lme4.
#   (3) lnCVR (late 2013-2018 vs early 1990-2006, the quartile periods of
#       update_descriptive_stats.R) with the s2.lnCVR sampling variance from
#       functions.R: (a) pooled per species (reproduces the published number),
#       (b) WITHIN SEX (species x sex cells), (c) WITHIN species x sex x
#       CONTRIBUTOR where both periods exist, at several minimum-n thresholds.
#       metafor::rma (REML) on each; rma.mv with a species random effect for the
#       multi-cell analyses. lnVR is reported alongside as a secondary quantity.
#   (4) MECHANICAL heterogeneity diagnostic: contributors per species-period
#       cell by period, a between-/within-contributor decomposition of the
#       within-cell variance by period, pooled vs within-contributor CV, and a
#       meta-regression of species lnCVR on the change in log(contributors).
#
# WHY: the referee (§3.6) argued that a rise in CV can be manufactured by
# pooling more measurers per species-period cell in the late period. The plan's
# triage found the within-contributor median CV rises only 0.036 -> 0.040 and
# the one contributor with paired cells gives a NEGATIVE lnCVR, so the variance
# headline must be re-tested with the source structure modelled explicitly.
#
# INPUTS:  data/derived/passer90.rda            (known sex, 73 spp, 12,571 rec)
#          data/derived/passer90_climate.rds    (P5 live-only build with ID_ABT + scaled_tmean;
#                                                optional - the script degrades gracefully)
#          data/derived/phylo_A_50trees.rds     (50 species correlation matrices of the published
#                                                brms wing model; built by _phylo_engine.R from
#                                                output/models/brm0_multiphylo.rda if absent)
#          data/raw/AvesDataLite-main/          (clootl tree cloud, offline; brms engine / fallback)
#          scripts/_phylo_engine.R, scripts/functions.R (s2.lnCVR), scripts/_sampling_config.R
# OUTPUTS: output/variance_results.rds          (every quoted fast-tier number; see $index;
#                                                $phylo = pointer + headline of the phylogenetic tier)
#          output/variance_phylo_results.rds    (glmmTMB engine: pooled, per-tree and variance
#                                                components of every phylogenetic tier)
#          output/models/variance_sigma_multiphylo.rda   (--engine brms --trees N run)
#          output/models/variance_sigma_SMOKE.rda        (--smoke run; NOT a result)
#          figures/variance_summary.png, figures/variance_lncvr_contributor.png
#
# RUN:  Rscript Analysis/scripts/atlantic_variance_sigma.R [--fast] [--trees N] [--engine glmmTMB|brms] [--smoke] [--no-tmean]
#         (no flags)  parts 2-4 + the phylogenetic glmmTMB tier on 50 trees (~25 min on a laptop)
#         --fast      parts 2-4 only (no phylogenetic tier)
#         --engine    glmmTMB (default) or brms (Stan; the Bayesian cross-check, Totoro)
#         --trees N   trees for the phylogenetic tier. Default 50 for glmmTMB; for brms
#                     50 on the server and 0 (= skip) locally.
#         --smoke     brms part 1 as a SMOKE TEST (implies --engine brms): 1 tree, chains = 2,
#                     iter = 400 (checks that the model compiles and runs; NOT a result).
#                     Locally wrap it in a ~21-min cap because there is no `timeout`
#                     on macOS:  perl -e 'alarm shift; exec @ARGV' 1260 Rscript ... --smoke
#         --no-tmean  drop scaled_tmean from the sigma sub-model even when the climate
#                     file is usable (fits the 8,282-record complete-case sample instead)
#       Works from the repo root, Analysis/, or Analysis/scripts/.
#
# Session used for development (parts 2-4 and the glmmTMB phylogenetic tier executed
# in full on 50 trees, 2026-09-09; the brms engine ran only as an aborted smoke test):
#   R 4.6.0; brms 2.23.0, rstan 2.32.7, glmmTMB 1.1.14, lme4 2.0.1,
#   metafor 5.0.1, data.table 1.18.4, dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3,
#   patchwork 1.3.2, posterior 1.7.0, ape 5.8.1, MCMCglmm 2.36, phytools 2.5.2,
#   prepR4pcm 0.5.0.9000, clootl (Aves 1.6 / Clements 2025).
# Non-obvious decisions are documented in REVISION_NOTES_P3.md.
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(tidyr); library(data.table)
  library(glmmTMB); library(lme4); library(metafor)
  library(ggplot2); library(patchwork)
})
set.seed(20260909)

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
ANALYSIS_DIR <- .find_analysis_dir()                 # = .../atlantic_birds/Analysis
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
source(script_path("_sampling_config.R"))            # SAMPLING, SAMPLING_CONTROL, .on_server
source(script_path("functions.R"))                   # s2.lnCVR()
source(script_path("_phylo_engine.R"))               # phylo_species_name, phylo_A_list, run_phylo_trees (uses the path helpers above)
dir.create(out_path("models"), showWarnings = FALSE, recursive = TRUE)
dir.create(fig_path(), showWarnings = FALSE, recursive = TRUE)

# --- command-line flags ------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
.flag <- function(name) any(args == name)
.opt  <- function(name, default) {
  i <- which(args == name); if (length(i) && length(args) > i[1]) args[i[1] + 1] else default
}
FAST     <- .flag("--fast")
SMOKE    <- .flag("--smoke")
NO_TMEAN <- .flag("--no-tmean")
ENGINE   <- .opt("--engine", "glmmTMB")
if (SMOKE && ENGINE != "brms") { message("[variance] --smoke is a Stan smoke test: engine set to brms"); ENGINE <- "brms" }
if (!ENGINE %in% c("glmmTMB", "brms")) stop("--engine must be 'glmmTMB' (default) or 'brms'")
N_TREES <- as.integer(.opt("--trees", if (ENGINE == "glmmTMB") 50L else if (.on_server) 50L else 0L))
if (SMOKE) N_TREES <- 1L
if (FAST)  N_TREES <- 0L
RUN_PHYLO_TMB <- ENGINE == "glmmTMB" && N_TREES > 0L    # part 1, glmmTMB propto engine (default)
RUN_BRMS      <- ENGINE == "brms"    && N_TREES > 0L    # part 1, Stan engine (cross-check)
MODE <- if (SMOKE) "smoke" else if (RUN_BRMS) "brms" else if (RUN_PHYLO_TMB) "glmmTMB-phylo" else "fast"
message(sprintf("[variance] mode: %s; engine = %s; trees = %d",
                c(smoke = "SMOKE TEST (brms)", brms = "full brms", `glmmTMB-phylo` = "fast tier + phylogenetic glmmTMB tier",
                  fast = "fast (glmmTMB/lme4/metafor only, no phylogeny)")[[MODE]], ENGINE, N_TREES))

res <- list(generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
            session   = list(R = R.version.string,
                             pkgs = sapply(c("brms","rstan","glmmTMB","lme4","metafor",
                                            "data.table","dplyr","ggplot2","posterior"),
                                          function(p) tryCatch(as.character(packageVersion(p)),
                                                               error = function(e) NA_character_))),
            mode = MODE, engine = ENGINE,
            definitions = c(
              period = "early = Year <= 2006, late = Year >= 2013 (quartile periods of update_descriptive_stats.R)",
              wing   = "conc.wing.length (coalesced right > left > generic), known-sex live adults, 73 species",
              lnCVR  = "log(CV_late/CV_early) + 1/(2(n_late-1)) - 1/(2(n_early-1)); s2 from functions.R::s2.lnCVR with cor(log mean, log sd) across cells per period",
              sigma_pct_decade = "100*(exp(b_sigma_scaled_yr * 10 / SD_year) - 1): % change in residual SD per decade",
              contributor = "Main_researcher; site = Municipality"))

# ===========================================================================
# 0. Data
# ===========================================================================
load(derived_path("passer90.rda"))
SD_YR <- as.numeric(attr(passer90$scaled_yr, "scaled:scale"))   # 5.0199 (pre-threshold known-sex sample)
if (!is.finite(SD_YR)) SD_YR <- readRDS(out_path("effect_scale.rds"))$sd_year_scaling
res$sd_year_scaling <- SD_YR
pct_decade <- function(b) 100 * (exp(b * 10 / SD_YR) - 1)

wing <- passer90 %>%
  filter(!is.na(conc.wing.length)) %>%
  mutate(scaled_yr  = as.numeric(scaled_yr),
         scaled_lat = as.numeric(scaled_lat),
         Sex        = factor(Sex, levels = c("Female", "Male")),
         season     = factor(season, levels = c("DJF", "MAM", "JJA", "SON")),
         wing_col   = factor(wing_col, levels = c("right", "left", "generic")),
         Main_researcher = factor(Main_researcher),
         Municipality    = factor(Municipality),
         spp = factor(spp))
res$n_wing <- nrow(wing)                                            # 8,478
# complete cases for the full sigma model (Altitude / season NA)
wing_cc <- wing %>% filter(!is.na(scaled_alt), !is.na(season))
res$n_wing_complete <- nrow(wing_cc)                                # 8,282
res$n_dropped_alt_na    <- sum(is.na(wing$scaled_alt))
res$n_dropped_season_na <- sum(is.na(wing$season) & !is.na(wing$scaled_alt))
res$n_contributors <- nlevels(droplevels(wing_cc$Main_researcher))
res$n_municipalities <- nlevels(droplevels(wing_cc$Municipality))
message(sprintf("wing records %d; complete-case for sigma models %d (%d contributors, %d municipalities)",
                res$n_wing, res$n_wing_complete, res$n_contributors, res$n_municipalities))

# Mundlak decomposition of year WITHIN vs BETWEEN contributor (for the sigma model)
wing_cc <- wing_cc %>%
  group_by(Main_researcher) %>%
  mutate(yr_src_mean   = mean(scaled_yr),
         yr_within_src = scaled_yr - yr_src_mean) %>%
  ungroup()

# --- optional climate covariate (scaled_tmean) -----------------------------
# P5's climate_extraction.R (2026-09-09 build) writes passer90_climate.rds with the
# ID_ABT key and rec_tmean / scaled_tmean (annual mean of (tmax+tmin)/2 in the
# capture year at the record's coordinates; WorldClim 2.1 historical monthly,
# CRU-TS downscaled). rec_tmean is NA where the coordinates fall outside the
# WorldClim land mask (coastal / island cells: Guarapari, Florianopolis,
# Mataraca ...), 772 of 12,571 records. The tmean-complete wing sample is kept as
# a SEPARATE object (wing_cc_t) so tiers S0-S6 stay on the 8,282-record
# complete-case sample and the tmean effect can be separated from the sample
# change (S3t vs S7). The brms model follows the plan (sigma ~ ... + scaled_tmean)
# on wing_cc_t unless --no-tmean.
clim_file <- derived_path("passer90_climate.rds")
USE_TMEAN <- FALSE
wing_cc_t <- NULL
res$scaled_tmean <- list(used = FALSE, reason = "file absent")
if (file.exists(clim_file)) {
  clim <- readRDS(clim_file)
  n_spp_clim <- length(unique(clim$Binomial))
  if (!"ID_ABT" %in% names(clim)) {
    res$scaled_tmean <- list(used = FALSE,
      reason = sprintf("passer90_climate.rds is a stale build (%d rows, %d species, no ID_ABT key) - cannot be joined to the live-only sample; P5 regenerates it with climate_extraction.R",
                       nrow(clim), n_spp_clim))
  } else {
    cl <- clim %>% select(ID_ABT, scaled_tmean, rec_tmean) %>% distinct()
    tmp <- wing_cc %>% left_join(cl, by = "ID_ABT")
    stopifnot(nrow(tmp) == nrow(wing_cc))                 # ID_ABT must be unique in the climate file
    cover <- mean(!is.na(tmp$scaled_tmean))
    na_by_period <- tmp %>%
      mutate(period = case_when(Year <= 2006 ~ "early (<=2006)", Year >= 2013 ~ "late (>=2013)", TRUE ~ "middle (2007-2012)")) %>%
      group_by(period) %>% summarise(n = n(), n_tmean_na = sum(is.na(scaled_tmean)), .groups = "drop")
    if (cover >= 0.90 && !NO_TMEAN) {
      wing_cc_t <- tmp %>% filter(!is.na(scaled_tmean)); USE_TMEAN <- TRUE
      res$scaled_tmean <- list(
        used = TRUE, coverage = cover, n = nrow(wing_cc_t), n_dropped = nrow(tmp) - nrow(wing_cc_t),
        climate_file_rows = nrow(clim), climate_file_species = n_spp_clim,
        na_by_period = na_by_period,
        cor_year_tmean = cor(wing_cc_t$scaled_yr, wing_cc_t$scaled_tmean),
        cor_lat_tmean  = cor(wing_cc_t$scaled_lat, wing_cc_t$scaled_tmean),
        na_top_municipalities = head(sort(table(tmp$Municipality[is.na(tmp$scaled_tmean)]), decreasing = TRUE), 5),
        definition = "rec_tmean = annual mean of (tmax + tmin)/2 in the capture year at the record coordinates (climate_extraction.R); scaled_tmean = scale() over the 12,571-record live-only sample",
        note = paste("Records without tmean are coastal/island coordinates outside the WorldClim land mask (P5).",
                     "Tiers S0-S6 use the full 8,282-record complete-case sample; S3t and S7 and the brms model use the tmean-complete subset."))
    } else if (NO_TMEAN) {
      res$scaled_tmean <- list(used = FALSE, coverage = cover,
        reason = "--no-tmean flag: scaled_tmean dropped by request; brms model on the 8,282-record complete-case sample",
        na_by_period = na_by_period)
    } else {
      res$scaled_tmean <- list(used = FALSE, coverage = cover,
        reason = sprintf("passer90_climate.rds covers only %.1f%% of the complete-case wing sample (< 90%%)", 100 * cover),
        na_by_period = na_by_period)
    }
  }
  rm(clim)
}
message("scaled_tmean: ", if (USE_TMEAN) sprintf("included (coverage %.1f%%, n = %d)", 100 * res$scaled_tmean$coverage, res$scaled_tmean$n)
                          else paste("DROPPED -", res$scaled_tmean$reason))

# ===========================================================================
# 2. FAST proxies: glmmTMB distributional (dispformula) models
# ===========================================================================
# Same complete-case sample for every tier so the year coefficient is comparable.
message("\n[2] glmmTMB sigma-model ladder ...")
fit_tmb <- function(mean_f, disp_f, data = wing_cc) {
  m <- glmmTMB(mean_f, dispformula = disp_f, data = data, REML = TRUE)
  m
}
mean_M0 <- conc.wing.length ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)
mean_M3 <- conc.wing.length ~ Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season +
  (1 + scaled_yr || spp) + (1 | Main_researcher) + (1 | Municipality)

tiers <- list(
  S0 = list(label = "S0: baseline mean (M0); sigma ~ year",
            mean = mean_M0, disp = ~ scaled_yr),
  S1 = list(label = "S1: baseline mean (M0); sigma ~ year + Sex",
            mean = mean_M0, disp = ~ scaled_yr + Sex),
  S2 = list(label = "S2: mean + contributor + site; sigma ~ year + Sex + (1|contributor)",
            mean = mean_M3, disp = ~ scaled_yr + Sex + (1 | Main_researcher)),
  S3 = list(label = "S3: full (M3 mean); sigma ~ year + Sex + (1|contributor) + (1|site)",
            mean = mean_M3, disp = ~ scaled_yr + Sex + (1 | Main_researcher) + (1 | Municipality)),
  S4 = list(label = "S4: S3 with sigma year split within/between contributor (Mundlak)",
            mean = mean_M3, disp = ~ yr_within_src + yr_src_mean + Sex + (1 | Main_researcher) + (1 | Municipality)),
  S5 = list(label = "S5: S3 + wing-column protocol proxy in sigma",
            mean = mean_M3, disp = ~ scaled_yr + Sex + wing_col + (1 | Main_researcher) + (1 | Municipality)),
  S6 = list(label = "S6: S3 with species intercepts in sigma",
            mean = mean_M3, disp = ~ scaled_yr + Sex + (1 | Main_researcher) + (1 | Municipality) + (1 | spp))
)
if (USE_TMEAN) {
  # tmean-complete subset: S3t isolates the sample change, S7 adds the covariate
  tiers$S3t <- list(label = "S3t: S3 refitted on the tmean-complete subset (no tmean term)",
                    mean = mean_M3, disp = ~ scaled_yr + Sex + (1 | Main_researcher) + (1 | Municipality),
                    data = "wing_cc_t")
  tiers$S7  <- list(label = "S7: S3t + scaled_tmean in sigma (the plan's sigma specification, minus phylogeny)",
                    mean = mean_M3, disp = ~ scaled_yr + scaled_tmean + Sex + (1 | Main_researcher) + (1 | Municipality),
                    data = "wing_cc_t")
}

tmb_fits <- list(); tmb_rows <- list()
for (nm in names(tiers)) {
  t0 <- Sys.time()
  tier_dat <- if (identical(tiers[[nm]]$data, "wing_cc_t")) wing_cc_t else wing_cc
  f <- tryCatch(fit_tmb(tiers[[nm]]$mean, tiers[[nm]]$disp, data = tier_dat), error = function(e) e)
  if (inherits(f, "error")) { message("  ", nm, " FAILED: ", conditionMessage(f)); next }
  tmb_fits[[nm]] <- f
  cf_d <- summary(f)$coefficients$disp
  cf_c <- summary(f)$coefficients$cond
  vc   <- VarCorr(f)
  sd_disp <- if (!is.null(vc$disp)) sapply(vc$disp, function(x) attr(x, "stddev")[[1]]) else numeric(0)
  rows <- data.frame(tier = nm, label = tiers[[nm]]$label, component = "sigma",
                     par = rownames(cf_d), estimate = cf_d[, "Estimate"], se = cf_d[, "Std. Error"],
                     z = cf_d[, "z value"], p = cf_d[, "Pr(>|z|)"], row.names = NULL)
  rows$lower <- rows$estimate - 1.96 * rows$se; rows$upper <- rows$estimate + 1.96 * rows$se
  rows$pct_decade <- ifelse(rows$par %in% c("scaled_yr", "yr_within_src", "yr_src_mean"),
                            pct_decade(rows$estimate), NA)
  rows$pct_decade_lo <- ifelse(is.na(rows$pct_decade), NA, pct_decade(rows$lower))
  rows$pct_decade_hi <- ifelse(is.na(rows$pct_decade), NA, pct_decade(rows$upper))
  # scaled_tmean is per SD of record-year temperature: % change in residual SD per SD tmean
  rows$pct_per_sd_tmean    <- ifelse(rows$par == "scaled_tmean", 100 * (exp(rows$estimate) - 1), NA)
  rows$pct_per_sd_tmean_lo <- ifelse(rows$par == "scaled_tmean", 100 * (exp(rows$lower) - 1), NA)
  rows$pct_per_sd_tmean_hi <- ifelse(rows$par == "scaled_tmean", 100 * (exp(rows$upper) - 1), NA)
  mrow <- data.frame(tier = nm, label = tiers[[nm]]$label, component = "mu",
                     par = "scaled_yr", estimate = cf_c["scaled_yr", "Estimate"], se = cf_c["scaled_yr", "Std. Error"],
                     z = cf_c["scaled_yr", "z value"], p = cf_c["scaled_yr", "Pr(>|z|)"], row.names = NULL)
  mrow$lower <- mrow$estimate - 1.96 * mrow$se; mrow$upper <- mrow$estimate + 1.96 * mrow$se
  mrow$pct_decade <- NA; mrow$pct_decade_lo <- NA; mrow$pct_decade_hi <- NA
  mrow$pct_per_sd_tmean <- NA; mrow$pct_per_sd_tmean_lo <- NA; mrow$pct_per_sd_tmean_hi <- NA
  tmb_rows[[nm]] <- rbind(rows, mrow)
  tmb_rows[[nm]]$n <- nobs(f); tmb_rows[[nm]]$logLik <- as.numeric(logLik(f))
  tmb_rows[[nm]]$converged <- isTRUE(f$fit$convergence == 0) && isTRUE(f$sdr$pdHess)
  tmb_rows[[nm]]$sd_sigma_re <- paste(names(sd_disp), round(sd_disp, 3), sep = "=", collapse = "; ")
  message(sprintf("  %s  sigma~year b = %+.4f (SE %.4f)  => %+.1f %%/decade   [%.1fs]", nm,
                  rows$estimate[rows$par %in% c("scaled_yr","yr_within_src")][1],
                  rows$se[rows$par %in% c("scaled_yr","yr_within_src")][1],
                  rows$pct_decade[rows$par %in% c("scaled_yr","yr_within_src")][1],
                  as.numeric(Sys.time() - t0, units = "secs")))
}
res$glmmTMB <- do.call(rbind, tmb_rows)
res$glmmTMB_note <- paste(
  "Tiers S0-S6 on the same 8,282-record complete-case sample (Altitude and season non-NA); S3t and S7 on the tmean-complete subset (see $scaled_tmean). REML.",
  "sigma coefficients are on the log-SD scale; pct_decade converts scaled_yr slopes to % change in residual SD per decade.",
  "S4 splits year into within-contributor deviation (yr_within_src) and contributor mean year (yr_src_mean).",
  "No phylogenetic term in these tiers; the phylogenetic tier (part 1: $phylo / variance_phylo_results.rds, glmmTMB propto engine, or $brms) refits every tier with it.")

# Random-effect SDs of the sigma sub-model in S3 (how much residual SD varies between contributors)
if (!is.null(tmb_fits$S3)) {
  vc <- VarCorr(tmb_fits$S3)
  res$glmmTMB_S3_sigma_re_sd <- sapply(vc$disp, function(x) attr(x, "stddev")[[1]])
  res$glmmTMB_S3_mu_re_sd    <- unlist(lapply(vc$cond, function(x) attr(x, "stddev")))
  # contributor-level residual SD multipliers exp(u_j): spread across contributors
  re_src <- ranef(tmb_fits$S3)$disp$Main_researcher[, 1]
  res$glmmTMB_S3_contributor_sd_multiplier_range <- round(range(exp(re_src)), 2)
  res$glmmTMB_S3_contributor_sd_multiplier_cv <- round(sd(exp(re_src)) / mean(exp(re_src)), 3)
}

# --- two-stage cross-check: log|residual| from an lme4 mean model -----------
message("[2b] two-stage log|residual| cross-check (lme4) ...")
lmm_mean <- lmer(mean_M3, data = wing_cc, REML = TRUE,
                 control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
wing_cc$resid_mean <- resid(lmm_mean)
# homoscedastic lme4 fit of the same M3 mean structure on the same sample: direct comparison
# for the mean-year slope with the heteroscedastic glmmTMB S3 fit (which re-weights records
# by their contributor/site-specific residual SD)
.cf <- summary(lmm_mean)$coefficients["scaled_yr", ]
res$lme4_mean_M3_year <- c(estimate = unname(.cf["Estimate"]), se = unname(.cf["Std. Error"]), t = unname(.cf["t value"]),
                           lower = unname(.cf["Estimate"] - 1.96 * .cf["Std. Error"]),
                           upper = unname(.cf["Estimate"] + 1.96 * .cf["Std. Error"]), n = nobs(lmm_mean))
res$mean_year_slope_comparison <- data.frame(
  model = c("lme4 homoscedastic, M3 mean structure (same complete-case sample)",
            "glmmTMB S3: M3 mean + sigma ~ year + Sex + (1|contributor) + (1|site)"),
  estimate = c(res$lme4_mean_M3_year[["estimate"]], res$glmmTMB$estimate[res$glmmTMB$tier == "S3" & res$glmmTMB$component == "mu"]),
  se       = c(res$lme4_mean_M3_year[["se"]],       res$glmmTMB$se[res$glmmTMB$tier == "S3" & res$glmmTMB$component == "mu"]),
  n = nobs(lmm_mean))
res$mean_year_slope_comparison$lower <- res$mean_year_slope_comparison$estimate - 1.96 * res$mean_year_slope_comparison$se
res$mean_year_slope_comparison$upper <- res$mean_year_slope_comparison$estimate + 1.96 * res$mean_year_slope_comparison$se
message("mean-model year slope, homoscedastic lme4 vs heteroscedastic glmmTMB S3:")
print(res$mean_year_slope_comparison, digits = 3)
res$two_stage_n_zero_resid <- sum(wing_cc$resid_mean == 0)
ts_dat <- wing_cc %>% filter(resid_mean != 0) %>% mutate(log_abs_r = log(abs(resid_mean)))
ts0 <- lmer(log_abs_r ~ scaled_yr + Sex + (1 | spp), data = ts_dat, REML = TRUE)
ts1 <- lmer(log_abs_r ~ scaled_yr + Sex + (1 | spp) + (1 | Main_researcher) + (1 | Municipality),
            data = ts_dat, REML = TRUE)
ts2 <- lmer(log_abs_r ~ yr_within_src + yr_src_mean + Sex + (1 | spp) + (1 | Main_researcher) + (1 | Municipality),
            data = ts_dat, REML = TRUE)
.ts_row <- function(m, tier, par) {
  cf <- summary(m)$coefficients
  data.frame(tier = tier, par = par, estimate = cf[par, "Estimate"], se = cf[par, "Std. Error"],
             t = cf[par, "t value"], n = nobs(m), row.names = NULL) %>%
    mutate(lower = estimate - 1.96 * se, upper = estimate + 1.96 * se,
           # E[log|r|] shifts by the same log-SD slope, so the % conversion applies directly
           pct_decade = pct_decade(estimate), pct_decade_lo = pct_decade(lower), pct_decade_hi = pct_decade(upper))
}
res$two_stage <- rbind(.ts_row(ts0, "T0: log|r| ~ year + Sex + (1|spp)", "scaled_yr"),
                       .ts_row(ts1, "T1: + (1|contributor) + (1|site)", "scaled_yr"),
                       .ts_row(ts2, "T2: Mundlak within-contributor year", "yr_within_src"),
                       .ts_row(ts2, "T2: Mundlak between-contributor year", "yr_src_mean"))
res$two_stage_note <- paste(
  "Stage 1: lmer mean model with the M3 structure on the complete-case sample; stage 2: log|residual| regressed on year etc.",
  "Slopes on log|r| estimate the log-SD slope (up to a constant), so pct_decade is comparable with the glmmTMB sigma rows.",
  sprintf("%d records with an exactly-zero residual were dropped from stage 2.", res$two_stage_n_zero_resid))
print(res$two_stage[, c("tier", "estimate", "se", "t", "pct_decade", "pct_decade_lo", "pct_decade_hi")])

# ===========================================================================
# 3. lnCVR late vs early: pooled, within sex, within species x sex x contributor
# ===========================================================================
message("\n[3] lnCVR meta-analyses ...")
quart <- wing %>%
  filter(Year <= 2006 | Year >= 2013) %>%
  mutate(period = ifelse(Year <= 2006, "early", "late"),
         Main_researcher = as.character(Main_researcher),
         Binomial = as.character(Binomial), Sex = as.character(Sex))
res$n_quartile_wing <- nrow(quart)
res$n_quartile_by_period <- table(quart$period)

# cell summaries at any grouping; returns wide table with early/late columns
cell_summary <- function(dat, ...) {
  dat %>% group_by(..., period) %>%
    summarise(n = n(), m = mean(conc.wing.length), s = sd(conc.wing.length), .groups = "drop")
}
# lnCVR + lnVR per unit (group vars) where both periods exist and both SDs are defined
lncvr_table <- function(cells, group_vars, min_n = 2L) {
  w <- cells %>% filter(n >= min_n, !is.na(s), s > 0) %>%
    pivot_wider(id_cols = all_of(group_vars), names_from = period, values_from = c(n, m, s)) %>%
    filter(!is.na(n_early), !is.na(n_late))
  if (nrow(w) == 0) return(w)
  # mean-SD correlation across cells within each period (Nakagawa et al. 2015 correction)
  ce <- suppressWarnings(cor(log(w$m_early), log(w$s_early), use = "complete.obs"))
  cl <- suppressWarnings(cor(log(w$m_late),  log(w$s_late),  use = "complete.obs"))
  if (!is.finite(ce)) ce <- 0; if (!is.finite(cl)) cl <- 0
  w$lnCVR <- log((w$s_late / w$m_late) / (w$s_early / w$m_early)) +
    1 / (2 * (w$n_late - 1)) - 1 / (2 * (w$n_early - 1))
  w$s2_lnCVR <- s2.lnCVR(w$m_early, w$s_early, w$n_early, ce, w$m_late, w$s_late, w$n_late, cl)
  w$lnVR <- log(w$s_late / w$s_early) + 1 / (2 * (w$n_late - 1)) - 1 / (2 * (w$n_early - 1))
  w$s2_lnVR <- 1 / (2 * (w$n_late - 1)) + 1 / (2 * (w$n_early - 1))
  attr(w, "cor_early") <- ce; attr(w, "cor_late") <- cl
  w
}
rma_row <- function(w, yi = "lnCVR", vi = "s2_lnCVR", label, cluster = NULL) {
  if (nrow(w) < 2) return(data.frame(analysis = label, k = nrow(w), estimate = NA, se = NA, lower = NA, upper = NA,
                                     pval = NA, I2 = NA, tau2 = NA, Q = NA, Qp = NA, model = "rma", row.names = NULL))
  m <- rma(yi = w[[yi]], vi = w[[vi]], method = "REML", test = "z")
  out <- data.frame(analysis = label, k = m$k, estimate = as.numeric(m$b), se = m$se, lower = m$ci.lb, upper = m$ci.ub,
                    pval = m$pval, I2 = m$I2, tau2 = m$tau2, Q = m$QE, Qp = m$QEp, model = "rma", row.names = NULL)
  if (!is.null(cluster) && length(unique(w[[cluster]])) > 1 && nrow(w) > length(unique(w[[cluster]]))) {
    w$.cell <- seq_len(nrow(w))
    mv <- tryCatch(rma.mv(yi = w[[yi]], V = w[[vi]], random = list(~ 1 | .cell, ~ 1 | .clu),
                          data = data.frame(.cell = w$.cell, .clu = w[[cluster]]), method = "REML"),
                   error = function(e) NULL)
    if (!is.null(mv)) {
      out <- rbind(out, data.frame(analysis = label, k = mv$k, estimate = as.numeric(mv$b), se = mv$se,
                                   lower = mv$ci.lb, upper = mv$ci.ub, pval = mv$pval, I2 = NA,
                                   tau2 = sum(mv$sigma2), Q = mv$QE, Qp = mv$QEp,
                                   model = paste0("rma.mv (+", cluster, ")"), row.names = NULL))
    }
  }
  out
}

# (a0) EXACT reproduction of the update_descriptive_stats.R code path. NOTE: that
# script groups ALL quartile records (wing NA included) and uses n() as the wing
# cell size, so n is inflated by records without a wing measurement (and the same
# wing n is reused for the bill lnCVR). Kept here only to validate against
# descriptive_summary.rds$lncvr_wing; the rows below use the non-NA wing count.
q_all <- passer90 %>% filter(Year <= 2006 | Year >= 2013) %>%
  mutate(period = ifelse(Year <= 2006, "early", "late"))
cells_orig <- q_all %>% group_by(Binomial, period) %>%
  summarise(n = n(), m = mean(conc.wing.length, na.rm = TRUE), s = sd(conc.wing.length, na.rm = TRUE), .groups = "drop")
w_orig <- lncvr_table(cells_orig, "Binomial", min_n = 1L)
lncvr_rows <- rma_row(w_orig, label = "pooled per species, EXACT update_descriptive_stats.R code path (n counts NA-wing records)")
res$lncvr_published_check <- list(
  reproduced = lncvr_rows[1, c("k", "estimate", "lower", "upper", "I2")],
  descriptive_summary = tryCatch({ ds <- readRDS(out_path("descriptive_summary.rds"))
    c(k = ds$lncvr_wing_k, estimate = ds$lncvr_wing, lower = ds$lncvr_wing_lo, upper = ds$lncvr_wing_hi, I2 = ds$lncvr_wing_I2) },
    error = function(e) NA),
  n_inflation = sum(cells_orig$n) - sum(cell_summary(quart, Binomial)$n),
  note = "update_descriptive_stats.R sets sample.n.wing = n() over all quartile records, including those with NA wing; the corrected rows below count non-NA wing records only")

# (a) pooled per species, corrected n (non-NA wing records; cells with >= 2 so an SD exists)
cells_spp <- cell_summary(quart, Binomial)
w_spp <- lncvr_table(cells_spp, "Binomial", min_n = 2L)
lncvr_rows <- rbind(lncvr_rows, rma_row(w_spp, label = "pooled per species (corrected n = non-NA wing records)"))
lncvr_rows <- rbind(lncvr_rows, rma_row(w_spp, yi = "lnVR", vi = "s2_lnVR", label = "pooled per species, lnVR"))
res$lncvr_species_table <- w_spp

# (b) within sex: species x sex cells
cells_sx <- cell_summary(quart, Binomial, Sex)
for (mn in c(2L, 5L, 10L)) {
  w <- lncvr_table(cells_sx, c("Binomial", "Sex"), min_n = mn)
  lncvr_rows <- rbind(lncvr_rows,
                      rma_row(w, label = sprintf("within sex (species x sex cells), n >= %d", mn), cluster = "Binomial"))
  if (mn == 5L) {
    lncvr_rows <- rbind(lncvr_rows, rma_row(w, yi = "lnVR", vi = "s2_lnVR",
                                            label = "within sex, n >= 5, lnVR", cluster = "Binomial"))
    for (sx in c("Female", "Male"))
      lncvr_rows <- rbind(lncvr_rows, rma_row(w[w$Sex == sx, ], label = sprintf("within sex, %s only, n >= 5", sx)))
    res$lncvr_sex_table <- w
  }
  if (mn == 2L) res$lncvr_sex_table_min2 <- w
}

# (c) within species x sex x contributor
cells_src <- cell_summary(quart, Binomial, Sex, Main_researcher)
src_tables <- list()
for (mn in c(2L, 3L, 5L, 10L)) {
  w <- lncvr_table(cells_src, c("Binomial", "Sex", "Main_researcher"), min_n = mn)
  src_tables[[as.character(mn)]] <- w
  lncvr_rows <- rbind(lncvr_rows,
                      rma_row(w, label = sprintf("within species x sex x contributor, n >= %d", mn), cluster = "Main_researcher"))
  if (mn == 5L)
    lncvr_rows <- rbind(lncvr_rows, rma_row(w, yi = "lnVR", vi = "s2_lnVR",
                                            label = "within species x sex x contributor, n >= 5, lnVR", cluster = "Main_researcher"))
}
res$lncvr_contributor_table <- src_tables[["5"]]
res$lncvr_contributor_table_min2 <- src_tables[["2"]]
# per-contributor summaries (k, mean, rma) at n >= 5 and n >= 2
per_src <- function(w, mn) {
  w %>% group_by(Main_researcher) %>%
    summarise(k = n(), mean_lnCVR = mean(lnCVR), median_lnCVR = median(lnCVR),
              n_negative = sum(lnCVR < 0), .groups = "drop") %>%
    mutate(min_n = mn)
}
res$lncvr_by_contributor <- bind_rows(per_src(src_tables[["2"]], 2L), per_src(src_tables[["3"]], 3L),
                                      per_src(src_tables[["5"]], 5L), per_src(src_tables[["10"]], 10L))
for (mn in c("5", "10")) {
  wc <- src_tables[[mn]] %>% filter(Main_researcher == "E.Carrano")
  if (nrow(wc) >= 2)
    lncvr_rows <- rbind(lncvr_rows, rma_row(wc, label = sprintf("E.Carrano only, n >= %s", mn)))
}
lncvr_rows$estimate_pct <- 100 * (exp(lncvr_rows$estimate) - 1)   # % change in CV (or SD for lnVR)
res$lncvr <- lncvr_rows
res$lncvr_cor_meansd <- c(species_early = attr(w_spp, "cor_early"), species_late = attr(w_spp, "cor_late"),
                          sex5_early = attr(res$lncvr_sex_table, "cor_early"), sex5_late = attr(res$lncvr_sex_table, "cor_late"),
                          src5_early = attr(src_tables[["5"]], "cor_early"), src5_late = attr(src_tables[["5"]], "cor_late"))
res$lncvr_note <- paste(
  "The pooled-per-species row uses exactly the update_descriptive_stats.R definition (cells with >= 2 records so an SD exists)",
  "and should match descriptive_summary.rds$lncvr_wing. 'rma.mv (+cluster)' rows add a random intercept for the cluster",
  "(species for the within-sex analysis; contributor for the within-contributor analysis) because cells within a cluster are not independent.",
  "The plan's 'E. Carrano k = 6, mean -0.46' does not reproduce at any single threshold: k = 7 (n >= 5) or 4 (n >= 10); see lncvr_by_contributor.")
print(res$lncvr[, c("analysis", "model", "k", "estimate", "lower", "upper", "I2")], digits = 3)

# ===========================================================================
# 4. Mechanical heterogeneity: contributors per cell, variance decomposition
# ===========================================================================
message("\n[4] mechanical heterogeneity diagnostics ...")
# 4a. contributors per species-period cell (wing records), by period
cpc <- quart %>% group_by(Binomial, period) %>%
  summarise(n_contrib = n_distinct(Main_researcher), n_site = n_distinct(Municipality), n = n(), .groups = "drop")
res$contributors_per_cell <- cpc
res$contributors_per_cell_summary <- cpc %>% group_by(period) %>%
  summarise(cells = n(), mean_contrib = mean(n_contrib), median_contrib = median(n_contrib),
            share_single_contrib = mean(n_contrib == 1),
            mean_sites = mean(n_site), median_sites = median(n_site), .groups = "drop")
both <- cpc %>% group_by(Binomial) %>% filter(n_distinct(period) == 2) %>% ungroup()
res$contributors_per_cell_summary_species_in_both <- both %>% group_by(period) %>%
  summarise(cells = n(), mean_contrib = mean(n_contrib), median_contrib = median(n_contrib), .groups = "drop")
# paired change in contributors per species (late - early), Wilcoxon signed-rank
pc <- both %>% select(Binomial, period, n_contrib) %>% pivot_wider(names_from = period, values_from = n_contrib)
res$contributors_paired <- list(k = nrow(pc), mean_early = mean(pc$early), mean_late = mean(pc$late),
                                n_increase = sum(pc$late > pc$early), n_decrease = sum(pc$late < pc$early),
                                wilcoxon_p = suppressWarnings(wilcox.test(pc$late, pc$early, paired = TRUE))$p.value)
# species x sex x period cells (matched to the within-sex lnCVR units)
cpc_sx <- quart %>% group_by(Binomial, Sex, period) %>%
  summarise(n_contrib = n_distinct(Main_researcher), n = n(), .groups = "drop")
res$contributors_per_sex_cell_summary <- cpc_sx %>% group_by(period) %>%
  summarise(cells = n(), mean_contrib = mean(n_contrib), median_contrib = median(n_contrib), .groups = "drop")
print(res$contributors_per_cell_summary)

# 4b. within-cell variance decomposition: between- vs within-contributor share, by period
decomp <- quart %>% group_by(Binomial, Sex, period) %>% filter(n() >= 5) %>%
  group_by(Binomial, Sex, period, Main_researcher) %>% mutate(src_mean = mean(conc.wing.length)) %>%
  group_by(Binomial, Sex, period) %>%
  summarise(n = n(), n_contrib = n_distinct(Main_researcher),
            var_total   = var(conc.wing.length),
            ss_within   = sum((conc.wing.length - src_mean)^2),
            ss_between  = sum((src_mean - mean(conc.wing.length))^2),
            share_between = ss_between / (ss_within + ss_between),
            cv_pooled   = sd(conc.wing.length) / mean(conc.wing.length),
            .groups = "drop")
res$variance_decomposition <- decomp
res$variance_decomposition_summary <- decomp %>% group_by(period) %>%
  summarise(cells = n(), mean_share_between = mean(share_between), median_share_between = median(share_between),
            median_cv_pooled = median(cv_pooled), mean_cv_pooled = mean(cv_pooled), .groups = "drop")
# within-contributor CV: species x sex x contributor x period cells with n >= 5 (and n >= 10)
for (mn in c(5L, 10L)) {
  cvw <- cells_src %>% filter(n >= mn, !is.na(s)) %>% mutate(cv = s / m) %>% group_by(period) %>%
    summarise(cells = n(), median_cv = median(cv), mean_cv = mean(cv), .groups = "drop")
  cvp <- cells_sx %>% filter(n >= mn, !is.na(s)) %>% mutate(cv = s / m) %>% group_by(period) %>%
    summarise(cells = n(), median_cv = median(cv), mean_cv = mean(cv), .groups = "drop")
  res[[sprintf("cv_by_period_min%d", mn)]] <- bind_rows(mutate(cvp, level = "species x sex (pooled across contributors)"),
                                                        mutate(cvw, level = "species x sex x contributor"))
}
print(res$cv_by_period_min5)
print(res$variance_decomposition_summary)

# 4c. meta-regression: species lnCVR vs change in log(contributors) and log(sites)
mr_dat <- w_spp %>% left_join(pc %>% rename(contrib_early = early, contrib_late = late), by = "Binomial") %>%
  left_join(cpc %>% select(Binomial, period, n_site) %>% pivot_wider(names_from = period, values_from = n_site) %>%
              rename(site_early = early, site_late = late), by = "Binomial") %>%
  mutate(d_log_contrib = log(contrib_late) - log(contrib_early),
         d_log_site    = log(site_late) - log(site_early))
mr1 <- rma(yi = lnCVR, vi = s2_lnCVR, mods = ~ d_log_contrib, data = mr_dat, method = "REML")
mr2 <- rma(yi = lnCVR, vi = s2_lnCVR, mods = ~ d_log_site, data = mr_dat, method = "REML")
mr3 <- rma(yi = lnCVR, vi = s2_lnCVR, mods = ~ d_log_contrib + d_log_site, data = mr_dat, method = "REML")
.mr_row <- function(m, label, mod) data.frame(
  analysis = label, moderator = mod, k = m$k, slope = m$b[mod, 1], se = m$se[which(rownames(m$b) == mod)],
  lower = m$ci.lb[which(rownames(m$b) == mod)], upper = m$ci.ub[which(rownames(m$b) == mod)],
  pval = m$pval[which(rownames(m$b) == mod)], intercept = m$b["intrcpt", 1],
  intercept_lower = m$ci.lb[1], intercept_upper = m$ci.ub[1], R2 = m$R2, QM_p = m$QMp, row.names = NULL)
res$lncvr_metaregression <- rbind(.mr_row(mr1, "species lnCVR ~ d log(contributors)", "d_log_contrib"),
                                  .mr_row(mr2, "species lnCVR ~ d log(sites)", "d_log_site"),
                                  .mr_row(mr3, "species lnCVR ~ d log(contributors) + d log(sites)", "d_log_contrib"),
                                  .mr_row(mr3, "species lnCVR ~ d log(contributors) + d log(sites)", "d_log_site"))
res$lncvr_metaregression_note <- paste(
  "Meta-regression of the per-species lnCVR (published definition) on the late-minus-early change in log number of",
  "contributors / municipalities. 'intercept' is the predicted lnCVR for a species whose contributor (site) count did not change.",
  "Slope > 0 with a shrunken intercept indicates the mechanical inflation the referee predicted.")
res$lncvr_metaregression_data <- mr_dat
print(res$lncvr_metaregression[, c("analysis", "moderator", "k", "slope", "lower", "upper", "pval", "intercept", "R2")], digits = 3)

# ===========================================================================
# 1a. PHYLOGENETIC distributional models — glmmTMB propto engine (default)
# ===========================================================================
# Every tier of the sigma ladder above is refitted with a phylogenetic species
# intercept (propto(0 + species_name | g, A), _phylo_engine.R) across the SAME 50
# trees as the published brms wing model, sigma sub-model via dispformula, REML,
# Rubin-pooled over trees. Same records as the non-phylogenetic tiers (8,282 /
# 7,651), so each tier's phylogenetic vs non-phylogenetic estimate is directly
# comparable. Headline models (REVISION_PLAN.md Phase 3): S7 = plan specification
# (sigma ~ scaled_yr + scaled_tmean + Sex + (1|contributor) + (1|site), 7,651
# records) and S3 = the --no-tmean variant on the 8,282-record sample.
if (RUN_PHYLO_TMB) {
  message(sprintf("\n[1a] phylogenetic glmmTMB distributional models across %d trees (propto + dispformula) ...", N_TREES))
  t_phylo0 <- Sys.time()
  # tree tip names (eBird synonyms as in atlantic_parallel.R); Herpsilochmus_sellowi is
  # absent from the trees and must be dropped (0 wing records in this sample: checked below)
  add_phylo_cols <- function(d) {
    d <- d %>% mutate(species_name = phylo_species_name(Binomial), spp = species_name,
                      scaled_lon = as.numeric(scaled_lon), scaled_alt = as.numeric(scaled_alt))
    as.data.frame(d[d$species_name != "Herpsilochmus_sellowi", , drop = FALSE])
  }
  ph_cc   <- add_phylo_cols(wing_cc)
  ph_cc_t <- if (USE_TMEAN) add_phylo_cols(wing_cc_t) else NULL
  n_dropped_sellowi <- nrow(wing_cc) - nrow(ph_cc)
  ph_species <- sort(unique(ph_cc$species_name))
  A_list <- phylo_A_list(ph_species, n_trees = N_TREES)
  A_list <- A_list[seq_len(min(N_TREES, length(A_list)))]
  A_cache <- derived_path("phylo_A_50trees.rds")
  message(sprintf("  %d records (%d dropped as H. sellowi), %d species, %d trees from %s",
                  nrow(ph_cc), n_dropped_sellowi, length(ph_species), length(A_list),
                  if (file.exists(A_cache)) basename(A_cache) else "clootl fallback"))

  # baseline residual SD for the phylogenetic-share denominator: sigma() is NA under a
  # dispformula, so use exp(b_sigma_Intercept) (female, reference season, covariates at
  # 0, sigma random effects at 0) per tree.
  phylo_share <- function(run) {
    b0 <- run$per_tree_fixed %>% filter(component == "disp", par == "(Intercept)") %>% select(tree, b_sigma_int = estimate)
    vc <- run$varcomp %>% left_join(b0, by = "tree")
    other_cols <- setdiff(names(vc), c("tree", "sd_phylo", "sigma", "phylo_prop", "converged", "b_sigma_int"))
    other_ss <- rowSums(as.matrix(vc[, other_cols, drop = FALSE])^2, na.rm = TRUE)
    vc$sigma_baseline <- exp(vc$b_sigma_int)
    vc$phylo_share_incl_resid <- vc$sd_phylo^2 / (vc$sd_phylo^2 + other_ss + vc$sigma_baseline^2)
    # tidy_phylo_fit() passes the SD names through data.frame(), so "spp:(Intercept)" arrives as "spp..Intercept."
    spp_int <- grep("^spp.*Intercept", names(vc), value = TRUE)[1]
    vc$phylo_share_between_species <- vc$sd_phylo^2 / (vc$sd_phylo^2 + vc[[spp_int]]^2)
    names(vc)[names(vc) == spp_int] <- "sd_spp_intercept"
    names(vc)[names(vc) == "spp.scaled_yr"] <- "sd_spp_year_slope"
    names(vc) <- sub("^Main_researcher.*Intercept.*$", "sd_contributor", names(vc))
    names(vc) <- sub("^Municipality.*Intercept.*$", "sd_site", names(vc))
    vc
  }

  phylo_runs <- list(); phylo_rows <- list(); phylo_pooled_all <- list(); phylo_vc <- list(); phylo_timing <- list()
  for (nm in names(tiers)) {
    dat_i <- if (identical(tiers[[nm]]$data, "wing_cc_t")) ph_cc_t else ph_cc
    r <- tryCatch(run_phylo_trees(tiers[[nm]]$mean, dat_i, A_list, dispformula = tiers[[nm]]$disp, verbose = FALSE),
                  error = function(e) e)
    if (inherits(r, "error")) { message("  ", nm, " FAILED: ", conditionMessage(r)); next }
    r$fits <- NULL
    r$varcomp <- phylo_share(r)
    r$formula <- c(mean = deparse1(tiers[[nm]]$mean), sigma = deparse1(tiers[[nm]]$disp),
                   phylo = "+ propto(0 + species_name | g, A)  [appended by fit_phylo_glmmtmb]")
    r$label <- tiers[[nm]]$label; r$n <- nrow(dat_i); r$n_species <- length(unique(dat_i$species_name))
    r$n_converged <- sum(r$varcomp$converged)
    # sensitivity: Rubin pool over the trees with a positive-definite Hessian only (the engine pools all trees).
    # Diagnosis (2026-09-09, S0 on 5/50 trees): a non-PD Hessian here means the optimizer settled in the
    # alternative mode where the non-phylogenetic species intercept SD (~16) absorbs the between-species
    # variance and the phylogenetic SD is ~0.1 (instead of ~0.005 and ~24): the two species-level terms are
    # weakly identified against each other on those trees. Fixed effects, including the sigma year term, are
    # identical to 3-4 decimals in either mode, so the *_conv pools differ from the all-tree pools only in
    # the 6th decimal; the per-tree varcomp table (sd_phylo vs sd_spp_intercept) shows which mode each tree hit.
    conv_trees <- r$varcomp$tree[r$varcomp$converged]
    r$pooled_converged <- if (length(conv_trees)) pool_rubin_df(r$per_tree_fixed %>% filter(tree %in% conv_trees)) else NULL
    phylo_runs[[nm]] <- r
    pooled <- r$pooled %>% mutate(component = recode(component, disp = "sigma", cond = "mu"), tier = nm, label = tiers[[nm]]$label)
    if (!is.null(r$pooled_converged))
      pooled <- pooled %>% left_join(r$pooled_converged %>% mutate(component = recode(component, disp = "sigma", cond = "mu")) %>%
                                       select(component, par, estimate_conv = estimate, se_conv = se, lower_conv = lower, upper_conv = upper, m_conv = m),
                                     by = c("component", "par"))
    phylo_pooled_all[[nm]] <- pooled
    rows <- pooled %>% filter(component == "sigma" | (component == "mu" & par == "scaled_yr")) %>%
      mutate(pct_decade    = ifelse(component == "sigma" & par %in% c("scaled_yr", "yr_within_src", "yr_src_mean"), pct_decade(estimate), NA),
             pct_decade_lo = ifelse(is.na(pct_decade), NA, pct_decade(lower)),
             pct_decade_hi = ifelse(is.na(pct_decade), NA, pct_decade(upper)),
             pct_per_sd_tmean    = ifelse(par == "scaled_tmean", 100 * (exp(estimate) - 1), NA),
             pct_per_sd_tmean_lo = ifelse(par == "scaled_tmean", 100 * (exp(lower) - 1), NA),
             pct_per_sd_tmean_hi = ifelse(par == "scaled_tmean", 100 * (exp(upper) - 1), NA),
             n = nrow(dat_i), n_species = r$n_species, n_trees = r$n_trees, n_converged = r$n_converged,
             secs = r$secs, engine = "glmmTMB propto + dispformula, REML, Rubin-pooled")
    phylo_rows[[nm]] <- rows
    phylo_vc[[nm]] <- r$varcomp %>% mutate(tier = nm, .before = 1)
    phylo_timing[[nm]] <- data.frame(tier = nm, n = nrow(dat_i), n_trees = r$n_trees, n_converged = r$n_converged,
                                     secs_total = r$secs, secs_per_tree = r$secs / r$n_trees)
    yr <- rows %>% filter(component == "sigma", par %in% c("scaled_yr", "yr_within_src")) %>% slice(1)
    message(sprintf("  %-3s sigma~year b = %+.4f (SE %.4f; between-tree var %.1e) => %+.1f %%/decade [%+.1f, %+.1f]; mu year %+.3f (SE %.3f); %d/%d trees pdHess; %.0f s",
                    nm, yr$estimate, yr$se, yr$between_tree_var, yr$pct_decade, yr$pct_decade_lo, yr$pct_decade_hi,
                    rows$estimate[rows$component == "mu"], rows$se[rows$component == "mu"], r$n_converged, r$n_trees, r$secs))
  }
  phylo_table  <- do.call(rbind, phylo_rows);  rownames(phylo_table) <- NULL
  phylo_pooled <- do.call(rbind, phylo_pooled_all); rownames(phylo_pooled) <- NULL
  phylo_vc_all <- bind_rows(phylo_vc); rownames(phylo_vc_all) <- NULL   # S0/S1 lack contributor/site columns
  phylo_timing <- do.call(rbind, phylo_timing)

  # --- like-for-like comparison with the non-phylogenetic tiers (same records) ------
  phylo_cmp <- phylo_table %>%
    select(tier, component, par, est_phylo = estimate, se_phylo = se, lo_phylo = lower, hi_phylo = upper,
           pct_decade_phylo = pct_decade, pct_decade_lo_phylo = pct_decade_lo, pct_decade_hi_phylo = pct_decade_hi, n_trees, n_converged) %>%
    inner_join(res$glmmTMB %>% select(tier, component, par, est_nophylo = estimate, se_nophylo = se, lo_nophylo = lower, hi_nophylo = upper,
                                      pct_decade_nophylo = pct_decade, pct_decade_lo_nophylo = pct_decade_lo, pct_decade_hi_nophylo = pct_decade_hi, n),
               by = c("tier", "component", "par")) %>%
    mutate(diff_est = est_phylo - est_nophylo, se_ratio = se_phylo / se_nophylo,
           excludes_zero_phylo = lo_phylo > 0 | hi_phylo < 0, excludes_zero_nophylo = lo_nophylo > 0 | hi_nophylo < 0)

  # --- headline models: S7 (plan spec, with tmean) and S3 (--no-tmean variant) -------
  headline_row <- function(tier, par) phylo_table %>% filter(tier == !!tier, component == "sigma", par == !!par)
  phylo_headline <- bind_rows(
    if (!is.null(phylo_runs$S7)) headline_row("S7", "scaled_yr") %>% mutate(role = "plan specification (with scaled_tmean), tmean-complete sample"),
    if (!is.null(phylo_runs$S7)) headline_row("S7", "scaled_tmean") %>% mutate(role = "plan specification: temperature term on sigma"),
    if (!is.null(phylo_runs$S3)) headline_row("S3", "scaled_yr") %>% mutate(role = "--no-tmean variant (no scaled_tmean), 8,282-record sample"))

  # --- single-tree (tree 1) fits of the headline models for quantities the pooled
  #     tables do not carry: sigma random-effect SDs, contributor SD multipliers, species slopes
  headline_tree1 <- list()
  for (nm in intersect(c("S3", "S7"), names(phylo_runs))) {
    dat_i <- if (identical(tiers[[nm]]$data, "wing_cc_t")) ph_cc_t else ph_cc
    f1 <- fit_phylo_glmmtmb(tiers[[nm]]$mean, dat_i, A_list[[1]], dispformula = tiers[[nm]]$disp)
    vc1 <- VarCorr(f1)
    re_src <- ranef(f1)$disp$Main_researcher[, 1]
    headline_tree1[[nm]] <- list(
      tree = 1L, converged = isTRUE(f1$sdr$pdHess), n = nobs(f1),
      sigma_re_sd = sapply(vc1$disp, function(x) attr(x, "stddev")[[1]]),
      mu_re_sd    = unlist(lapply(vc1$cond[setdiff(names(vc1$cond), "g")], function(x) attr(x, "stddev"))),
      sd_phylo    = attr(vc1$cond$g, "stddev")[[1]],
      contributor_sd_multiplier_range = round(range(exp(re_src)), 2),
      species_year_slopes = tryCatch(species_slopes_phylo(f1, "scaled_yr", "spp"), error = function(e) conditionMessage(e)),
      disp_coef = summary(f1)$coefficients$disp, cond_coef = summary(f1)$coefficients$cond)
  }

  phylo_secs_total <- as.numeric(Sys.time() - t_phylo0, units = "secs")
  phylo_out <- list(
    generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
    engine = "glmmTMB propto (Williams et al. 2025) + dispformula; REML; Rubin's rules over trees (_phylo_engine.R)",
    glmmTMB_version = as.character(packageVersion("glmmTMB")), R = R.version.string,
    n_trees = length(A_list), tree_source = if (file.exists(A_cache)) "phylo_A_50trees.rds (the 50 trees of the published brms wing model, brm0_multiphylo.rda)" else "clootl fallback (_phylo_engine.R)",
    n_records = c(complete_case = nrow(ph_cc), tmean_complete = if (USE_TMEAN) nrow(ph_cc_t) else NA_integer_),
    n_species = length(ph_species), species_dropped = c(Herpsilochmus_sellowi = n_dropped_sellowi),
    scaled_tmean_used = USE_TMEAN, sd_year_scaling = SD_YR,
    tiers = lapply(phylo_runs, function(r) r[c("label", "formula", "n", "n_species", "n_trees", "n_converged", "secs")]),
    table = phylo_table,                 # sigma pars + mu year per tier (mirrors variance_results.rds$glmmTMB)
    pooled_all = phylo_pooled,           # every pooled fixed effect (mu and sigma) per tier; *_conv = pooled over pdHess trees only
    per_tree_fixed = do.call(rbind, lapply(names(phylo_runs), function(nm) cbind(tier = nm, phylo_runs[[nm]]$per_tree_fixed))),
    varcomp = phylo_vc_all,              # per tree: sd_phylo, spp/contributor/site SDs, exp(b_sigma_int), phylo shares, pdHess
    varcomp_summary = phylo_vc_all %>% group_by(tier) %>%
      summarise(across(c(sd_phylo, sigma_baseline, phylo_share_incl_resid, phylo_share_between_species, sd_spp_intercept, sd_spp_year_slope),
                       list(mean = ~mean(.x, na.rm = TRUE), lo = ~quantile(.x, .025, na.rm = TRUE), hi = ~quantile(.x, .975, na.rm = TRUE))),
                n_converged = sum(converged), n_trees = n(), .groups = "drop"),
    comparison = phylo_cmp,              # phylogenetic vs non-phylogenetic, same tier, same records
    headline = phylo_headline,
    headline_tree1 = headline_tree1,
    timing = phylo_timing, secs_total = phylo_secs_total,
    definitions = c(
      sigma_pct_decade = res$definitions[["sigma_pct_decade"]],
      pct_per_sd_tmean = "100*(exp(b_sigma_scaled_tmean) - 1): % change in residual SD per SD of record-year mean temperature",
      phylo_share_incl_resid = "sd_phylo^2 / (sd_phylo^2 + sd_spp_int^2 + sd_spp_yr^2 + sd_contributor^2 + sd_site^2 + exp(b_sigma_Intercept)^2), per tree (comparable with atlantic_parallel.R's phylogenetic signal; the residual SD is the sigma-model baseline)",
      phylo_share_between_species = "sd_phylo^2 / (sd_phylo^2 + sd_spp_int^2): share of the between-species intercept variance that is phylogenetic",
      converged = "glmmTMB positive-definite Hessian (sdr$pdHess) per tree; Rubin pooling (estimate, se, lower, upper) uses all trees; *_conv columns pool the pdHess trees only (m_conv of them)",
      intervals = "estimate +/- 1.96 * Rubin SE (within-tree variance + (1 + 1/m) * between-tree variance)"),
    note = paste("The species random slope (1 + scaled_yr || spp) stays in the mean model alongside the phylogenetic intercept,",
                 "as in atlantic_parallel.R; the non-phylogenetic species intercept SD collapses towards 0 because the phylogenetic",
                 "term absorbs it (phylogenetic signal ~0.95 in the published wing model)."))
  saveRDS(phylo_out, out_path("variance_phylo_results.rds"))
  message(sprintf("wrote %s (%d tiers, %d trees, %.1f min)", out_path("variance_phylo_results.rds"), length(phylo_runs), length(A_list), phylo_secs_total / 60))
  # pointer + headline in the fast-tier file (fast-tier slots untouched)
  res$phylo <- list(file = "variance_phylo_results.rds", generated = phylo_out$generated, engine = phylo_out$engine,
                    n_trees = phylo_out$n_trees, tree_source = phylo_out$tree_source, n_records = phylo_out$n_records,
                    headline = phylo_headline, table = phylo_table, comparison = phylo_cmp, timing = phylo_timing, secs_total = phylo_secs_total)
  message("phylogenetic vs non-phylogenetic sigma year effect (%/decade), same records:")
  print(phylo_cmp %>% filter(component == "sigma", par %in% c("scaled_yr", "yr_within_src", "yr_src_mean", "scaled_tmean")) %>%
          select(tier, par, est_nophylo, est_phylo, se_nophylo, se_phylo, pct_decade_nophylo, pct_decade_phylo, n_converged), digits = 3)
}

# ===========================================================================
# Save fast results (before any Stan work) and draw the figures
# ===========================================================================
res$index <- c(
  glmmTMB = "sigma-model ladder S0-S6 on 8,282 records; S3t/S7 on the tmean-complete subset: tier, component (mu/sigma), par, estimate, se, pct_decade, pct_per_sd_tmean",
  scaled_tmean = "whether scaled_tmean was used, coverage, n, NA by period, cor(year, tmean)",
  two_stage = "lme4 log|residual| cross-check",
  lncvr = "all lnCVR / lnVR meta-analyses (analysis, model, k, estimate, lower, upper, I2, estimate_pct)",
  lncvr_by_contributor = "per-contributor k / mean lnCVR at each min-n threshold",
  contributors_per_cell_summary = "mean/median contributors per species-period cell (wing records)",
  variance_decomposition_summary = "between-contributor share of within-cell variance by period",
  cv_by_period_min5 = "median CV pooled vs within-contributor, n >= 5",
  lncvr_metaregression = "species lnCVR vs change in log contributors / sites",
  phylo = "glmmTMB propto engine (default --engine): pointer to variance_phylo_results.rds + headline (S7 with tmean, S3 without) + phylo-vs-no-phylo comparison of every tier",
  brms = "Rubin-pooled brms distributional model (only after --engine brms --trees N)",
  brms_smoke = "brms smoke-test summary (NOT a result)")
saveRDS(res, out_path("variance_results.rds"))
message("wrote ", out_path("variance_results.rds"), " (fast parts)")

# ---- figures ---------------------------------------------------------------
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))
# A. lnCVR forest across analysis levels (rma rows, lnCVR only)
fa <- res$lncvr %>% filter(!grepl("lnVR|EXACT", analysis), model == "rma", k >= 2) %>%
  mutate(analysis = factor(analysis, levels = rev(unique(analysis))),
         lab = sprintf("k = %d", k),
         x_lab = max(upper, na.rm = TRUE) + 0.05)
pA <- ggplot(fa, aes(x = estimate, y = analysis)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
  geom_point(size = 2.2) +
  geom_text(aes(label = lab, x = x_lab), hjust = 0, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = "lnCVR, late (2013-2018) vs early (1990-2006)  [95% CI]", y = NULL,
       title = "A. Wing-length CV change at three levels of aggregation")
# B. sigma-model year effect as % change in residual SD per decade
fb <- res$glmmTMB %>% filter(component == "sigma", par %in% c("scaled_yr", "yr_within_src", "yr_src_mean")) %>%
  mutate(lab = ifelse(par == "scaled_yr", tier, paste0(tier, ifelse(par == "yr_within_src", " within-contributor", " between-contributor"))),
         src = "glmmTMB dispformula") %>%
  select(lab, pct_decade, pct_decade_lo, pct_decade_hi, src)
fb2 <- res$two_stage %>% transmute(lab = tier, pct_decade, pct_decade_lo, pct_decade_hi, src = "lme4 log|residual|")
fb <- bind_rows(fb, fb2)
if (!is.null(res$phylo)) {   # phylogenetic tier (glmmTMB propto, Rubin-pooled over trees): same tiers, "+ phylo" rows
  fb3 <- res$phylo$table %>% filter(component == "sigma", par %in% c("scaled_yr", "yr_within_src", "yr_src_mean")) %>%
    mutate(lab = paste0(ifelse(par == "scaled_yr", tier, paste0(tier, ifelse(par == "yr_within_src", " within-contributor", " between-contributor"))),
                        " + phylogeny"),
           src = sprintf("glmmTMB dispformula + phylogeny (%d trees)", res$phylo$n_trees)) %>%
    select(lab, pct_decade, pct_decade_lo, pct_decade_hi, src)
  # interleave: each phylogenetic row directly under its non-phylogenetic tier
  ord <- unlist(lapply(fb$lab[fb$src == "glmmTMB dispformula"], function(l) c(l, paste0(l, " + phylogeny"))))
  fb <- bind_rows(fb, fb3) %>% mutate(lab = factor(lab, levels = rev(c(ord, fb2$lab))))
} else {
  fb <- fb %>% mutate(lab = factor(lab, levels = rev(unique(lab))))
}
pB <- ggplot(fb, aes(x = pct_decade, y = lab, colour = src)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = pct_decade_lo, xmax = pct_decade_hi), width = 0.2, orientation = "y") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = setNames(c("#1b6ca8", "#c0504d", "#2e8b57"),
                                        c("glmmTMB dispformula", "lme4 log|residual|",
                                          sprintf("glmmTMB dispformula + phylogeny (%d trees)", if (!is.null(res$phylo)) res$phylo$n_trees else N_TREES))),
                      name = NULL) +
  guides(colour = guide_legend(nrow = if (!is.null(res$phylo)) 2 else 1)) +
  labs(x = "% change in residual SD of wing length per decade  [95% CI]", y = NULL,
       title = "B. Continuous sigma models: year effect on residual SD") +
  theme(legend.position = "bottom")
# C. contributors per species-period cell
pC <- ggplot(cpc, aes(x = n_contrib, fill = period)) +
  geom_histogram(binwidth = 1, position = "dodge", colour = "white") +
  scale_fill_manual(values = c(early = "#7fa8c9", late = "#d9822b"), name = NULL) +
  labs(x = "contributors per species x period cell (wing records)", y = "species",
       title = sprintf("C. Contributors per cell: mean %.1f early -> %.1f late",
                       res$contributors_per_cell_summary$mean_contrib[1], res$contributors_per_cell_summary$mean_contrib[2])) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85))
# D. between-contributor share of within-cell variance by period
pD <- ggplot(decomp, aes(x = period, y = share_between, fill = period)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8) +
  scale_fill_manual(values = c(early = "#7fa8c9", late = "#d9822b"), guide = "none") +
  labs(x = NULL, y = "between-contributor share of within-cell variance",
       title = "D. Between-contributor share of within-cell variance")
fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(if (!is.null(res$phylo)) 1.6 else 1.2, 1))
ggsave(fig_path("variance_summary.png"), fig, width = 15, height = if (!is.null(res$phylo)) 12 else 10, dpi = 200, bg = "white")

# per-cell within-contributor lnCVR forest (n >= 5)
wc5 <- src_tables[["5"]] %>% mutate(cell = paste(gsub("_", " ", Binomial), Sex, Main_researcher, sep = " / "),
                                    se = sqrt(s2_lnCVR)) %>% arrange(Main_researcher, lnCVR) %>%
  mutate(cell = factor(cell, levels = unique(cell)))
pE <- ggplot(wc5, aes(x = lnCVR, y = cell, colour = Main_researcher)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lnCVR - 1.96 * se, xmax = lnCVR + 1.96 * se), width = 0.2, orientation = "y") +
  geom_point(size = 2) +
  labs(x = "lnCVR late vs early within species x sex x contributor (n >= 5 in both periods)", y = NULL, colour = "contributor",
       title = sprintf("Within-contributor lnCVR, k = %d cells, %d contributors", nrow(wc5), n_distinct(wc5$Main_researcher)))
ggsave(fig_path("variance_lncvr_contributor.png"), pE, width = 10, height = 6, dpi = 200, bg = "white")
message("wrote figures/variance_summary.png, figures/variance_lncvr_contributor.png")

# ===========================================================================
# 1b. brms distributional model with phylogeny — Stan engine, kept as the Bayesian
#     cross-check (--engine brms --trees N on Totoro; --smoke locally)
# ===========================================================================
if (!RUN_BRMS) {
  message("\n[1b] brms distributional model SKIPPED (engine = ", ENGINE,
          if (RUN_PHYLO_TMB) "; the phylogenetic tier ran in glmmTMB above" else "; use --engine brms --smoke locally or --engine brms --trees N on the server", ").")
  quit(save = "no", status = 0)
}

suppressMessages({ library(brms); library(ape); library(MCMCglmm); library(prepR4pcm); library(phytools)
  library(future.apply); library(posterior) })
message("\n[1] brms distributional model: ", if (SMOKE) "SMOKE TEST (1 tree, 2 chains, 400 iter)" else sprintf("%d trees", N_TREES))

# --- phylogeny: same recipe as atlantic_parallel.R -----------------------
# tmean-complete subset when scaled_tmean enters sigma (plan specification); else the full complete-case sample
dat <- as.data.frame(if (USE_TMEAN) wing_cc_t else wing_cc)
message(sprintf("brms data: %d records (%s)", nrow(dat), if (USE_TMEAN) "tmean-complete subset" else "complete-case sample, no tmean"))
dat$species_name <- gsub("_", " ", dat$Binomial)
ebird_synonyms <- c(
  "Antilophia galeata"       = "Chiroxiphia galeata",
  "Tachyphonus cristatus"    = "Loriotus cristatus",
  "Pyrrhocoma ruficeps"      = "Thlypopsis pyrrhocoma",
  "Pyriglena pernambucensis" = "Pyriglena leuconota",
  "Tangara sayaca"           = "Thraupis sayaca",
  "Tangara cayana"           = "Stilpnia cayana",
  "Tangara palmarum"         = "Thraupis palmarum",
  "Tangara peruviana"        = "Stilpnia peruviana",
  "Dixiphia pipra"           = "Pseudopipra pipra",
  "Tiaris fuliginosus"       = "Asemospiza fuliginosa")
hit <- dat$species_name %in% names(ebird_synonyms)
dat$species_name[hit] <- ebird_synonyms[dat$species_name[hit]]
spp_data <- unique(dat$species_name)
# offline-first: point clootl at the repo copy of AvesData if the env var is unset
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
  if (dir.exists(raw_path("AvesDataLite-main"))) Sys.setenv(AVESDATA_PATH = raw_path("AvesDataLite-main"))
  else clootl::get_avesdata_repo(path = raw_path())
}
.check <- pr_get_tree(spp_data, source = "clootl", n_tree = 1)
if (length(.check$unmatched) > 0) {
  message("Removing ", length(.check$unmatched), " species absent from eBird taxonomy: ", paste(.check$unmatched, collapse = ", "))
  spp_data <- setdiff(spp_data, .check$unmatched)
}
got <- pr_get_tree(spp_data, source = "clootl", n_tree = 100, cache = TRUE)
trees <- got$tree
clootl_version <- pr_cite_tree(got, format = "text")
rec <- reconcile_tree(x = dat, tree = trees[[1]], x_species = "species_name", fuzzy = TRUE, resolve = "flag")
print(reconcile_summary(rec))
aligned <- reconcile_apply(rec, data = dat, tree = trees[[1]], species_col = "species_name", drop_unresolved = TRUE)
dat <- aligned$data
dat$species_name <- gsub(" ", "_", dat$species_name)
dat$spp <- dat$species_name
norm_us <- function(x) gsub(" ", "_", x)
keep_tips <- norm_us(aligned$tree$tip.label)
trees_pruned <- lapply(trees, function(t) { t$tip.label <- norm_us(t$tip.label); ape::keep.tip(t, intersect(keep_tips, t$tip.label)) })
class(trees_pruned) <- "multiPhylo"
tree_samp <- if (N_TREES >= length(trees_pruned)) trees_pruned else sample(trees_pruned, N_TREES)
res$brms_n <- nrow(dat); res$brms_n_species <- length(unique(dat$spp))
message(sprintf("brms sample: %d records, %d species", nrow(dat), length(unique(dat$spp))))

make_A <- function(tree) {
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

# --- model ---------------------------------------------------------------
sigma_rhs <- if (USE_TMEAN) {
  "scaled_yr + scaled_tmean + Sex + (1 | Main_researcher) + (1 | Municipality)"
} else {
  "scaled_yr + Sex + (1 | Main_researcher) + (1 | Municipality)"
}
bform <- bf(as.formula(paste(
  "conc.wing.length ~ 1 + Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season +",
  "(1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)) + (1 | Main_researcher) + (1 | Municipality)")),
  as.formula(paste("sigma ~", sigma_rhs)))
res$brms_formula <- c(mu = deparse1(bform$formula), sigma = deparse1(bform$pforms$sigma))
# Priors: mean part as in atlantic_parallel.R; sigma part weakly informative on the log-SD scale
# (residual SD of wing ~ 3-4 mm -> log ~ 1.3; the glmmTMB S3 intercept is ~1.4).
priors <- c(prior(normal(71, 15), class = Intercept),
            prior(normal(0, 10),  class = b),
            prior(cauchy(0, 1),   class = sd),
            prior(normal(1.4, 1), class = Intercept, dpar = sigma),
            prior(normal(0, 1),   class = b,  dpar = sigma),
            prior(cauchy(0, 1),   class = sd, dpar = sigma))

if (SMOKE) {
  SAMPLING$chains <- 2L; SAMPLING$iter <- 400L; SAMPLING$warmup <- 200L; SAMPLING$cores <- 2L
  SAMPLING$workers <- 1L
  SAMPLING_CONTROL <- list(adapt_delta = 0.9, max_treedepth = 10)   # smoke test only: speed, not inference
}
fit_one <- function(tree) {
  A <- make_A(tree)
  brm(bform, data = dat, data2 = list(A = A), family = gaussian(), prior = priors,
      iter = SAMPLING$iter, warmup = SAMPLING$warmup, chains = SAMPLING$chains,
      cores = SAMPLING$cores, backend = SAMPLING$backend, control = SAMPLING_CONTROL,
      seed = 20260909, refresh = if (SMOKE) 50 else 0)
}

t_start <- Sys.time()
if (SMOKE || SAMPLING$workers == 1L) {
  fits <- lapply(tree_samp, fit_one)
} else {
  plan(multisession, workers = SAMPLING$workers)
  fits <- future_lapply(tree_samp, fit_one, future.seed = TRUE)
}
elapsed_min <- as.numeric(Sys.time() - t_start, units = "mins")
message(sprintf("sampling finished in %.1f min", elapsed_min))

# --- Rubin's rules across trees (as in atlantic_parallel.R), all b_ and sd_ pars ----
pool_rubin <- function(fits, pars) {
  m <- length(fits)
  draws <- lapply(fits, function(f) as_draws_df(f)[, pars, drop = FALSE])
  means <- t(sapply(draws, function(d) sapply(d, mean)))
  vars  <- t(sapply(draws, function(d) sapply(d, var)))
  if (m == 1) { means <- matrix(means, nrow = 1, dimnames = list(NULL, pars)); vars <- matrix(vars, nrow = 1, dimnames = list(NULL, pars)) }
  qbar <- colMeans(means); ubar <- colMeans(vars)
  b    <- if (m > 1) apply(means, 2, var) else 0 * qbar
  tot  <- ubar + (1 + 1/m) * b
  se   <- sqrt(tot)
  data.frame(par = pars, estimate = qbar, se = se, lower = qbar - 1.96*se, upper = qbar + 1.96*se,
             between_tree_var = b, within_tree_var = ubar, row.names = NULL)
}
all_pars <- variables(fits[[1]])
pars <- all_pars[grepl("^b_|^sd_", all_pars)]
rubin_summary <- pool_rubin(fits, pars)
rubin_summary$pct_decade <- ifelse(rubin_summary$par == "b_sigma_scaled_yr", pct_decade(rubin_summary$estimate), NA)
rubin_summary$pct_decade_lo <- ifelse(rubin_summary$par == "b_sigma_scaled_yr", pct_decade(rubin_summary$lower), NA)
rubin_summary$pct_decade_hi <- ifelse(rubin_summary$par == "b_sigma_scaled_yr", pct_decade(rubin_summary$upper), NA)
rubin_summary$pct_per_sd_tmean <- ifelse(rubin_summary$par == "b_sigma_scaled_tmean", 100 * (exp(rubin_summary$estimate) - 1), NA)
print(rubin_summary[grepl("scaled_yr|Sex|tmean", rubin_summary$par), ], digits = 3)

# phylogenetic signal (share of the between-species + residual variance at the sigma intercept)
phylo_sig <- sapply(fits, function(f) {
  d <- as_draws_df(f)
  num <- d$sd_species_name__Intercept^2
  den <- num + d$sd_spp__Intercept^2 + d$sd_spp__scaled_yr^2 + d$sd_Main_researcher__Intercept^2 +
    d$sd_Municipality__Intercept^2 + exp(d$b_sigma_Intercept)^2
  mean(num / den)
})
# convergence diagnostics per fit
diag <- do.call(rbind, lapply(seq_along(fits), function(i) {
  s <- summary(fits[[i]])$fixed
  data.frame(tree = i, max_rhat = max(c(s$Rhat, summary(fits[[i]])$spec_pars$Rhat), na.rm = TRUE),
             min_bulk_ess = min(s$Bulk_ESS, na.rm = TRUE), divergent = sum(rstan::get_divergent_iterations(fits[[i]]$fit)))
}))
brms_out <- list(rubin_summary = rubin_summary, phylo_signal = c(mean = mean(phylo_sig),
                 lwr = as.numeric(quantile(phylo_sig, .025)), upr = as.numeric(quantile(phylo_sig, .975))),
                 diagnostics = diag, n = nrow(dat), n_species = length(unique(dat$spp)), n_trees = length(fits),
                 sampling = SAMPLING[c("chains", "iter", "warmup", "backend")], control = SAMPLING_CONTROL,
                 elapsed_min = elapsed_min, clootl_version = clootl_version, scaled_tmean_used = USE_TMEAN,
                 formula = res$brms_formula)

if (SMOKE) {
  brms_out$WARNING <- "SMOKE TEST: 1 tree, 2 chains x 400 iterations, adapt_delta 0.9. Checks that the model compiles and samples. NOT A RESULT - do not quote."
  res$brms_smoke <- brms_out
  save(fits, rubin_summary, clootl_version, file = out_path("models", "variance_sigma_SMOKE.rda"))
  message("SMOKE TEST complete; saved output/models/variance_sigma_SMOKE.rda (not a result)")
} else {
  res$brms <- brms_out
  save(fits, rubin_summary, clootl_version, file = out_path("models", "variance_sigma_multiphylo.rda"))
  message("saved output/models/variance_sigma_multiphylo.rda")
}
saveRDS(res, out_path("variance_results.rds"))
message("wrote ", out_path("variance_results.rds"), " (with brms slot)")
