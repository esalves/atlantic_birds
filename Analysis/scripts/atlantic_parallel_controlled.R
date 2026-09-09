# atlantic_parallel_controlled.R
# ---------------------------------------------------------------------------
# WHAT: Re-fits the primary wing-length model (atlantic_parallel.R) with the
# artefact controls the referee asked for (REVISION_PLAN.md §3 Phase 1, "THE
# DECISION GATE"): contributor and site intercepts, the wing-column protocol
# proxy, longitude/altitude, season and moult, an individual (ring) intercept,
# contributor-specific year slopes, a Mundlak within/between decomposition of
# year (within municipality; also within contributor) and latitude (within
# species), an unknown-sex replication and a first-captures-only sensitivity.
# Writes a tidy before/after table of the scaled-year coefficient, a
# leave-one-contributor-out jackknife, per-contributor year slopes for the six
# contributors that span the record, and the species-specific slopes for the
# caterpillar figure.
#
# WHY: only 6 of 42 contributors and 15 of 99 municipalities span both the
# early (<= 2006) and late (>= 2013) periods, and the wing column used flips
# 87 % right-wing -> 67 % unspecified (REVISION_NOTES_P0.md). A year slope that
# does not survive these controls is a sampling-provenance contrast, not a
# phenotypic trend. The outcome of M3/M5 chooses between Scenario A and B in
# REVISION_PLAN.md §2.
#
# >>> THREE-TIER COMPUTATIONAL DIRECTIVE (REVISION_PLAN.md §3 Phase 1) <<<
#   Tier 1  --fast-lme4 (default)  lme4 REML, no phylogeny, all models (M0..M7
#           plus the finer ladder used in REVISION_PLAN.md §0), on BOTH samples
#           (see SAMPLES). Under two minutes. Sets the decision gate provisionally.
#   PHYLOGENETIC TIER  --trees N [--engine glmmTMB] (default engine for --trees)
#           glmmTMB REML with the phylogenetic species intercept as a `propto`
#           covariance structure (Williams, McGillycuddy, Drobniak, Bolker, Warton
#           & Nakagawa 2025, bioRxiv 10.64898/2025.12.20.695312), fitted on the
#           IDENTICAL 50 trees as the published brms analysis (correlation matrices
#           cached from brm0_multiphylo.rda by _phylo_engine.R) and Rubin-pooled.
#           glmmtmb_validation_wing.R shows this reproduces the published brms M0
#           to 2-3 decimals (year -0.922 [-1.400, -0.444] vs -0.922 [-1.403, -0.440])
#           in ~0.5 s per fit, so EVERY Tier-1 specification (39 specs, both
#           samples) is carried across the trees, not just M0/M3/M5/M6.
#           --models all (default) or a comma list, e.g. --models M0,M3,M5,M6.
#           Writes controlled_wing_phylo_results.rds / _species_slopes.rds; the
#           decision gate in that file is evaluated on the phylogenetic estimates.
#   Tier 2  --trees 1 --engine brms   brms with the phylogenetic term on ONE tree,
#           full chains/iterations from _sampling_config.R, to check that the
#           Bayesian posteriors agree with REML (R-hat, bulk-ESS). ~1-2 h/model
#           on the server. Kept as the Bayesian cross-check (Totoro).
#   Tier 3  --trees 50 --engine brms  Rubin-pooled 50-tree brms fits ONLY for the
#           definitive models M0, M3, M5, M6 (the default brms model set;
#           override with --models M0,M3). 8 models x 50 trees = 400 Stan fits
#           would take hundreds of CPU hours and is NOT to be attempted.
#   --sample full|cc               which sample the brms tiers use (default full;
#           ignored by --engine glmmTMB, which fits both samples).
#   --smoke                        LOCAL SMOKE TEST of the brms code path only:
#           1 tree, chains = 2, iter = 400, warmup = 200, model M3 only,
#           written to *_smoke.* files. Its numbers are NOT results.
#
# SAMPLES. Altitude is NA for 190 wing records and the capture date for 6.
#   "cc"   complete cases (8,282), the sample REVISION_PLAN.md §0 used.
#   "full" all 8,478 wing records: altitude imputed from the nearest coordinate
#          site with a recorded altitude (see DECISIONS), season "unrecorded" for
#          the 6 undated records.
#   The 190 altitude-NA records are not random: 103 are E. Carrano's 1996-1999
#   Ilha Rasa (Guaraqueçaba) records, i.e. the early anchor of the one
#   contributor who spans 1996-2017 with 1,040 wing records. Dropping them
#   removes within-contributor temporal contrast, so "full" is the primary
#   sample and "cc" the plan-comparable one; both are reported.
#
# MODELS (fixed part shown; every model keeps (1 + scaled_yr || spp), i.e.
# uncorrelated species intercepts and year slopes, matching the brms structure;
# brms models add (1 | gr(species_name, cov = A)) for the phylogeny):
#   M0  wing ~ Sex + scaled_yr + scaled_lat                       [baseline]
#   M1  M0 + (1 | src) + (1 | site)                    src = Main_researcher,
#                                                      site = Municipality
#   M2  M1 + wing_col + scaled_lon + scaled_alt        wing_col = right/left/generic
#   M3  M2 + season + molt3 + (1 | ind)                ind = Ring x Binomial
#   M4  M3 with (1 + scaled_yr || src)                 contributor year slopes
#   M5  M3 with Mundlak: yr_within_site + yr_site_mean,
#                        lat_within_spp + lat_spp_mean; species random slope on
#                        yr_within_site (REWB, Bell & Jones 2015 - the plan's spec)
#   M6  M3 on UNKNOWN-sex records (passer90_allsex.rda, !known_sex), no Sex
#   M7  M3 on first captures (first WING record per Ring x Binomial); the ring
#       term is dropped because every individual then has one record
# Finer ladder for the comparison with REVISION_PLAN.md §0 and for attribution:
#   M1a_src (M0 + src), M1b_site (M0 + site), M2a_wingcol (M1 + wing_col),
#   M3a_season (M2 + season), M3b_ring (M3a + ind; = M3 without moult),
#   M5_rsTotal / M5src_rsTotal (Mundlak with the species random slope on TOTAL
#   scaled_yr instead of the within component - the estimate is sensitive to this),
#   M5src (year decomposed within/between CONTRIBUTOR, REWB),
#   M6_m0 / M6_m1 (unknown-sex records with the M0 / M1 structure),
#   M3_longrun (M3 on the 6 long-running contributors: >= 8 yr, >= 200 wing
#   records), M_carrano (M3 without src on E. Carrano's records alone),
#   M0_cc / M1_cc (M0 / M1 on the complete-case sample: the sample-restriction
#   effect on its own). "_cc" suffix = complete-case sample; no suffix = full.
#   *_noanom: M1 / M3 / M5src refitted without contributor x municipality x year
#   blocks (n >= 10) whose species x sex-centred mean wing deviates by more than
#   +/- 6 mm from the species means (a sign-blind protocol / data-entry screen;
#   the flagged blocks are listed in $anomalous_blocks).
#   Also written: a leave-one-contributor-out jackknife of M1 / M3 / M3_cc and
#   per-contributor year slopes for every contributor with >= 8 sampling years.
#
# INPUTS:  data/derived/passer90.rda, data/derived/passer90_allsex.rda,
#          scripts/_phylo_engine.R + data/derived/phylo_A_50trees.rds (glmmTMB tier;
#            the cache is built from output/models/brm0_multiphylo.rda on first use),
#          scripts/_sampling_config.R (brms tiers), AvesDataLite (brms tiers)
# OUTPUTS: output/controlled_wing_results.rds          (Tier 1; list, see $table,
#            $before_after, $jackknife, $spanning_contributor_slopes, $decision_gate)
#          output/controlled_wing_species_slopes.rds   (Tier 1; M3 and M0 slopes)
#          output/controlled_wing_phylo_results.rds    (glmmTMB phylogenetic tier;
#            $before_after, $table, $varcomp, $decision_gate, $comparison_tier1,
#            $comparison_published_M0, $engine, $timings)
#          output/controlled_wing_phylo_species_slopes.rds (glmmTMB tier; $M3, $M3_cc,
#            $M0 pooled over trees, same columns as the Tier-1 file + tree_name)
#          output/controlled_wing_brms_results.rds     (brms Tiers 2-3; Rubin-pooled)
#          output/models/controlled_<model>.rda        (brms Tiers 2-3; per-tree fits)
#          output/controlled_wing_brms_smoke.rds, output/models/controlled_smoke_M3.rda
#                                                      (--smoke only)
# RUN:  Rscript Analysis/scripts/atlantic_parallel_controlled.R --fast-lme4
#       Rscript Analysis/scripts/atlantic_parallel_controlled.R --trees 50            # glmmTMB, all specs (~1 h)
#       Rscript Analysis/scripts/atlantic_parallel_controlled.R --trees 50 --models M0,M3,M5,M6
#       Rscript Analysis/scripts/atlantic_parallel_controlled.R --trees 1 --engine brms
#       Rscript Analysis/scripts/atlantic_parallel_controlled.R --trees 50 --engine brms --models M0,M3 --sample cc
#       Rscript Analysis/scripts/atlantic_parallel_controlled.R --smoke
#       (works from the repo root, Analysis/, or Analysis/scripts/)
#
# NON-OBVIOUS DECISIONS (also in REVISION_NOTES_P1.md):
#   * Altitude imputation ("full" sample): for each record without Altitude the
#     median altitude of the NEAREST coordinate site (great-circle distance)
#     that has one, using every record in passer90_allsex. The municipality
#     median was rejected: E. Carrano's Guaraqueçaba records are from Ilha Rasa
#     (an island) while the only Guaraqueçaba altitude on file (563 m) is
#     R. Bobato's inland farm. Distances are stored in $sample_sizes$altitude_imputation.
#     Altitude has a negligible coefficient and municipality is also a random
#     intercept, so the imputation cannot drive the year effect; keeping the
#     records can.
#   * Moult is a 3-level factor (no / Yes / unrecorded) rather than a filter:
#     Molt is NA for 43 % of early vs 13 % of late wing records, so dropping
#     NA would delete a sixth of the sample non-randomly. M3b_ring (= M3
#     without moult) isolates the moult term's contribution.
#   * Individual = interaction(Ring, Binomial) (36 ring strings recur on > 1
#     species); unringed birds get a unique singleton level (their own record
#     id), so they do not collapse into one spurious "NA" individual.
#   * First captures (M7) = first WING-measured record per Ring x Binomial
#     (ordered Year, month, day), unringed kept: 7,791 records.
#     REVISION_NOTES_P0.md's 7,734 takes the first record of ANY kind and then
#     requires a wing value; for a wing model the first wing measurement is the
#     relevant first capture.
#   * Mundlak / REWB: the species random slope is on the same within-cluster
#     year variable whose fixed effect is being estimated (the plan's spec).
#     With the random slope on total scaled_yr instead (the _rsTotal variants)
#     the within-contributor slope moves from ~0 to ~-0.27 and its SE halves;
#     both are reported because the choice matters.
#   * mm_per_decade uses SD(year) = 5.0199 (the SD scaled_yr was actually
#     standardised with, attr(passer90, "scaling")$year), not the 5.05 of the
#     final sample (REVISION_NOTES_P0.md §2 effect_scale).
#   * lme4 intervals are Wald (estimate +/- 1.96 SE) and REML. lme4 has no
#     phylogenetic term: the species intercept absorbs it (the plan's baseline
#     reproduced the brms -0.92). The glmmTMB phylogenetic tier adds the
#     propto(0 + species_name | g, A) intercept on each of the 50 trees; its
#     intervals are Rubin-pooled Wald intervals (pooled SE = sqrt(mean within-tree
#     SE^2 + (1 + 1/m) between-tree variance)); `z` = estimate / pooled SE (the
#     `t` column repeats `z` so Tier-1 consumers of $before_after keep working).
#     Species names in that tier are the tree tip names (phylo_species_name();
#     e.g. Tiaris_fuliginosus -> Asemospiza_fuliginosa); `spp` (the species
#     random-slope grouping) uses the same names, a 1:1 relabelling of Binomial.
#     Species absent from the trees (none in passer90 / passer90_allsex as of
#     2026-09-09; Herpsilochmus_sellowi would be) are dropped with a message and
#     listed in $engine$dropped_species. Mundlak means are re-derived after the
#     relabelling (identical here because no species is dropped).
#   * Species slopes in the glmmTMB tier: per tree slope_j = fixed + BLUP (SE in
#     quadrature, species_slopes_phylo()); across trees slope = mean, SE =
#     sqrt(mean(SE^2) + (1 + 1/m) var(slope)) (Rubin), same for ranef and its SE.
#   * Informed restarts (glmmTMB tier): every spec has both an i.i.d. species
#     intercept and the phylogenetic one (as the published brms model). From
#     glmmTMB's default start a tree-dependent subset of fits stops on the i.i.d.
#     side of that ridge (phylo SD ~ 0.1, non-PD Hessian, REML objective ~24 units
#     WORSE than the phylogenetic side; year estimate identical to 3 decimals).
#     Non-PD fits are restarted from the tree's own parameter layout with the
#     propto scale at SD 20 and accepted only if PD and not worse in objective;
#     $before_after carries converged (final), converged_default_start, n_retried,
#     n_retry_accepted, and $per_model[[m]]$retry the per-tree objectives.
# Session (Tier 1 run, 2026-09-09): R 4.6.0, lme4 2.0.1, dplyr 1.2.1.
# glmmTMB tier (2026-09-09): glmmTMB 1.1.14 (CRAN; `propto` available), TMB per session.
# brms tiers: brms 2.23.0, rstan 2.32.7 (local) / cmdstanr (server),
# prepR4pcm 0.5.0.9000, clootl 0.1.4, phytools 2.5.2, MCMCglmm 2.36, ape 5.8.1.
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(lme4) })
set.seed(20240101)

# --- command-line flags ------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
MODE <- if (!is.null(get_flag("--smoke"))) "smoke" else
        if (!is.null(get_flag("--trees"))) "trees" else "fast"
N_TREES <- if (MODE == "trees") as.integer(get_flag("--trees")) else if (MODE == "smoke") 1L else 0L
if (MODE == "trees" && (is.na(N_TREES) || N_TREES < 1L)) stop("--trees needs a positive integer")
# --engine: glmmTMB (default for --trees; propto phylogenetic tier) or brms (Tiers 2-3; --smoke)
ENGINE <- if (MODE == "smoke") "brms" else if (MODE == "trees") match.arg(get_flag("--engine", "glmmTMB"), c("glmmTMB", "brms")) else "lme4"
BRMS_MODELS <- if (MODE == "smoke") "M3" else
               strsplit(get_flag("--models", if (ENGINE == "glmmTMB") "all" else "M0,M3,M5,M6"), ",")[[1]]
BRMS_SAMPLE <- match.arg(get_flag("--sample", "full"), c("full", "cc"))
message(sprintf("[controlled] mode = %s%s", MODE,
                if (MODE != "fast") sprintf(", engine = %s, trees = %d, models = %s, sample = %s", ENGINE, N_TREES,
                                            paste(BRMS_MODELS, collapse = "/"),
                                            if (ENGINE == "glmmTMB") "full + cc (all specs)" else BRMS_SAMPLE) else ""))

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
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
dir.create(out_path("models"), showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# 1. Data: analysis frames shared by every tier
# ===========================================================================
load(derived_path("passer90.rda"))          # known-sex, 73 spp, 12,571 records
load(derived_path("passer90_allsex.rda"))   # + unknown-sex records, same centring
scaling <- attr(passer90, "scaling")
SD_YEAR <- unname(scaling$year["scale"])    # 5.0199: the SD scaled_yr is per
stopifnot(abs(SD_YEAR - 5.02) < 0.01)
ALT_CENTRE <- unname(scaling$alt["center"]); ALT_SCALE <- unname(scaling$alt["scale"])

# --- altitude imputation table: nearest coordinate site with a recorded altitude ---
haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(sqrt(pmin(1, a)))
}
alt_sites <- passer90_allsex %>% filter(!is.na(Altitude), !is.na(Latitude_decimal_degrees)) %>%
  group_by(Latitude_decimal_degrees, Longitude_decimal_degrees) %>%
  summarise(alt_site = median(Altitude), .groups = "drop")
na_sites <- passer90_allsex %>% filter(is.na(Altitude), !is.na(Latitude_decimal_degrees)) %>%
  distinct(Latitude_decimal_degrees, Longitude_decimal_degrees, Municipality, Locality, Main_researcher)
na_sites$alt_imp <- NA_real_; na_sites$dist_km <- NA_real_
for (i in seq_len(nrow(na_sites))) {
  dk <- haversine_km(na_sites$Latitude_decimal_degrees[i], na_sites$Longitude_decimal_degrees[i],
                     alt_sites$Latitude_decimal_degrees, alt_sites$Longitude_decimal_degrees)
  k <- which.min(dk); na_sites$alt_imp[i] <- alt_sites$alt_site[k]; na_sites$dist_km[i] <- dk[k]
}
alt_lookup <- na_sites %>% select(Latitude_decimal_degrees, Longitude_decimal_degrees, alt_imp, dist_km) %>% distinct()

prep_frame <- function(d) {
  d %>%
    left_join(alt_lookup, by = c("Latitude_decimal_degrees", "Longitude_decimal_degrees")) %>%
    mutate(scaled_yr  = as.numeric(scaled_yr),
           scaled_lat = as.numeric(scaled_lat),
           scaled_lon = as.numeric(scaled_lon),
           alt_imputed = is.na(Altitude),
           alt_filled  = ifelse(alt_imputed, alt_imp, Altitude),
           scaled_alt_cc  = as.numeric(scaled_alt),                    # NA where Altitude NA
           scaled_alt     = ifelse(is.na(alt_filled), 0, (alt_filled - ALT_CENTRE) / ALT_SCALE),
           Sex   = factor(Sex, levels = c("Female", "Male")),   # explicit factor (reviewer H1)
           spp   = factor(Binomial),
           src   = factor(Main_researcher),
           site  = factor(Municipality),
           # individual = Ring x species; unringed birds are their own singleton
           ind   = factor(ifelse(is.na(Ring), paste0("unringed_", ID_ABT),
                                 paste(Ring, Binomial, sep = "|"))),
           wing_col = factor(wing_col, levels = c("right", "left", "generic")),
           season_cc = factor(as.character(season), levels = c("DJF", "MAM", "JJA", "SON")),
           season   = factor(ifelse(is.na(season), "unrecorded", as.character(season)),
                             levels = c("DJF", "MAM", "JJA", "SON", "unrecorded")),
           molt3    = factor(ifelse(is.na(Molt), "unrecorded", Molt),
                             levels = c("no", "Yes", "unrecorded")),
           has_wing = !is.na(conc.wing.length))
}
# first WING record per individual (unringed kept); call on a wing-only frame
flag_first_capture <- function(w) {
  w %>% arrange(Year, month, day) %>% group_by(Ring, Binomial) %>%
    mutate(first_capture = is.na(Ring) | row_number() == 1L) %>% ungroup()
}
add_mundlak <- function(d) {
  d %>% group_by(site) %>% mutate(yr_site_mean = mean(scaled_yr)) %>% ungroup() %>%
    group_by(src)  %>% mutate(yr_src_mean  = mean(scaled_yr)) %>% ungroup() %>%
    group_by(spp)  %>% mutate(lat_spp_mean = mean(scaled_lat)) %>% ungroup() %>%
    mutate(yr_within_site = scaled_yr - yr_site_mean,
           yr_within_src  = scaled_yr - yr_src_mean,
           lat_within_spp = scaled_lat - lat_spp_mean)
}
finish <- function(d) d %>% droplevels() %>% add_mundlak()

known   <- prep_frame(passer90)
w_all   <- known %>% filter(has_wing) %>% flag_first_capture() %>% finish()      # "full", 8,478
w_cc    <- w_all %>% filter(!alt_imputed, !is.na(season_cc)) %>% finish()          # "cc",   8,282
w_first <- w_all %>% filter(first_capture) %>% finish()
w_first_cc <- w_cc %>% filter(first_capture) %>% finish()
long_running <- w_all %>% group_by(src) %>%
  summarise(n_years = n_distinct(Year), n_wing = n(), .groups = "drop") %>%
  filter(n_years >= 8, n_wing >= 200) %>% pull(src) %>% as.character()
w_long    <- w_all %>% filter(as.character(src) %in% long_running) %>% finish()
w_long_cc <- w_cc  %>% filter(as.character(src) %in% long_running) %>% finish()
w_carrano <- w_all %>% filter(as.character(src) == "E.Carrano") %>% finish()

# --- sign-blind screen for anomalous contributor x municipality x year blocks -------
# Species x sex-centred wing deviations (centred on the whole wing sample) averaged per
# contributor x municipality x year block with >= 10 records. Blocks whose mean deviation
# exceeds +/- ANOM_MM (6 mm ~ 1.5 residual SD of the M3 fit, ~8 % of a wing) are flagged
# irrespective of sign: within one contributor such shifts are protocol / data-entry
# changes, not phenotypes. The screen exists because two contributors' late blocks sit
# 9-12 mm below species means and generate most of the within-contributor "decline".
ANOM_MM <- 6
block_screen <- w_all %>% group_by(spp, Sex) %>% mutate(dev = conc.wing.length - mean(conc.wing.length)) %>% ungroup() %>%
  group_by(src, site, Year) %>%
  summarise(n = n(), n_spp = n_distinct(spp), mean_dev = mean(dev), wing_col = paste(sort(unique(as.character(wing_col))), collapse = "/"), .groups = "drop") %>%
  mutate(flagged = n >= 10 & abs(mean_dev) > ANOM_MM) %>% arrange(desc(abs(mean_dev)))
anom_blocks <- block_screen %>% filter(flagged)
w_noanom <- w_all %>% anti_join(anom_blocks %>% select(src, site, Year), by = c("src", "site", "Year")) %>% finish()
message(sprintf("[screen] %d of %d blocks (n >= 10) flagged at |mean dev| > %g mm, %d records; sample without them = %d",
                nrow(anom_blocks), sum(block_screen$n >= 10), ANOM_MM, sum(anom_blocks$n), nrow(w_noanom)))
if (nrow(anom_blocks)) print(as.data.frame(anom_blocks), digits = 3)

unknown <- prep_frame(passer90_allsex %>% filter(!known_sex))
u_all <- unknown %>% filter(has_wing) %>% flag_first_capture() %>% finish()        # 3,532
u_cc  <- u_all %>% filter(!alt_imputed, !is.na(season_cc)) %>% finish()

stopifnot(nrow(w_all) == 8478, nrow(passer90) == 12571, nlevels(droplevels(known$spp)) == 73,
          !anyNA(w_all$scaled_alt), !anyNA(u_all$scaled_alt))
# (one of the 73 species has no wing record, so the wing models see 72 species)
message(sprintf("[data] wing records: full %d | cc %d | first captures %d/%d | long-running (%d src) %d/%d | E.Carrano %d | unknown-sex %d/%d",
                nrow(w_all), nrow(w_cc), nrow(w_first), nrow(w_first_cc), length(long_running),
                nrow(w_long), nrow(w_long_cc), nrow(w_carrano), nrow(u_all), nrow(u_cc)))

alt_imp_summary <- w_all %>% filter(alt_imputed) %>%
  group_by(src, site, Locality) %>%
  summarise(n = n(), yr_min = min(Year), yr_max = max(Year), alt_imputed_m = first(alt_filled),
            dist_km = first(dist_km), .groups = "drop") %>% arrange(desc(n))
sample_sizes <- list(
  wing_all = nrow(w_all), wing_complete_cases = nrow(w_cc),
  altitude_na_wing = sum(w_all$alt_imputed), season_na_wing = sum(is.na(w_all$season_cc)),
  altitude_imputation = as.data.frame(alt_imp_summary),
  altitude_imputation_note = "nearest coordinate site with a recorded altitude (all passer90_allsex records); municipality median rejected (Ilha Rasa vs inland farm)",
  first_captures_full = nrow(w_first), first_captures_cc = nrow(w_first_cc),
  first_captures_P0_definition = 7734L,
  long_running_contributors = long_running, long_running_wing_full = nrow(w_long), long_running_wing_cc = nrow(w_long_cc),
  carrano_wing = nrow(w_carrano), carrano_years = range(w_carrano$Year),
  unknown_sex_wing_full = nrow(u_all), unknown_sex_wing_cc = nrow(u_cc),
  unknown_sex_label_Unknown = sum(u_all$Sex == "Unknown", na.rm = TRUE),
  molt_unrecorded_share = mean(w_all$molt3 == "unrecorded"),
  anomalous_block_threshold_mm = ANOM_MM, anomalous_blocks = as.data.frame(anom_blocks),
  anomalous_records = sum(anom_blocks$n), wing_without_anomalous_blocks = nrow(w_noanom),
  individuals_full = nlevels(w_all$ind), individuals_with_repeats_full = sum(table(w_all$ind) > 1),
  sd_year = SD_YEAR, mean_wing_mm = mean(w_all$conc.wing.length))

# ===========================================================================
# 2. Model specifications (shared vocabulary for lme4 and brms)
# ===========================================================================
FX_BASE  <- "Sex + scaled_yr + scaled_lat"
FX_NOSEX <- "scaled_yr + scaled_lat"
FX_M2    <- "wing_col + scaled_lon + scaled_alt"
FX_M3    <- "season + molt3"
FX_MUND  <- "Sex + yr_within_site + yr_site_mean + lat_within_spp + lat_spp_mean"
FX_MUNDS <- "Sex + yr_within_src + yr_src_mean + lat_within_spp + lat_spp_mean"
RE_SPP   <- "(1 + scaled_yr || spp)"
RE_SPP_WSITE <- "(1 + yr_within_site || spp)"     # REWB: random slope on the within variable
RE_SPP_WSRC  <- "(1 + yr_within_src || spp)"
RE_SRC   <- "(1 | src)"; RE_SITE <- "(1 | site)"; RE_IND <- "(1 | ind)"
RE_SRCSL <- "(1 + scaled_yr || src)"
j <- function(...) paste(c(...), collapse = " + ")
spec <- function(fx, re, data, label) list(fx = fx, re = re, data = data, label = label)

# Ladder on one sample; `d` = frame name for the main models, `d_first` / `d_long` for the subsets
ladder <- function(sfx, d, d_first, d_long) {
  L <- list(
    M2a_wingcol = spec(j(FX_BASE, "wing_col"), j(RE_SPP, RE_SRC, RE_SITE), d, "M1 + wing-column protocol proxy"),
    M2          = spec(j(FX_BASE, FX_M2), j(RE_SPP, RE_SRC, RE_SITE), d, "M2 = M1 + wing_col + scaled_lon + scaled_alt"),
    M3a_season  = spec(j(FX_BASE, FX_M2, "season"), j(RE_SPP, RE_SRC, RE_SITE), d, "M2 + season"),
    M3b_ring    = spec(j(FX_BASE, FX_M2, "season"), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), d, "M2 + season + (1|individual)  [= M3 without moult]"),
    M3          = spec(j(FX_BASE, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), d, "M3 = M2 + season + moult + (1|individual)"),
    M4          = spec(j(FX_BASE, FX_M2, FX_M3), j(RE_SPP, RE_SRCSL, RE_SITE, RE_IND), d, "M4 = M3 with (1 + scaled_yr || contributor)"),
    M5          = spec(j(FX_MUND, FX_M2, FX_M3), j(RE_SPP_WSITE, RE_SRC, RE_SITE, RE_IND), d,
                       "M5 = M3 with Mundlak year within/between municipality, latitude within/among species (REWB: species slope on yr_within_site)"),
    M5_rsTotal  = spec(j(FX_MUND, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), d, "M5 with the species random slope on total scaled_yr"),
    M5src       = spec(j(FX_MUNDS, FX_M2, FX_M3), j(RE_SPP_WSRC, RE_SRC, RE_SITE, RE_IND), d,
                       "M3 with Mundlak year within/between CONTRIBUTOR (REWB: species slope on yr_within_src)"),
    M5src_rsTotal = spec(j(FX_MUNDS, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), d, "M5src with the species random slope on total scaled_yr"),
    M7          = spec(j(FX_BASE, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE), d_first, "M7 = M3 on first (wing) captures only; no ring term (one record per individual)"),
    M3_longrun  = spec(j(FX_BASE, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), d_long, "M3 on the 6 long-running contributors (>= 8 yr, >= 200 wing records)")
  )
  names(L) <- paste0(names(L), sfx); L
}
SPECS <- c(
  list(M0       = spec(FX_BASE, RE_SPP, "w_all", "M0 baseline (= published brms structure)"),
       M1a_src  = spec(FX_BASE, j(RE_SPP, RE_SRC), "w_all", "M0 + (1|contributor)"),
       M1b_site = spec(FX_BASE, j(RE_SPP, RE_SITE), "w_all", "M0 + (1|municipality)"),
       M1       = spec(FX_BASE, j(RE_SPP, RE_SRC, RE_SITE), "w_all", "M1 = M0 + (1|contributor) + (1|municipality)"),
       M0_cc    = spec(FX_BASE, RE_SPP, "w_cc", "M0 on the complete-case sample"),
       M1_cc    = spec(FX_BASE, j(RE_SPP, RE_SRC, RE_SITE), "w_cc", "M1 on the complete-case sample")),
  ladder("",    "w_all", "w_first",    "w_long"),
  ladder("_cc", "w_cc",  "w_first_cc", "w_long_cc"),
  list(M_carrano = spec(j(FX_BASE, FX_M2, FX_M3), j(RE_SPP, RE_SITE, RE_IND), "w_carrano",
                        "M3 without contributor term on E. Carrano's records alone (1996-2017)"),
       M1_noanom    = spec(FX_BASE, j(RE_SPP, RE_SRC, RE_SITE), "w_noanom", "M1 without the flagged anomalous blocks"),
       M3_noanom    = spec(j(FX_BASE, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), "w_noanom", "M3 without the flagged anomalous blocks"),
       M5src_noanom = spec(j(FX_MUNDS, FX_M2, FX_M3), j(RE_SPP_WSRC, RE_SRC, RE_SITE, RE_IND), "w_noanom",
                           "M5src (within/between contributor, REWB) without the flagged anomalous blocks"),
       M5src_rsTotal_noanom = spec(j(FX_MUNDS, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), "w_noanom",
                                   "M5src_rsTotal without the flagged anomalous blocks"),
       M6_m0 = spec(FX_NOSEX, RE_SPP, "u_all", "unknown-sex records, M0 structure (no Sex)"),
       M6_m1 = spec(FX_NOSEX, j(RE_SPP, RE_SRC, RE_SITE), "u_all", "unknown-sex records, M1 structure (no Sex)"),
       M6    = spec(j(FX_NOSEX, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), "u_all", "M6 = M3 on unknown-sex records (no Sex)"),
       M6_cc = spec(j(FX_NOSEX, FX_M2, FX_M3), j(RE_SPP, RE_SRC, RE_SITE, RE_IND), "u_cc", "M6 on the complete-case unknown-sex sample"))
)
lme4_formula <- function(s) as.formula(paste("conc.wing.length ~ 1 +", s$fx, "+", s$re))
frames <- list(w_all = w_all, w_cc = w_cc, w_first = w_first, w_first_cc = w_first_cc,
               w_long = w_long, w_long_cc = w_long_cc, w_carrano = w_carrano, w_noanom = w_noanom,
               u_all = u_all, u_cc = u_cc)
# Spec guard shared by the lme4 and glmmTMB tiers. Subset frames (E. Carrano alone):
# (1 | ind) is not identifiable when every individual has one record, and a factor
# with a single level (he used only the right-wing column) cannot enter the fixed
# part. Drop them and say so in the label. The single-level test must work for
# character AND factor columns: nlevels() of a character vector is 0, which silently
# dropped Sex from every known-sex model (reviewer H1).
guard_spec <- function(s, dat) {
  if (grepl("ind", s$re) && nlevels(dat$ind) == nrow(dat)) {
    s$re <- gsub("\\s*\\+\\s*\\(1 \\| ind\\)", "", s$re)
    s$label <- paste(s$label, "[ring term dropped: no repeated individuals]")
  }
  for (fct in c("Sex", "wing_col", "season", "molt3")) {
    if (grepl(paste0("\\b", fct, "\\b"), s$fx) && length(unique(na.omit(dat[[fct]]))) < 2) {
      s$fx <- gsub(paste0("\\s*\\+\\s*\\b", fct, "\\b|\\b", fct, "\\b\\s*\\+\\s*"), "", s$fx)
      s$label <- paste0(s$label, " [", fct, " dropped: single level]")
    }
  }
  s
}

# ===========================================================================
# 3. Tier 1: lme4 REML screening
# ===========================================================================
YEAR_TERMS <- c("scaled_yr", "yr_within_site", "yr_site_mean", "yr_within_src", "yr_src_mean")
ctrl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
tidy_lmer <- function(fit, name, s, dat) {
  cf <- summary(fit)$coefficients
  msgs <- unlist(fit@optinfo$conv$lme4$messages)
  data.frame(model = name, label = s$label, term = rownames(cf),
             estimate = cf[, "Estimate"], se = cf[, "Std. Error"],
             ci_lo = cf[, "Estimate"] - 1.96 * cf[, "Std. Error"],
             ci_hi = cf[, "Estimate"] + 1.96 * cf[, "Std. Error"],
             t = cf[, "t value"],
             N = nobs(fit), n_spp = nlevels(dat$spp), n_src = nlevels(dat$src), n_site = nlevels(dat$site),
             n_ind = if (grepl("ind", s$re)) nlevels(dat$ind) else NA_integer_,
             sample = s$data, fixed = s$fx, random = s$re,
             singular = isSingular(fit),
             convergence_msg = if (length(msgs)) paste(msgs, collapse = "; ") else "",
             row.names = NULL, stringsAsFactors = FALSE) %>%
    mutate(year_term = term %in% YEAR_TERMS,
           mm_per_decade    = ifelse(year_term, estimate / SD_YEAR * 10, NA_real_),
           mm_per_decade_lo = ifelse(year_term, ci_lo / SD_YEAR * 10, NA_real_),
           mm_per_decade_hi = ifelse(year_term, ci_hi / SD_YEAR * 10, NA_real_),
           pct_per_decade   = mm_per_decade / sample_sizes$mean_wing_mm * 100)
}

if (MODE == "fast") {
  fits <- list(); rows <- list(); timings <- c()
  for (nm in names(SPECS)) {
    s <- guard_spec(SPECS[[nm]], frames[[SPECS[[nm]]$data]]); dat <- frames[[s$data]]
    SPECS[[nm]] <- s
    t0 <- Sys.time()
    fit <- lmer(lme4_formula(s), data = dat, REML = TRUE, control = ctrl)
    timings[nm] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fits[[nm]] <- fit
    rows[[nm]] <- tidy_lmer(fit, nm, s, dat)
    yr <- rows[[nm]] %>% filter(year_term)
    message(sprintf("[lme4] %-14s N=%5d  %s  (%.1fs)%s", nm, nobs(fit),
                    paste(sprintf("%s=%.3f [%.3f, %.3f] t=%.2f", yr$term, yr$estimate, yr$ci_lo, yr$ci_hi, yr$t), collapse = " | "),
                    timings[nm], if (isSingular(fit)) "  SINGULAR" else ""))
  }
  results <- bind_rows(rows)

  # --- variance components ---------------------------------------------------
  varcomp <- bind_rows(lapply(names(fits), function(nm) {
    v <- as.data.frame(VarCorr(fits[[nm]]))
    data.frame(model = nm, grp = v$grp, var1 = v$var1, sd = v$sdcor, stringsAsFactors = FALSE)
  }))

  # --- leave-one-contributor-out jackknife (M1 on full; M3 on full and cc) --------
  jackknife <- function(spec_name, dat) {
    s <- SPECS[[spec_name]]; f <- lme4_formula(s)
    bind_rows(lapply(levels(dat$src), function(k) {
      d <- dat %>% filter(as.character(src) != k) %>% finish()
      fit <- lmer(f, data = d, REML = TRUE, control = ctrl)
      cf <- summary(fit)$coefficients["scaled_yr", ]
      data.frame(model = spec_name, dropped = k, n_dropped = sum(as.character(dat$src) == k),
                 estimate = cf[1], se = cf[2], ci_lo = cf[1] - 1.96 * cf[2], ci_hi = cf[1] + 1.96 * cf[2],
                 t = cf[3], N = nobs(fit), row.names = NULL)
    })) %>% mutate(mm_per_decade = estimate / SD_YEAR * 10)
  }
  message("[jackknife] leave-one-contributor-out: M1 (full), M3 (full), M3_cc ...")
  jk <- bind_rows(jackknife("M1", w_all), jackknife("M3", w_all), jackknife("M3_cc", w_cc))
  full_est <- results %>% filter(term == "scaled_yr", model %in% c("M1", "M3", "M3_cc")) %>% select(model, full_estimate = estimate)
  jk <- jk %>% left_join(full_est, by = "model") %>% mutate(shift = estimate - full_estimate) %>% arrange(model, shift)
  jk_summary <- jk %>% group_by(model) %>%
    summarise(full_estimate = first(full_estimate), min = min(estimate), max = max(estimate),
              most_negative_when_dropping = dropped[which.min(estimate)],
              least_negative_when_dropping = dropped[which.max(estimate)],
              n_excl_zero = sum(ci_hi < 0), n_fits = n(), .groups = "drop")
  print(jk_summary)

  # --- per-contributor year slopes for the spanning / long-running contributors -----
  # wing ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) + (1 | site) within each
  # contributor with >= 8 distinct years (the within-contributor trend, one contributor at a time)
  span_src <- w_all %>% group_by(src) %>% summarise(n_years = n_distinct(Year), n = n(), yr_min = min(Year), yr_max = max(Year),
                                                    n_site = n_distinct(site), n_spp = n_distinct(spp), .groups = "drop") %>%
    filter(n_years >= 8) %>% arrange(desc(n))
  spanning_slopes <- bind_rows(lapply(seq_len(nrow(span_src)), function(i) {
    k <- as.character(span_src$src[i]); d <- w_all %>% filter(as.character(src) == k) %>% finish()
    re <- if (nlevels(d$site) > 1) "(1 + scaled_yr || spp) + (1 | site)" else "(1 + scaled_yr || spp)"
    fit <- tryCatch(lmer(as.formula(paste("conc.wing.length ~ Sex + scaled_yr + scaled_lat +", re)), data = d, REML = TRUE, control = ctrl),
                    error = function(e) lmer(conc.wing.length ~ Sex + scaled_yr + (1 + scaled_yr || spp), data = d, REML = TRUE, control = ctrl))
    cf <- summary(fit)$coefficients["scaled_yr", ]
    # descriptive: species x sex-centred wing (centred on the whole sample) by period within the contributor
    cen <- w_all %>% group_by(spp, Sex) %>% mutate(dev = conc.wing.length - mean(conc.wing.length)) %>% ungroup() %>%
      filter(as.character(src) == k)
    data.frame(src = k, n = nrow(d), yr_min = span_src$yr_min[i], yr_max = span_src$yr_max[i], n_years = span_src$n_years[i],
               n_site = span_src$n_site[i], n_spp = span_src$n_spp[i],
               wing_col = paste(names(which(table(d$wing_col) > 0)), collapse = "/"),
               estimate = cf[1], se = cf[2], ci_lo = cf[1] - 1.96 * cf[2], ci_hi = cf[1] + 1.96 * cf[2], t = cf[3],
               mm_per_decade = cf[1] / SD_YEAR * 10, singular = isSingular(fit),
               centred_dev_early = mean(cen$dev[cen$Year <= 2006]), n_early = sum(cen$Year <= 2006),
               centred_dev_late = mean(cen$dev[cen$Year >= 2013]), n_late = sum(cen$Year >= 2013), row.names = NULL)
  }))
  print(spanning_slopes %>% select(src, n, yr_min, yr_max, estimate, ci_lo, ci_hi, t, centred_dev_early, centred_dev_late), digits = 3)

  # --- species-specific slopes (fixed + random) from M3 (full) and M0 --------------
  species_slopes <- function(fit, dat, name) {
    fe <- fixef(fit)[["scaled_yr"]]; fe_se <- summary(fit)$coefficients["scaled_yr", "Std. Error"]
    re <- ranef(fit, condVar = TRUE)
    # lme4 expands (1 + scaled_yr || spp) into (1|spp) + (0 + scaled_yr|spp); locate the
    # component that carries the scaled_yr column whichever way ranef() names it.
    comp <- NULL
    for (k in names(re)) if ("scaled_yr" %in% colnames(re[[k]])) comp <- re[[k]]
    stopifnot(!is.null(comp))
    pv <- attr(comp, "postVar"); jcol <- match("scaled_yr", colnames(comp))
    # lme4 >= 2.0 merges the two || terms and returns postVar as a list of 1x1xn arrays
    # (one per term); older versions return a q x q x n array.
    re_var <- if (is.list(pv)) as.numeric(pv[["scaled_yr"]])
              else if (length(dim(pv)) == 3) pv[jcol, jcol, ] else as.numeric(pv)
    spp_info <- dat %>% group_by(spp) %>%
      summarise(n = n(), yr_min = min(Year), yr_max = max(Year), n_src = n_distinct(src),
                n_site = n_distinct(site), mean_wing = mean(conc.wing.length), .groups = "drop") %>%
      mutate(spp = as.character(spp))
    data.frame(model = name, spp = rownames(comp), ranef_yr = comp[["scaled_yr"]],
               se_ranef = sqrt(re_var), stringsAsFactors = FALSE) %>%
      mutate(slope = fe + ranef_yr, fixed = fe, se_fixed = fe_se,
             se_total = sqrt(se_ranef^2 + se_fixed^2),
             lo = slope - 1.96 * se_total, hi = slope + 1.96 * se_total,
             mm_per_decade = slope / SD_YEAR * 10,
             mm_per_decade_lo = lo / SD_YEAR * 10, mm_per_decade_hi = hi / SD_YEAR * 10) %>%
      left_join(spp_info, by = "spp") %>%
      mutate(pct_per_decade = mm_per_decade / mean_wing * 100) %>%
      arrange(slope)
  }
  sl_M3 <- species_slopes(fits$M3, w_all, "M3"); sl_M3cc <- species_slopes(fits$M3_cc, w_cc, "M3_cc")
  sl_M0 <- species_slopes(fits$M0, w_all, "M0")
  slope_summary <- function(sl) list(
    n_species = nrow(sl), n_negative = sum(sl$slope < 0),
    n_excl_zero_negative = sum(sl$hi < 0), n_excl_zero_positive = sum(sl$lo > 0),
    median_mm_per_decade = median(sl$mm_per_decade),
    fixed_slope = sl$fixed[1], sd_random_slope = sd(sl$ranef_yr))
  slopes <- list(generated = Sys.time(), M3 = sl_M3, M3_cc = sl_M3cc, M0 = sl_M0,
                 summary = list(M3 = slope_summary(sl_M3), M3_cc = slope_summary(sl_M3cc), M0 = slope_summary(sl_M0)),
                 sd_year = SD_YEAR,
                 note = paste("slope = fixed scaled_yr + species random slope (lme4 REML; M3 = full controls on the full sample,",
                              "M3_cc = complete cases, M0 = baseline); se_total adds the fixed-effect SE and the conditional SD of the",
                              "BLUP in quadrature; 'excl. zero' uses slope +/- 1.96 se_total. The brms posteriors (Tiers 2-3) supersede these for the figure."))
  saveRDS(slopes, out_path("controlled_wing_species_slopes.rds"))

  # --- decision gate (provisional, Tier 1) -----------------------------------
  g <- function(m, term = "scaled_yr") results %>% filter(model == m, term == !!term)
  m0 <- g("M0"); m3 <- g("M3"); m3cc <- g("M3_cc")
  m5w <- g("M5", "yr_within_site"); m5wcc <- g("M5_cc", "yr_within_site"); m5b <- g("M5", "yr_site_mean")
  m5sw <- g("M5src", "yr_within_src"); m5sb <- g("M5src", "yr_src_mean")
  m5sw_alt <- g("M5src_rsTotal", "yr_within_src"); m5w_alt <- g("M5_rsTotal", "yr_within_site")
  m6 <- g("M6"); m7 <- g("M7"); mc <- g("M_carrano")
  excl0 <- function(r) r$ci_hi < 0
  scenario_full <- if (excl0(m3) && excl0(m5w)) "A" else "B"
  scenario_cc   <- if (excl0(m3cc) && excl0(m5wcc)) "A" else "B"
  decision_gate <- list(
    tier = "1 (lme4 REML, Wald CI; provisional until the brms tiers run)",
    rule = "Scenario A requires the controlled year effect (M3) AND the within-municipality year slope (M5, REWB) to exclude zero; otherwise B",
    scenario_full_sample = scenario_full, scenario_complete_cases = scenario_cc,
    scenario = if (scenario_full == scenario_cc) scenario_full else "B (samples disagree; see note)",
    M0 = m0[, c("estimate", "ci_lo", "ci_hi", "t", "N")],
    M3 = m3[, c("estimate", "ci_lo", "ci_hi", "t", "N", "mm_per_decade")],
    M3_cc = m3cc[, c("estimate", "ci_lo", "ci_hi", "t", "N", "mm_per_decade")],
    M5_within_site = m5w[, c("estimate", "ci_lo", "ci_hi", "t")], M5_within_site_cc = m5wcc[, c("estimate", "ci_lo", "ci_hi", "t")],
    M5_within_site_rsTotal = m5w_alt[, c("estimate", "ci_lo", "ci_hi", "t")],
    M5_between_site = m5b[, c("estimate", "ci_lo", "ci_hi", "t")],
    M5src_within = m5sw[, c("estimate", "ci_lo", "ci_hi", "t")], M5src_within_rsTotal = m5sw_alt[, c("estimate", "ci_lo", "ci_hi", "t")],
    M5src_between = m5sb[, c("estimate", "ci_lo", "ci_hi", "t")],
    M6 = m6[, c("estimate", "ci_lo", "ci_hi", "t", "N")], M7 = m7[, c("estimate", "ci_lo", "ci_hi", "t", "N")],
    M_carrano = mc[, c("estimate", "ci_lo", "ci_hi", "t", "N")],
    M3_noanom = g("M3_noanom")[, c("estimate", "ci_lo", "ci_hi", "t", "N")],
    M5src_within_noanom = g("M5src_noanom", "yr_within_src")[, c("estimate", "ci_lo", "ci_hi", "t")],
    M5src_within_rsTotal_noanom = g("M5src_rsTotal_noanom", "yr_within_src")[, c("estimate", "ci_lo", "ci_hi", "t")],
    share_of_M0_remaining_in_M3 = m3$estimate / m0$estimate,
    share_of_M0_remaining_in_M3_cc = m3cc$estimate / m0$estimate,
    jackknife = jk_summary)
  message(sprintf("[gate] Tier-1 scenario: full = %s, cc = %s | M3 = %.3f [%.3f, %.3f] (cc %.3f [%.3f, %.3f]) | within-site = %.3f [%.3f, %.3f] | within-contributor = %.3f [%.3f, %.3f]",
                  scenario_full, scenario_cc, m3$estimate, m3$ci_lo, m3$ci_hi, m3cc$estimate, m3cc$ci_lo, m3cc$ci_hi,
                  m5w$estimate, m5w$ci_lo, m5w$ci_hi, m5sw$estimate, m5sw$ci_lo, m5sw$ci_hi))

  before_after <- results %>% filter(year_term) %>%
    select(model, label, term, estimate, se, ci_lo, ci_hi, t, N, n_spp, n_src, n_site, n_ind, sample,
           mm_per_decade, mm_per_decade_lo, mm_per_decade_hi, pct_per_decade, singular, convergence_msg)

  out <- list(
    generated = Sys.time(), mode = "fast-lme4",
    session = list(R = R.version.string, lme4 = as.character(packageVersion("lme4")),
                   dplyr = as.character(packageVersion("dplyr")), optimizer = "bobyqa, calc.derivs = FALSE, REML"),
    sd_year = SD_YEAR, mean_wing_mm = sample_sizes$mean_wing_mm,
    sample_sizes = sample_sizes,
    specs = bind_rows(lapply(names(SPECS), function(n) data.frame(model = n, label = SPECS[[n]]$label,
                                                                   fixed = SPECS[[n]]$fx, random = SPECS[[n]]$re,
                                                                   sample = SPECS[[n]]$data))),
    table = results,             # all fixed effects, all models
    before_after = before_after, # year terms only
    varcomp = varcomp,
    jackknife = jk, jackknife_summary = jk_summary,
    spanning_contributor_slopes = spanning_slopes,
    block_screen = as.data.frame(block_screen), anomalous_blocks = as.data.frame(anom_blocks),
    decision_gate = decision_gate,
    timings_sec = timings,
    note = paste("lme4 REML re-fit of the primary wing model with artefact controls (REVISION_PLAN.md Phase 1, Tier 1).",
                 "CIs are Wald (+/- 1.96 SE). mm_per_decade = estimate / SD(year) * 10 with SD(year) = 5.0199 (the scaling SD).",
                 "Models without suffix use the full 8,478-record sample (altitude imputed from the nearest site, season 'unrecorded' for 6 undated records);",
                 "'_cc' models use the 8,282 complete cases of REVISION_PLAN.md §0.",
                 "No phylogenetic term in this tier; species intercepts absorb it. The brms tiers (--trees) are the inferential result."))
  saveRDS(out, out_path("controlled_wing_results.rds"))
  message("[write] ", out_path("controlled_wing_results.rds"), " and controlled_wing_species_slopes.rds")
  print(before_after %>% select(model, term, estimate, ci_lo, ci_hi, t, N, mm_per_decade), digits = 3, row.names = FALSE)
  quit(save = "no", status = 0)
}

# ===========================================================================
# 4. Phylogenetic tier, --engine glmmTMB (default for --trees N): every Tier-1
#    specification with propto(0 + species_name | g, A) on each tree, Rubin-pooled
# ===========================================================================
if (ENGINE == "glmmTMB") {
  # _phylo_engine.R uses out_path()/derived_path()/raw_path(), so it is sourced
  # after the path helpers. It loads glmmTMB and provides phylo_species_name(),
  # phylo_A_list(), fit_phylo_glmmtmb(), tidy_phylo_fit(), pool_rubin_df(),
  # run_phylo_trees(), species_slopes_phylo(). Nothing from it is re-implemented here.
  source(script_path("_phylo_engine.R"))
  T_START <- Sys.time()

  # --- species names on the tree; drop species the trees do not carry ----------------
  cache_file <- derived_path("phylo_A_50trees.rds")
  if (!file.exists(cache_file)) invisible(phylo_A_list(unique(phylo_species_name(passer90$Binomial)), n_trees = N_TREES))
  tree_species <- if (file.exists(cache_file)) readRDS(cache_file)$species else character(0)
  NOT_ON_TREE <- "Herpsilochmus_sellowi"           # known absentee; anything else missing is dropped too
  dropped_species <- character(0)
  attach_phylo_names <- function(d) {
    d$species_name <- phylo_species_name(d$Binomial)
    drop <- setdiff(unique(d$species_name), setdiff(tree_species, NOT_ON_TREE))
    if (length(tree_species) == 0) drop <- intersect(unique(d$species_name), NOT_ON_TREE)
    if (length(drop)) {
      message("[phylo] dropping ", length(drop), " species absent from the trees: ", paste(drop, collapse = ", "))
      dropped_species <<- union(dropped_species, drop)
    }
    d %>% filter(!species_name %in% drop) %>%
      mutate(spp = factor(species_name)) %>%                    # species random slope on the same (tree) names
      select(-yr_site_mean, -yr_src_mean, -lat_spp_mean, -yr_within_site, -yr_within_src, -lat_within_spp) %>%
      finish()                                                  # re-derive Mundlak terms on the matched species set
  }
  phylo_frames <- lapply(frames, attach_phylo_names)
  binomial_map <- bind_rows(lapply(phylo_frames, function(d) distinct(d, Binomial, species_name))) %>% distinct()
  stopifnot(!any(duplicated(binomial_map$species_name)))       # relabelling is 1:1
  stopifnot(nrow(phylo_frames$w_all) == nrow(w_all))           # no species lost on the primary frame

  # --- which specs ------------------------------------------------------------------
  PHYLO_MODELS <- if (identical(BRMS_MODELS, "all")) names(SPECS) else BRMS_MODELS
  unknown_models <- setdiff(PHYLO_MODELS, names(SPECS))
  if (length(unknown_models)) stop("unknown model(s): ", paste(unknown_models, collapse = ", "))
  SLOPE_MODELS <- intersect(c("M0", "M3", "M3_cc"), PHYLO_MODELS)   # species slopes are taken from these fits
  A_cache <- list()
  A_for <- function(species) {                                  # one A_list per species set (subset frames differ)
    key <- paste(sort(species), collapse = "|")
    if (is.null(A_cache[[key]])) A_cache[[key]] <<- phylo_A_list(sort(species), n_trees = N_TREES)
    A_cache[[key]]
  }

  # --- informed restart for trees that stop on the i.i.d. side of the species-intercept ridge ---
  # Every spec carries BOTH an i.i.d. species intercept (1 | spp) and the phylogenetic one
  # (propto), as the published brms model did. Their variances are nearly exchangeable, and
  # from glmmTMB's default start (all log-SDs = 0) nlminb stops, on a tree-dependent subset
  # of trees, at a point where the phylogenetic SD is ~0.1 and the Hessian is not positive
  # definite (a saddle / local optimum: on M1a_src its REML objective is 24730.8 against
  # 24707.2 on the phylogenetic side, so it is NOT the REML optimum). The fixed effects are
  # identical to 3 decimals on both sides (probed 2026-09-09), so nothing inferential hangs
  # on this, but the convergence flag and phylo_prop do. Remedy, using only engine functions:
  # for each non-PD tree take that tree's own parameter layout (fit_phylo_glmmtmb(doFit =
  # FALSE): the mapped theta entries encode the tree's Cholesky factor, so a start vector
  # must be built PER TREE - a constant one silently replaces A by the identity), restart
  # with the free propto scale at log(20) (~ the phylogenetic SD the converged trees reach),
  # accept the refit only if its Hessian is PD AND its REML objective is not worse, and
  # re-pool with pool_rubin_df(). Per-tree objectives on both sides are kept in $retry.
  PHYLO_START_SD <- 20
  retry_nonpd <- function(res, formula, data, A_list) {
    t0 <- Sys.time()
    obj0 <- vapply(res$fits, function(f) f$fit$objective, numeric(1))
    retry <- data.frame(tree = seq_along(A_list), pdHess_default = res$varcomp$converged, objective_default = obj0,
                        retried = FALSE, pdHess_final = res$varcomp$converged, objective_final = obj0, accepted = FALSE)
    for (i in which(!res$varcomp$converged)) {
      st <- fit_phylo_glmmtmb(formula, data, A_list[[i]], doFit = FALSE)
      th <- st$parameters$theta; free <- which(!is.na(st$mapArg$theta))   # last free theta = propto scale
      th[free[length(free)]] <- log(PHYLO_START_SD)
      fit2 <- tryCatch(suppressWarnings(fit_phylo_glmmtmb(formula, data, A_list[[i]], start = list(theta = th))), error = function(e) NULL)
      if (!is.null(fit2) && !isTRUE(fit2$sdr$pdHess)) {                 # second attempt: also shrink the i.i.d. species intercept start
        th[free[1]] <- log(0.5)
        fit2 <- tryCatch(suppressWarnings(fit_phylo_glmmtmb(formula, data, A_list[[i]], start = list(theta = th))), error = function(e) NULL)
      }
      retry$retried[i] <- TRUE
      if (!is.null(fit2)) { retry$pdHess_final[i] <- isTRUE(fit2$sdr$pdHess); retry$objective_final[i] <- fit2$fit$objective }
      if (!is.null(fit2) && isTRUE(fit2$sdr$pdHess) && fit2$fit$objective <= obj0[i] + 1e-6) {
        retry$accepted[i] <- TRUE
        td <- tidy_phylo_fit(fit2, tree = i)
        res$fits[[i]] <- fit2
        res$per_tree_fixed <- rbind(res$per_tree_fixed[res$per_tree_fixed$tree != i, ], td$fixed)
        res$varcomp[res$varcomp$tree == i, names(td$varcomp)] <- td$varcomp
      } else { retry$pdHess_final[i] <- res$varcomp$converged[i]; retry$objective_final[i] <- obj0[i] }
    }
    res$per_tree_fixed <- res$per_tree_fixed[order(res$per_tree_fixed$tree, res$per_tree_fixed$component), ]
    res$pooled <- pool_rubin_df(res$per_tree_fixed)
    res$varcomp_summary <- res$varcomp %>% summarise(across(where(is.numeric) & !tree, list(mean = mean, lo = ~quantile(.x, .025), hi = ~quantile(.x, .975))))
    res$retry <- retry; res$secs_retry <- as.numeric(Sys.time() - t0, units = "secs"); res$secs <- res$secs + res$secs_retry
    res
  }

  # --- fit every spec across the trees ---------------------------------------------
  phylo_runs <- list(); rows <- list(); table_rows <- list(); vc_rows <- list(); ridge_rows <- list(); timings <- c(); errors <- list()
  slope_per_tree <- list()
  for (nm in PHYLO_MODELS) {
    s <- guard_spec(SPECS[[nm]], phylo_frames[[SPECS[[nm]]$data]]); dat <- phylo_frames[[s$data]]
    SPECS[[nm]] <- s
    A_list <- A_for(levels(droplevels(dat$spp)))
    f <- lme4_formula(s)                                        # same syntax in glmmTMB; propto term appended by the engine
    # fits are kept transiently (objective values for the restart rule, species slopes), then dropped
    res <- tryCatch(suppressWarnings(run_phylo_trees(f, dat, A_list, keep_fits = TRUE, verbose = FALSE)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      errors[[nm]] <- conditionMessage(res)
      message(sprintf("[glmmTMB] %-20s FAILED: %s", nm, conditionMessage(res))); next
    }
    n_conv_default <- sum(res$varcomp$converged)
    res <- retry_nonpd(res, f, dat, A_list)
    if (nm %in% SLOPE_MODELS)                                   # per-tree species slopes now, while the fits exist
      slope_per_tree[[nm]] <- bind_rows(lapply(seq_along(res$fits), function(i) cbind(tree = i, species_slopes_phylo(res$fits[[i]], "scaled_yr", "spp"))))
    res$fits <- NULL
    timings[nm] <- res$secs
    phylo_runs[[nm]] <- res[c("pooled", "per_tree_fixed", "varcomp", "varcomp_summary", "n_trees", "secs", "secs_retry", "retry")]
    n_conv <- sum(res$varcomp$converged); phylo_prop <- mean(res$varcomp$phylo_prop, na.rm = TRUE)
    n_retried <- sum(res$retry$retried); n_accepted <- sum(res$retry$accepted)
    # Ridge diagnostic (after the restarts). Trees are classed by which species-intercept term
    # carries the variance: "phylo" (propto SD >= i.i.d. SD) or "iid". If any tree is still on
    # the i.i.d. side, phylo_prop is bimodal and its mean is not interpretable; the year
    # estimate by side is recorded so the reader can check it is the same. species_level_prop
    # (phylo + i.i.d. species intercept variance over the total) is invariant to the side.
    vc_t <- res$varcomp
    spp_int_col <- grep("^spp\\.\\.Intercept\\.$", names(vc_t), value = TRUE)
    sd_spp_int <- if (length(spp_int_col)) vc_t[[spp_int_col]] else 0
    other_cols <- setdiff(names(vc_t)[sapply(vc_t, is.numeric)], c("tree", "sd_phylo", "sigma", "phylo_prop"))
    tot_var <- vc_t$sd_phylo^2 + rowSums(as.matrix(vc_t[, other_cols, drop = FALSE])^2) + vc_t$sigma^2
    vc_t$species_level_prop <- (vc_t$sd_phylo^2 + sd_spp_int^2) / tot_var
    vc_t$mode <- ifelse(vc_t$sd_phylo^2 >= sd_spp_int^2, "phylo", "iid")
    yr_main <- intersect(c("scaled_yr", "yr_within_site", "yr_within_src"), res$per_tree_fixed$par[res$per_tree_fixed$component == "cond"])[1]
    yr_by_tree <- res$per_tree_fixed %>% filter(component == "cond", par == yr_main) %>% select(tree, estimate, se) %>%
      left_join(vc_t %>% select(tree, mode, converged), by = "tree")
    by_mode <- function(md, what) { x <- yr_by_tree[[what]][yr_by_tree$mode == md]; if (length(x)) mean(x) else NA_real_ }
    ridge <- data.frame(model = nm, year_term = yr_main,
                        n_trees = nrow(vc_t), n_phylo_mode = sum(vc_t$mode == "phylo"), n_iid_mode = sum(vc_t$mode == "iid"),
                        pdHess_phylo_mode = sum(vc_t$converged[vc_t$mode == "phylo"]), pdHess_iid_mode = sum(vc_t$converged[vc_t$mode == "iid"]),
                        yr_est_phylo_mode = by_mode("phylo", "estimate"), yr_est_iid_mode = by_mode("iid", "estimate"),
                        yr_se_phylo_mode = by_mode("phylo", "se"), yr_se_iid_mode = by_mode("iid", "se"),
                        yr_est_range_across_trees = diff(range(yr_by_tree$estimate)),
                        sd_phylo_phylo_mode = mean(vc_t$sd_phylo[vc_t$mode == "phylo"]), sd_spp_int_iid_mode = mean(sd_spp_int[vc_t$mode == "iid"]),
                        phylo_prop_phylo_mode = mean(vc_t$phylo_prop[vc_t$mode == "phylo"]),
                        species_level_prop = mean(vc_t$species_level_prop), stringsAsFactors = FALSE)
    ridge_rows[[nm]] <- ridge
    phylo_runs[[nm]]$varcomp <- vc_t                            # per-tree table now carries mode + species_level_prop
    cond <- res$pooled %>% filter(component == "cond")
    tab <- data.frame(model = nm, label = s$label, term = cond$par,
                      estimate = cond$estimate, se = cond$se, ci_lo = cond$lower, ci_hi = cond$upper,
                      t = cond$z, z = cond$z,
                      within_tree_var = cond$ubar, between_tree_var = cond$between_tree_var,
                      N = nrow(dat), n_spp = nlevels(droplevels(dat$spp)), n_src = nlevels(droplevels(dat$src)),
                      n_site = nlevels(droplevels(dat$site)),
                      n_ind = if (grepl("ind", s$re)) nlevels(droplevels(dat$ind)) else NA_integer_,
                      sample = s$data, fixed = s$fx, random = s$re,
                      n_trees = res$n_trees, converged = n_conv, converged_default_start = n_conv_default,
                      n_retried = n_retried, n_retry_accepted = n_accepted,
                      secs = res$secs, secs_retry = res$secs_retry, phylo_prop = phylo_prop,
                      species_level_prop = ridge$species_level_prop, n_phylo_mode = ridge$n_phylo_mode, n_iid_mode = ridge$n_iid_mode,
                      row.names = NULL, stringsAsFactors = FALSE) %>%
      mutate(year_term = term %in% YEAR_TERMS,
             mm_per_decade    = ifelse(year_term, estimate / SD_YEAR * 10, NA_real_),
             mm_per_decade_lo = ifelse(year_term, ci_lo / SD_YEAR * 10, NA_real_),
             mm_per_decade_hi = ifelse(year_term, ci_hi / SD_YEAR * 10, NA_real_),
             pct_per_decade   = mm_per_decade / sample_sizes$mean_wing_mm * 100)
    table_rows[[nm]] <- tab; rows[[nm]] <- tab %>% filter(year_term)
    # variance components: long format like the Tier-1 $varcomp, mean and 95 % range over trees
    vcs <- res$varcomp %>% select(-tree, -converged, -phylo_prop)
    vc_rows[[nm]] <- bind_rows(lapply(names(vcs), function(cn) data.frame(
      model = nm, component = cn, sd = mean(vcs[[cn]], na.rm = TRUE),
      sd_lo = unname(quantile(vcs[[cn]], .025, na.rm = TRUE)), sd_hi = unname(quantile(vcs[[cn]], .975, na.rm = TRUE)),
      stringsAsFactors = FALSE))) %>%
      bind_rows(data.frame(model = nm, component = "phylo_prop", sd = phylo_prop,
                           sd_lo = unname(quantile(res$varcomp$phylo_prop, .025, na.rm = TRUE)),
                           sd_hi = unname(quantile(res$varcomp$phylo_prop, .975, na.rm = TRUE))))
    yr <- rows[[nm]]
    message(sprintf("[glmmTMB] %-20s N=%5d  %s  | phylo %.2f | pdHess %d/%d (default start %d; %d restarted, %d accepted) | %.0fs (+%.0fs restarts)", nm, nrow(dat),
                    paste(sprintf("%s=%.3f [%.3f, %.3f] z=%.2f", yr$term, yr$estimate, yr$ci_lo, yr$ci_hi, yr$z), collapse = " | "),
                    phylo_prop, n_conv, res$n_trees, n_conv_default, n_retried, n_accepted, res$secs, res$secs_retry))
  }
  results <- bind_rows(table_rows); before_after <- bind_rows(rows); varcomp <- bind_rows(vc_rows); ridge <- bind_rows(ridge_rows)
  retry_all <- bind_rows(lapply(names(phylo_runs), function(nm) cbind(model = nm, phylo_runs[[nm]]$retry)))
  message(sprintf("[restart] %d of %d fits were non-PD from the default start; %d restarted, %d accepted (PD and REML objective not worse); %d fits remain non-PD. Mean objective gain on accepted restarts: %.2f",
                  sum(!retry_all$pdHess_default), nrow(retry_all), sum(retry_all$retried), sum(retry_all$accepted), sum(!retry_all$pdHess_final),
                  if (any(retry_all$accepted)) mean((retry_all$objective_default - retry_all$objective_final)[retry_all$accepted]) else NA_real_))
  message(sprintf("[ridge] trees still on the i.i.d. side of the species-intercept ridge after restarts: %d of %d fits; max |year estimate, phylo side - iid side| = %.4f",
                  sum(ridge$n_iid_mode), sum(ridge$n_trees), suppressWarnings(max(abs(ridge$yr_est_phylo_mode - ridge$yr_est_iid_mode), na.rm = TRUE))))

  # --- species slopes (computed per tree inside the loop), pooled across trees ----------
  pool_species_slopes <- function(nm) {
    per_tree <- slope_per_tree[[nm]]; dat <- phylo_frames[[SPECS[[nm]]$data]]
    if (is.null(per_tree)) return(NULL)
    m <- length(unique(per_tree$tree))
    fe <- before_after %>% filter(model == nm, term == "scaled_yr")
    spp_info <- dat %>% group_by(species_name) %>%
      summarise(Binomial = first(Binomial), n = n(), yr_min = min(Year), yr_max = max(Year), n_src = n_distinct(src),
                n_site = n_distinct(site), mean_wing = mean(conc.wing.length), .groups = "drop")
    # (summarise() evaluates sequentially: compute the pooled SEs before overwriting slope/ranef with their means)
    per_tree %>% group_by(tree_name = spp) %>%
      summarise(se_ranef = sqrt(mean(se_ranef^2) + (1 + 1/m) * var(ranef)),
                se_total = sqrt(mean(se_total^2) + (1 + 1/m) * var(slope)),
                slope_tree_min = min(slope), slope_tree_max = max(slope), n_trees = n(),
                ranef_yr = mean(ranef), slope = mean(slope), .groups = "drop") %>%
      mutate(model = nm, fixed = fe$estimate, se_fixed = fe$se,
             lo = slope - 1.96 * se_total, hi = slope + 1.96 * se_total,
             mm_per_decade = slope / SD_YEAR * 10, mm_per_decade_lo = lo / SD_YEAR * 10, mm_per_decade_hi = hi / SD_YEAR * 10) %>%
      left_join(spp_info, by = c("tree_name" = "species_name")) %>%
      mutate(spp = Binomial, pct_per_decade = mm_per_decade / mean_wing * 100) %>%   # `spp` = ABT Binomial as in the Tier-1 file
      select(model, spp, ranef_yr, se_ranef, slope, fixed, se_fixed, se_total, lo, hi,
             mm_per_decade, mm_per_decade_lo, mm_per_decade_hi, n, yr_min, yr_max, n_src, n_site, mean_wing, pct_per_decade,
             tree_name, slope_tree_min, slope_tree_max, n_trees) %>%
      arrange(slope) %>% as.data.frame()
  }
  slope_summary <- function(sl) if (is.null(sl)) NULL else list(
    n_species = nrow(sl), n_negative = sum(sl$slope < 0),
    n_excl_zero_negative = sum(sl$hi < 0), n_excl_zero_positive = sum(sl$lo > 0),
    median_mm_per_decade = median(sl$mm_per_decade),
    fixed_slope = sl$fixed[1], sd_random_slope = sd(sl$ranef_yr))
  sl_M3 <- pool_species_slopes("M3"); sl_M3cc <- pool_species_slopes("M3_cc"); sl_M0 <- pool_species_slopes("M0")
  if (!is.null(sl_M3)) message(sprintf("[glmmTMB M3] species slopes: %d of %d negative, %d intervals < 0, %d > 0, median %.2f mm/decade",
                                       sum(sl_M3$slope < 0), nrow(sl_M3), sum(sl_M3$hi < 0), sum(sl_M3$lo > 0), median(sl_M3$mm_per_decade)))
  slopes <- list(generated = Sys.time(), engine = sprintf("glmmTMB propto, %d trees, Rubin-pooled", N_TREES), n_trees = N_TREES,
                 M3 = sl_M3, M3_cc = sl_M3cc, M0 = sl_M0,
                 summary = list(M3 = slope_summary(sl_M3), M3_cc = slope_summary(sl_M3cc), M0 = slope_summary(sl_M0)),
                 sd_year = SD_YEAR,
                 note = paste("slope = pooled fixed scaled_yr + species random slope (glmmTMB REML with the phylogenetic propto",
                              "species intercept, one fit per tree). Per tree: se_total = sqrt(SE_fixed^2 + conditional SD_BLUP^2);",
                              "across trees: slope = mean, se = sqrt(mean(se^2) + (1 + 1/m) var(slope)) (Rubin). `spp` = ABT Binomial",
                              "(as in controlled_wing_species_slopes.rds), `tree_name` = tip name used in the fit. M3 = full controls on",
                              "the full sample, M3_cc = complete cases, M0 = baseline. Intervals are +/- 1.96 se_total."))
  saveRDS(slopes, out_path("controlled_wing_phylo_species_slopes.rds"))

  # --- decision gate (same rule as Tier 1, on the phylogenetic estimates) -------------
  g <- function(m, term = "scaled_yr") { r <- results %>% filter(model == m, term == !!term); if (nrow(r)) r else NULL }
  pick <- function(r, cols = c("estimate", "ci_lo", "ci_hi", "z", "N")) if (is.null(r)) NULL else r[, intersect(cols, names(r))]
  excl0 <- function(r) !is.null(r) && r$ci_hi < 0
  m0 <- g("M0"); m3 <- g("M3"); m3cc <- g("M3_cc")
  m5w <- g("M5", "yr_within_site"); m5wcc <- g("M5_cc", "yr_within_site"); m5b <- g("M5", "yr_site_mean")
  m5sw <- g("M5src", "yr_within_src"); m5sb <- g("M5src", "yr_src_mean")
  m5sw_alt <- g("M5src_rsTotal", "yr_within_src"); m5w_alt <- g("M5_rsTotal", "yr_within_site")
  gate_ok <- !is.null(m3) && !is.null(m5w)
  scenario_full <- if (!gate_ok) "not evaluable (M3 or M5 not fitted)" else if (excl0(m3) && excl0(m5w)) "A" else "B"
  scenario_cc   <- if (is.null(m3cc) || is.null(m5wcc)) "not evaluable" else if (excl0(m3cc) && excl0(m5wcc)) "A" else "B"
  decision_gate <- list(
    tier = sprintf("phylogenetic (glmmTMB propto, %d trees, Rubin-pooled Wald CI, REML)", N_TREES),
    rule = "Scenario A requires the controlled year effect (M3) AND the within-municipality year slope (M5, REWB) to exclude zero; otherwise B",
    scenario_full_sample = scenario_full, scenario_complete_cases = scenario_cc,
    scenario = if (identical(scenario_full, scenario_cc)) scenario_full else paste0(scenario_full, " (full) / ", scenario_cc, " (cc)"),
    M0 = pick(m0), M3 = pick(m3, c("estimate", "ci_lo", "ci_hi", "z", "N", "mm_per_decade")),
    M3_cc = pick(m3cc, c("estimate", "ci_lo", "ci_hi", "z", "N", "mm_per_decade")),
    M5_within_site = pick(m5w), M5_within_site_cc = pick(m5wcc), M5_within_site_rsTotal = pick(m5w_alt),
    M5_between_site = pick(m5b),
    M5src_within = pick(m5sw), M5src_within_rsTotal = pick(m5sw_alt), M5src_between = pick(m5sb),
    M6 = pick(g("M6")), M7 = pick(g("M7")), M_carrano = pick(g("M_carrano")),
    M3_noanom = pick(g("M3_noanom")),
    M5src_within_noanom = pick(g("M5src_noanom", "yr_within_src")),
    M5src_within_rsTotal_noanom = pick(g("M5src_rsTotal_noanom", "yr_within_src")),
    share_of_M0_remaining_in_M3 = if (!is.null(m0) && !is.null(m3)) m3$estimate / m0$estimate else NA_real_,
    share_of_M0_remaining_in_M3_cc = if (!is.null(m0) && !is.null(m3cc)) m3cc$estimate / m0$estimate else NA_real_)
  if (gate_ok) message(sprintf("[gate] phylogenetic scenario: full = %s, cc = %s | M3 = %.3f [%.3f, %.3f] | within-site = %.3f [%.3f, %.3f]%s",
                               scenario_full, scenario_cc, m3$estimate, m3$ci_lo, m3$ci_hi, m5w$estimate, m5w$ci_lo, m5w$ci_hi,
                               if (!is.null(m5sw)) sprintf(" | within-contributor = %.3f [%.3f, %.3f]", m5sw$estimate, m5sw$ci_lo, m5sw$ci_hi) else ""))

  # --- comparison with the Tier-1 lme4 estimates and the published brms M0 ------------
  comparison_tier1 <- NULL; tier1_file <- out_path("controlled_wing_results.rds")
  if (file.exists(tier1_file)) {
    t1 <- readRDS(tier1_file)
    comparison_tier1 <- before_after %>%
      select(model, term, N, phylo_est = estimate, phylo_se = se, phylo_lo = ci_lo, phylo_hi = ci_hi, phylo_prop) %>%
      inner_join(t1$before_after %>% select(model, term, N_lme4 = N, lme4_est = estimate, lme4_se = se, lme4_lo = ci_lo, lme4_hi = ci_hi),
                 by = c("model", "term")) %>%
      mutate(delta_est = phylo_est - lme4_est, se_ratio = phylo_se / lme4_se,
             excl0_phylo = phylo_hi < 0 | phylo_lo > 0, excl0_lme4 = lme4_hi < 0 | lme4_lo > 0,
             same_conclusion = excl0_phylo == excl0_lme4) %>% as.data.frame()
    message(sprintf("[compare] vs Tier 1 (lme4): max |delta estimate| = %.3f (%s %s); SE ratio range %.2f-%.2f; %d/%d rows same exclude-zero verdict",
                    max(abs(comparison_tier1$delta_est)),
                    comparison_tier1$model[which.max(abs(comparison_tier1$delta_est))], comparison_tier1$term[which.max(abs(comparison_tier1$delta_est))],
                    min(comparison_tier1$se_ratio), max(comparison_tier1$se_ratio), sum(comparison_tier1$same_conclusion), nrow(comparison_tier1)))
  } else message("[compare] controlled_wing_results.rds not found: run --fast-lme4 first for the Tier-1 comparison")
  comparison_published <- NULL; val_file <- out_path("glmmtmb_validation_wing.rds")
  if (!is.null(m0) && file.exists(val_file)) {
    v <- readRDS(val_file); br <- as.data.frame(v$brms_rubin)
    par_col <- intersect(c("par", "term", "parameter"), names(br))[1]
    br_yr <- br[grepl("scaled_yr", br[[par_col]]), , drop = FALSE]
    est_col <- intersect(c("estimate", "Estimate", "mean"), names(br))[1]
    lo_col <- intersect(c("lower", "l-95% CI", "Q2.5", "ci_lo"), names(br))[1]; hi_col <- intersect(c("upper", "u-95% CI", "Q97.5", "ci_hi"), names(br))[1]
    if (nrow(br_yr) == 1 && !is.na(est_col)) {
      comparison_published <- data.frame(
        source = c("published brms M0 (brm0_multiphylo.rda, Rubin over 50 trees)", "glmmTMB validation refit of the brms data (glmmtmb_validation_wing.R)",
                   sprintf("this script, M0, glmmTMB propto %d trees", N_TREES)),
        estimate = c(br_yr[[est_col]], v$pooled$estimate[v$pooled$par == "scaled_yr"], m0$estimate),
        ci_lo = c(if (!is.na(lo_col)) br_yr[[lo_col]] else NA, v$pooled$lower[v$pooled$par == "scaled_yr"], m0$ci_lo),
        ci_hi = c(if (!is.na(hi_col)) br_yr[[hi_col]] else NA, v$pooled$upper[v$pooled$par == "scaled_yr"], m0$ci_hi),
        phylo_prop = c(0.95, mean(v$varcomp$phylo_prop), m0$phylo_prop), stringsAsFactors = FALSE)
      print(comparison_published, digits = 4)
    }
  }

  total_secs <- as.numeric(Sys.time() - T_START, units = "secs")
  out <- list(
    generated = Sys.time(), mode = "trees", engine = list(
      name = "glmmTMB", version = as.character(packageVersion("glmmTMB")), TMB = as.character(packageVersion("TMB")),
      covstruct = "propto(0 + species_name | g, A): phylogenetic species intercept, A = per-tree species correlation matrix",
      trees = sprintf("%d trees; correlation matrices cached from brm0_multiphylo.rda (data/derived/phylo_A_50trees.rds), i.e. the published brms trees", N_TREES),
      pooling = "Rubin's rules over trees (Nakagawa & de Villemereuil 2019): estimate = mean; se^2 = mean(se^2) + (1 + 1/m) var(estimate)",
      estimation = "REML; Wald intervals +/- 1.96 pooled SE; z = estimate / pooled SE (column t repeats z)",
      restarts = sprintf("fits with a non-PD Hessian from glmmTMB's default start (i.i.d. side of the species-intercept ridge) are refitted from the tree's own parameter layout with the propto scale started at SD %g; accepted only if PD and the REML objective is not worse (per-tree record: $per_model[[model]]$retry, all models: $retry)", PHYLO_START_SD),
      reference = "Williams, McGillycuddy, Drobniak, Bolker, Warton & Nakagawa 2025, bioRxiv 10.64898/2025.12.20.695312; validated on this data by glmmtmb_validation_wing.R",
      engine_script = "scripts/_phylo_engine.R", dropped_species = dropped_species, binomial_to_tree_name = as.data.frame(binomial_map),
      R = R.version.string, dplyr = as.character(packageVersion("dplyr")), host = Sys.info()[["nodename"]]),
    n_trees = N_TREES, models = PHYLO_MODELS, errors = errors,
    sd_year = SD_YEAR, mean_wing_mm = sample_sizes$mean_wing_mm, sample_sizes = sample_sizes,
    specs = bind_rows(lapply(PHYLO_MODELS, function(n) data.frame(model = n, label = SPECS[[n]]$label, fixed = SPECS[[n]]$fx,
                                                                   random = SPECS[[n]]$re, sample = SPECS[[n]]$data))),
    table = results,             # all cond fixed effects, all models, pooled
    before_after = before_after, # year terms only (+ n_trees, converged, secs, phylo_prop)
    varcomp = varcomp,           # long: model, component, sd (mean over trees), sd_lo, sd_hi
    ridge = ridge,               # per model: trees on the phylo vs i.i.d. side of the species-intercept ridge, year estimate by side
    retry = retry_all,           # per model x tree: pdHess and REML objective from the default start and after the informed restart
    per_model = phylo_runs,      # pooled (cond + disp), per_tree_fixed, varcomp per tree (+ mode, species_level_prop), varcomp_summary, retry, n_trees, secs
    decision_gate = decision_gate,
    comparison_tier1 = comparison_tier1, comparison_published_M0 = comparison_published,
    timings_sec = timings, total_secs = total_secs,
    note = paste("glmmTMB phylogenetic tier of the artefact-controlled wing model (REVISION_PLAN.md Phase 1).",
                 "Every Tier-1 specification refitted with the propto phylogenetic species intercept on each of the",
                 N_TREES, "trees of the published analysis and Rubin-pooled. CIs are pooled Wald (+/- 1.96 SE).",
                 "mm_per_decade = estimate / SD(year) * 10 with SD(year) = 5.0199. Species names are tree tip names",
                 "(phylo_species_name); no species dropped on 2026-09-09. Same samples and guards as Tier 1 (controlled_wing_results.rds)."))
  saveRDS(out, out_path("controlled_wing_phylo_results.rds"))
  message(sprintf("[write] %s and controlled_wing_phylo_species_slopes.rds (total %.1f min, %d models, %d fits)",
                  out_path("controlled_wing_phylo_results.rds"), total_secs / 60, length(timings), sum(before_after$n_trees[!duplicated(before_after$model)])))
  print(before_after %>% select(model, term, estimate, ci_lo, ci_hi, z, N, mm_per_decade, phylo_prop, converged, converged_default_start, secs), digits = 3, row.names = FALSE)
  quit(save = "no", status = 0)
}

# ===========================================================================
# 5. Tiers 2-3 (and --smoke), --engine brms: brms with the phylogenetic term, Rubin-pooled
# ===========================================================================
suppressMessages({
  library(brms); library(ape); library(MCMCglmm); library(prepR4pcm)
  library(phytools); library(future.apply); library(posterior)
})
source(script_path("_sampling_config.R"))   # SAMPLING, SAMPLING_CONTROL (host-aware)
if (MODE == "smoke") {
  SAMPLING$chains <- 2L; SAMPLING$iter <- 400L; SAMPLING$warmup <- 200L
  SAMPLING$cores <- 2L; SAMPLING$workers <- 1L
  message("[smoke] chains = 2, iter = 400, warmup = 200, 1 tree, model M3 only. NOT A RESULT.")
}
# brms model set draws on the requested sample
brms_sfx <- if (BRMS_SAMPLE == "cc") "_cc" else ""
brms_spec_name <- function(nm) {
  if (nm %in% c("M0", "M1")) return(if (BRMS_SAMPLE == "cc") paste0(nm, "_cc") else nm)
  cand <- paste0(nm, brms_sfx); if (cand %in% names(SPECS)) cand else nm
}

# --- tree retrieval and species reconciliation (copied from atlantic_parallel.R) ---
passer90$species_name <- gsub("_", " ", passer90$Binomial)
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
  "Tiaris fuliginosus"       = "Asemospiza fuliginosa"
)
hit <- passer90$species_name %in% names(ebird_synonyms)
passer90$species_name[hit] <- ebird_synonyms[passer90$species_name[hit]]
spp_data <- unique(passer90$species_name)

# AvesData repo: the 100-tree dated sample sets. Use the copy already in
# data/raw/AvesDataLite-main when AVESDATA_PATH is unset (avoids a 700 MB
# download); otherwise fall back to clootl's downloader exactly as atlantic_parallel.R.
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
  local_aves <- raw_path("AvesDataLite-main")
  if (dir.exists(local_aves)) clootl::set_avesdata_repo_path(local_aves, overwrite = TRUE)
  else clootl::get_avesdata_repo(path = raw_path())
}
.check <- pr_get_tree(spp_data, source = "clootl", n_tree = 1)
if (length(.check$unmatched) > 0) {
  message("Removing ", length(.check$unmatched), " species absent from eBird taxonomy: ",
          paste(.check$unmatched, collapse = ", "))
  spp_data <- setdiff(spp_data, .check$unmatched)
}
stopifnot(length(.check$unmatched) <= 1)
got   <- pr_get_tree(spp_data, source = "clootl", n_tree = 100, cache = TRUE)
trees <- got$tree
clootl_version <- pr_cite_tree(got, format = "text")

rec <- reconcile_tree(x = passer90, tree = trees[[1]], x_species = "species_name",
                      fuzzy = TRUE, resolve = "flag")
print(reconcile_summary(rec))
aligned <- reconcile_apply(rec, data = passer90, tree = trees[[1]],
                           species_col = "species_name", drop_unresolved = TRUE)
# Binomial -> tip label map, applied to every analysis frame (incl. unknown-sex)
name_map <- aligned$data %>% transmute(Binomial, species_name = gsub(" ", "_", species_name)) %>% distinct()
attach_tree_names <- function(d) {
  d %>% select(-any_of("species_name")) %>% inner_join(name_map, by = "Binomial") %>%
    mutate(spp = factor(species_name)) %>%                # both grouping terms use matched names
    select(-yr_site_mean, -yr_src_mean, -lat_spp_mean, -yr_within_site, -yr_within_src, -lat_within_spp) %>%
    finish()                                              # re-derive Mundlak terms on the matched species set
}
norm_us   <- function(x) gsub(" ", "_", x)
keep_tips <- norm_us(aligned$tree$tip.label)
trees_pruned <- lapply(trees, function(t) {
  t$tip.label <- norm_us(t$tip.label)
  ape::keep.tip(t, intersect(keep_tips, t$tip.label))
})
class(trees_pruned) <- "multiPhylo"
tree_samp <- sample(trees_pruned, N_TREES)        # seed set above; 50 on the server

make_A <- function(tree) {
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A   <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

priors <- c(prior(normal(71, 15), class = Intercept),
            prior(normal(0, 10),  class = b),
            prior(cauchy(0, 1),   class = sd),
            prior(cauchy(0, 1),   class = sigma))

brms_frames <- lapply(frames, attach_tree_names)
brms_formula <- function(s) {
  as.formula(paste("conc.wing.length ~ 1 +", s$fx, "+", s$re, "+ (1 | gr(species_name, cov = A))"))
}
fit_one <- function(tree, s) {
  A <- make_A(tree)
  brm(brms_formula(s), data = brms_frames[[s$data]], data2 = list(A = A),
      family = gaussian(), prior = priors,
      iter = SAMPLING$iter, warmup = SAMPLING$warmup, chains = SAMPLING$chains,
      cores = SAMPLING$cores, backend = SAMPLING$backend, control = SAMPLING_CONTROL,
      seed = 20240101)
}

# --- Rubin's rules pooling across trees (copied from atlantic_parallel.R, generalised
# to every population-level parameter) ---
# For each fixed effect: pooled estimate = mean of per-tree posterior means;
# total variance = within-tree var (mean of per-tree posterior variances)
# + between-tree var (variance of per-tree means) inflated by (1 + 1/m).
# (Nakagawa & de Villemereuil 2019, Syst. Biol. 68:632-641.)
pool_rubin <- function(fits, pars = NULL) {
  m <- length(fits)
  if (is.null(pars)) pars <- grep("^b_", variables(fits[[1]]), value = TRUE)
  draws <- lapply(fits, function(f) as_draws_df(f)[, pars, drop = FALSE])
  means <- t(sapply(draws, function(d) sapply(d, mean)))
  vars  <- t(sapply(draws, function(d) sapply(d, var)))
  if (m == 1) { means <- matrix(means, nrow = 1, dimnames = list(NULL, pars)); vars <- matrix(vars, nrow = 1, dimnames = list(NULL, pars)) }
  qbar  <- colMeans(means)                 # pooled point estimate
  ubar  <- colMeans(vars)                  # within-imputation variance
  b     <- if (m > 1) apply(means, 2, var) else 0 * qbar   # between-tree variance (0 for one tree)
  tot   <- ubar + (1 + 1/m) * b            # total variance
  se    <- sqrt(tot)
  data.frame(par = pars, estimate = qbar, se = se,
             lower = qbar - 1.96*se, upper = qbar + 1.96*se, m_trees = m, row.names = NULL)
}
diag_one <- function(f) {
  s <- summary(f)$fixed
  data.frame(par = paste0("b_", rownames(s)), rhat = s$Rhat, bulk_ess = s$Bulk_ESS, tail_ess = s$Tail_ESS,
             divergent = sum(nuts_params(f, pars = "divergent__")$Value), row.names = NULL)
}

plan(multisession, workers = SAMPLING$workers)
brms_results <- list(); timings <- c(); model_files <- c()
for (nm in BRMS_MODELS) {
  sn <- brms_spec_name(nm); s <- SPECS[[sn]]; if (is.null(s)) stop("unknown model ", nm)
  message(sprintf("[brms] %s (spec %s, N = %d) on %d tree(s): %s", nm, sn, nrow(brms_frames[[s$data]]), N_TREES, deparse1(brms_formula(s))))
  t0 <- Sys.time()
  fits <- future_lapply(tree_samp, fit_one, s = s, future.seed = TRUE)
  timings[nm] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  pooled <- pool_rubin(fits) %>%
    mutate(model = nm, spec = sn, label = s$label, N = nrow(brms_frames[[s$data]]),
           year_term = sub("^b_", "", par) %in% YEAR_TERMS,
           mm_per_decade = ifelse(year_term, estimate / SD_YEAR * 10, NA_real_),
           mm_per_decade_lo = ifelse(year_term, lower / SD_YEAR * 10, NA_real_),
           mm_per_decade_hi = ifelse(year_term, upper / SD_YEAR * 10, NA_real_))
  diags <- bind_rows(lapply(seq_along(fits), function(i) cbind(tree = i, diag_one(fits[[i]]))))
  brms_results[[nm]] <- list(pooled = pooled, diagnostics = diags, minutes = timings[nm])
  print(pooled %>% select(par, estimate, lower, upper, mm_per_decade), digits = 3)
  model_files[nm] <- if (MODE == "smoke") out_path("models", sprintf("controlled_smoke_%s.rda", nm))
                     else out_path("models", sprintf("controlled_%s%s.rda", nm, brms_sfx))
  rubin_summary <- pooled
  save(fits, rubin_summary, diags, clootl_version, file = model_files[nm])
  message("[write] ", model_files[nm], sprintf(" (%.1f min)", timings[nm]))
}

# --- species-slope posteriors from M3 (caterpillar figure) -------------------------
species_posteriors <- NULL
if ("M3" %in% BRMS_MODELS) {
  load(model_files["M3"])
  # slope_j = b_scaled_yr + r_spp[j, scaled_yr]; species-level effects are summarised on the
  # posterior draws stacked across trees (they are not Rubin-pooled)
  sl <- lapply(fits, function(fit) {
    d <- as_draws_df(fit); rs <- grep("^r_spp__scaled_yr\\[", names(d), value = TRUE)
    m <- as.matrix(d[, rs]) + d$b_scaled_yr
    colnames(m) <- sub("^r_spp__scaled_yr\\[(.*),Intercept\\]$", "\\1", rs); m
  })
  sl <- do.call(rbind, sl)
  species_posteriors <- data.frame(spp = colnames(sl), slope = colMeans(sl),
                                   lo = apply(sl, 2, quantile, 0.025), hi = apply(sl, 2, quantile, 0.975),
                                   p_negative = colMeans(sl < 0), row.names = NULL) %>%
    mutate(mm_per_decade = slope / SD_YEAR * 10) %>% arrange(slope)
  message(sprintf("[brms M3] species slopes: %d negative of %d, %d with 95%% interval < 0",
                  sum(species_posteriors$slope < 0), nrow(species_posteriors), sum(species_posteriors$hi < 0)))
}

out <- list(generated = Sys.time(), mode = MODE, n_trees = N_TREES, models = BRMS_MODELS, sample = BRMS_SAMPLE,
            sampling = SAMPLING, clootl_version = clootl_version,
            n_species_tree = length(keep_tips),
            results = brms_results,
            pooled = bind_rows(lapply(brms_results, `[[`, "pooled")),
            species_slopes_M3 = species_posteriors,
            sd_year = SD_YEAR, mean_wing_mm = sample_sizes$mean_wing_mm, sample_sizes = sample_sizes,
            session = list(R = R.version.string, brms = as.character(packageVersion("brms")),
                           prepR4pcm = as.character(packageVersion("prepR4pcm")), clootl = as.character(packageVersion("clootl"))),
            note = if (MODE == "smoke") "SMOKE TEST (1 tree, 2 chains, 400 iterations): code-path check only, NOT a result."
                   else "Rubin-pooled across trees (Nakagawa & de Villemereuil 2019); intervals are pooled-SE normal intervals.")
rds <- if (MODE == "smoke") out_path("controlled_wing_brms_smoke.rds") else out_path(sprintf("controlled_wing_brms_results%s.rds", brms_sfx))
saveRDS(out, rds)
message("[write] ", rds)
