# atlantic_drmsem.R
# ---------------------------------------------------------------------------
# WHAT: Distributional piecewise SEM (drmSEM) of wing length in live-measured
# adult Atlantic Forest passerines (73 species, 1990-2018). Two endogenous
# nodes: record-level temperature (WorldClim Tmean at the capture coordinates
# and year <- year + latitude) and wing length, whose MEAN is modelled on
# temperature, sex and latitude with species AND contributor (Main_researcher)
# random intercepts, and whose residual SD (sigma) is modelled on year and
# temperature. Section 10 refits the wing node with a phylogenetic species
# effect (Pagel's lambda, 50 trees, Rubin pooling).
#
# STATUS - SUPPLEMENTARY, EXPLORATORY (2026-09, revision Phase 4f;
# REVISION_PLAN.md S1 "S2.7 drmSEM" and S3 "Phase 4f"). This script is NOT the
# primary evidence for any claim in the manuscript:
#   * the primary wing trend is the artefact-controlled brms model
#     (atlantic_parallel_controlled.R, Phase 1);
#   * the primary variance evidence is the brms distributional sigma model
#     with contributor / site / sex terms (atlantic_variance_sigma.R, Phase 3).
# drmSEM is reported in the Supplement only, as a check that the sign of the
# sigma-channel paths agrees under an explicit graph, with intervals (no
# P-values in the main text: one inferential standard throughout).
#
# AUTHORSHIP / SOFTWARE DISCLOSURE (must appear in the Supplement): drmSEM
# (and drmTMB, its fitting engine) and prepR4pcm (tree retrieval and name
# reconciliation, section 10) are co-authored by E.S.A.S. and, per
# REVISION_PLAN S2.7, are to be disclosed as unpublished software distributed
# as development versions (r-universe / GitHub). Versions used for this edit:
# R 4.6.0, drmTMB 0.1.4, drmSEM 0.5.0, prepR4pcm 0.5.0.9000, brms 2.23.0.
# A short validation note (drmSEM i.i.d. paths vs the brms fits on the same
# data, same signs / magnitudes) accompanies the Supplement table.
#
# D-SEPARATION CLAIMS TO BE REPORTED (Shipley's test, Fisher's C):
#   Full-record graph (the ONLY graph reported):
#     year -> tmean <- lat ;  tmean, Sex, lat -> mean(wing) ;
#     year, tmean -> sigma(wing) ;  random intercepts (1|spp), (1|Main_researcher)
#     [and relmat() in section 10] are stripped from the causal edge set by
#     drmSEM and add NO claims (verified 2026-09-09 on a 1,500-record subset:
#     the claim set is identical with and without the contributor intercept).
#     Exactly ONE missing edge => one claim on Fisher's C with df = 2:
#        Sex _||_ scaled_tmean | {scaled_yr, scaled_lat}
#     Report: C, df = 2, P, the LR of the single claim, and N / species.
#     (Last full run before this edit, 2026-06-15, N = 7,810, 72 spp, no
#      contributor intercept: C = 0.51, df = 2, P = 0.78 - superseded; refit.)
#   Retired arthropod graph (NOT reported, see MODEL SCOPE below): two claims
#     (Sex _||_ tmean | yr, lat ; Sex _||_ arthro | yr, tmean, lat) on df = 4;
#     rejected in June 2026 (C = 13.5, P = 0.009). No coefficient from a
#     rejected graph is quoted anywhere in the revision.
#
# SIGMA-CHANNEL WORDING (referee S2.7): the sigma submodel gives TWO separable
# effects - later years -> larger residual SD, warmer records -> smaller
# residual SD. Year and record temperature are essentially uncorrelated in
# these data (year -> tmean path ~ 0 in the June run), so the variance rise
# is NOT "linked to thermal extremes"; the manuscript text must present the
# two sigma paths as independent associations.
#
# WHY the wing node carries (1 | Main_researcher): the revision's central
# finding is that observer turnover generates apparent trends (REVISION_PLAN
# S0). Any model whose sigma submodel contains year must therefore absorb
# between-contributor mean differences, otherwise a change in who measured
# masquerades as a change in residual variance. drmTMB accepts a second
# i.i.d. random intercept (tested 2026-09-09, also together with relmat()).
# A contributor term in the sigma submodel itself is not attempted here (the
# brms sigma model in Phase 3 carries it).
#
# INPUTS:
#   data/derived/passer90.rda          live-only, known-sex, 73 spp, with
#                                      Main_researcher / Municipality / ID_ABT
#                                      (rebuild_passer90_live.R, Phase 0)
#   data/derived/passer90_climate.rds  the same records + record-level WorldClim
#                                      Tmean (climate_extraction.R). MUST carry
#                                      ID_ABT: a file without it predates the
#                                      Phase 0 rebuild (89-species, museum
#                                      specimens included) and is refused.
#   [retired] archive/predicts/dat_no_grassland.rds - only with the explicit
#                                      --arthro-override flag (see MODEL SCOPE)
# OUTPUTS (Analysis/output, mode-suffixed "_noarthro"):
#   drmsem_results_<mode>.rds / .md    Fisher's C, d-sep claims, path tables,
#                                      effect decomposition, pooled phylo paths
#   drmsem_effects_cache_<mode>.rds, drmsem_phylo_prep_<mode>.rds,
#   drmsem_phylo_results_<mode>.rds    caches (git-ignored)
#   figures/drmsem_{dag,coef_forest,effects_forest,sigma_curves}_<mode>.png
#
# RUN:
#   Rscript Analysis/scripts/atlantic_drmsem.R                # full supplement run
#   Rscript Analysis/scripts/atlantic_drmsem.R --no-phylo     # skip section 10
#   Rscript Analysis/scripts/atlantic_drmsem.R --smoke        # 1,500-record subset,
#        no phylo, tiny B; writes to a temp dir, touches NOTHING in the repo
#   Rscript Analysis/scripts/atlantic_drmsem.R --allow-stale-climate
#        fall back to passer90.rda + STATIC WorldClim climatology when
#        passer90_climate.rds is stale (spatial, not temporal, temperature)
#   Rscript Analysis/scripts/atlantic_drmsem.R --arthro-override=PREDICTS_RETIRED
#        re-run the RETIRED 3-node PREDICTS graph for provenance only
#
# REQUIREMENTS:
#   pak::pak(c("drmTMB", "drmSEM"), repos = "https://itchyshin.r-universe.dev")
#   Section 10 additionally: pak::pak("itchyshin/prepR4pcm");
#   install.packages(c("clootl", "phytools", "ape")) - the same tree-retrieval
#   / name-reconciliation stack as atlantic_parallel.R.
#
# HISTORY: v3 (2026-06) fitted an endogenous PREDICTS arthropod node in a
# 3-node graph; that graph was rejected by d-separation and PREDICTS itself
# was retired in Phase 4e (archive/predicts/README.md). v4 (2026-09): arthropod
# node off with an explicit override, contributor intercept on the wing node,
# climate-file staleness guard, --smoke mode, results carry the disclosure.
# ---------------------------------------------------------------------------

# --- Installation (run once) ------------------------------------------------
# install.packages("pak")
# pak::pak("itchyshin/drmTMB")
# pak::pak("itchyshin/drmSEM")

library(drmTMB)
library(drmSEM)
library(brms)
library(dplyr)
library(ggplot2)

# ============================================================================
# 0. PATH HELPERS  (mirrors pattern from atlantic_parallel.R)
# ============================================================================

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

# ============================================================================
# CONFIG: simulation cost + result caching
# ----------------------------------------------------------------------------
# The simulate-based effect calls in section 7 cost ~ B * EFFECT_NSIM * N_rows.
#   EFFECT_NSIM = inner distributional draws; feeds ONLY the distribution-
#                 mediated leg. 2000 was overkill (default is 50) — it sharpens
#                 a quantity that is ~0 here. 200 is ample.
#   EFFECT_B    = coefficient draws; controls CI width across all legs.
# Lower both for quick iteration, raise EFFECT_B for the final/publication run.
# SEM_SEED makes the Monte-Carlo effects reproducible (they were not before).
# Effect tables are cached to CACHE_FILE and reused unless the SEM data changes
# (or RECOMPUTE = TRUE), so re-sourcing the script does not re-simulate.
# ============================================================================
# --- COMMAND-LINE FLAGS -------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.flag <- function(f) any(.args == f)
.flag_value <- function(prefix) {
  hit <- grep(paste0("^", prefix, "="), .args, value = TRUE)
  if (length(hit)) sub(paste0("^", prefix, "="), "", hit[1]) else ""
}
SMOKE               <- .flag("--smoke")
NO_PHYLO            <- .flag("--no-phylo")
ALLOW_STALE_CLIMATE <- .flag("--allow-stale-climate") || SMOKE
ARTHRO_OVERRIDE     <- identical(.flag_value("--arthro-override"), "PREDICTS_RETIRED")

# --- MODEL SCOPE ------------------------------------------------------------
# The PREDICTS arthropod node is RETIRED (revision Phase 4e, referee S2.6):
# 53 of 56 PREDICTS studies are single-year snapshots and the 17 Atlantic
# Forest studies span 1998-2009 only, so the "arthropod decline" was a contrast
# between independent studies, not a time series; the 3-node graph was also
# rejected by d-separation. The inputs now live in archive/predicts/ (README
# there). INCLUDE_ARTHRO is therefore FALSE and CANNOT be switched on by
# editing this file: the guard below stops the script if the constant is set
# TRUE without the explicit command-line override
#     --arthro-override=PREDICTS_RETIRED
# which exists only to reproduce the archived June-2026 numbers for the
# response letter. Nothing from that model is reported.
#
# INCLUDE_ARTHRO = FALSE -> full-record temperature / variance graph (2 nodes),
#                          ALL live-measured records 1990-2018. THE model.
# INCLUDE_ARTHRO = TRUE  -> retired 3-node graph restricted to years with
#                          PREDICTS coverage (~1998-2008). Override only.
#
# FUTURE HOOK (Phase 4a/4b): if a locality x year litter-ant or GBIF occupancy
# index becomes available (atlantic_ants_index.R / atlantic_gbif_occupancy.R,
# Analysis/data/derived/ants_*.rds), it should enter as a NEW node
#   ants_index ~ scaled_yr + scaled_lat (+ contributor / method terms)
# joined on Municipality x Year, NOT through this PREDICTS switch. Not
# implemented here because the P4a outputs are not final.
INCLUDE_ARTHRO <- FALSE
if (isTRUE(INCLUDE_ARTHRO) && !ARTHRO_OVERRIDE) {
  stop("INCLUDE_ARTHRO is retired (Phase 4e). Re-run with ",
       "--arthro-override=PREDICTS_RETIRED if you really need the archived ",
       "3-node PREDICTS graph; nothing from it may be reported.", call. = FALSE)
}
if (ARTHRO_OVERRIDE) {
  INCLUDE_ARTHRO <- TRUE
  warning("RETIRED PREDICTS arthropod node enabled by explicit override. ",
          "This graph was rejected by d-separation (June 2026) and PREDICTS ",
          "carries no within-study temporal information; for provenance only.",
          call. = FALSE)
}
.mode <- if (INCLUDE_ARTHRO) "arthro" else "noarthro"

# --- SMOKE MODE ---------------------------------------------------------------
# --smoke: 1,500-record subsample, no phylogenetic section, tiny Monte-Carlo
# budgets, and EVERY output (caches, results, figures) redirected to a temp
# directory so a smoke run can never overwrite the files the manuscript reads.
if (SMOKE) {
  SMOKE_DIR <- file.path(tempdir(), "drmsem_smoke")
  dir.create(SMOKE_DIR, showWarnings = FALSE, recursive = TRUE)
  out_path <- function(...) file.path(SMOKE_DIR, ...)
  fig_path <- function(...) file.path(SMOKE_DIR, "figures", ...)
  dir.create(fig_path(), showWarnings = FALSE, recursive = TRUE)
  message("SMOKE MODE: subset fit, no phylo, outputs -> ", SMOKE_DIR,
          " (nothing in the repo is written). Numbers are NOT results.")
}

# --- Spatial arthropod aggregation ------------------------------------------
# Instead of a single Brazil-wide year mean, we aggregate PREDICTS sites within
# ARTHRO_RADIUS_KM of each bird locality and compute an inverse-distance-weighted
# mean for that locality × year. This gives each bird a more spatially relevant
# food-availability index and introduces within-year variation (birds at different
# localities get different values), partially addressing the level-mismatch.
# ARTHRO_MIN_SITES: minimum number of PREDICTS sites needed for a locality × year
# cell; cells below this threshold fall back to the AF-region year mean.
ARTHRO_RADIUS_KM <- 200L   # search radius for nearby PREDICTS sites
ARTHRO_MIN_SITES <-   2L   # min sites required before using local estimate

SEM_SEED    <- 20250611
EFFECT_B    <- if (SMOKE) 20L else 200L   # coefficient draws (CI precision); 50 for dev
EFFECT_NSIM <- if (SMOKE) 20L else 200L   # inner sims, distribution-mediated leg
RECOMPUTE   <- SMOKE    # TRUE = ignore the cache and recompute
CACHE_FILE  <- out_path(paste0("drmsem_effects_cache_", .mode, ".rds"))
# MODEL_TAG (cache key for node FORMULAS) is defined at the end of section 3,
# once it is known whether the contributor intercept is available (HAS_SRC).

# --- Section 10 (phylogenetic SEM, Route A) ---------------------------------
RUN_PHYLO        <- !NO_PHYLO && !SMOKE   # --no-phylo / --smoke skip the phylo variant
N_PHYLO_TREES    <- 50L                   # trees sampled for Rubin pooling
LAMBDA_GRID      <- c(0, 0.25, 0.5, 0.75, 0.9, 1)  # Pagel's lambda grid (AIC-selected)
PHYLO_PREP_CACHE <- out_path(paste0("drmsem_phylo_prep_", .mode, ".rds"))   # reconciled trees + matched data
PHYLO_CACHE      <- out_path(paste0("drmsem_phylo_results_", .mode, ".rds")) # pooled paths + per-tree lambda
RECOMPUTE_PHYLO  <- FALSE                 # TRUE = ignore phylo caches and refit

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Bird data. passer90.rda (Phase 0 rebuild) is ALWAYS loaded: it defines the
# live-only analytical sample and carries the provenance columns. The
# time-resolved temperature comes from passer90_climate.rds when that file is
# current; a stale file (built before the Phase 0 rebuild, hence without
# ID_ABT: 89 species, museum specimens included) is refused unless
# --allow-stale-climate, in which case the STATIC WorldClim climatology in
# passer90.rda is used instead (spatial, not temporal, temperature).
load(derived_path("passer90.rda"))   # -> passer90
.prov_cols <- intersect(c("ID_ABT", "Main_researcher", "Municipality", "Locality",
                          "Ring", "Recapture", "season", "wing_col", "Status"),
                        names(passer90))

.static_fallback <- function() {
  warning(
    "Using STATIC Annual_mean_temperature from passer90.rda.\n",
    "This is a per-locality climatology: it captures spatial but NOT temporal\n",
    "temperature variation, so the year -> temperature -> wing path reflects\n",
    "spatial Bergmann's rule, not the temporal trend. Run climate_extraction.R\n",
    "(live-only sample, Phase 5) to regenerate passer90_climate.rds.",
    call. = FALSE)
  b <- passer90
  b$scaled_tmean <- as.numeric(scale(b$Annual_mean_temperature))
  b
}

climate_path <- derived_path("passer90_climate.rds")
TEMP_SOURCE  <- NA_character_
if (file.exists(climate_path)) {
  birds <- readRDS(climate_path)
  if (!"ID_ABT" %in% names(birds)) {
    msg <- sprintf(paste0(
      "passer90_climate.rds is STALE: %d rows / %d species and no ID_ABT column, ",
      "i.e. built before the Phase 0 live-only rebuild (REVISION_PLAN S1 ",
      "'Housekeeping'). Re-run climate_extraction.R, or pass ",
      "--allow-stale-climate to fall back to the static climatology."),
      nrow(birds), length(unique(birds$spp)))
    if (!ALLOW_STALE_CLIMATE) stop(msg, call. = FALSE)
    message(msg)
    birds <- .static_fallback(); TEMP_SOURCE <- "static_worldclim_climatology"
  } else {
    # Guarantee the live-only sample and attach any provenance column the
    # climate script did not carry (it inherits passer90's columns, so this is
    # normally a no-op).
    birds <- birds %>% semi_join(passer90 %>% select(ID_ABT), by = "ID_ABT")
    .missing <- setdiff(.prov_cols, names(birds))
    if (length(.missing)) {
      birds <- birds %>%
        left_join(passer90 %>% select(all_of(c("ID_ABT", .missing))), by = "ID_ABT")
    }
    TEMP_SOURCE <- "worldclim_record_level"
    message("Using time-resolved temperature (scaled_tmean from WorldClim), ",
            nrow(birds), " live-only records.")
  }
} else {
  message("passer90_climate.rds not found.")
  birds <- .static_fallback(); TEMP_SOURCE <- "static_worldclim_climatology"
}

# PREDICTS arthropod records: RETIRED (Phase 4e). Loaded ONLY under the explicit
# override, from archive/predicts/ (or the old derived/ path if a copy remains).
REPO_DIR <- dirname(ANALYSIS_DIR)
if (INCLUDE_ARTHRO) {
  .cand <- c(file.path(REPO_DIR, "archive", "predicts", "dat_no_grassland.rds"),
             derived_path("dat_no_grassland.rds"))
  .hit  <- .cand[file.exists(.cand)]
  if (!length(.hit)) {
    stop("dat_no_grassland.rds not found in archive/predicts/ or data/derived/. ",
         "The PREDICTS inputs were archived in Phase 4e; see ",
         "archive/predicts/README.md.", call. = FALSE)
  }
  load(.hit[1])   # -> dat_no_grassland
  message("Loaded RETIRED PREDICTS records from ", .hit[1])
} else {
  message("PREDICTS arthropod data not loaded: the node is retired and its ",
          "inputs are archived under archive/predicts/ (Phase 4e).")
}

# ============================================================================
# 2. BUILD SPATIALLY-RESOLVED ARTHROPOD AGGREGATE (AF-region, locality × year)
# ============================================================================
# Strategy: instead of a single Brazil-wide year mean, we build a two-level
# arthropod index.
#
#   Level 1 — per PREDICTS site × year mean log(abundance / effort).
#   Level 2 — for each bird locality × year, inverse-distance-weighted mean of
#             all PREDICTS sites within ARTHRO_RADIUS_KM. Bird localities with
#             fewer than ARTHRO_MIN_SITES nearby sites fall back to the AF-region
#             year mean (all sites within the AF bounding box).
#
# This introduces within-year spatial variation in the arthropod index (birds at
# different localities get different values) and restricts the PREDICTS data to
# the Atlantic Forest region, partially addressing the level-mismatch described
# in the design notes.
#
# Sample_midpoint is a Date stored as integer days since 1970-01-01; extract the
# calendar year with format(), NOT as.integer().

if (INCLUDE_ARTHRO) {

  # Haversine distance in km (vectorised over pred sites for a single bird loc)
  .haversine_km <- function(lon1, lat1, lon2, lat2) {
    to_rad <- pi / 180
    dlat <- (lat2 - lat1) * to_rad
    dlon <- (lon2 - lon1) * to_rad
    a <- sin(dlat / 2)^2 +
         cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
    6371 * 2 * asin(pmin(1, sqrt(a)))
  }

  # --- Level 1: per site × year mean ----------------------------------------
  site_year_arthro <- dat_no_grassland %>%
    mutate(
      Year          = as.integer(format(Sample_midpoint, "%Y")),
      log_abund_eff = log((Measurement / Sampling_effort) + 0.01)
    ) %>%
    filter(Year >= 1990, Year <= 2018) %>%
    group_by(SSBS, Longitude, Latitude, Year) %>%
    summarise(site_mean = mean(log_abund_eff, na.rm = TRUE),
              n_recs    = dplyr::n(),
              .groups   = "drop")

  # AF-region bounding box (same as bird locality cloud; used for fallback)
  af_lon <- range(birds$Longitude_decimal_degrees, na.rm = TRUE)
  af_lat <- range(birds$Latitude_decimal_degrees,  na.rm = TRUE)
  af_sites <- site_year_arthro %>%
    filter(Longitude >= af_lon[1] - 0.5, Longitude <= af_lon[2] + 0.5,
           Latitude  >= af_lat[1] - 0.5, Latitude  <= af_lat[2] + 0.5)

  # AF-region fallback: year mean across all AF sites (≥5 records)
  arthro_af_year <- af_sites %>%
    group_by(Year) %>%
    summarise(arthro_obs_mean = mean(site_mean, na.rm = TRUE),
              n_arthro        = dplyr::n(),
              .groups         = "drop") %>%
    filter(n_arthro >= 5)

  # --- Level 2: locality × year IDW mean ------------------------------------
  bird_locs <- birds %>%
    filter(!is.na(Longitude_decimal_degrees), !is.na(Latitude_decimal_degrees)) %>%
    distinct(Longitude_decimal_degrees, Latitude_decimal_degrees)

  arthro_by_loc_year <- do.call(rbind, lapply(seq_len(nrow(bird_locs)), function(i) {
    blon <- bird_locs$Longitude_decimal_degrees[i]
    blat <- bird_locs$Latitude_decimal_degrees[i]
    d_km <- .haversine_km(blon, blat, af_sites$Longitude, af_sites$Latitude)
    nearby <- af_sites[d_km <= ARTHRO_RADIUS_KM, ]
    if (nrow(nearby) == 0) return(NULL)
    d_near  <- d_km[d_km <= ARTHRO_RADIUS_KM]
    weights <- 1 / (d_near + 1)   # +1 avoids division by zero at d=0
    nearby %>%
      mutate(w = weights) %>%
      group_by(Year) %>%
      summarise(
        arthro_obs_mean = sum(site_mean * w) / sum(w),
        n_arthro        = dplyr::n(),
        .groups         = "drop"
      ) %>%
      filter(n_arthro >= ARTHRO_MIN_SITES) %>%
      mutate(Longitude_decimal_degrees = blon,
             Latitude_decimal_degrees  = blat,
             arthro_source             = "local")
  }))

  message(sprintf(
    "Arthropod spatial index: %d locality × year cells with local estimate (radius %d km); %d AF-year fallback cells.",
    nrow(arthro_by_loc_year),
    ARTHRO_RADIUS_KM,
    nrow(arthro_af_year)
  ))
}

# ============================================================================
# 3. PREPARE COMBINED DATA FRAME
# ============================================================================
# Join the spatially-resolved arthropod index onto individual bird records.
# Priority: (1) local IDW estimate (locality × year); (2) AF-region year mean.
# Records in years with no arthropod coverage at either level are dropped.
# arthropod variable is standardised at the record level after joining.

birds <- birds %>% mutate(wing_length = conc.wing.length)  # clean name for formulas

if (INCLUDE_ARTHRO) {
  birds <- birds %>%
    # Local estimate (locality × year)
    left_join(arthro_by_loc_year %>%
                select(Longitude_decimal_degrees, Latitude_decimal_degrees,
                       Year, arthro_obs_mean, arthro_source),
              by = c("Longitude_decimal_degrees", "Latitude_decimal_degrees", "Year")) %>%
    # Fallback: AF-region year mean for localities with no local estimate
    left_join(arthro_af_year %>% select(Year, arthro_af_mean = arthro_obs_mean),
              by = "Year") %>%
    mutate(
      arthro_source   = if_else(is.na(arthro_obs_mean), "af_region", arthro_source),
      arthro_obs_mean = if_else(is.na(arthro_obs_mean), arthro_af_mean, arthro_obs_mean)
    ) %>%
    select(-arthro_af_mean) %>%
    filter(!is.na(arthro_obs_mean)) %>%
    mutate(arthro_obs_std = as.numeric(scale(arthro_obs_mean)))

  message(sprintf("Arthropod source breakdown: %s",
    paste(names(table(birds$arthro_source)),
          table(birds$arthro_source), sep = "=", collapse = ", ")))
}

birds <- birds %>%
  mutate(spp        = as.factor(spp),
         scaled_yr  = as.numeric(scaled_yr),    # scale() matrices -> numeric
         scaled_lat = as.numeric(scaled_lat)) %>%  # (before filter(): dplyr >= 1.1
  filter(!is.na(wing_length),                     #  deprecates 1-col matrices there)
         !is.na(scaled_tmean),
         !is.na(scaled_lat),
         !is.na(Sex))

# Contributor intercept on the wing node (Phase 4f). Main_researcher is complete
# (0 % NA, 42 contributors among wing records; audit_sources.rds). If the
# column is somehow absent the model falls back to species-only intercepts and
# says so - the results file records which was fitted (src_intercept).
HAS_SRC <- "Main_researcher" %in% names(birds) && !all(is.na(birds$Main_researcher))
if (HAS_SRC) {
  n_na_src <- sum(is.na(birds$Main_researcher))
  if (n_na_src > 0) message("Dropping ", n_na_src, " records with NA Main_researcher.")
  birds <- birds %>% filter(!is.na(Main_researcher)) %>%
    mutate(Main_researcher = factor(Main_researcher))
} else {
  message("Main_researcher not available: fitting WITHOUT the contributor intercept.")
}

if (SMOKE) {
  set.seed(SEM_SEED)
  birds <- birds[sample(nrow(birds), min(1500L, nrow(birds))), ] %>% droplevels()
  message("SMOKE MODE: subsampled to ", nrow(birds), " records.")
}

# Cache key for node FORMULAS: bump whenever a node formula changes (the data
# signature alone cannot see a structural change). v4 = contributor intercept.
MODEL_TAG <- paste0("v4-", .mode,
                    if (HAS_SRC) "-src" else "-nosrc",
                    if (INCLUDE_ARTHRO) paste0("-spatArthro", ARTHRO_RADIUS_KM, "km") else "",
                    if (SMOKE) "-SMOKE" else "")

# Wing-node MEAN formula, shared by section 4 (i.i.d. species intercept) and
# section 10 (relmat() phylogenetic species effect). `env` must see K for relmat.
wing_mu_formula <- function(species_term, env = parent.frame()) {
  rhs <- c("scaled_tmean",
           if (INCLUDE_ARTHRO) "arthro_obs_std",
           "Sex", "scaled_lat", species_term,
           if (HAS_SRC) "(1 | Main_researcher)")
  stats::as.formula(paste("wing_length ~", paste(rhs, collapse = " + ")), env = env)
}
# drmTMB::bf() captures its arguments with substitute() and insists on literal
# formula calls, so a formula held in a variable must be spliced in with
# bquote() and evaluated in `env` (which must see K for relmat()).
wing_bf <- function(species_term, env = parent.frame()) {
  f <- wing_mu_formula(species_term, env)
  eval(bquote(drmTMB::bf(.(f), sigma ~ scaled_yr + scaled_tmean)), envir = env)
}

if (INCLUDE_ARTHRO) {
  cat(sprintf(
    "SEM dataset [arthropod model]: %d records, %d species, %d years with arthropod data (range %d–%d)\n",
    nrow(birds), length(unique(birds$spp)), length(unique(birds$Year)),
    min(birds$Year), max(birds$Year)))
} else {
  cat(sprintf(
    "SEM dataset [full-data temp/variance model]: %d records, %d species, years %d–%d\n",
    nrow(birds), length(unique(birds$spp)), min(birds$Year), max(birds$Year)))
}
if (HAS_SRC) cat(sprintf("Contributors (Main_researcher): %d\n", nlevels(birds$Main_researcher)))
cat("Wing-node mean formula:", deparse1(wing_mu_formula("(1 | spp)")), "\n")

# ============================================================================
# 4. FIT drmSEM
# ============================================================================
#
# INCLUDE_ARTHRO = TRUE → three endogenous nodes:
#   Node 1 — scaled_tmean:   temperature, driven by year and latitude
#   Node 2 — arthro_obs_std: arthropod abundance, driven by year, temperature,
#                            and latitude (RETIRED, override only; the
#                            temperature / latitude edges are explained in
#                            HISTORY in the header and in section 9)
#   Node 3 — wing_length:    the primary response, with a sigma submodel
#
# INCLUDE_ARTHRO = FALSE → two endogenous nodes (no arthropod node, arthropod
#   term dropped from wing_length): the full-data temperature/variance model.
#
# In both cases the sigma submodel sigma ~ scaled_yr + scaled_tmean tests whether
# morphological variability changes over time and with temperature (parallel to
# the lnCVR analysis). Remaining missing edges define the Fisher's C d-sep claims.

# The wing node's mean formula is built by wing_mu_formula(): species intercept
# (1 | spp), plus (1 | Main_researcher) when available, plus arthro_obs_std only
# under the retired override. Random-effect terms are not causal edges: drmSEM
# strips them from the DAG, so the d-sep claim set is unchanged by HAS_SRC.
wing_node <- drm_node(wing_bf("(1 | spp)"), family = gaussian())

if (INCLUDE_ARTHRO) {
  sem_fit <- drm_sem(
    scaled_tmean = drm_node(
      drmTMB::bf(scaled_tmean ~ scaled_yr + scaled_lat + (1 | spp)),
      family = gaussian()
    ),
    # Observed arthropod abundance ~ year + temperature + latitude. The
    # scaled_lat edge is a statistical adjustment for the year-level aggregate's
    # sampling structure (LR ~ 276 claim), NOT a causal claim - see HISTORY.
    arthro_obs_std = drm_node(
      drmTMB::bf(arthro_obs_std ~ scaled_yr + scaled_tmean + scaled_lat),
      family = gaussian()
    ),
    wing_length = wing_node,
    data = birds
  )
} else {
  sem_fit <- drm_sem(
    scaled_tmean = drm_node(
      drmTMB::bf(scaled_tmean ~ scaled_yr + scaled_lat + (1 | spp)),
      family = gaussian()
    ),
    wing_length = wing_node,
    data = birds
  )
}

# ============================================================================
# 5. DIAGNOSTICS AND OVERALL FIT
# ============================================================================

# Node-level convergence, available samplers for d-sep and simulation
check_sem(sem_fit)

# Fisher's C statistic: global test of DAG fit
# p > 0.05 → data are consistent with the DAG (no missing paths required)
fc <- fisher_c(sem_fit)
cat("\n--- Fisher's C ---\n"); print(fc)

# Individual d-separation claims (missing-edge independence tests)
dsep_tab <- dsep(sem_fit)
cat("\n--- D-separation tests ---\n"); print(dsep_tab)

# ============================================================================
# 6. PATH COEFFICIENTS
# ============================================================================

# Raw (unstandardised) paths, labelled by distributional component (mu / sigma)
cat("\n--- Path table ---\n")
path_table <- paths(sem_fit)
print(path_table)

# Standardised paths (SD of X scaling, preserves response unit for mu paths)
cat("\n--- Standardised paths (sd_x) ---\n")
path_std <- standardize(sem_fit, method = "sd_x")
print(path_std)

# Fully standardised (latent method; both X and Y rescaled)
cat("\n--- Fully standardised paths (latent) ---\n")
print(standardize(sem_fit, method = "latent"))

# ============================================================================
# 7. DIRECT, INDIRECT, AND DISTRIBUTION-MEDIATED EFFECTS
# ============================================================================
# Effects propagate through ALL distributional components of mediating nodes
# (not just means), so the `distribution_mediated` column captures the sigma-
# channel contribution. With TWO mediators (scaled_tmean, arthro_obs_std) the
# year → wing indirect effect now has two routes, so we also report each channel
# separately via `through =`.
#
# Interpretation key:
#   total_path        = direct + indirect (mean + distribution mediated)
#   mean_mediated     = indirect effect via mediator MEANS
#   distribution_mediated = indirect effect via mediator sigma/nu
#   through = "scaled_tmean"   → the temperature channel only
#   through = "arthro_obs_std" → the arthropod channel only

# Cache key: invalidate whenever the SEM input data OR model structure changes
# (MODEL_TAG captures formula changes the data summary alone would miss).
.data_sig <- function(d) paste(
  MODEL_TAG, nrow(d), length(unique(d$spp)),
  round(sum(d$wing_length), 3), round(sum(d$scaled_tmean), 3),
  if ("arthro_obs_std" %in% names(d)) round(sum(d$arthro_obs_std), 3) else "noarth",
  round(sum(as.integer(factor(d$Sex))), 0),
  if (HAS_SRC) paste0("src", nlevels(d$Main_researcher)) else "nosrc",
  sep = "-"
)
sig <- .data_sig(birds)

if (!RECOMPUTE && file.exists(CACHE_FILE) &&
    identical(tryCatch(readRDS(CACHE_FILE)$sig, error = function(e) NA), sig)) {
  message("Loading cached SEM effects from ", CACHE_FILE,
          " (set RECOMPUTE <- TRUE to refit).")
  .cache         <- readRDS(CACHE_FILE)
  effects_yr     <- .cache$effects_yr
  total_yr       <- .cache$total_yr
  effects_yr_tmp <- .cache$effects_yr_tmp
  effects_yr_art <- .cache$effects_yr_art
  effects_temp   <- .cache$effects_temp
} else {
  message("Computing SEM effects (B = ", EFFECT_B, ", nsim = ", EFFECT_NSIM,
          ") and caching to ", CACHE_FILE, " ...")
  # NOTE: indirect_effects() has no `method` arg (it always reports both the
  # mean- and distribution-mediated legs); only total_effects() takes `method`.
  effects_yr <- indirect_effects(            # year → wing, BOTH mediators
    sem_fit, from = "scaled_yr", to = "wing_length",
    B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)

  total_yr <- total_effects(
    sem_fit, from = "scaled_yr", to = "wing_length",
    method = "simulate", B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)

  # Per-channel decomposition only makes sense with >1 mediator (arthropod mode).
  if (INCLUDE_ARTHRO) {
    effects_yr_tmp <- indirect_effects(      # year → wing, temperature channel
      sem_fit, from = "scaled_yr", to = "wing_length", through = "scaled_tmean",
      B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)
    effects_yr_art <- indirect_effects(      # year → wing, arthropod channel
      sem_fit, from = "scaled_yr", to = "wing_length", through = "arthro_obs_std",
      B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)
  } else {
    effects_yr_tmp <- NULL; effects_yr_art <- NULL
  }

  effects_temp <- indirect_effects(          # temperature → wing (via arthropod if present)
    sem_fit, from = "scaled_tmean", to = "wing_length",
    B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)

  saveRDS(list(sig = sig, EFFECT_B = EFFECT_B, EFFECT_NSIM = EFFECT_NSIM,
               seed = SEM_SEED, effects_yr = effects_yr, total_yr = total_yr,
               effects_yr_tmp = effects_yr_tmp, effects_yr_art = effects_yr_art,
               effects_temp = effects_temp), CACHE_FILE)
}

cat("\n--- Effects of scaled_yr → wing_length (all mediators) ---\n"); print(effects_yr)
cat("\n--- Total effects: scaled_yr → wing_length ---\n");             print(total_yr)
if (INCLUDE_ARTHRO) {
  cat("\n--- scaled_yr → wing_length via temperature channel ---\n");  print(effects_yr_tmp)
  cat("\n--- scaled_yr → wing_length via arthropod channel ---\n");    print(effects_yr_art)
}
cat("\n--- Effects of scaled_tmean → wing_length ---\n");              print(effects_temp)

# ============================================================================
# 8. VISUALISE THE DAG
# ============================================================================
# Component-labelled arrows: solid = mu paths, styled differently for sigma.

# plot() draws by side effect, so under Rscript it would leave a stray Rplots.pdf;
# the publication DAG is built explicitly in section 12a.
if (interactive()) { dag_plot <- plot(sem_fit); print(dag_plot) }

# Optional: save for the manuscript
# ggsave(file.path("..", "Manuscript", "images", "fig-drmsem-dag.png"),
#        dag_plot, width = 7, height = 5, dpi = 300)

# ============================================================================
# 9. SENSITIVITY: ARTHROPOD AS EXOGENOUS PREDICTED INDEX (retired main model)
# ============================================================================
# The endogenous observed-arthropod node is now the MAIN model (sections 2–8).
# This section documents the RETIRED alternative — the arthropod abundance index
# predicted by the fitted brm_no_grass model and treated as EXOGENOUS —
# for sensitivity comparison.
#
# WHY IT WAS RETIRED: the predicted index is a near-linear function of year, so
# once temperature is time-resolved it shares temperature's temporal trend and
# the d-separation claim (arthro ⟂ tmean | year, lat) is strongly rejected
# (Fisher's C p ≈ 1e-26). It is unbiased as an exogenous covariate but cannot
# sit as an endogenous node (perfect collinearity with year), so it cannot carry
# the scaled_tmean → arthropod edge the data require.
#
# TRADE-OFF vs the main model: the predicted index has correct (year-level)
# replication and a smoother signal, but forces the misspecified DAG; the
# observed aggregate fixes the DAG at the cost of anticonservative SEs on the
# arthropod node. Compare directions/magnitudes across the two.
#
# Requires archive/predicts/brm_no_grass.rda (archived in Phase 4e; never loaded by default).
# Uncomment to run.

# # --- 9a. Predicted exogenous index (population-level, per year) ----------
# arthro_env <- new.env()
# load(file.path(REPO_DIR, "archive", "predicts", "brm_no_grass.rda"), envir = arthro_env)  # archived Phase 4e
# brm_arthro  <- arthro_env$brm_no_grass
# # Recover the date-scaling brm_arthro used (scale() on Sample_midpoint Dates).
# arthro_sp_num  <- as.numeric(dat_no_grassland$Sample_midpoint)
# arthro_yr_mean <- mean(arthro_sp_num, na.rm = TRUE)
# arthro_yr_sd   <- sd(arthro_sp_num,   na.rm = TRUE)
# new_dates_num  <- as.numeric(as.Date(paste0(1990:2018, "-07-01")))
# newdata_arthro <- data.frame(
#   Sample_midpoint = as.Date(paste0(1990:2018, "-07-01")),
#   scaled_yr       = (new_dates_num - arthro_yr_mean) / arthro_yr_sd,
#   log_sampling    = 0)
# epred <- brms::posterior_epred(brm_arthro, newdata = newdata_arthro, re_formula = NA)
# arthro_idx <- data.frame(
#   Year             = 1990:2018,
#   arthro_index_std = as.numeric(scale(log(apply(epred, 2, median)))))
#
# birds_s9 <- birds %>%
#   dplyr::select(-dplyr::any_of(c("arthro_obs_mean", "arthro_obs_std", "n_arthro"))) %>%
#   left_join(arthro_idx, by = "Year") %>%
#   filter(!is.na(arthro_index_std))
#
# # --- 9b. Fit with arthropod EXOGENOUS (no arthropod node) ----------------
# sem_s9 <- drm_sem(
#   scaled_tmean = drm_node(
#     drmTMB::bf(scaled_tmean ~ scaled_yr + scaled_lat + (1 | spp)),
#     family = gaussian()
#   ),
#   wing_length = drm_node(
#     drmTMB::bf(wing_length ~ scaled_tmean + arthro_index_std +
#                  Sex + scaled_lat + (1 | spp),
#                sigma ~ scaled_yr + scaled_tmean),
#     family = gaussian()
#   ),
#   data = birds_s9
# )
#
# check_sem(sem_s9)
# print(fisher_c(sem_s9))                 # expect strong rejection (see above)
# print(paths(sem_s9))
# print(standardize(sem_s9, method = "sd_x"))
# print(indirect_effects(sem_s9, from = "scaled_yr", to = "wing_length",
#                         B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED))

# ============================================================================
# 10. PHYLOGENETIC SEM (Route A: relmat() + drm_phylo_cov)
# ============================================================================
# Phylogeny-corrected refit of the main DAG. The wing_length node's species
# intercept becomes a PHYLOGENETIC random effect via relmat(1 | sp_tip, K = K),
# where K is an evolutionary correlation matrix built by drm_phylo_cov() under
# Pagel's lambda. drmSEM recognises relmat() as a structured-effect marker and
# strips it from the causal edge set, so paths(), dsep() and fisher_c() still
# describe the fixed-effect DAG — the phylogeny only enters each node's
# likelihood. Two uncertainty layers are propagated:
#   - lambda is SELECTED PER TREE by AIC over LAMBDA_GRID (focal node = wing);
#   - coefficients are POOLED ACROSS N_PHYLO_TREES trees by Rubin's rules,
#     mirroring atlantic_parallel.R::pool_rubin().
# The temperature node keeps its i.i.d. (1 | spp) intercept (it indexes the
# habitat a species was sampled in, not an evolving trait); swap in relmat()
# there too if you want it phylogenetic.
#
# This block is self-contained but reuses the SAME tree-retrieval / name-
# reconciliation stack as atlantic_parallel.R (prepR4pcm + clootl). It is gated
# on RUN_PHYLO and on those packages being installed, and caches both the
# reconciled trees and the pooled results so re-sourcing does not refit.

if (isTRUE(RUN_PHYLO) &&
    all(vapply(c("ape", "phytools", "prepR4pcm", "clootl"),
               requireNamespace, logical(1), quietly = TRUE))) {

  suppressPackageStartupMessages({
    library(ape); library(phytools); library(prepR4pcm)
  })
  set.seed(SEM_SEED)

  # --- 10a. Tree retrieval + name reconciliation (cached) -------------------
  # ABT (2018) -> current eBird/Clements genus changes, identical to
  # atlantic_parallel.R (clootl errors on names absent from the eBird taxonomy).
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
    # NB: Herpsilochmus sellowi (Caatinga Antwren) is intentionally NOT mapped —
    # it is absent from the clootl eBird/Clements 2025 taxonomy under both that
    # name and its current Radinopsyche sellowi reclassification, so it is dropped
    # below (same decision as atlantic_parallel.R).
  )

  if (!RECOMPUTE_PHYLO && file.exists(PHYLO_PREP_CACHE) &&
      identical(tryCatch(readRDS(PHYLO_PREP_CACHE)$sig, error = function(e) NA), sig)) {
    message("Loading cached phylo prep from ", PHYLO_PREP_CACHE)
    .prep        <- readRDS(PHYLO_PREP_CACHE)
    birds_phylo  <- .prep$birds_phylo
    tree_samp    <- .prep$tree_samp
    clootl_ver_p <- .prep$clootl_version
  } else {
    # Build a clean scientific-name column to match against tree tips.
    base_name <- if ("Binomial" %in% names(birds)) {
      gsub("_", " ", as.character(birds$Binomial))
    } else if ("species_name" %in% names(birds)) {
      gsub("_", " ", as.character(birds$species_name))
    } else {
      gsub("_", " ", as.character(birds$spp))
    }
    birds$sp_tip <- base_name
    hit <- birds$sp_tip %in% names(ebird_synonyms)
    birds$sp_tip[hit] <- ebird_synonyms[birds$sp_tip[hit]]
    spp_data <- unique(birds$sp_tip)

    # clootl needs the AvesData repo for the dated 100-tree sample sets; only
    # download when AVESDATA_PATH is unset/missing (the download dominates cost).
    if (!nzchar(Sys.getenv("AVESDATA_PATH")) ||
        !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
      clootl::get_avesdata_repo(path = raw_path())
    }
    # Sanity check before the expensive 100-tree pull: n_tree = 1 REPORTS the
    # names absent from the clootl eBird/Clements taxonomy instead of erroring.
    # Strip them from spp_data so the 100-tree call doesn't error hard, and let
    # reconcile_apply(drop_unresolved = TRUE) below also drop any tree-absent
    # species from the model. Herpsilochmus sellowi is the one known, expected
    # miss (see the synonym note above); stopifnot(<= 1) keeps that tolerance
    # but halts if a NEW unmatched name appears, so additions aren't dropped
    # silently — add a synonym above and re-run. Mirrors atlantic_parallel.R.
    .check <- pr_get_tree(spp_data, source = "clootl", n_tree = 1)
    if (length(.check$unmatched) > 0) {
      message("Removing ", length(.check$unmatched),
              " species absent from eBird taxonomy: ",
              paste(.check$unmatched, collapse = ", "))
      spp_data <- setdiff(spp_data, .check$unmatched)
    }
    stopifnot(length(.check$unmatched) <= 1)
    got    <- pr_get_tree(spp_data, source = "clootl", n_tree = 100, cache = TRUE)
    trees  <- got$tree
    clootl_ver_p <- pr_cite_tree(got, format = "text")

    rec <- reconcile_tree(x = birds, tree = trees[[1]], x_species = "sp_tip",
                          fuzzy = TRUE, resolve = "flag")
    print(reconcile_summary(rec))           # inspect; reconcile_override() for misses
    aligned <- reconcile_apply(rec, data = birds, tree = trees[[1]],
                               species_col = "sp_tip", drop_unresolved = TRUE)
    birds_phylo <- aligned$data
    birds_phylo$sp_tip <- factor(gsub(" ", "_", birds_phylo$sp_tip))

    # Prune all 100 trees to the resolved set (underscore labels), subsample N.
    norm_us   <- function(x) gsub(" ", "_", x)
    keep_tips <- norm_us(aligned$tree$tip.label)
    trees_pruned <- lapply(trees, function(t) {
      t$tip.label <- norm_us(t$tip.label)
      ape::keep.tip(t, intersect(keep_tips, t$tip.label))
    })
    class(trees_pruned) <- "multiPhylo"
    tree_samp <- sample(trees_pruned, min(N_PHYLO_TREES, length(trees_pruned)))

    saveRDS(list(birds_phylo = birds_phylo, tree_samp = tree_samp,
                 clootl_version = clootl_ver_p, sig = sig), PHYLO_PREP_CACHE)
  }

  cat(sprintf(
    "\nPhylo SEM dataset: %d records, %d species, %d trees\n",
    nrow(birds_phylo), nlevels(birds_phylo$sp_tip), length(tree_samp)))

  # --- 10b. Helpers ---------------------------------------------------------
  # Underscore tip labels + coerce to ultrametric (pruning leaves clootl trees
  # very slightly non-ultrametric, which the covariance build dislikes).
  prep_tree <- function(tree) {
    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    if (!ape::is.ultrametric(tree)) {
      tree <- phytools::force.ultrametric(tree, method = "nnls")
    }
    tree
  }
  # Pagel's-lambda evolutionary correlation matrix for relmat().
  make_K <- function(tree, lambda) {
    drmSEM::drm_phylo_cov(prep_tree(tree), model = "lambda", lambda = lambda)
  }
  # Restrict data to species the matrix covers (defensive: every relmat level
  # must be a row of K).
  sub_to_K <- function(K) birds_phylo[birds_phylo$sp_tip %in% rownames(K), , drop = FALSE]

  # Focal (wing_length) node, fitted alone for AIC-based lambda selection.
  # The arthropod term is included only in INCLUDE_ARTHRO mode (matches main DAG).
  fit_wing_node <- function(K) {
    bf_wing <- wing_bf("relmat(1 | sp_tip, K = K)", env = environment())
    drmTMB::drmTMB(bf_wing, data = sub_to_K(K), family = gaussian())
  }
  # Full phylogenetic SEM at a given K (same node set as the main model; the
  # wing node keeps the (1 | Main_researcher) intercept alongside relmat() -
  # drmTMB accepts both in one node, tested 2026-09-09 on a subset).
  fit_phylo_sem <- function(K) {
    d <- sub_to_K(K)
    wing_node_K <- drm_node(wing_bf("relmat(1 | sp_tip, K = K)", env = environment()),
                            family = gaussian())
    if (INCLUDE_ARTHRO) {
      drm_sem(
        scaled_tmean = drm_node(
          drmTMB::bf(scaled_tmean ~ scaled_yr + scaled_lat + (1 | spp)),
          family = gaussian()),
        arthro_obs_std = drm_node(
          drmTMB::bf(arthro_obs_std ~ scaled_yr + scaled_tmean + scaled_lat),
          family = gaussian()),
        wing_length = wing_node_K,
        data = d)
    } else {
      drm_sem(
        scaled_tmean = drm_node(
          drmTMB::bf(scaled_tmean ~ scaled_yr + scaled_lat + (1 | spp)),
          family = gaussian()),
        wing_length = wing_node_K,
        data = d)
    }
  }
  # Pick lambda minimising AIC of the focal node on a single tree.
  select_lambda <- function(tree) {
    aics <- vapply(LAMBDA_GRID, function(l) {
      out <- tryCatch(stats::AIC(fit_wing_node(make_K(tree, l))),
                      error = function(e) NA_real_)
      out
    }, numeric(1))
    if (all(is.na(aics))) return(NA_real_)
    LAMBDA_GRID[which.min(aics)]
  }
  # Rubin's-rules pooling of a list of raw path tables (estimate + std.error).
  pool_rubin_paths <- function(path_list) {
    key <- function(d) paste(d$from, d$to, d$component, d$term, sep = " | ")
    ids <- key(path_list[[1]])
    m   <- length(path_list)
    E <- sapply(path_list, function(d) d$estimate[match(ids, key(d))])
    S <- sapply(path_list, function(d) d$std.error[match(ids, key(d))])
    if (is.null(dim(E))) { E <- matrix(E, nrow = 1); S <- matrix(S, nrow = 1) }
    qbar <- rowMeans(E, na.rm = TRUE)
    ubar <- rowMeans(S^2, na.rm = TRUE)            # within-tree variance
    b    <- apply(E, 1L, stats::var, na.rm = TRUE) # between-tree variance
    se   <- sqrt(ubar + (1 + 1 / m) * b)           # total (Rubin) SE
    z    <- qbar / se
    data.frame(path = ids, estimate = qbar, se = se,
               lower = qbar - 1.96 * se, upper = qbar + 1.96 * se,
               z = z, p.value = 2 * stats::pnorm(-abs(z)),
               stringsAsFactors = FALSE)
  }

  # --- 10c. Fit across trees + pool (cached) --------------------------------
  if (!RECOMPUTE_PHYLO && file.exists(PHYLO_CACHE) &&
      identical(readRDS(PHYLO_CACHE)$n_trees, length(tree_samp)) &&
      identical(readRDS(PHYLO_CACHE)$model_tag, MODEL_TAG) &&
      identical(readRDS(PHYLO_CACHE)$sig, sig)) {
    message("Loading cached phylo results from ", PHYLO_CACHE)
    .phy        <- readRDS(PHYLO_CACHE)
    phylo_pool  <- .phy$phylo_pool
    lambda_hat  <- .phy$lambda_hat
  } else {
    message("Fitting phylogenetic SEM across ", length(tree_samp),
            " trees (AIC lambda grid: ",
            paste(LAMBDA_GRID, collapse = ", "), ") ...")
    per_tree <- lapply(seq_along(tree_samp), function(i) {
      tree <- tree_samp[[i]]
      l_star <- select_lambda(tree)
      if (is.na(l_star)) return(NULL)
      sem_i <- tryCatch(fit_phylo_sem(make_K(tree, l_star)),
                        error = function(e) NULL)
      if (is.null(sem_i)) return(NULL)
      list(lambda = l_star, paths = as.data.frame(paths(sem_i)))
    })
    per_tree <- Filter(Negate(is.null), per_tree)
    if (!length(per_tree)) {
      warning("No phylogenetic SEM fit succeeded; skipping pooled output.")
    } else {
      lambda_hat <- vapply(per_tree, `[[`, numeric(1), "lambda")
      phylo_pool <- pool_rubin_paths(lapply(per_tree, `[[`, "paths"))
      saveRDS(list(phylo_pool = phylo_pool, lambda_hat = lambda_hat,
                   n_trees = length(tree_samp), lambda_grid = LAMBDA_GRID,
                   model_tag = MODEL_TAG, sig = sig),
              PHYLO_CACHE)
    }
  }

  # --- 10d. Report ----------------------------------------------------------
  if (exists("phylo_pool")) {
    cat("\n--- Pagel's lambda selected per tree (AIC) ---\n")
    print(summary(lambda_hat)); cat("table:\n"); print(table(lambda_hat))
    cat("\n--- Phylogeny-corrected path table (Rubin-pooled over trees) ---\n")
    print(phylo_pool, row.names = FALSE, digits = 4)
    cat("\nCompare these directions/magnitudes with the i.i.d. paths() in",
        "section 6 and with brm0_multiphylo.rda (the formal authority).\n")
  }

} else if (isTRUE(RUN_PHYLO)) {
  message("Section 10 skipped: install ape, phytools, prepR4pcm, clootl to run ",
          "the phylogenetic SEM (same stack as atlantic_parallel.R).")
}

# ============================================================================
# 11. EXPORT RESULTS  (persist findings to the repo — runs every time)
# ============================================================================
# Writes a machine-readable .rds and a human-readable .md so Fisher's C, the
# d-separation claims, the path tables, the effect decomposition, and (if
# fitted) the pooled phylogenetic paths are recorded in the codebase rather than
# living only on the console. Ungated by the effect cache, so it always reflects
# the current run. Files are mode-specific (arthro vs noarthro).
RESULTS_RDS <- out_path(paste0("drmsem_results_", .mode, ".rds"))
RESULTS_MD  <- out_path(paste0("drmsem_results_", .mode, ".md"))

results <- list(
  model_tag  = MODEL_TAG,
  mode       = .mode,
  generated  = as.character(Sys.time()),
  smoke      = SMOKE,                       # TRUE = subset smoke test, NOT a result
  status     = "supplementary, exploratory (REVISION_PLAN S2.7 / Phase 4f)",
  n_records  = nrow(birds),
  n_species  = length(unique(birds$spp)),
  n_contributors = if (HAS_SRC) nlevels(birds$Main_researcher) else NA_integer_,
  src_intercept  = HAS_SRC,                 # (1 | Main_researcher) on the wing node
  temperature_source = TEMP_SOURCE,
  wing_mu_formula = deparse1(wing_mu_formula("(1 | spp)")),
  year_range = range(birds$Year),
  dsep_claims_to_report = if (INCLUDE_ARTHRO) c(
    "Sex _||_ scaled_tmean | {scaled_yr, scaled_lat}",
    "Sex _||_ arthro_obs_std | {scaled_yr, scaled_tmean, scaled_lat}") else
    "Sex _||_ scaled_tmean | {scaled_yr, scaled_lat}   (Fisher's C on df = 2)",
  software = list(
    R = R.version.string,
    drmTMB    = as.character(utils::packageVersion("drmTMB")),
    drmSEM    = as.character(utils::packageVersion("drmSEM")),
    prepR4pcm = if (requireNamespace("prepR4pcm", quietly = TRUE))
                  as.character(utils::packageVersion("prepR4pcm")) else NA_character_,
    disclosure = paste("drmSEM/drmTMB and prepR4pcm are co-authored by E.S.A.S.;",
                       "development versions (r-universe / GitHub), to be disclosed",
                       "as unpublished software in the Supplement (REVISION_PLAN S2.7).")),
  fisher_c   = if (exists("fc"))         as.data.frame(fc)         else NULL,
  dsep       = if (exists("dsep_tab"))   as.data.frame(dsep_tab)   else NULL,
  paths_raw  = if (exists("path_table")) as.data.frame(path_table) else NULL,
  paths_sdx  = if (exists("path_std"))   as.data.frame(path_std)   else NULL,
  effects    = list(
    yr            = if (exists("effects_yr"))     effects_yr     else NULL,
    yr_via_tmean  = if (exists("effects_yr_tmp")) effects_yr_tmp else NULL,
    yr_via_arthro = if (exists("effects_yr_art")) effects_yr_art else NULL,
    tmean         = if (exists("effects_temp"))   effects_temp   else NULL,
    total_yr      = if (exists("total_yr"))       total_yr       else NULL
  ),
  phylo = if (exists("phylo_pool")) list(
    pooled_paths = phylo_pool,
    lambda_hat   = if (exists("lambda_hat")) lambda_hat else NULL,
    n_trees      = if (exists("tree_samp")) length(tree_samp) else NA_integer_
  ) else NULL
)
saveRDS(results, RESULTS_RDS)

# Human-readable .md (code-fenced table dumps via capture.output).
.fence <- function(title, obj) {
  if (is.null(obj)) character(0)
  else c("", paste0("## ", title), "```", capture.output(print(obj)), "```")
}

.lines <- c(
  paste0("# drmSEM results — ", MODEL_TAG, if (SMOKE) "  [SMOKE TEST — NOT A RESULT]" else ""),
  paste0("_generated ", results$generated, " · N = ", results$n_records,
         " records, ", results$n_species, " species, years ",
         results$year_range[1], "–", results$year_range[2],
         if (HAS_SRC) paste0(", ", results$n_contributors, " contributors") else "",
         "; temperature: ", TEMP_SOURCE, "_"),
  "",
  "_Supplementary, exploratory analysis (REVISION_PLAN §2.7 / Phase 4f). ",
  "Wing-node mean formula: `", results$wing_mu_formula, "`. ",
  "d-separation claim(s) to report: ", paste(results$dsep_claims_to_report, collapse = "; "), ". ",
  results$software$disclosure, "_",
  .fence("Fisher's C (global DAG fit)",            results$fisher_c),
  .fence("d-separation claims",                    results$dsep),
  .fence("Path coefficients (raw)",                results$paths_raw),
  .fence("Path coefficients (standardised, sd_x)", results$paths_sdx)
)
if (!is.null(results$phylo)) .lines <- c(.lines,
  "", paste0("## Phylogeny-corrected paths (Rubin-pooled across ",
             results$phylo$n_trees, " trees)"),
  paste0("Pagel's lambda (AIC-selected per tree), range: ",
         paste(range(results$phylo$lambda_hat), collapse = "–")),
  "```", capture.output(print(results$phylo$pooled_paths)), "```")
.lines <- c(.lines,
  .fence("Effect: scaled_yr -> wing_length (all mediators)", results$effects$yr),
  .fence("Effect: scaled_yr -> wing_length via temperature", results$effects$yr_via_tmean),
  .fence("Effect: scaled_yr -> wing_length via arthropod",   results$effects$yr_via_arthro),
  .fence("Effect: scaled_tmean -> wing_length",              results$effects$tmean))
writeLines(.lines, RESULTS_MD)
message("Wrote results to ", RESULTS_RDS, " and ", RESULTS_MD)

# ============================================================================
# 12. FIGURES  (export result visualisations — runs every time)
# ============================================================================
# Four core figures, saved as mode-tagged PNGs in Analysis/figures/. Each is
# built from objects already in the session and wrapped in tryCatch so a single
# plotting failure never aborts the run.
FIG_DIR <- fig_path()
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
fpath     <- function(name) file.path(FIG_DIR, paste0(name, "_", .mode, ".png"))
.try_save <- function(name, plot_obj, w = 7, h = 5) {
  tryCatch({
    ggplot2::ggsave(fpath(name), plot_obj, width = w, height = h, dpi = 300)
    message("  saved ", fpath(name))
  }, error = function(e) message("  figure '", name, "' failed: ", conditionMessage(e)))
}

# --- 12a. Annotated DAG (self-contained ggplot path diagram) ---------------
# Built from standardize(sem_fit, "sd_x") so it always renders (the package's
# plot(sem_fit) draws by side-effect / needs a graph-layout backend and ggsave'd
# blank). The wing-length response is drawn as TWO co-equal response nodes —
# mean(wing) and sd(wing) — reflecting that drmSEM models the mean and the
# residual SD as separate distributional responses, each with its own submodel.
# Mean-model paths point into mean(wing); the variance (sigma) paths point into
# sd(wing). Both receive ordinary directed arrows coloured by standardized
# coefficient (solid = p<0.05, dashed = n.s.), so the residual-SD channel reads
# as a response on equal footing with the mean rather than a dotted annotation.
.try_save("drmsem_dag", local({
  ps  <- as.data.frame(standardize(sem_fit, method = "sd_x"))
  pos <- data.frame(
    node  = c("scaled_yr","scaled_lat","Sex","scaled_tmean","arthro_obs_std","mean_wing","sd_wing"),
    x     = c(0, 0, 0, 1.2, 1.2, 2.6, 2.6),
    y     = c(3.2, 2.0, 0.6, 3.2, 1.3, 2.7, 0.9),
    role  = c("exogenous","exogenous","exogenous","mediator","mediator","response","response"),
    label = c("scaled_yr","scaled_lat","Sex","scaled_tmean","arthro_obs_std",
              "mean(wing)","sd(wing)"),
    stringsAsFactors = FALSE)
  xy <- function(n, col) pos[[col]][match(n, pos$node)]
  # Every path is a directed edge into a node. Retarget the two distributional
  # components of wing_length onto the two response nodes.
  e <- ps
  e$fromn <- ifelse(e$term == "SexMale", "Sex", e$from)
  e$ton   <- e$to
  e$ton[e$to == "wing_length" & e$component == "mu"]    <- "mean_wing"
  e$ton[e$to == "wing_length" & e$component == "sigma"] <- "sd_wing"
  # Draw only nodes that actually appear in this model (so the no-arthropod
  # variant doesn't show a phantom arthro_obs_std node).
  e <- e[e$fromn %in% pos$node & e$ton %in% pos$node, ]
  e$x0 <- xy(e$fromn,"x"); e$y0 <- xy(e$fromn,"y")
  e$x1 <- xy(e$ton,"x");   e$y1 <- xy(e$ton,"y")
  e$linef <- ifelse(e$p.value < 0.05, "significant", "n.s.")
  # Stagger label position along the edge by target node so the crossing
  # Sex -> mean(wing) and year -> sd(wing) labels do not collide mid-plot.
  e$lf <- ifelse(e$ton == "sd_wing", 0.74, 0.56)
  present  <- unique(c(e$fromn, e$ton))
  pos_draw <- pos[pos$node %in% present, ]
  # Colour-scale limits exclude the Sex outlier (std coef ~2.5) so the gradient
  # stays informative for the climate/arthropod paths; Sex saturates.
  .clim <- max(0.6, stats::quantile(abs(e$std.estimate[e$fromn != "Sex"]), 0.95, na.rm = TRUE))
  ggplot() +
    geom_segment(data = e,
      aes(x = x0, y = y0, xend = x1, yend = y1, colour = std.estimate, linetype = linef),
      arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"), linewidth = 0.7) +
    geom_text(data = e,
      aes(x = x0 + lf*(x1-x0), y = y0 + lf*(y1-y0),
          label = sprintf("%.2f", std.estimate), colour = std.estimate),
      size = 2.7, fontface = "bold") +
    geom_label(data = pos_draw, aes(x = x, y = y, label = label, fill = role),
      colour = "black", size = 3.1) +
    scale_colour_gradient2(low = "#b2182b", mid = "grey75", high = "#1b7837",
      midpoint = 0, limits = c(-.clim, .clim), oob = scales::squish,
      name = "std. coef") +
    scale_fill_manual(values = c(exogenous = "#ECECEC", mediator = "#CDE7DD",
      response = "#FCE3C8"), name = NULL) +
    scale_linetype_manual(values = c(significant = "solid", `n.s.` = "dashed"), name = NULL) +
    coord_cartesian(xlim = c(-0.3, 3.1), ylim = c(0.2, 3.7)) +
    theme_void() + theme(legend.position = "bottom") +
    labs(title = paste0("drmSEM path diagram [", .mode,
                        "] — standardized coefficients; wing mean and SD are separate responses"))
}), w = 9, h = 6.5)

# --- 12b. Coefficient forest (raw paths ± 95% CI, faceted mu vs sigma) ------
# Raw estimates keep units honest (mu ≈ mm per SD; sigma on the log scale), so
# the two components are faceted with free axes rather than forced onto one scale.
.try_save("drmsem_coef_forest", local({
  pd <- as.data.frame(path_table)
  pd$label <- paste(pd$from, "→", pd$to, "(", pd$term, ")")
  pd$lo  <- pd$estimate - 1.96 * pd$std.error
  pd$hi  <- pd$estimate + 1.96 * pd$std.error
  pd$sig <- ifelse(pd$p.value < 0.05, "p < 0.05", "n.s.")
  ggplot(pd, aes(estimate, reorder(label, estimate), colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi)) +
    facet_wrap(~ component, scales = "free", ncol = 1) +
    scale_colour_manual(values = c("p < 0.05" = "#1b7837", "n.s." = "grey60")) +
    labs(x = "coefficient (mu: mm per SD of predictor; sigma: log scale)",
         y = NULL, colour = NULL,
         title = paste0("drmSEM path coefficients [", .mode, "]")) +
    theme_minimal()
}), h = 6)

# --- 12c. Effect-decomposition forest (year and temperature → wing) --------
if (exists("effects_yr") && exists("effects_temp")) .try_save("drmsem_effects_forest", local({
  mk <- function(e, lbl) { d <- as.data.frame(e); d$panel <- lbl; d }
  ed <- rbind(mk(effects_yr, "scaled_yr → wing"),
              mk(effects_temp, "scaled_tmean → wing"))
  ed$quantity <- factor(ed$quantity,
    levels = c("total_path","direct","indirect","mean_mediated","distribution_mediated"))
  ggplot(ed, aes(estimate, quantity)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
    facet_wrap(~ panel, scales = "free_x") +
    labs(x = "effect on wing length (mm)", y = NULL,
         title = paste0("Effect decomposition [", .mode, "]")) +
    theme_minimal()
}), w = 9, h = 4)

# --- 12d. Distributional (sigma) channel: fold-change in residual SD --------
# exp(slope * x) over +/- 2 SD of each sigma predictor, with a pointwise band
# from the slope SE. Shows the novel variance result (e.g. warmer -> less
# variable wings; the year trend whose sign is under scrutiny).
.try_save("drmsem_sigma_curves", local({
  sr <- as.data.frame(path_table)
  sr <- sr[sr$component == "sigma", , drop = FALSE]
  lab <- c(scaled_yr = "Year (SD)", scaled_tmean = "Temperature (SD)")
  xx <- seq(-2, 2, length.out = 60)
  cur <- do.call(rbind, lapply(seq_len(nrow(sr)), function(i) {
    b <- sr$estimate[i]; se <- sr$std.error[i]
    f_lo <- exp((b - 1.96 * se) * xx); f_hi <- exp((b + 1.96 * se) * xx)
    data.frame(
      predictor = ifelse(sr$term[i] %in% names(lab), lab[sr$term[i]], sr$term[i]),
      x = xx, fold = exp(b * xx),
      lo = pmin(f_lo, f_hi), hi = pmax(f_lo, f_hi))
  }))
  ggplot(cur, aes(x, fold)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_line(linewidth = 1, colour = "#762a83") +
    facet_wrap(~ predictor, scales = "free_x") +
    labs(x = "predictor (SD units)",
         y = "fold-change in residual SD of wing length",
         title = paste0("Distributional (sigma) channel [", .mode, "]")) +
    theme_minimal()
}), w = 8, h = 4)

