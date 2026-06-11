# atlantic_drmsem.R
# ---------------------------------------------------------------------------
# Distributional piecewise SEM (drmSEM) for Atlantic Forest passerines.
#
# RESEARCH QUESTION: Is the bird body size decline over time (1990–2018)
# mediated by rising temperatures (Bergmann's rule) and/or declining
# arthropod abundance (food limitation), and do these drivers also affect
# the VARIANCE of body size — not just the mean?
#
# CAUSAL GRAPH (DAG):
#
#   scaled_yr ──→ scaled_tmean ─────────────────→ wing_length
#        │             │                              ↑
#        │             └──────→ arthro_obs_std ───────┘
#        └────────────────────────────┘
#
#   sigma(wing_length) ~ scaled_yr + scaled_tmean
#   [Tests whether morphological variability changes over time/temperature,
#    complementing the lnCVR meta-analysis in the primary manuscript.]
#
# KEY DESIGN DECISIONS:
#   (1) Arthropod node is ENDOGENOUS (arthro_obs_std ~ scaled_yr + scaled_tmean).
#       It uses the OBSERVED annual mean log(abundance / effort) from the raw
#       PREDICTS data (dat_no_grassland), NOT model predictions. EARLIER this
#       node was exogenous and built from the fitted brm_no_ants_no_grass model;
#       that index is a near-linear function of year, so with TIME-RESOLVED
#       temperature it shared temperature's temporal trend and the d-separation
#       claim (arthro ⟂ tmean | year, lat) was strongly rejected (Fisher's C
#       p ≈ 1e-26). Making arthropods endogenous and adding the scaled_tmean →
#       arthropod edge models that dependence instead of assuming it away, so
#       the DAG is no longer misspecified on that edge. The retired exogenous
#       predicted-index variant is kept (commented) in section 9 as a sensitivity.
#       CAVEAT: the observed aggregate is a YEAR-LEVEL value broadcast to every
#       bird in that year, so the arthropod node is fitted on replicated values
#       and its SEs are anticonservative; per-bird scaled_tmean as a predictor of
#       a year-level outcome is a level mismatch. Treat the arthropod node's
#       coefficients/SEs as approximate (the wing_length and temperature nodes
#       are unaffected). Years with < 5 arthropod records are dropped, so the
#       analysis is restricted to years with arthropod coverage.
#
#   (2) The MAIN SEM (sections 4–8) uses i.i.d. species intercepts (1 | spp).
#       drmTMB/drmSEM DO support phylogenetic effects (this corrects an earlier
#       note): section 10 fits a phylogeny-corrected variant via the structured-
#       effect marker relmat(1 | species_name, K), where K is an evolutionary
#       covariance built by drmSEM::drm_phylo_cov(). drmSEM strips relmat() from
#       the causal edge set, so paths(), dsep(), fisher_c() and the effect
#       calculus keep operating on the fixed-effect DAG while each node's
#       likelihood carries the phylogenetic random effect. Pagel's lambda is
#       selected per tree by AIC and coefficients are pooled across 50 trees
#       (Rubin's rules). The primary brms analysis (brm0_multiphylo.rda) remains
#       the formal authority; treat the drmSEM phylo paths as corroboration and
#       compare directions/magnitudes.
#
#   (3) Single complete-case dataset (from passer90_climate.rds). Multiple
#       imputation is not propagated here (50-tree phylogenetic uncertainty IS,
#       in section 10). Run the SEM on each imputed dataset and pool coefficients
#       manually (Rubin's rules) if MI uncertainty propagation is also needed.
#
# REQUIREMENTS:
#   pak::pak(c("drmTMB", "drmSEM"),
#            repos = "https://itchyshin.r-universe.dev")
#   passer90_climate.rds — run climate_extraction.R (WorldClim) first
#   dat_no_grassland.rds — from atlantic_birds_ms.Rmd (raw arthropod records;
#     the main model no longer needs brm_no_ants_no_grass.rda — only the
#     commented exogenous-index sensitivity in section 9 does)
#   Section 10 (phylogenetic variant) additionally needs:
#     pak::pak("itchyshin/prepR4pcm"); install.packages(c("clootl","phytools","ape"))
#     — same tree-retrieval / name-reconciliation stack as atlantic_parallel.R.
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

.find_file <- function(fname, subdir = "Analysis") {
  dirs <- getwd(); d <- getwd()
  for (i in 1:8) { d <- dirname(d); dirs <- c(dirs, d) }
  cand <- unique(c(file.path(dirs, fname), file.path(dirs, subdir, fname)))
  hit  <- cand[file.exists(cand)]
  if (!length(hit)) stop("Could not locate '", fname,
      "'. Run from inside the atlantic_birds repo.")
  normalizePath(hit[1], winslash = "/")
}
ANALYSIS_DIR <- dirname(.find_file("passer90.rda"))
apath <- function(...) file.path(ANALYSIS_DIR, ...)

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
SEM_SEED    <- 20250611
EFFECT_B    <- 200L     # coefficient draws  (CI precision)        — 50 for dev
EFFECT_NSIM <- 200L     # inner sims, distribution-mediated leg    — 50 for dev
RECOMPUTE   <- FALSE    # TRUE = ignore the cache and recompute
CACHE_FILE  <- apath("drmsem_effects_cache.rds")

# --- Section 10 (phylogenetic SEM, Route A) ---------------------------------
RUN_PHYLO        <- TRUE                  # FALSE skips the whole phylo variant
N_PHYLO_TREES    <- 50L                   # trees sampled for Rubin pooling
LAMBDA_GRID      <- c(0, 0.25, 0.5, 0.75, 0.9, 1)  # Pagel's lambda grid (AIC-selected)
PHYLO_PREP_CACHE <- apath("drmsem_phylo_prep.rds")   # reconciled trees + matched data
PHYLO_CACHE      <- apath("drmsem_phylo_results.rds") # pooled paths + per-tree lambda
RECOMPUTE_PHYLO  <- FALSE                 # TRUE = ignore phylo caches and refit

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Bird data: prefer time-resolved WorldClim temperature (passer90_climate.rds);
# fall back to static Annual_mean_temperature from passer90.rda if not yet generated.
climate_path <- apath("passer90_climate.rds")
if (file.exists(climate_path)) {
  birds <- readRDS(climate_path)
  message("Using time-resolved temperature (scaled_tmean from WorldClim).")
} else {
  warning(
    "passer90_climate.rds not found — falling back to Annual_mean_temperature.\n",
    "This is a STATIC per-locality climatology: it captures spatial but NOT\n",
    "temporal temperature variation. The year → temperature → body size\n",
    "mediation test will reflect spatial Bergmann's rule, not the temporal trend.\n",
    "Run climate_extraction.R (WorldClim) first to generate passer90_climate.rds.",
    call. = FALSE
  )
  load(apath("passer90.rda"))   # → passer90
  birds <- passer90
  birds$scaled_tmean <- as.numeric(scale(birds$Annual_mean_temperature))
}

# Raw arthropod records (PREDICTS, no ants, forest biomes only). The MAIN model
# aggregates these directly (section 2); the fitted brm_no_ants_no_grass model is
# only needed for the commented exogenous-index sensitivity in section 9.
load(apath("dat_no_grassland.rds"))   # → dat_no_grassland

# ============================================================================
# 2. BUILD OBSERVED YEAR-LEVEL ARTHROPOD AGGREGATE (1990–2018)
# ============================================================================
# Observed annual mean of log(abundance / sampling effort) from the raw PREDICTS
# records — the ENDOGENOUS arthropod node's data (see design decision (1)).
# Sample_midpoint is a Date; extract the calendar YEAR with format(), NOT
# as.integer() (which would return days since 1970-01-01).
arthro_by_year <- dat_no_grassland %>%
  mutate(
    Year          = as.integer(format(Sample_midpoint, "%Y")),
    # small pseudocount before logging to handle zero counts
    log_abund_eff = log((Measurement / Sampling_effort) + 0.01)
  ) %>%
  group_by(Year) %>%
  summarise(
    arthro_obs_mean = mean(log_abund_eff, na.rm = TRUE),
    n_arthro        = dplyr::n(),
    .groups         = "drop"
  ) %>%
  filter(n_arthro >= 5,                    # require ≥5 arthropod records/yr
         Year >= 1990, Year <= 2018)

# ============================================================================
# 3. PREPARE COMBINED DATA FRAME
# ============================================================================
# Join the observed year-level arthropod aggregate onto individual bird records,
# keep complete cases, then standardise the arthropod variable at the RECORD
# level (so its SD reflects the analysed sample, matching scaled_tmean etc.).
# Records in years without arthropod coverage (n_arthro < 5) are dropped.

birds <- birds %>%
  left_join(arthro_by_year, by = "Year") %>%
  mutate(wing_length = conc.wing.length) %>%      # clean name for formula parsing
  filter(!is.na(wing_length),
         !is.na(scaled_tmean),
         !is.na(arthro_obs_mean),
         !is.na(scaled_lat),
         !is.na(Sex)) %>%
  mutate(spp            = as.factor(spp),
         arthro_obs_std = as.numeric(scale(arthro_obs_mean)))

cat(sprintf(
  "SEM dataset: %d individual records, %d species, %d years with arthropod data (range %d–%d)\n",
  nrow(birds), length(unique(birds$spp)), length(unique(birds$Year)),
  min(birds$Year), max(birds$Year)
))

# ============================================================================
# 4. FIT drmSEM
# ============================================================================
#
# Three endogenous nodes (scaled_yr and scaled_lat are exogenous):
#   Node 1 — scaled_tmean:   temperature, driven by year and latitude
#   Node 2 — arthro_obs_std: arthropod abundance, driven by year AND temperature
#                            (the edge that resolves the temp⟂arthro misspec.)
#   Node 3 — wing_length:    the primary response, with a sigma submodel
#
# Remaining missing edges define the d-separation claims Fisher's C tests, e.g.
# scaled_yr ⟂ wing_length | {scaled_tmean, arthro_obs_std, scaled_lat, Sex}:
# once both mediators are conditioned on, year should have no residual DIRECT
# effect on the wing-length mean.

sem_fit <- drm_sem(

  # Node 1: temperature driven by temporal trend and latitude
  scaled_tmean = drm_node(
    drmTMB::bf(scaled_tmean ~ scaled_yr + scaled_lat + (1 | spp)),
    family = gaussian()
  ),

  # Node 2: observed arthropod abundance, driven by year and temperature.
  # (Year-level outcome broadcast to bird-records — SEs anticonservative; see
  #  design decision (1). No species term: it is a year-level covariate.)
  arthro_obs_std = drm_node(
    drmTMB::bf(arthro_obs_std ~ scaled_yr + scaled_tmean),
    family = gaussian()
  ),

  # Node 3: wing length with a sigma submodel
  # sigma ~ scaled_yr + scaled_tmean tests whether morphological variability
  # changes over time and with temperature (parallel to the lnCVR analysis).
  wing_length = drm_node(
    drmTMB::bf(wing_length ~ scaled_tmean + arthro_obs_std +
                 Sex + scaled_lat + (1 | spp),
               sigma ~ scaled_yr + scaled_tmean),
    family = gaussian()
  ),

  data = birds
)

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
cat("\n--- D-separation tests ---\n"); print(dsep(sem_fit))

# ============================================================================
# 6. PATH COEFFICIENTS
# ============================================================================

# Raw (unstandardised) paths, labelled by distributional component (mu / sigma)
cat("\n--- Path table ---\n")
path_table <- paths(sem_fit)
print(path_table)

# Standardised paths (SD of X scaling, preserves response unit for mu paths)
cat("\n--- Standardised paths (sd_x) ---\n")
print(standardize(sem_fit, method = "sd_x"))

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

# Cache key: invalidate whenever the SEM input data changes.
.data_sig <- function(d) paste(
  nrow(d), length(unique(d$spp)),
  round(sum(d$wing_length),    3), round(sum(d$scaled_tmean),   3),
  round(sum(d$arthro_obs_std), 3), round(sum(as.integer(factor(d$Sex))), 0),
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

  effects_yr_tmp <- indirect_effects(        # year → wing, temperature channel
    sem_fit, from = "scaled_yr", to = "wing_length", through = "scaled_tmean",
    B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)

  effects_yr_art <- indirect_effects(        # year → wing, arthropod channel
    sem_fit, from = "scaled_yr", to = "wing_length", through = "arthro_obs_std",
    B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)

  effects_temp <- indirect_effects(          # temperature → wing (now also via arthropod)
    sem_fit, from = "scaled_tmean", to = "wing_length",
    B = EFFECT_B, nsim = EFFECT_NSIM, seed = SEM_SEED)

  saveRDS(list(sig = sig, EFFECT_B = EFFECT_B, EFFECT_NSIM = EFFECT_NSIM,
               seed = SEM_SEED, effects_yr = effects_yr, total_yr = total_yr,
               effects_yr_tmp = effects_yr_tmp, effects_yr_art = effects_yr_art,
               effects_temp = effects_temp), CACHE_FILE)
}

cat("\n--- Effects of scaled_yr → wing_length (all mediators) ---\n"); print(effects_yr)
cat("\n--- Total effects: scaled_yr → wing_length ---\n");             print(total_yr)
cat("\n--- scaled_yr → wing_length via temperature channel ---\n");    print(effects_yr_tmp)
cat("\n--- scaled_yr → wing_length via arthropod channel ---\n");      print(effects_yr_art)
cat("\n--- Effects of scaled_tmean → wing_length ---\n");              print(effects_temp)

# ============================================================================
# 8. VISUALISE THE DAG
# ============================================================================
# Component-labelled arrows: solid = mu paths, styled differently for sigma.

dag_plot <- plot(sem_fit)
print(dag_plot)

# Optional: save for the manuscript
# ggsave(file.path("..", "Manuscript", "images", "fig-drmsem-dag.png"),
#        dag_plot, width = 7, height = 5, dpi = 300)

# ============================================================================
# 9. SENSITIVITY: ARTHROPOD AS EXOGENOUS PREDICTED INDEX (retired main model)
# ============================================================================
# The endogenous observed-arthropod node is now the MAIN model (sections 2–8).
# This section documents the RETIRED alternative — the arthropod abundance index
# predicted by the fitted brm_no_ants_no_grass model and treated as EXOGENOUS —
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
# Requires brm_no_ants_no_grass.rda (not loaded by default — see section 1).
# Uncomment to run.

# # --- 9a. Predicted exogenous index (population-level, per year) ----------
# arthro_env <- new.env()
# load(apath("brm_no_ants_no_grass.rda"), envir = arthro_env)
# brm_arthro  <- arthro_env$brm_no_ants_no_grass
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
    "Tiaris fuliginosus"       = "Asemospiza fuliginosa"
  )

  if (!RECOMPUTE_PHYLO && file.exists(PHYLO_PREP_CACHE)) {
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
      clootl::get_avesdata_repo(path = ANALYSIS_DIR)
    }
    stopifnot(length(pr_get_tree(spp_data, source = "clootl",
                                 n_tree = 1)$unmatched) == 0)
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
                 clootl_version = clootl_ver_p), PHYLO_PREP_CACHE)
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
  fit_wing_node <- function(K) {
    drmTMB::drmTMB(
      drmTMB::bf(wing_length ~ scaled_tmean + arthro_obs_std + Sex +
                   scaled_lat + relmat(1 | sp_tip, K = K),
                 sigma ~ scaled_yr + scaled_tmean),
      data = sub_to_K(K), family = gaussian())
  }
  # Full phylogenetic SEM at a given K (same 3-node DAG as the main model).
  fit_phylo_sem <- function(K) {
    d <- sub_to_K(K)
    drm_sem(
      scaled_tmean = drm_node(
        drmTMB::bf(scaled_tmean ~ scaled_yr + scaled_lat + (1 | spp)),
        family = gaussian()),
      arthro_obs_std = drm_node(
        drmTMB::bf(arthro_obs_std ~ scaled_yr + scaled_tmean),
        family = gaussian()),
      wing_length = drm_node(
        drmTMB::bf(wing_length ~ scaled_tmean + arthro_obs_std + Sex +
                     scaled_lat + relmat(1 | sp_tip, K = K),
                   sigma ~ scaled_yr + scaled_tmean),
        family = gaussian()),
      data = d)
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
      identical(readRDS(PHYLO_CACHE)$n_trees, length(tree_samp))) {
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
                   n_trees = length(tree_samp), lambda_grid = LAMBDA_GRID),
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

