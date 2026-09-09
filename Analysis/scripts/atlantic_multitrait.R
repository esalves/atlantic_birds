# atlantic_multitrait.R
# ---------------------------------------------------------------------------
# WHAT: Fast (lme4, REML) multi-trait screen of temporal trends in the SIX
# morphological traits carried by ATLANTIC BIRD TRAITS for the live-only,
# known-sex analytical sample (73 species): log body mass, bill width, bill
# length, tail length, tarsus length, and — as a cross-reference only — the
# coalesced wing length (the authoritative controlled wing analysis is Phase 1,
# atlantic_parallel_controlled.R). For every trait the same cumulative sequence
# of provenance / site / geography / season / individual controls is fitted,
# followed by a Mundlak within/between decomposition of the year effect by
# contributor and by municipality. For body mass a capture-hour (diurnal
# fattening) model is fitted on the records with a parseable clock time.
#
# WHY (REVISION_PLAN.md §0 "Multi-trait scrutiny", §1 "§2.5 Mass and
# allometry", "Bill width and other traits", Phase 2 items 2-4): the referee's
# provenance critique applies to every trait, not just wing length. The plan's
# §0 table was produced from session scratch scripts (§6 of the plan, "planning
# aids, not results"); this script is the reproducible version of that table,
# with sample sizes and complete-case handling made explicit, so the manuscript
# can read every number from Analysis/output/multitrait_results.rds.
#
# MODELS (per trait; y = trait; REML; Wald 95 % CI = estimate ± 1.96 SE):
#   M0   y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)      [= brms
#        fixed/random structure without the phylogenetic term]
#   M1a  M0 + (1 | src)                    src  = Main_researcher (contributor)
#   M1b  M0 + (1 | site)                   site = Municipality
#   M1   M0 + (1 | src) + (1 | site)
#   M1cc M1 refitted on the M2 complete-case sample (separates covariate
#        adjustment from the change of sample when season / altitude are NA)
#   M2   M1 + scaled_lon + scaled_alt + season
#   M3   M2 + (1 | ring)                   ring = Ring x Binomial; unringed
#        records get a unique level each (they carry no information about the
#        individual variance but are not dropped)
#   MWsrc  Mundlak, M3 set: scaled_yr replaced by yr_within_src + yr_mean_src
#          (within-contributor and between-contributor year slopes)
#   MWsite Mundlak, M3 set: scaled_yr replaced by yr_within_site + yr_mean_site
#   HOUR (mass only) ln_body_mass ~ Sex + scaled_yr + scaled_lat + hour_num +
#        season + (1 + scaled_yr || spp) + (1 | src) + (1 | site), on records
#        with 5 <= hour_num <= 19
#
# MODES (2026-09-09):
#   (no flag)     lme4 fast tier as above (the Tier-1 screen). Writes
#                 output/multitrait_results.rds and output/multitrait_table.md
#                 (the markdown gains a phylogenetic-tier section whenever
#                 output/multitrait_phylo_results.rds exists on disk).
#   --trees N     PHYLOGENETIC TIER: the same models refitted with glmmTMB's
#                 `propto` covariance structure (Williams, McGillycuddy, Drobniak,
#                 Bolker, Warton & Nakagawa 2025, bioRxiv 10.64898/2025.12.20.695312)
#                 through the shared engine scripts/_phylo_engine.R, on the
#                 IDENTICAL N (default 50) trees of the published brms analysis,
#                 Rubin-pooled across trees. Per trait (log body mass, bill width,
#                 bill length, tail, tarsus; wing is left to Phase 1): M0, M1,
#                 M3 and the within-contributor Mundlak decomposition (MWsrc);
#                 plus the body-mass capture-hour pair (HOUR0 / HOUR1). Reads the
#                 lme4 tier from output/multitrait_results.rds (run the default
#                 mode first) for the comparison table and writes
#                 output/multitrait_phylo_results.rds, then regenerates
#                 output/multitrait_table.md with both tiers. Minutes.
#   --engine glmmTMB|brms   engine of the phylogenetic tier. Only glmmTMB is
#                 implemented in this script (there never was a brms path here;
#                 the brms body-mass model is atlantic_parallel_mass.R and the
#                 bivariate brms cross-check is atlantic_bivariate_wing_mass.R
#                 --engine brms). --engine brms stops with that message.
#
# UNITS: scaled_yr is (Year - 2009.487) / 5.0199 (attr "scaled:scale" of
# passer90$scaled_yr, i.e. the SD of the pre-threshold known-sex sample; see
# REVISION_NOTES_P0.md). Per-decade change = beta * 10 / 5.0199. Percent per
# decade = 100 * per_decade / mean(trait) for mm traits and
# 100 * (exp(per_decade) - 1) for log mass.
#
# INPUTS:  data/derived/passer90.rda            (known sex, 73 spp, 42 columns)
#          output/effect_scale.rds              (optional; P0's published-model
#                                               conversions, read for the
#                                               isometry comparison)
# OUTPUTS: output/multitrait_results.rds        list: $table (tidy, one row per
#          trait x model x term), $coverage (contributors / municipalities /
#          spanning contributors per trait), $hour (mass diurnal model),
#          $isometry, $plan_comparison, $meta
#          output/multitrait_phylo_results.rds  (--trees) list: $table (pooled,
#          one row per trait x model x term, with n, trees converged, per-decade
#          units), $varcomp (phylogenetic SD / proportion per model), $hour,
#          $comparison_lme4 (phylogenetic vs lme4 estimate, SE ratio, whether
#          the zero-exclusion verdict changes), $timing, $meta
#          output/multitrait_table.md           markdown tables pasted into
#          REVISION_NOTES_P2.md (lme4 tier + phylogenetic tier)
#
# RUN:  Rscript Analysis/scripts/atlantic_multitrait.R              # lme4 tier, ~30 s
#       Rscript Analysis/scripts/atlantic_multitrait.R --trees 50   # phylogenetic tier, minutes
#       (from the repo root, Analysis/, or Analysis/scripts/)
# Session used for the numbers in REVISION_NOTES_P2.md: R 4.6.0, lme4 2.0.1,
# Matrix, dplyr 1.1.x. Nothing stochastic; set.seed recorded by convention.
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(lme4) })
set.seed(20260909)

# --- command line: default = lme4 tier; --trees N = phylogenetic tier (glmmTMB) ---
args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag) {
  i <- grep(paste0("^", flag, "(=|$)"), args)[1]
  if (is.na(i)) return(NULL)
  v <- sub(paste0("^", flag, "=?"), "", args[i]); if (!nzchar(v)) v <- args[i + 1]
  v
}
MODE    <- if (!is.null(get_flag("--trees"))) "trees" else "lme4"
N_TREES <- if (MODE == "trees") as.integer(get_flag("--trees")) else 0L
if (MODE == "trees" && (is.na(N_TREES) || N_TREES < 1L)) stop("--trees needs a positive integer")
ENGINE  <- if (!is.null(get_flag("--engine"))) get_flag("--engine") else "glmmTMB"
if (MODE == "trees" && tolower(ENGINE) != "glmmtmb")
  stop("Only --engine glmmTMB is implemented in atlantic_multitrait.R (no brms path exists here; ",
       "see atlantic_parallel_mass.R for the brms body-mass model and atlantic_bivariate_wing_mass.R --engine brms).")
ENGINE <- if (MODE == "trees") "glmmTMB" else "lme4"

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
# Shared phylogenetic engine (glmmTMB propto; uses the path helpers above):
# phylo_species_name(), phylo_A_list(), fit_phylo_glmmtmb(), run_phylo_trees(), pool_rubin_df(), tidy_phylo_fit()
if (MODE == "trees") source(file.path(ANALYSIS_DIR, "scripts", "_phylo_engine.R"))

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
load(derived_path("passer90.rda"))
stopifnot(all(passer90$Status == "live"), all(passer90$known_sex))

SD_YR    <- as.numeric(attr(passer90$scaled_yr, "scaled:scale"))   # 5.0199
CTR_YR   <- as.numeric(attr(passer90$scaled_yr, "scaled:center"))  # 2009.487
DEC      <- 10 / SD_YR                                             # SD-years per decade
EARLY_MAX <- 2006; LATE_MIN <- 2013                                # quartile periods (P0)

d <- passer90 %>%
  mutate(
    scaled_yr  = as.numeric(scaled_yr),      # 1-column matrices -> numeric
    scaled_lat = as.numeric(scaled_lat),
    Sex        = factor(Sex, levels = c("Female", "Male")),
    src        = factor(Main_researcher),
    site       = factor(Municipality),
    # individual = Ring x Binomial (36 ring strings recur across species, P0);
    # unringed records get a unique level so they are kept but uninformative
    ring       = factor(ifelse(is.na(Ring), paste0("unringed_", ID_ABT),
                               paste(Ring, Binomial, sep = "|"))),
    season     = factor(season, levels = c("DJF", "MAM", "JJA", "SON")),
    # wing kept for cross-reference only (Phase 1 is authoritative)
    wing       = conc.wing.length
  )
stopifnot(!anyNA(d$src), !anyNA(d$scaled_lat))
n_site_na <- sum(is.na(d$site))

traits <- c(
  ln_body_mass        = "log body mass (ln g)",
  `Bill_width.mm.`    = "bill width (mm)",
  `Bill_length.mm.`   = "bill length (mm)",
  `Tail_length.mm.`   = "tail length (mm)",
  `Tarsus_length.mm.` = "tarsus length (mm)",
  wing                = "wing length, coalesced (mm) [cross-reference; Phase 1 is authoritative]"
)
is_log <- c(ln_body_mass = TRUE, `Bill_width.mm.` = FALSE, `Bill_length.mm.` = FALSE,
            `Tail_length.mm.` = FALSE, `Tarsus_length.mm.` = FALSE, wing = FALSE)

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------
# Fit with lmer(REML), capture convergence / singularity, return a tidy row set
fit_row <- function(formula, data, trait, model, terms_keep, note = "") {
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  msgs <- character(0)
  fit <- withCallingHandlers(
    tryCatch(lmer(formula, data = data, REML = TRUE, control = ctrl),
             error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL }),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  if (is.null(fit)) {
    return(data.frame(trait = trait, model = model, term = terms_keep, estimate = NA_real_,
                      se = NA_real_, t = NA_real_, lower = NA_real_, upper = NA_real_,
                      n = nrow(data), n_spp = NA_integer_, n_src = NA_integer_, n_site = NA_integer_,
                      singular = NA, messages = paste(msgs, collapse = " | "),
                      formula = deparse1(formula), stringsAsFactors = FALSE))
  }
  cf <- summary(fit)$coefficients
  ng <- ngrps(fit)
  rows <- lapply(terms_keep, function(tm) {
    if (!tm %in% rownames(cf)) return(NULL)
    est <- cf[tm, "Estimate"]; se <- cf[tm, "Std. Error"]
    data.frame(trait = trait, model = model, term = tm, estimate = est, se = se,
               t = cf[tm, "t value"], lower = est - 1.96 * se, upper = est + 1.96 * se,
               n = nobs(fit), n_spp = as.integer(ng[["spp"]]),
               n_src = as.integer(if ("src" %in% names(ng)) ng[["src"]] else NA),
               n_site = as.integer(if ("site" %in% names(ng)) ng[["site"]] else NA),
               singular = isSingular(fit),
               messages = paste(unique(msgs), collapse = " | "),
               formula = deparse1(formula), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$note <- note
  out
}

# Convert a beta on scaled_yr to per-decade and percent-per-decade units
add_units <- function(tab, trait_mean, log_scale) {
  yr_terms <- c("scaled_yr", "yr_within_src", "yr_mean_src", "yr_within_site", "yr_mean_site")
  is_yr <- tab$term %in% yr_terms
  tab$per_decade <- ifelse(is_yr, tab$estimate * DEC, NA_real_)
  tab$per_decade_lower <- ifelse(is_yr, tab$lower * DEC, NA_real_)
  tab$per_decade_upper <- ifelse(is_yr, tab$upper * DEC, NA_real_)
  pct <- function(x) if (log_scale) 100 * (exp(x) - 1) else 100 * x / trait_mean
  tab$pct_per_decade       <- ifelse(is_yr, pct(tab$per_decade), NA_real_)
  tab$pct_per_decade_lower <- ifelse(is_yr, pct(tab$per_decade_lower), NA_real_)
  tab$pct_per_decade_upper <- ifelse(is_yr, pct(tab$per_decade_upper), NA_real_)
  tab$trait_mean <- trait_mean
  tab
}

# Row selector on a tidy table (lme4 tier by default; the phylogenetic tier passes its own)
get <- function(tr, model, term = "scaled_yr", tab = table_all) tab[tab$trait == tr & tab$model == model & tab$term == term, ]

if (MODE == "lme4") {   # ======================= lme4 fast tier (sections 3-7) =======================
# ---------------------------------------------------------------------------
# 3. Per-trait model sequence
# ---------------------------------------------------------------------------
results  <- list()
coverage <- list()

for (tr in names(traits)) {
  dt <- d[!is.na(d[[tr]]), ]
  dt$y <- dt[[tr]]
  trait_mean <- if (is_log[[tr]]) mean(exp(dt$y)) else mean(dt$y)   # raw-scale mean for % units
  # --- coverage / provenance of this trait --------------------------------
  per_src <- dt %>% group_by(src) %>%
    summarise(n = n(), early = any(Year <= EARLY_MAX), late = any(Year >= LATE_MIN),
              yr_min = min(Year), yr_max = max(Year), .groups = "drop")
  spanning <- per_src %>% filter(early & late)
  coverage[[tr]] <- data.frame(
    trait = tr, label = traits[[tr]], n = nrow(dt), pct_of_sample = 100 * nrow(dt) / nrow(d),
    n_spp = n_distinct(dt$spp), n_src = n_distinct(dt$src), n_site = n_distinct(dt$site),
    n_src_both_periods = nrow(spanning), n_rec_src_both_periods = sum(spanning$n),
    src_both_periods = paste(spanning$src, collapse = "; "),
    n_early = sum(dt$Year <= EARLY_MAX), n_late = sum(dt$Year >= LATE_MIN),
    n_site_na = sum(is.na(dt$site)), n_season_na = sum(is.na(dt$season)),
    n_alt_na = sum(is.na(dt$scaled_alt)), n_ring_na = sum(is.na(dt$Ring)),
    trait_mean_raw = trait_mean, stringsAsFactors = FALSE)

  # --- model sequence ------------------------------------------------------
  base_terms <- c("scaled_yr", "scaled_lat", "SexMale")
  f0  <- y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)
  f1a <- update(f0, . ~ . + (1 | src))
  f1b <- update(f0, . ~ . + (1 | site))
  f1  <- update(f0, . ~ . + (1 | src) + (1 | site))
  f2  <- update(f1, . ~ . + scaled_lon + scaled_alt + season)
  f3  <- update(f2, . ~ . + (1 | ring))
  dcc <- dt %>% filter(!is.na(site), !is.na(scaled_lon), !is.na(scaled_alt), !is.na(season))
  # Mundlak decompositions on the M3 (complete-case) sample
  dcc <- dcc %>% group_by(src)  %>% mutate(yr_mean_src  = mean(scaled_yr), yr_within_src  = scaled_yr - yr_mean_src) %>%
                 group_by(site) %>% mutate(yr_mean_site = mean(scaled_yr), yr_within_site = scaled_yr - yr_mean_site) %>% ungroup()
  fWsrc  <- y ~ Sex + yr_within_src  + yr_mean_src  + scaled_lat + scaled_lon + scaled_alt + season +
                (1 + yr_within_src  || spp) + (1 | src) + (1 | site) + (1 | ring)
  fWsite <- y ~ Sex + yr_within_site + yr_mean_site + scaled_lat + scaled_lon + scaled_alt + season +
                (1 + yr_within_site || spp) + (1 | src) + (1 | site) + (1 | ring)

  cat(sprintf("[%s] n = %d, complete-case (M2/M3) n = %d ...\n", tr, nrow(dt), nrow(dcc)))
  res <- rbind(
    fit_row(f0,  dt,                           tr, "M0_baseline",       base_terms),
    fit_row(f1a, dt[!is.na(dt$site), ],        tr, "M1a_src",           base_terms),
    fit_row(f1b, dt[!is.na(dt$site), ],        tr, "M1b_site",          base_terms),
    fit_row(f1,  dt[!is.na(dt$site), ],        tr, "M1_src_site",       base_terms),
    fit_row(f1,  dcc,                          tr, "M1cc_src_site_on_M2_sample", base_terms,
            note = "M1 structure on the M2/M3 complete-case sample"),
    fit_row(f2,  dcc,                          tr, "M2_geo_season",     c(base_terms, "scaled_lon", "scaled_alt", "seasonMAM", "seasonJJA", "seasonSON")),
    fit_row(f3,  dcc,                          tr, "M3_ring",           c(base_terms, "scaled_lon", "scaled_alt", "seasonMAM", "seasonJJA", "seasonSON")),
    fit_row(fWsrc,  dcc, tr, "MWsrc_mundlak_contributor",  c("yr_within_src", "yr_mean_src", "scaled_lat", "SexMale")),
    fit_row(fWsite, dcc, tr, "MWsite_mundlak_municipality", c("yr_within_site", "yr_mean_site", "scaled_lat", "SexMale"))
  )
  results[[tr]] <- add_units(res, trait_mean, is_log[[tr]])
}
table_all <- do.call(rbind, results); rownames(table_all) <- NULL
coverage  <- do.call(rbind, coverage); rownames(coverage) <- NULL

# ---------------------------------------------------------------------------
# 4. Body mass: diurnal (capture-hour) model
# ---------------------------------------------------------------------------
dh <- d %>% filter(!is.na(ln_body_mass), !is.na(hour_num), hour_num >= 5, hour_num <= 19,
                   !is.na(site), !is.na(season))
dh$y <- dh$ln_body_mass
dh_all_hours <- d %>% filter(!is.na(ln_body_mass), !is.na(hour_num))
fH0 <- y ~ Sex + scaled_yr + scaled_lat + season + (1 + scaled_yr || spp) + (1 | src) + (1 | site)
fH1 <- update(fH0, . ~ . + hour_num)
hour_tab <- rbind(
  fit_row(fH0, dh, "ln_body_mass", "HOUR0_no_hour",  c("scaled_yr", "scaled_lat", "SexMale"),
          note = "records with parseable Hour in [05:00, 19:00], same sample as HOUR1"),
  fit_row(fH1, dh, "ln_body_mass", "HOUR1_with_hour", c("scaled_yr", "scaled_lat", "SexMale", "hour_num"),
          note = "hour_num in decimal hours; 100*estimate = % mass per hour")
)
hour_tab <- add_units(hour_tab, mean(exp(dh$y)), TRUE)
hr <- hour_tab[hour_tab$model == "HOUR1_with_hour" & hour_tab$term == "hour_num", ]
hour <- list(
  table = hour_tab,
  n_mass_with_hour_any = nrow(dh_all_hours),
  n_mass_hour_05_19_complete = nrow(dh),
  hour_range_used = c(5, 19),
  pct_per_hour = 100 * hr$estimate, pct_per_hour_lower = 100 * hr$lower, pct_per_hour_upper = 100 * hr$upper,
  t_hour = hr$t,
  year_beta_with_hour = hour_tab$estimate[hour_tab$model == "HOUR1_with_hour" & hour_tab$term == "scaled_yr"],
  year_t_with_hour    = hour_tab$t[hour_tab$model == "HOUR1_with_hour" & hour_tab$term == "scaled_yr"],
  hour_distribution   = quantile(dh$hour_num, c(0, .1, .25, .5, .75, .9, 1)),
  note = "Diurnal fattening check (REVISION_PLAN.md §2.5). hour_num parsed from 'H:MM'; 'morning'/'afternoon' are NA."
)

# ---------------------------------------------------------------------------
# 5. Isometry: expected mass change for the observed wing change (mass ∝ L^3)
# ---------------------------------------------------------------------------
# For a fractional wing change w, isometry predicts a fractional mass change of
# (1 + w)^3 - 1. Computed here for the lme4 M0 and M1 estimates of this script;
# the published brms (50-tree) version is in P0's effect_scale.rds.
iso_row <- function(model, years = 23) {
  w <- get("wing", model); m <- get("ln_body_mass", model)
  w_frac <- w$estimate * (years / SD_YR) / w$trait_mean          # fractional wing change over `years`
  w_lo   <- w$lower * (years / SD_YR) / w$trait_mean
  w_hi   <- w$upper * (years / SD_YR) / w$trait_mean
  m_frac <- exp(m$estimate * (years / SD_YR)) - 1
  m_lo   <- exp(m$lower * (years / SD_YR)) - 1
  m_hi   <- exp(m$upper * (years / SD_YR)) - 1
  data.frame(model = model, years = years,
             wing_pct = 100 * w_frac, wing_pct_lower = 100 * w_lo, wing_pct_upper = 100 * w_hi,
             isometric_mass_pct_expected = 100 * ((1 + w_frac)^3 - 1),
             isometric_mass_pct_expected_lower = 100 * ((1 + w_lo)^3 - 1),
             isometric_mass_pct_expected_upper = 100 * ((1 + w_hi)^3 - 1),
             observed_mass_pct = 100 * m_frac, observed_mass_pct_lower = 100 * m_lo, observed_mass_pct_upper = 100 * m_hi,
             n_wing = w$n, n_mass = m$n, stringsAsFactors = FALSE)
}
isometry <- list(
  lme4 = rbind(iso_row("M0_baseline"), iso_row("M1_src_site"), iso_row("M3_ring")),
  note = paste("Isometric expectation: mass ∝ wing^3, so expected fractional mass change = (1+w)^3-1.",
               "Wing and mass rows come from DIFFERENT samples here (all wing records vs all mass records);",
               "the shared-record contrast is in bivariate_results.rds (atlantic_bivariate_wing_mass.R).",
               "Over-record change uses 23 years (1995 -> 2018), as effect_scale.rds does.")
)
if (file.exists(out_path("effect_scale.rds"))) {
  es <- readRDS(out_path("effect_scale.rds"))
  isometry$published_brms_from_effect_scale <- es[intersect(names(es), c("isometry", "wing", "mass", "using_scaling_sd"))]
}

# ---------------------------------------------------------------------------
# 6. Comparison with the REVISION_PLAN.md §0 multi-trait table
# ---------------------------------------------------------------------------
plan <- data.frame(
  trait = c("wing", "ln_body_mass", "Bill_width.mm.", "Bill_length.mm.", "Tail_length.mm.", "Tarsus_length.mm."),
  plan_n = c(8478, 11256, 3209, 7697, 8872, 4361),
  plan_M0_beta = c(-0.92, -0.0070, 0.23, 0.02, 0.29, 0.42),
  plan_M0_t    = c(-3.8, -1.5, 2.5, 0.4, 1.0, 5.0),
  plan_M1_beta = c(-0.53, -0.0057, -0.04, 0.16, 0.31, 0.24),
  plan_M1_t    = c(-2.6, -1.0, -0.5, 1.7, 1.2, 1.5),
  plan_full_beta = c(-0.33, -0.0038, -0.05, 0.16, 0.28, 0.22),
  plan_full_t    = c(-1.9, -0.8, -0.6, 1.7, 1.1, 1.4),
  plan_within_src_beta = c(0.02, -0.0035, -0.05, 0.09, 0.12, 0.15),
  plan_within_src_t    = c(0.08, -0.7, -0.7, 0.8, 0.5, 0.9),
  stringsAsFactors = FALSE)
pick_val <- function(tr, model, term, col) { r <- get(tr, model, term); if (nrow(r)) r[[col]] else NA_real_ }
cmp <- lapply(seq_len(nrow(plan)), function(i) {
  tr <- plan$trait[i]
  data.frame(
    here_n_M0    = pick_val(tr, "M0_baseline", "scaled_yr", "n"),
    here_M0_beta = pick_val(tr, "M0_baseline", "scaled_yr", "estimate"),
    here_M0_t    = pick_val(tr, "M0_baseline", "scaled_yr", "t"),
    here_M1_beta = pick_val(tr, "M1_src_site", "scaled_yr", "estimate"),
    here_M1_t    = pick_val(tr, "M1_src_site", "scaled_yr", "t"),
    here_n_M3    = pick_val(tr, "M3_ring", "scaled_yr", "n"),
    here_M3_beta = pick_val(tr, "M3_ring", "scaled_yr", "estimate"),
    here_M3_t    = pick_val(tr, "M3_ring", "scaled_yr", "t"),
    here_within_src_beta  = pick_val(tr, "MWsrc_mundlak_contributor", "yr_within_src", "estimate"),
    here_within_src_t     = pick_val(tr, "MWsrc_mundlak_contributor", "yr_within_src", "t"),
    here_between_src_beta = pick_val(tr, "MWsrc_mundlak_contributor", "yr_mean_src", "estimate"),
    here_between_src_t    = pick_val(tr, "MWsrc_mundlak_contributor", "yr_mean_src", "t"),
    here_within_site_beta = pick_val(tr, "MWsite_mundlak_municipality", "yr_within_site", "estimate"),
    here_within_site_t    = pick_val(tr, "MWsite_mundlak_municipality", "yr_within_site", "t"))
})
plan_comparison <- cbind(plan, do.call(rbind, cmp))
# flag rows where the plan's coefficient differs from this run by more than 10 % or the t-value by > 0.3
rel <- function(a, b) ifelse(is.na(a) | is.na(b), NA, abs(a - b) / pmax(abs(a), abs(b), 1e-8))
plan_comparison$flag_M0   <- rel(plan_comparison$plan_M0_beta, plan_comparison$here_M0_beta) > 0.10 | abs(plan_comparison$plan_M0_t - plan_comparison$here_M0_t) > 0.3
plan_comparison$flag_M1   <- rel(plan_comparison$plan_M1_beta, plan_comparison$here_M1_beta) > 0.10 | abs(plan_comparison$plan_M1_t - plan_comparison$here_M1_t) > 0.3
plan_comparison$flag_full <- rel(plan_comparison$plan_full_beta, plan_comparison$here_M3_beta) > 0.10 | abs(plan_comparison$plan_full_t - plan_comparison$here_M3_t) > 0.3
plan_comparison$flag_within_src <- rel(plan_comparison$plan_within_src_beta, plan_comparison$here_within_src_beta) > 0.10 | abs(plan_comparison$plan_within_src_t - plan_comparison$here_within_src_t) > 0.3

# ---------------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------------
meta <- list(
  generated = as.character(Sys.time()),
  r_version = R.version.string,
  lme4_version = as.character(packageVersion("lme4")),
  sd_year = SD_YR, centre_year = CTR_YR, sd_years_per_decade = DEC,
  early_max = EARLY_MAX, late_min = LATE_MIN,
  n_sample = nrow(d), n_species = n_distinct(d$spp), n_site_na_in_sample = n_site_na,
  estimation = "lme4::lmer, REML = TRUE, bobyqa; Wald 95% CI = estimate +/- 1.96 SE; no phylogenetic term (fast screen)",
  ring_definition = "interaction(Ring, Binomial); unringed records given a unique level each",
  complete_case = "M2/M3/Mundlak models drop records with NA Municipality, longitude, Altitude or season; M1cc reports M1 on that same sample",
  wing_row = "wing is included only as a cross-reference to Phase 1 (atlantic_parallel_controlled.R); do not quote the wing row from this file in the manuscript"
)
out <- list(table = table_all, coverage = coverage, hour = hour, isometry = isometry,
            plan_comparison = plan_comparison, meta = meta)
saveRDS(out, out_path("multitrait_results.rds"))
phylo <- if (file.exists(out_path("multitrait_phylo_results.rds"))) readRDS(out_path("multitrait_phylo_results.rds")) else NULL

} else {   # ======================= phylogenetic tier (section 10) =======================
# ---------------------------------------------------------------------------
# 10. Phylogenetic tier: glmmTMB propto across N_TREES trees (engine: _phylo_engine.R)
# ---------------------------------------------------------------------------
# The lme4 tier is read from disk (comparison + coverage); nothing of it is refitted here.
if (!file.exists(out_path("multitrait_results.rds")))
  stop("output/multitrait_results.rds not found: run `Rscript atlantic_multitrait.R` (lme4 tier) first.")
lme4_res  <- readRDS(out_path("multitrait_results.rds"))
table_all <- lme4_res$table; coverage <- lme4_res$coverage; hour <- lme4_res$hour
t_phylo <- Sys.time()

# Species -> tree tips (eBird synonyms); the species term (spp) and the phylogenetic
# term (species_name) share the tip name, as in the published brms model.
# Herpsilochmus_sellowi is absent from the tree (engine note; 0 records in this sample).
d$species_name <- phylo_species_name(d$Binomial)
d <- d[d$species_name != "Herpsilochmus_sellowi", , drop = FALSE]
d$spp <- d$species_name
A_all <- phylo_A_list(sort(unique(d$species_name)), n_trees = N_TREES)   # the published 50 trees (cached)
N_TREES <- length(A_all)
A_for <- function(data) { s <- sort(unique(as.character(data$species_name))); lapply(A_all, function(A) A[s, s, drop = FALSE]) }
cat(sprintf("Phylogenetic tier: glmmTMB %s propto, %d trees, %d species in the sample\n",
            as.character(packageVersion("glmmTMB")), N_TREES, n_distinct(d$species_name)))

# Convergence repair around the engine loop (identical to atlantic_bivariate_wing_mass.R;
# candidate addition to _phylo_engine.R). With both (1 | spp) and the phylogenetic
# term, the REML surface has two modes; glmmTMB's default all-zero start reaches
# the worse one (phylogenetic SD ~ 0, iid species SD absorbing it, objective
# higher by ~24 units, non-positive-definite Hessian) on some trees. Trees without
# a positive-definite Hessian are refitted from the free variance parameters of a
# converged tree of the SAME model (or a generic seed), inserted into that tree's
# own theta vector (propto entries encode A and are tree-specific); the refit is
# kept only if it converged and its objective is not worse. Both counts are kept.
# A usable tree = positive-definite Hessian AND finite fixed-effect SEs. Trees still
# unusable after two seeded restarts (seed from a converged tree, then a generic
# seed) are EXCLUDED from the Rubin pooling — their SEs are undefined and pooling
# them returns NA — and listed in $excluded_trees; n_converged counts the trees
# actually pooled. Point estimates of excluded trees stay in $per_tree_fixed.
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

# One model across trees -> tidy pooled rows (same columns as the lme4 table where they exist)
pfit <- function(formula, data, trait, model, terms_keep, note = "") {
  t1 <- Sys.time()
  r <- run_phylo_robust(formula, data, A_for(data))
  secs <- as.numeric(Sys.time() - t1, units = "secs")
  po <- r$pooled %>% filter(component == "cond", par %in% terms_keep)
  # lme4-"singular" analogue: the species year-slope SD at the zero boundary (mean over
  # trees < 1e-3 x sigma). The non-phylogenetic species INTERCEPT SD is ~0 by design
  # whenever the phylogenetic term is present (absorbed) and is not flagged.
  vs <- r$varcomp_summary
  slope_col <- grep("^spp\\.(?!.*Intercept).*_mean$", names(vs), perl = TRUE, value = TRUE)[1]
  sd_slope <- if (is.na(slope_col)) NA_real_ else vs[[slope_col]]
  tab <- data.frame(trait = trait, model = model, term = po$par, estimate = po$estimate, se = po$se, t = po$z,
                    lower = po$lower, upper = po$upper, n = nrow(data), n_spp = n_distinct(data$species_name),
                    n_src = if ("src" %in% all.vars(formula)) nlevels(droplevels(data$src)) else NA_integer_,
                    n_site = if ("site" %in% all.vars(formula)) nlevels(droplevels(data$site)) else NA_integer_,
                    n_trees = r$n_trees, n_converged = r$n_converged, n_converged_first_pass = r$n_converged_first_pass,
                    n_repaired = length(r$repaired_trees), between_tree_var = po$between_tree_var,
                    sd_phylo_mean = vs$sd_phylo_mean, phylo_prop_mean = vs$phylo_prop_mean, sd_spp_slope_mean = sd_slope,
                    boundary_spp_slope = isTRUE(sd_slope < 1e-3 * vs$sigma_mean),
                    singular = NA, messages = "", formula = deparse1(formula), note = note, stringsAsFactors = FALSE)
  vc <- cbind(trait = trait, model = model, n = nrow(data), n_trees = r$n_trees, n_converged = r$n_converged,
              as.data.frame(r$varcomp_summary), stringsAsFactors = FALSE)
  yr <- po[po$par %in% c("scaled_yr", "yr_within_src"), ][1, ]
  cat(sprintf("  [%s | %s] %s = %s [%s, %s] z = %.2f | phylo SD %.4g, phylo prop %.3f | %d/%d pdHess (first pass %d) | %.0f s\n",
              trait, model, yr$par, signif(yr$estimate, 4), signif(yr$lower, 4), signif(yr$upper, 4), yr$z,
              r$varcomp_summary$sd_phylo_mean, r$varcomp_summary$phylo_prop_mean, r$n_converged, r$n_trees, r$n_converged_first_pass, secs))
  list(table = tab, varcomp = vc, per_tree = list(fixed = r$per_tree_fixed, varcomp = r$varcomp),
       timing = data.frame(trait = trait, model = model, n = nrow(data), n_trees = r$n_trees, n_converged = r$n_converged,
                           n_converged_first_pass = r$n_converged_first_pass, n_retried = length(r$retried_trees),
                           n_repaired = length(r$repaired_trees), secs = secs, secs_per_tree = secs / r$n_trees, stringsAsFactors = FALSE))
}

phylo_traits <- setdiff(names(traits), "wing")     # wing: Phase 1 (atlantic_parallel_controlled.R) is authoritative
p_tab <- list(); p_vc <- list(); p_tree <- list(); p_time <- list()
base_terms <- c("scaled_yr", "scaled_lat", "SexMale")
for (tr in phylo_traits) {
  dt <- d[!is.na(d[[tr]]), ]; dt$y <- dt[[tr]]
  trait_mean <- if (is_log[[tr]]) mean(exp(dt$y)) else mean(dt$y)
  dcc <- dt %>% filter(!is.na(site), !is.na(scaled_lon), !is.na(scaled_alt), !is.na(season)) %>%
    group_by(src) %>% mutate(yr_mean_src = mean(scaled_yr), yr_within_src = scaled_yr - yr_mean_src) %>% ungroup() %>% as.data.frame()
  f0 <- y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)
  f1 <- update(f0, . ~ . + (1 | src) + (1 | site))
  f3 <- update(f1, . ~ . + scaled_lon + scaled_alt + season + (1 | ring))
  fWsrc <- y ~ Sex + yr_within_src + yr_mean_src + scaled_lat + scaled_lon + scaled_alt + season +
    (1 + yr_within_src || spp) + (1 | src) + (1 | site) + (1 | ring)
  cat(sprintf("[%s] n = %d, complete-case (M3) n = %d, %d species ...\n", tr, nrow(dt), nrow(dcc), n_distinct(dt$species_name)))
  fits <- list(
    pfit(f0, dt,                    tr, "M0_baseline", base_terms),
    pfit(f1, dt[!is.na(dt$site), ], tr, "M1_src_site", base_terms),
    pfit(f3, dcc,                   tr, "M3_ring",     c(base_terms, "scaled_lon", "scaled_alt", "seasonMAM", "seasonJJA", "seasonSON")),
    pfit(fWsrc, dcc, tr, "MWsrc_mundlak_contributor", c("yr_within_src", "yr_mean_src", "scaled_lat", "SexMale")))
  p_tab[[tr]]  <- add_units(bind_rows(lapply(fits, `[[`, "table")), trait_mean, is_log[[tr]])
  p_vc[[tr]]   <- bind_rows(lapply(fits, `[[`, "varcomp"))
  p_time[[tr]] <- bind_rows(lapply(fits, `[[`, "timing"))
  p_tree[[tr]] <- setNames(lapply(fits, `[[`, "per_tree"), sapply(fits, function(x) x$table$model[1]))
}
p_table <- as.data.frame(bind_rows(p_tab)); p_varcomp <- as.data.frame(bind_rows(p_vc)); p_timing <- as.data.frame(bind_rows(p_time))

# Body mass, capture hour (same sample rule as the lme4 tier)
dh <- d %>% filter(!is.na(ln_body_mass), !is.na(hour_num), hour_num >= 5, hour_num <= 19, !is.na(site), !is.na(season)) %>% as.data.frame()
dh$y <- dh$ln_body_mass
fH0 <- y ~ Sex + scaled_yr + scaled_lat + season + (1 + scaled_yr || spp) + (1 | src) + (1 | site)
fH1 <- update(fH0, . ~ . + hour_num)
cat(sprintf("[ln_body_mass HOUR] n = %d records with 05:00-19:00 clock time ...\n", nrow(dh)))
hf <- list(pfit(fH0, dh, "ln_body_mass", "HOUR0_no_hour",   c("scaled_yr", "scaled_lat", "SexMale"), note = "records with parseable Hour in [05:00, 19:00]"),
           pfit(fH1, dh, "ln_body_mass", "HOUR1_with_hour", c("scaled_yr", "scaled_lat", "SexMale", "hour_num"), note = "hour_num in decimal hours; 100*estimate = % mass per hour"))
hour_tab_p <- add_units(bind_rows(lapply(hf, `[[`, "table")), mean(exp(dh$y)), TRUE)
hr <- hour_tab_p[hour_tab_p$model == "HOUR1_with_hour" & hour_tab_p$term == "hour_num", ]
hour_p <- list(table = hour_tab_p, varcomp = bind_rows(lapply(hf, `[[`, "varcomp")), n_mass_hour_05_19_complete = nrow(dh),
               pct_per_hour = 100 * hr$estimate, pct_per_hour_lower = 100 * hr$lower, pct_per_hour_upper = 100 * hr$upper, z_hour = hr$t,
               year_beta_with_hour = hour_tab_p$estimate[hour_tab_p$model == "HOUR1_with_hour" & hour_tab_p$term == "scaled_yr"],
               year_z_with_hour    = hour_tab_p$t[hour_tab_p$model == "HOUR1_with_hour" & hour_tab_p$term == "scaled_yr"],
               year_beta_without_hour = hour_tab_p$estimate[hour_tab_p$model == "HOUR0_no_hour" & hour_tab_p$term == "scaled_yr"])
p_timing <- bind_rows(p_timing, bind_rows(lapply(hf, `[[`, "timing")))
p_tree$ln_body_mass_hour <- setNames(lapply(hf, `[[`, "per_tree"), c("HOUR0_no_hour", "HOUR1_with_hour"))

# Comparison with the lme4 tier (same trait, model, term, sample)
yr_terms <- c("scaled_yr", "yr_within_src", "yr_mean_src", "hour_num")
l <- bind_rows(table_all, hour$table) %>% filter(term %in% yr_terms) %>%
  transmute(trait, model, term, lme4_estimate = estimate, lme4_se = se, lme4_lower = lower, lme4_upper = upper, lme4_n = n, lme4_singular = singular)
p <- bind_rows(p_table, hour_tab_p) %>% filter(term %in% yr_terms) %>%
  transmute(trait, model, term, phylo_estimate = estimate, phylo_se = se, phylo_lower = lower, phylo_upper = upper, phylo_n = n,
            n_trees, n_converged, phylo_pct_per_decade = pct_per_decade, phylo_pct_lower = pct_per_decade_lower, phylo_pct_upper = pct_per_decade_upper)
comparison <- inner_join(l, p, by = c("trait", "model", "term")) %>%
  mutate(same_sample = lme4_n == phylo_n, delta_estimate = phylo_estimate - lme4_estimate,
         delta_in_lme4_se = delta_estimate / lme4_se, se_ratio_phylo_over_lme4 = phylo_se / lme4_se,
         lme4_excludes_zero = lme4_lower > 0 | lme4_upper < 0, phylo_excludes_zero = phylo_lower > 0 | phylo_upper < 0,
         verdict_changes = lme4_excludes_zero != phylo_excludes_zero) %>% as.data.frame()

phylo <- list(
  table = p_table, varcomp = p_varcomp, hour = hour_p, comparison_lme4 = comparison, timing = p_timing, per_tree = p_tree,
  meta = list(generated = as.character(Sys.time()), engine = "glmmTMB propto (Williams et al. 2025) via _phylo_engine.R; REML; Rubin's rules across trees",
              n_trees = N_TREES, trees = "the 50 clootl trees of the published brms analysis (phylo_A_list(); cache data/derived/phylo_A_50trees.rds from models/brm0_multiphylo.rda)",
              glmmTMB_version = as.character(packageVersion("glmmTMB")), r_version = R.version.string,
              sd_year = SD_YR, centre_year = CTR_YR, sd_years_per_decade = DEC, n_sample = nrow(d), n_species = n_distinct(d$species_name),
              models = c(M0_baseline = "y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) + phylo",
                         M1_src_site = "M0 + (1 | src) + (1 | site)",
                         M3_ring = "M1 + scaled_lon + scaled_alt + season + (1 | ring), complete cases",
                         MWsrc_mundlak_contributor = "M3 with scaled_yr -> yr_within_src + yr_mean_src",
                         HOUR0_no_hour = "ln mass ~ Sex + scaled_yr + scaled_lat + season + (1 + scaled_yr || spp) + (1 | src) + (1 | site) + phylo, 05-19 h",
                         HOUR1_with_hour = "HOUR0 + hour_num",
                         phylo = "propto(0 + species_name | g, A) appended by fit_phylo_glmmtmb()"),
              intervals = "estimate +/- 1.96 * Rubin SE (within-tree Wald variance + (1 + 1/m) between-tree variance); column t holds the pooled z",
              converged = "n_converged = trees with a positive-definite Hessian after the seeded-restart repair; n_converged_first_pass = before it",
              lme4_tier_generated = lme4_res$meta$generated,
              wing = "wing is not fitted here; the phylogenetic controlled wing models are Phase 1 (controlled_wing_phylo_results.rds)",
              elapsed_min = as.numeric(difftime(Sys.time(), t_phylo, units = "mins"))))
saveRDS(phylo, out_path("multitrait_phylo_results.rds"))
}   # ======================= end of the two tiers =======================

# ---------------------------------------------------------------------------
# 8. Markdown tables (pasted into REVISION_NOTES_P2.md): lme4 tier, then the
#    phylogenetic tier when multitrait_phylo_results.rds is available
# ---------------------------------------------------------------------------
fmt <- function(x, d = 2) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
cell <- function(tr, model, term = "scaled_yr") {
  r <- get(tr, model, term)
  if (!nrow(r)) return("—")
  dg <- if (is_log[[tr]]) 4 else 2
  sprintf("%s [%s, %s] t = %s (n = %s)%s", fmt(r$estimate, dg), fmt(r$lower, dg), fmt(r$upper, dg),
          fmt(r$t, 1), format(r$n, big.mark = ","), if (isTRUE(r$singular)) " (singular)" else "")
}
md <- c(
  sprintf("| Trait | n (%% of 12,571) | contributors (spanning) | M0 baseline β_yr | M1 +src +site | M1 on M3 sample | M2 +lon +alt +season | M3 +ring | Within-contributor (Mundlak) | Between-contributor | Within-municipality |"),
  "|---|---:|---:|---|---|---|---|---|---|---|---|")
for (tr in names(traits)) {
  cv <- coverage[coverage$trait == tr, ]
  md <- c(md, sprintf("| %s | %s (%.1f %%) | %d (%d) | %s | %s | %s | %s | %s | %s | %s | %s |",
    traits[[tr]], format(cv$n, big.mark = ","), cv$pct_of_sample, cv$n_src, cv$n_src_both_periods,
    cell(tr, "M0_baseline"), cell(tr, "M1_src_site"), cell(tr, "M1cc_src_site_on_M2_sample"),
    cell(tr, "M2_geo_season"), cell(tr, "M3_ring"),
    cell(tr, "MWsrc_mundlak_contributor", "yr_within_src"), cell(tr, "MWsrc_mundlak_contributor", "yr_mean_src"),
    cell(tr, "MWsite_mundlak_municipality", "yr_within_site")))
}
md <- c(md, "",
  sprintf("β_yr is per SD-year (SD = %.4f yr); multiply by %.3f for per-decade units. Wald 95 %% CI in brackets; lme4 REML, no phylogenetic term.", SD_YR, DEC),
  sprintf("Body mass diurnal model (n = %s records with 05:00–19:00 clock time): +%.3f %% mass per hour [%.3f, %.3f], t = %.2f; year slope with hour in the model %.4f (t = %.2f).",
          format(hour$n_mass_hour_05_19_complete, big.mark = ","), hour$pct_per_hour, hour$pct_per_hour_lower,
          hour$pct_per_hour_upper, hour$t_hour, hour$year_beta_with_hour, hour$year_t_with_hour))
if (!is.null(phylo)) {
  pcell <- function(tr, model, term = "scaled_yr") {
    r <- get(tr, model, term, tab = phylo$table)
    if (!nrow(r)) return("—")
    dg <- if (is_log[[tr]]) 4 else 2
    sprintf("%s [%s, %s] z = %s (n = %s; %d/%d trees)%s", fmt(r$estimate, dg), fmt(r$lower, dg), fmt(r$upper, dg),
            fmt(r$t, 1), format(r$n, big.mark = ","), r$n_converged, r$n_trees,
            if (isTRUE(r$boundary_spp_slope)) " (species-slope SD at zero)" else "")
  }
  pprop <- function(tr, model) { v <- phylo$varcomp[phylo$varcomp$trait == tr & phylo$varcomp$model == model, ]
    if (!nrow(v)) "—" else sprintf("%.3f (SD %s)", v$phylo_prop_mean, signif(v$sd_phylo_mean, 3)) }
  hp <- phylo$hour
  md <- c(md, "", sprintf("### Phylogenetic tier (glmmTMB propto, %d trees, Rubin-pooled; generated %s)", phylo$meta$n_trees, phylo$meta$generated), "",
    "| Trait | M0 baseline β_yr | M1 +src +site | M3 +lon +alt +season +ring | Within-contributor (Mundlak) | Between-contributor | Phylogenetic proportion, M1 (phylo SD) |",
    "|---|---|---|---|---|---|---|")
  for (tr in setdiff(names(traits), "wing")) {
    md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s |", traits[[tr]],
      pcell(tr, "M0_baseline"), pcell(tr, "M1_src_site"), pcell(tr, "M3_ring"),
      pcell(tr, "MWsrc_mundlak_contributor", "yr_within_src"), pcell(tr, "MWsrc_mundlak_contributor", "yr_mean_src"), pprop(tr, "M1_src_site")))
  }
  cmp <- phylo$comparison_lme4
  md <- c(md, "",
    sprintf("Same fixed and random structure as the lme4 rows plus the phylogenetic species term propto(0 + species_name | g, A) on the %d published trees; interval = estimate ± 1.96 × Rubin SE; \"k/%d trees\" = trees with a positive-definite Hessian after the seeded-restart repair; \"(species-slope SD at zero)\" = the species year-slope variance sits at the boundary (the lme4 \"singular\" analogue; fixed effects unaffected). Wing is Phase 1's.",
            phylo$meta$n_trees, phylo$meta$n_trees),
    sprintf("Phylogenetic vs lme4 (%d matched estimates): median |Δβ| / lme4 SE = %.3f (max %.3f); median SE ratio %.3f (range %.3f–%.3f); zero-exclusion verdict changes in %d of %d.",
            nrow(cmp), median(abs(cmp$delta_in_lme4_se)), max(abs(cmp$delta_in_lme4_se)), median(cmp$se_ratio_phylo_over_lme4),
            min(cmp$se_ratio_phylo_over_lme4), max(cmp$se_ratio_phylo_over_lme4), sum(cmp$verdict_changes), nrow(cmp)),
    sprintf("Body mass diurnal model, phylogenetic (n = %s): +%.3f %% mass per hour [%.3f, %.3f], z = %.2f; year slope with hour %.4f (z = %.2f), without hour %.4f.",
            format(hp$n_mass_hour_05_19_complete, big.mark = ","), hp$pct_per_hour, hp$pct_per_hour_lower, hp$pct_per_hour_upper, hp$z_hour,
            hp$year_beta_with_hour, hp$year_z_with_hour, hp$year_beta_without_hour))
}
writeLines(md, out_path("multitrait_table.md"))

# ---------------------------------------------------------------------------
# 9. Console summary
# ---------------------------------------------------------------------------
if (MODE == "trees") {
  cat("\n=== Phylogenetic tier: year slopes (Rubin-pooled; t column = z) ===\n")
  print(phylo$table %>% filter(term %in% c("scaled_yr", "yr_within_src", "yr_mean_src")) %>%
          transmute(trait, model, term, estimate = signif(estimate, 4), se = signif(se, 3), z = round(t, 2), lower = signif(lower, 4), upper = signif(upper, 4),
                    n, conv = paste0(n_converged, "/", n_trees), pct_per_decade = round(pct_per_decade, 2)) %>% as.data.frame(), row.names = FALSE)
  cat("\n=== Phylogenetic vs lme4 ===\n")
  print(phylo$comparison_lme4 %>% transmute(trait, model, term, lme4 = signif(lme4_estimate, 4), phylo = signif(phylo_estimate, 4),
                                            d_in_se = round(delta_in_lme4_se, 3), se_ratio = round(se_ratio_phylo_over_lme4, 3), verdict_changes) %>% as.data.frame(), row.names = FALSE)
  cat("\n=== Timing (seconds per model, all trees) ===\n"); print(phylo$timing, row.names = FALSE)
  cat(sprintf("\nWrote %s and %s (%.1f min)\n", out_path("multitrait_phylo_results.rds"), out_path("multitrait_table.md"), phylo$meta$elapsed_min))
  quit(save = "no", status = 0)
}
cat("\n=== Coverage ===\n"); print(coverage[, c("trait", "n", "pct_of_sample", "n_src", "n_site", "n_src_both_periods", "n_rec_src_both_periods", "src_both_periods")])
cat("\n=== Year slopes (scaled_yr and Mundlak components) ===\n")
print(table_all %>% filter(term %in% c("scaled_yr", "yr_within_src", "yr_mean_src", "yr_within_site", "yr_mean_site")) %>%
        transmute(trait, model, term, estimate = signif(estimate, 4), se = signif(se, 3), t = round(t, 2),
                  lower = signif(lower, 4), upper = signif(upper, 4), n, pct_per_decade = round(pct_per_decade, 2), singular) %>%
        as.data.frame(), row.names = FALSE)
cat("\n=== Body mass diurnal model ===\n"); print(hour_tab[, c("model", "term", "estimate", "se", "t", "n")])
cat(sprintf("hour effect: %+.3f %% per hour [%.3f, %.3f], t = %.2f; n = %d (05-19 h), %d with any parsed hour\n",
            hour$pct_per_hour, hour$pct_per_hour_lower, hour$pct_per_hour_upper, hour$t_hour,
            hour$n_mass_hour_05_19_complete, hour$n_mass_with_hour_any))
cat("\n=== Isometry (lme4, different samples for wing and mass) ===\n"); print(isometry$lme4)
cat("\n=== Plan comparison (REVISION_PLAN.md §0 table vs this script) ===\n")
print(plan_comparison %>% transmute(trait, plan_n, here_n_M0, plan_M0_beta, here_M0_beta = signif(here_M0_beta, 3),
                                    plan_M0_t, here_M0_t = round(here_M0_t, 2), plan_M1_beta, here_M1_beta = signif(here_M1_beta, 3),
                                    plan_M1_t, here_M1_t = round(here_M1_t, 2), plan_full_beta, here_M3_beta = signif(here_M3_beta, 3),
                                    plan_full_t, here_M3_t = round(here_M3_t, 2), plan_within_src_beta,
                                    here_within_src_beta = signif(here_within_src_beta, 3), plan_within_src_t,
                                    here_within_src_t = round(here_within_src_t, 2)))
cat("\nWrote", out_path("multitrait_results.rds"), "and", out_path("multitrait_table.md"), "\n")
