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
#          output/multitrait_table.md           markdown table pasted into
#          REVISION_NOTES_P2.md
#
# RUN:  Rscript Analysis/scripts/atlantic_multitrait.R
#       (from the repo root, Analysis/, or Analysis/scripts/; ~1-2 min)
# Session used for the numbers in REVISION_NOTES_P2.md: R 4.6.0, lme4 2.0.1,
# Matrix, dplyr 1.1.x. Nothing stochastic; set.seed recorded by convention.
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(lme4) })
set.seed(20260909)

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
get <- function(tr, model, term = "scaled_yr") table_all[table_all$trait == tr & table_all$model == model & table_all$term == term, ]
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

# ---------------------------------------------------------------------------
# 8. Markdown table (pasted into REVISION_NOTES_P2.md)
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
writeLines(md, out_path("multitrait_table.md"))

# ---------------------------------------------------------------------------
# 9. Console summary
# ---------------------------------------------------------------------------
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
