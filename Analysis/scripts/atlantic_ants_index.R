# atlantic_ants_index.R
# ---------------------------------------------------------------------------
# WHAT: Builds a region-matched, provenance-controlled index of litter-ant
# community richness for the Brazilian Atlantic Forest, 1990-2019, from the
# ATLANTIC ANTS compilation (Silva et al. 2022, Ecology, doi 10.1002/ecy.3580,
# CC-BY), and estimates its temporal trend with contributor / site / method /
# effort / habitat controls. It replaces the PREDICTS South-America "arthropod
# abundance trend" retired in Phase 4e (REVISION_PLAN.md section 3, Phase 4a;
# section 7 documents why this is the best available substitute).
#
# WHAT THE INDEX IS: the number of ant species (distinct Genus + Species, with
# contributor morphospecies codes counted as species within a contributor) found
# in one standardised sampling event, adjusted for the number of Winkler samples
# or pitfall traps in that event. It is a litter/ground-ant COMMUNITY RICHNESS
# index per standardised sample.
# WHAT IT IS NOT: prey biomass, prey availability for birds, or arthropod
# abundance. Only ~10 % of standardised records carry true abundance; those are
# analysed separately as a secondary model and reported with the same caveats.
#
# SAMPLING EVENT (unit of analysis). Primary = CAMPAIGN: contributor file
# (fl.nm) x Method (Winkler | Pitfall) x Start.year/month/day x coordinates
# rounded to 0.01 degree (~1 km) x Habitat.Type, all ID.codLoc codes pooled.
# This is the level at which the metadata define the effort fields ("number of
# Winkler / pitfall samples" for the sampled site) and the only level that is
# coherent across contributors: for several large contributor files ID.codLoc
# is a per-record code (hundreds of "plot-dates" with exactly 1 species and
# 20-600 traps), for PERD 2005-06 it is 80 plot codes sharing one
# Winkler.Number = 360, and in the Servio Ribeiro PERD/PEIT pitfall files the
# effort field itself varies among records of the same plot and year (1..30,
# the same species repeated) - i.e. it is a trap identifier there. Effort is
# therefore built in two steps: (i) per plot-date (campaign key + ID.codLoc)
# effort = max(effort field) - the trap count where the field is an
# identifier, the reported constant otherwise; (ii) per campaign effort = the
# shared value if all plot-dates report the same number (a campaign total),
# else the sum of plot-date efforts. Sensitivity = PLOT-DATE level (step (i)
# only), reported with its richness-equals-1 share as a warning. Neither level
# can be verified for every contributor without the data owners; a
# contributor-constant convention is absorbed by the contributor intercept, so
# only within-contributor changes in reporting convention bias the within
# slope - which is exactly what happens inside PERD (see notes). A
# constant-protocol check therefore keeps, within each contributor x method
# series, only campaigns at the series' modal effort.
# Events whose campaign spans >= 2 calendar years are dropped.
#
# MODELS (glmmTMB, negative binomial 2 unless Poisson is preferred by AIC):
#   (a) pooled:   richness ~ scaled_year + Method * log(effort) + habitat3
#                            + (1|fl.nm) + (1|locality) + (1|site)
#                 plus offset(log(effort)) and no-random-effect variants;
#   (b) Mundlak:  scaled_year split into the contributor mean (between) and the
#                 deviation from it (within); the within slope is the temporal
#                 signal that survives contributor turnover; a within-SITE
#                 version uses only resampled sites;
#   (c) (a)+(b) restricted to monitoring programmes with >= 6 distinct sampling
#                 years (literature-compilation files excluded), per-programme
#                 GLM slopes and a random-slope model;
#   (d) abundance: total ants per event (Measurement.Type == "abundance"),
#                 offset(log(effort)), same random effects;
#   (e) year-as-factor model -> an annual adjusted-richness index (Phase 4f);
#   (f) sensitivity: leave-one-contributor-out on the within slope; pooled
#                 model without literature files; constant-protocol subset;
#                 plot-date event definition.
# Year is standardised with the BIRD constants (centre 2009.487, SD 5.020 yr,
# from output/effect_scale.rds) so ant slopes are on the same "per SD of year"
# scale as the wing model; per-decade and % per-decade conversions are stored.
#
# INPUTS:  data/raw/atlantic_ants/Atlantic_Ants_dataset_v1_2021-05-26d.zip
#          (unzipped to a temp dir at run time; data/raw is never written)
#          output/effect_scale.rds (P0; falls back to hard-coded constants)
# OUTPUTS: data/derived/ants_events.rds            campaign events (primary)
#          data/derived/ants_events_plotdate.rds   plot-date events (sensitivity)
#          data/derived/ants_records_std.rds       standardised records (for 4f)
#          output/ants_results.rds                 all numbers quoted in the paper
#          output/ants_results.md                  human-readable table of the same
#          output/ants_annual_index.rds            year-factor adjusted index
#          figures/ants_events_by_year.png         events per year by contributor
#          figures/ants_year_effect.png            year effects + within-programme fits
#
# RUN:  Rscript Analysis/scripts/atlantic_ants_index.R
#       (works from the repo root, Analysis/, or Analysis/scripts/; ~3-5 min
#        including the leave-one-out loop; no Stan - the full version runs on
#        a laptop, nothing here is a smoke test)
# Session: R 4.6.0, data.table 1.18.4, glmmTMB 1.1.14, ggplot2 4.0.3,
#          patchwork 1.3.2.
# ---------------------------------------------------------------------------

suppressMessages({
  library(data.table); library(glmmTMB); library(ggplot2); library(patchwork)
})
set.seed(20260909)   # glmmTMB is deterministic; recorded for the convention

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
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
dir.create(fig_path(), showWarnings = FALSE, recursive = TRUE)

res <- list(generated = Sys.time(),
            session = c(R = R.version.string,
                        sapply(c("data.table", "glmmTMB", "ggplot2", "patchwork"),
                               function(p) as.character(packageVersion(p)))))
short_name <- function(x) sub("_validchar.*$", "", sub("^ATLANTIC_ANTS_", "", x))

# ---------------------------------------------------------------------------
# 1. Read the raw compilation (tab-delimited, Latin-1, unquoted free text)
# ---------------------------------------------------------------------------
zip_file <- raw_path("atlantic_ants", "Atlantic_Ants_dataset_v1_2021-05-26d.zip")
txt_name <- "ATLANTIC_ANTS_dataset.txt"
txt_file <- raw_path("atlantic_ants", "extracted", txt_name)   # convenience copy, optional
if (!file.exists(txt_file)) {
  stopifnot(file.exists(zip_file))
  tmpd <- tempfile("aants_"); dir.create(tmpd)
  unzip(zip_file, files = txt_name, exdir = tmpd)
  txt_file <- file.path(tmpd, txt_name)
}
raw <- fread(txt_file, sep = "\t", quote = "", encoding = "Latin-1", showProgress = FALSE)
res$n_raw_records <- nrow(raw)
message(sprintf("[ants] raw records: %d, fields: %d, contributor files: %d",
                nrow(raw), ncol(raw), uniqueN(raw$fl.nm)))

# ---------------------------------------------------------------------------
# 2. Region, quality and period filters
# ---------------------------------------------------------------------------
# ma.lmt == 1 is the data authors' own Atlantic Forest limit flag, preferable to
# a bounding box; the bbox is applied as well as a guard and its effect recorded.
d <- raw[Country == "BRAZIL"]
res$n_brazil <- nrow(d)
d <- d[ma.lmt == 1]
res$n_brazil_af_limit <- nrow(d)
in_bbox <- d$Latitude.y >= -31 & d$Latitude.y <= -5 & d$Longitude.x >= -57 & d$Longitude.x <= -34
res$n_outside_bbox_dropped <- sum(!in_bbox, na.rm = TRUE) + sum(is.na(in_bbox))
d <- d[which(in_bbox)]
d <- d[Exclude != 1]
res$n_after_exclude <- nrow(d)
d <- d[!is.na(Start.year) & Start.year >= 1990 & Start.year <= 2019]
res$n_1990_2019 <- nrow(d)

# Effort fields are character in the file (a "5to15" token in Pitfall.Number;
# thousands separators in Total.Ant.Abundance). Non-numeric -> NA.
d[, winkler_n := suppressWarnings(as.integer(Winkler.Number))]
d[, pitfall_n := suppressWarnings(as.integer(Pitfall.Number))]
d[, abundance := suppressWarnings(as.numeric(gsub(",", "", Total.Ant.Abundance)))]
res$pitfall_number_nonnumeric <- d[Method == "Pitfall" & !is.na(Pitfall.Number) &
                                    Pitfall.Number != "NA" & is.na(pitfall_n), .N]

# Standardised records: pure Winkler with known sample number, or pure pitfall
# with known trap number. Mixed-method records ("Pitfall and Winkler", baited
# pitfalls, vegetation pitfalls, ...) are not comparable and are dropped.
d[, std := (Method %in% "Winkler" & !is.na(winkler_n) & winkler_n > 0) |
           (Method %in% "Pitfall" & !is.na(pitfall_n) & pitfall_n > 0)]
res$n_winkler_all      <- d[Method %in% "Winkler", .N]
res$n_winkler_with_n   <- d[Method %in% "Winkler" & std, .N]
res$n_pitfall_all      <- d[Method %in% "Pitfall", .N]
res$n_pitfall_with_n   <- d[Method %in% "Pitfall" & std, .N]
s <- d[std == TRUE]
s[, effort := ifelse(Method == "Winkler", winkler_n, pitfall_n)]

# Campaign span: drop events lumped over >= 2 calendar years (no year resolution)
s[, span := End.year - Start.year]
res$n_span_ge2_dropped <- s[!is.na(span) & (span >= 2 | span < 0), .N]
s <- s[is.na(span) | (span >= 0 & span <= 1)]
res$n_std_records <- nrow(s)
res$share_abundance_records <- mean(s$Measurement.Type == "abundance")

# Taxon: Genus + Species. Morphospecies codes (sp01, ...) are contributor-
# specific, but every event is nested in one contributor file, so they count as
# distinct species within an event without cross-file matching.
s[, taxon := paste(Genus, Species)]

# Literature-compilation files (several source papers per file) are kept in the
# pooled models as contributors, but are not "monitoring programmes".
s[, literature_file := grepl("LITERATURE|PUBLISHED_DATA", fl.nm)]

# Coarse habitat class from the free-text Habitat.Type (200+ levels):
forest_kw <- "forest|floresta|mata|ombroph|arboreal restinga|riparian|cipo|cabruca|araucaria|mangrove|paludosa"
open_kw   <- "pasture|grassland|agro|agri|plantation|eucalyptus|pine|urban|cerrado|vineyard|road|crop|cane|coffee|cocoa|corn|rehabilitation|restoration|dune|brush|campo|savanna|caatinga|carrasco|monoculture|forestry|reflorest"
ht <- tolower(s$Habitat.Type)
is_forest <- grepl(forest_kw, ht); is_open <- grepl(open_kw, ht)
s[, habitat3 := fifelse(is.na(ht), "unknown",
                  fifelse(is_forest & !is_open, "forest",
                  fifelse(is_open & !is_forest, "open",
                  fifelse(is_forest & is_open, "mixed", "other"))))]
res$habitat3_records <- s[, .N, by = habitat3][order(-N)]

# ---------------------------------------------------------------------------
# 3. Standardised sampling events (campaign = primary; plot-date = sensitivity)
# ---------------------------------------------------------------------------
s[, lat2 := round(Latitude.y, 2)][, lon2 := round(Longitude.x, 2)]
s[, locality := sprintf("%.1f_%.1f", round(Latitude.y, 1), round(Longitude.x, 1))]  # ~11 km cell
s[, loc_code := fifelse(is.na(ID.codLoc) | ID.codLoc == "NA", "", ID.codLoc)]
s[, site_campaign := paste(fl.nm, Method, lat2, lon2, sep = "|")]
s[, site_plot     := paste(site_campaign, loc_code, sep = "|")]
s[, event_campaign := paste(site_campaign, Start.year, Start.month, Start.day, Habitat.Type, sep = "|")]
s[, event_plot     := paste(site_plot,     Start.year, Start.month, Start.day, Habitat.Type, sep = "|")]
# Plot-date effort: max of the effort field; flag plot-dates where the field
# varies among records (a sample identifier, not a design count)
s[, effort_raw := effort]
s[, `:=`(effort = max(effort_raw), effort_is_id = uniqueN(effort_raw) > 1L), by = event_plot]
res$effort_is_id <- list(
  n_records = s[effort_is_id == TRUE, .N],
  n_plot_dates = s[, .(id = effort_is_id[1]), by = event_plot][id == TRUE, .N],
  by_contributor = s[, .(plot_dates = uniqueN(event_plot),
                         share_plot_dates_varying = uniqueN(event_plot[effort_is_id]) / uniqueN(event_plot),
                         n_records = .N), by = .(contributor = short_name(fl.nm), Method)][share_plot_dates_varying > 0][order(-n_records)])

# Year standardised with the BIRD constants so slopes are per SD of the wing model
yr_center <- 2009.487; yr_scale <- 5.020
es_file <- out_path("effect_scale.rds")
if (file.exists(es_file)) {
  es <- readRDS(es_file)
  if (!is.null(es$scaling_reference$year)) {
    yr_center <- unname(es$scaling_reference$year["center"])
    yr_scale  <- unname(es$scaling_reference$year["scale"])
  }
}
res$year_scaling <- c(center = yr_center, scale = yr_scale)
per_decade <- 10 / yr_scale   # multiply a per-SD slope by this to get per decade

build_events <- function(s, level = c("campaign", "plot")) {
  level <- match.arg(level)
  s <- copy(s)
  s[, event := if (level == "campaign") event_campaign else event_plot]
  s[, site  := if (level == "campaign") site_campaign  else site_plot]
  ev <- s[, .(fl.nm = fl.nm[1], contributor = short_name(fl.nm[1]), literature_file = literature_file[1],
              Method = Method[1], year = Start.year[1], month = Start.month[1],
              lat = mean(Latitude.y), lon = mean(Longitude.x), lat2 = lat2[1], lon2 = lon2[1],
              locality = locality[1], site = site[1], n_plot_codes = uniqueN(loc_code),
              habitat_type = Habitat.Type[1], habitat3 = habitat3[1], disturbance = Disturbance[1],
              # plot level: the plot-date effort; campaign level: shared value if all plots
              # report the same number (campaign total), else the sum of plot efforts
              effort = { pe <- unique(data.table(loc_code, effort))$effort
                         as.numeric(if (uniqueN(pe) == 1L) pe[1] else sum(pe)) },
              effort_is_id = any(effort_is_id),
              n_records = .N, richness = uniqueN(taxon), n_genera = uniqueN(Genus),
              all_abundance = all(Measurement.Type == "abundance"),
              abundance = if (all(Measurement.Type == "abundance") && !anyNA(abundance)) sum(abundance) else NA_real_,
              source = Source.Data[1]),
          by = event]
  ev[, scaled_year := (year - yr_center) / yr_scale]
  ev[, log_effort := log(effort)]
  ev[, Method := factor(Method, levels = c("Winkler", "Pitfall"))]
  ev[, habitat3 := droplevels(factor(habitat3, levels = c("forest", "open", "mixed", "other", "unknown")))]
  ev[, year_between := mean(scaled_year), by = fl.nm]
  ev[, year_within  := scaled_year - year_between]
  ev[, year_site_mean := mean(scaled_year), by = site]
  ev[, year_within_site := scaled_year - year_site_mean]
  setattr(ev, "event_level", level)
  ev
}
ev  <- build_events(s, "campaign")   # primary
evP <- build_events(s, "plot")       # sensitivity

summarise_events <- function(ev) list(
  n_events = nrow(ev), n_events_by_method = ev[, .N, by = Method],
  n_contributors = uniqueN(ev$fl.nm), n_localities_11km = uniqueN(ev$locality),
  n_sites = uniqueN(ev$site), n_sites_multi_year = ev[, uniqueN(year), by = site][V1 > 1, .N],
  n_events_sites_multi_year = ev[site %in% ev[, uniqueN(year), by = site][V1 > 1, site], .N],
  events_per_year = ev[, .N, by = year][order(year)], years_covered = sort(unique(ev$year)),
  share_events_multi_plot = mean(ev$n_plot_codes > 1),
  richness_summary_by_method = ev[, as.list(summary(richness)), by = Method],
  effort_quantiles_by_method = ev[, as.list(quantile(effort, c(0, .1, .5, .9, 1))), by = Method],
  n_events_abundance = ev[!is.na(abundance) & abundance > 0, .N],
  share_richness_1 = mean(ev$richness == 1),          # per-record codes masquerading as plots
  share_multi_code = mean(ev$n_plot_codes > 1))
res$events_campaign <- summarise_events(ev)
res$events_plot     <- summarise_events(evP)
res$richness1_by_contributor_plotlevel <- evP[, .(n_plot_dates = .N, share_richness_1 = mean(richness == 1),
                                                 median_effort = as.numeric(median(effort))), by = contributor][share_richness_1 > 0.5 & n_plot_dates >= 20][order(-n_plot_dates)]
res$n_events_alt <- c(contributor_method_coord_year = s[, uniqueN(paste(fl.nm, Method, lat2, lon2, Start.year))],
                      campaign_used = nrow(ev), plot_date_level = nrow(evP))

contrib <- ev[, .(contributor = contributor[1], literature_file = literature_file[1], n_events = .N,
                  n_years = uniqueN(year), year_min = min(year), year_max = max(year),
                  methods = paste(sort(unique(as.character(Method))), collapse = "/"),
                  n_sites = uniqueN(site), n_abund = sum(!is.na(abundance)),
                  mean_richness = mean(richness), median_effort = as.numeric(median(effort))),
              by = fl.nm][order(-n_years, -n_events)]
res$contributors <- contrib
long_files <- contrib[n_years >= 6 & !literature_file, fl.nm]
res$long_programmes <- contrib[fl.nm %in% long_files]
res$long_programmes_excluded_literature <- contrib[n_years >= 6 & literature_file, .(fl.nm, n_events, n_years)]
ev[, long_programme := fl.nm %in% long_files]

saveRDS(ev,  derived_path("ants_events.rds"))
saveRDS(evP, derived_path("ants_events_plotdate.rds"))
saveRDS(s[, .(event_campaign, event_plot, fl.nm, literature_file, Method, Start.year, Start.month, Start.day,
              Latitude.y, Longitude.x, locality, site_campaign, site_plot, loc_code, Habitat.Type, habitat3,
              effort_raw, effort, effort_is_id, Genus, Species, taxon, Morphospecies, Measurement.Type, abundance, Source.Data)],
        derived_path("ants_records_std.rds"))
message(sprintf("[ants] standardised records %d -> campaign events %d (%d Winkler, %d pitfall; plot-date level %d), %d contributors, %d 11-km localities, %d sites (%d sampled in >1 year), years %d-%d",
                nrow(s), nrow(ev), ev[Method == "Winkler", .N], ev[Method == "Pitfall", .N], nrow(evP),
                res$events_campaign$n_contributors, res$events_campaign$n_localities_11km,
                res$events_campaign$n_sites, res$events_campaign$n_sites_multi_year, min(ev$year), max(ev$year)))

# ---------------------------------------------------------------------------
# 4. Model helpers
# ---------------------------------------------------------------------------
fit_ok <- function(fit) {
  (is.null(fit$fit$convergence) || fit$fit$convergence == 0) && !isTRUE(fit$sdr$pdHess == FALSE)
}
report_fit <- function(fit, label) {
  message(sprintf("[ants] %-34s converged=%s pdHess=%s AIC=%.1f", label, fit_ok(fit),
                  isTRUE(fit$sdr$pdHess), AIC(fit)))
  invisible(fit)
}
tidy_year <- function(fit, terms, label, n, level = "campaign") {
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients$cond
  ci <- tryCatch(confint(fit, parm = "beta_", method = "wald", level = 0.95), error = function(e) NULL)
  rbindlist(lapply(terms, function(t) {
    if (!t %in% rownames(cf)) return(NULL)
    lo <- if (!is.null(ci) && t %in% rownames(ci)) ci[t, 1] else cf[t, 1] - 1.96 * cf[t, 2]
    hi <- if (!is.null(ci) && t %in% rownames(ci)) ci[t, 2] else cf[t, 1] + 1.96 * cf[t, 2]
    data.table(event_level = level, model = label, term = t, n_events = n, converged = fit_ok(fit),
               est = cf[t, 1], se = cf[t, 2], lo = lo, hi = hi, z = cf[t, 3], p = cf[t, 4],
               est_decade = cf[t, 1] * per_decade, lo_decade = lo * per_decade, hi_decade = hi * per_decade,
               pct_decade = 100 * (exp(cf[t, 1] * per_decade) - 1),
               pct_lo_decade = 100 * (exp(lo * per_decade) - 1),
               pct_hi_decade = 100 * (exp(hi * per_decade) - 1))
  }))
}
RE3 <- "(1 | fl.nm) + (1 | locality) + (1 | site)"
f_pool    <- as.formula(paste("richness ~ scaled_year + Method * log_effort + habitat3 +", RE3))
f_offset  <- as.formula(paste("richness ~ scaled_year + Method + habitat3 + offset(log_effort) +", RE3))
f_mundlak <- as.formula(paste("richness ~ year_within + year_between + Method * log_effort + habitat3 +", RE3))
f_site    <- as.formula(paste("richness ~ year_within_site + year_site_mean + Method * log_effort + habitat3 +", RE3))

# Core model set, run on any event table (primary and sensitivity levels)
fit_core <- function(ev, level, fam = NULL) {
  m_pois <- glmmTMB(f_pool, data = ev, family = poisson)
  m_nb   <- glmmTMB(f_pool, data = ev, family = nbinom2)
  pr <- residuals(m_pois, type = "pearson")
  use_nb <- if (is.null(fam)) AIC(m_nb) < AIC(m_pois) else identical(fam, "nbinom2")
  fam_fun <- if (use_nb) nbinom2 else poisson
  m_a <- if (use_nb) m_nb else m_pois
  m_a_off <- glmmTMB(f_offset, data = ev, family = fam_fun)
  m_naive <- glmmTMB(richness ~ scaled_year + Method * log_effort, data = ev, family = fam_fun)
  m_b     <- glmmTMB(f_mundlak, data = ev, family = fam_fun)
  m_bs    <- glmmTMB(f_site, data = ev, family = fam_fun)
  list(models = list(pooled = m_a, pooled_offset = m_a_off, naive = m_naive, mundlak = m_b, mundlak_site = m_bs),
       family = if (use_nb) "nbinom2" else "poisson", aic = c(poisson = AIC(m_pois), nbinom2 = AIC(m_nb)),
       poisson_dispersion = sum(pr^2) / df.residual(m_pois),
       nb_theta = if (use_nb) sigma(m_nb) else NA_real_,
       year_effects = rbindlist(list(
         tidy_year(m_naive, "scaled_year", "naive: no random effects", nrow(ev), level),
         tidy_year(m_a,     "scaled_year", "pooled: contributor + locality + site RE, effort covariate", nrow(ev), level),
         tidy_year(m_a_off, "scaled_year", "pooled: effort as offset", nrow(ev), level),
         tidy_year(m_b,     c("year_within", "year_between"), "Mundlak by contributor", nrow(ev), level),
         tidy_year(m_bs,    c("year_within_site", "year_site_mean"), "Mundlak by site", nrow(ev), level))))
}

# ---------------------------------------------------------------------------
# 5. Primary analysis (campaign events)
# ---------------------------------------------------------------------------
core <- fit_core(ev, "campaign")
fam <- if (core$family == "nbinom2") nbinom2 else poisson
res$family_used <- core$family; res$aic_pooled <- core$aic
res$poisson_dispersion <- core$poisson_dispersion; res$nb_theta <- core$nb_theta
m_a <- core$models$pooled; m_b <- core$models$mundlak
res$fixed_effects_pooled  <- summary(m_a)$coefficients$cond
res$fixed_effects_mundlak <- summary(m_b)$coefficients$cond
res$random_sd_pooled <- sapply(VarCorr(m_a)$cond, function(x) attr(x, "stddev"))
V <- vcov(m_b)$cond
dif <- fixef(m_b)$cond["year_within"] - fixef(m_b)$cond["year_between"]
se_dif <- sqrt(V["year_within", "year_within"] + V["year_between", "year_between"] - 2 * V["year_within", "year_between"])
res$mundlak_within_minus_between <- c(est = unname(dif), se = se_dif, z = unname(dif / se_dif),
                                      p = 2 * pnorm(-abs(unname(dif / se_dif))))
res$cor_year_logeffort <- c(overall = cor(ev$year, ev$log_effort),
                            within_contributor = cor(ev$year_within, ev$log_effort - ev[, mean(log_effort), by = fl.nm][ev, on = "fl.nm", V1]))

# (c) long-running monitoring programmes (>= 6 distinct years, not literature files)
evL <- ev[long_programme == TRUE]
evL[, year_between := mean(scaled_year), by = fl.nm][, year_within := scaled_year - year_between]
evL[, habitat3 := droplevels(habitat3)]
m_c  <- glmmTMB(f_pool, data = evL, family = fam)
m_cb <- glmmTMB(f_mundlak, data = evL, family = fam)
res$fixed_effects_long_mundlak <- summary(m_cb)$coefficients$cond
# random slopes per programme: try the full covariance, then a diagonal one, then drop the site RE
rs_forms <- list(
  full = richness ~ scaled_year + Method * log_effort + habitat3 + (1 + scaled_year | fl.nm) + (1 | locality) + (1 | site),
  diag = richness ~ scaled_year + Method * log_effort + habitat3 + diag(1 + scaled_year | fl.nm) + (1 | locality) + (1 | site),
  diag_nosite = richness ~ scaled_year + Method * log_effort + habitat3 + diag(1 + scaled_year | fl.nm) + (1 | locality))
m_cs <- NULL; res$random_slope_spec <- NA_character_
for (nm in names(rs_forms)) {
  fit <- tryCatch(suppressWarnings(glmmTMB(rs_forms[[nm]], data = evL, family = fam)), error = function(e) NULL)
  if (!is.null(fit) && fit_ok(fit)) { m_cs <- fit; res$random_slope_spec <- nm; break }
}
prog_slopes <- if (!is.null(m_cs)) {
  re_cs <- ranef(m_cs)$cond$fl.nm
  data.table(fl.nm = rownames(re_cs), slope_conditional = fixef(m_cs)$cond["scaled_year"] + re_cs[, "scaled_year"])
} else data.table(fl.nm = long_files, slope_conditional = NA_real_)
res$random_slope_sd <- if (!is.null(m_cs)) attr(VarCorr(m_cs)$cond$fl.nm, "stddev")["scaled_year"] else NA_real_
# independent per-programme GLMs (transparent, no shrinkage)
prog_glm <- rbindlist(lapply(long_files, function(f) {
  dd <- evL[fl.nm == f]
  eff_terms <- if (uniqueN(dd$Method) > 1 && dd[, uniqueN(effort), by = Method][, all(V1 > 1)]) "Method * log_effort" else
               if (uniqueN(dd$Method) > 1) "Method" else if (uniqueN(dd$effort) > 1) "log_effort" else NULL
  form <- as.formula(paste("richness ~ scaled_year", if (!is.null(eff_terms)) paste("+", eff_terms) else ""))
  fit <- tryCatch(glmmTMB(form, data = dd, family = fam), error = function(e) NULL)
  out <- data.table(fl.nm = f, contributor = short_name(f), n_events = nrow(dd), n_years = uniqueN(dd$year),
                    year_min = min(dd$year), year_max = max(dd$year),
                    methods = paste(sort(unique(as.character(dd$Method))), collapse = "/"),
                    est = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_)
  if (!is.null(fit)) {
    cf <- summary(fit)$coefficients$cond
    out[, `:=`(est = cf["scaled_year", 1], se = cf["scaled_year", 2],
               lo = cf["scaled_year", 1] - 1.96 * cf["scaled_year", 2], hi = cf["scaled_year", 1] + 1.96 * cf["scaled_year", 2])]
  }
  out
}))
prog_glm <- merge(prog_glm, prog_slopes, by = "fl.nm", all.x = TRUE)
prog_glm[, `:=`(pct_decade = 100 * (exp(est * per_decade) - 1),
                pct_lo_decade = 100 * (exp(lo * per_decade) - 1), pct_hi_decade = 100 * (exp(hi * per_decade) - 1))]
res$programme_slopes <- prog_glm[order(-n_years)]

# (d) abundance (secondary)
evA <- ev[!is.na(abundance) & abundance > 0]
evA[, habitat3 := droplevels(habitat3)]
evA[, year_between := mean(scaled_year), by = fl.nm][, year_within := scaled_year - year_between]
res$abundance_contributors <- evA[, .(contributor = contributor[1], n_events = .N, n_years = uniqueN(year), year_min = min(year),
                                      year_max = max(year), methods = paste(unique(as.character(Method)), collapse = "/")),
                                  by = fl.nm][order(-n_years)]
m_d <- m_db <- NULL
if (nrow(evA) >= 30 && uniqueN(evA$fl.nm) >= 2) {
  m_d  <- glmmTMB(as.formula(paste("abundance ~ scaled_year + Method + habitat3 + offset(log_effort) +", RE3)), data = evA, family = nbinom2)
  m_db <- glmmTMB(as.formula(paste("abundance ~ year_within + year_between + Method + habitat3 + offset(log_effort) +", RE3)), data = evA, family = nbinom2)
}

# (e) year-as-factor adjusted annual index; reference = year with most contributors
ref_year <- ev[, uniqueN(fl.nm), by = year][order(-V1, year)][1, year]
ev[, year_f := relevel(factor(year), ref = as.character(ref_year))]
m_e <- glmmTMB(as.formula(paste("richness ~ year_f + Method * log_effort + habitat3 +", RE3)), data = ev, family = fam)
cf_e <- summary(m_e)$coefficients$cond
yr_rows <- grep("^year_f", rownames(cf_e), value = TRUE)
annual <- rbind(data.table(year = ref_year, est = 0, se = NA_real_),
                data.table(year = as.integer(sub("year_f", "", yr_rows)), est = cf_e[yr_rows, 1], se = cf_e[yr_rows, 2]))
annual <- merge(annual, ev[, .(n_events = .N, n_contributors = uniqueN(fl.nm)), by = year], by = "year")
annual[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
setattr(annual, "reference_year", ref_year)
setattr(annual, "note", "log adjusted richness relative to the reference year, from a year-factor GLMM on campaign events with contributor/locality/site random intercepts, Method x log(effort) and habitat3")
saveRDS(annual, out_path("ants_annual_index.rds"))
res$annual_index <- annual

# (f) sensitivity: leave-one-contributor-out on the Mundlak within slope
loo_files <- contrib[n_events >= 20, fl.nm]
loo <- rbindlist(lapply(loo_files, function(f) {
  dd <- ev[fl.nm != f]
  dd[, year_between := mean(scaled_year), by = fl.nm][, year_within := scaled_year - year_between]
  dd[, habitat3 := droplevels(habitat3)]
  fit <- tryCatch(suppressWarnings(glmmTMB(f_mundlak, data = dd, family = fam)), error = function(e) NULL)
  if (is.null(fit)) return(data.table(dropped = short_name(f), n_events = nrow(dd)))
  cf <- summary(fit)$coefficients$cond
  data.table(dropped = short_name(f), n_events_dropped = nrow(ev) - nrow(dd), n_events = nrow(dd), converged = fit_ok(fit),
             within = cf["year_within", 1], within_se = cf["year_within", 2],
             within_pct_decade = 100 * (exp(cf["year_within", 1] * per_decade) - 1))
}), fill = TRUE)
res$loo_within <- loo[order(within)]
# sensitivity: without literature-compilation files
evNL <- ev[literature_file == FALSE]
evNL[, year_between := mean(scaled_year), by = fl.nm][, year_within := scaled_year - year_between]
evNL[, habitat3 := droplevels(habitat3)]
m_nl  <- glmmTMB(f_pool,    data = evNL, family = fam)
m_nlb <- glmmTMB(f_mundlak, data = evNL, family = fam)
# sensitivity: constant-protocol subset - within each contributor x method series keep
# only campaigns at the series' modal effort, so the within slope compares like with like
ev[, modal_effort := as.numeric(names(which.max(table(effort)))), by = .(fl.nm, Method)]
evK <- ev[effort == modal_effort]
evK[, year_between := mean(scaled_year), by = fl.nm][, year_within := scaled_year - year_between]
evK[, habitat3 := droplevels(habitat3)]
m_k  <- glmmTMB(f_pool,    data = evK, family = fam)
m_kb <- glmmTMB(f_mundlak, data = evK, family = fam)
res$constant_protocol <- list(n_events = nrow(evK), n_contributors = uniqueN(evK$fl.nm),
                              n_series_multi_year = evK[, uniqueN(year), by = .(fl.nm, Method)][V1 > 1, .N],
                              series = evK[, .(n_events = .N, n_years = uniqueN(year), year_min = min(year), year_max = max(year),
                                               effort = effort[1]), by = .(contributor, Method)][n_years > 1][order(-n_years)])
# sensitivity: plot-date event definition (full core set; ID.codLoc is a per-record code for some files)
coreP <- fit_core(evP, "plot", fam = core$family)
res$plot_level <- list(family = coreP$family, aic = coreP$aic, poisson_dispersion = coreP$poisson_dispersion,
                       fixed_effects_pooled = summary(coreP$models$pooled)$coefficients$cond)

# ---------------------------------------------------------------------------
# 6. Collect the year effects
# ---------------------------------------------------------------------------
year_effects <- rbindlist(list(
  core$year_effects,
  tidy_year(m_c,   "scaled_year", "monitoring programmes (>= 6 yr): pooled", nrow(evL)),
  tidy_year(m_cb,  c("year_within", "year_between"), "monitoring programmes (>= 6 yr): Mundlak by contributor", nrow(evL)),
  tidy_year(m_nl,  "scaled_year", "without literature files: pooled", nrow(evNL)),
  tidy_year(m_nlb, c("year_within", "year_between"), "without literature files: Mundlak by contributor", nrow(evNL)),
  tidy_year(m_k,   "scaled_year", "constant protocol (modal effort per series): pooled", nrow(evK)),
  tidy_year(m_kb,  c("year_within", "year_between"), "constant protocol (modal effort per series): Mundlak by contributor", nrow(evK)),
  tidy_year(m_d,   "scaled_year", "abundance (secondary): pooled, offset", nrow(evA)),
  tidy_year(m_db,  c("year_within", "year_between"), "abundance (secondary): Mundlak by contributor", nrow(evA)),
  coreP$year_effects
), fill = TRUE)
res$year_effects <- year_effects
res$models <- c(core$models,
                list(long_pooled = m_c, long_mundlak = m_cb, long_random_slope = m_cs,
                     no_literature_pooled = m_nl, no_literature_mundlak = m_nlb,
                     constant_protocol_pooled = m_k, constant_protocol_mundlak = m_kb,
                     abundance = m_d, abundance_mundlak = m_db, year_factor = m_e),
                setNames(coreP$models, paste0("plotdate_", names(coreP$models))))
for (nm in names(res$models)) if (!is.null(res$models[[nm]])) report_fit(res$models[[nm]], nm)
res$converged <- lapply(res$models, function(m) if (is.null(m)) NA else fit_ok(m))
saveRDS(res, out_path("ants_results.rds"))

# Markdown table of the quoted numbers
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
md <- c("# ATLANTIC ANTS litter-ant richness index: year effects",
        "",
        sprintf("Generated %s. Campaign events = %d (%d Winkler, %d pitfall) from %d contributor files, %d ~11-km localities, %d sites (%d resampled in >1 year); years %d-%d. Family: %s (Poisson dispersion %.2f). Year is per SD of the bird model (%.3f yr); per-decade columns multiply by %.3f. Plot-date sensitivity: %d events (%.0f %% with richness = 1).",
                format(res$generated), nrow(ev), ev[Method == "Winkler", .N], ev[Method == "Pitfall", .N],
                res$events_campaign$n_contributors, res$events_campaign$n_localities_11km, res$events_campaign$n_sites,
                res$events_campaign$n_sites_multi_year, min(ev$year), max(ev$year),
                res$family_used, res$poisson_dispersion, yr_scale, per_decade, nrow(evP), 100 * res$events_plot$share_richness_1),
        "",
        "| Level | Model | Term | N events | Estimate (per SD yr) | 95% CI | z | % per decade [95% CI] |",
        "|---|---|---|---|---|---|---|---|",
        year_effects[, sprintf("| %s | %s | %s | %d | %s | [%s, %s] | %s | %s [%s, %s] |", event_level, model, term, n_events,
                               fmt(est), fmt(lo), fmt(hi), fmt(z, 2), fmt(pct_decade, 1), fmt(pct_lo_decade, 1), fmt(pct_hi_decade, 1))],
        "",
        sprintf("Mundlak within - between (contributor, campaign level): %s (SE %s, z = %s, p = %s).",
                fmt(res$mundlak_within_minus_between["est"]), fmt(res$mundlak_within_minus_between["se"]),
                fmt(res$mundlak_within_minus_between["z"], 2), fmt(res$mundlak_within_minus_between["p"], 3)),
        "",
        "## Monitoring programmes (>= 6 sampling years, literature files excluded): independent per-programme slopes",
        "",
        "| Contributor file | N events | Years | Methods | Slope per SD yr (SE) | % per decade [95% CI] | Random-slope conditional |",
        "|---|---|---|---|---|---|---|",
        prog_glm[order(-n_years), sprintf("| %s | %d | %d (%d-%d) | %s | %s (%s) | %s [%s, %s] | %s |", contributor, n_events, n_years, year_min, year_max,
                                          methods, fmt(est), fmt(se), fmt(pct_decade, 1), fmt(pct_lo_decade, 1), fmt(pct_hi_decade, 1), fmt(slope_conditional))],
        "",
        "## Leave-one-contributor-out (contributors with >= 20 events): Mundlak within slope",
        "",
        "| Dropped | Events dropped | Within slope (SE) | % per decade |",
        "|---|---|---|---|",
        loo[order(within), sprintf("| %s | %d | %s (%s) | %s |", dropped, n_events_dropped, fmt(within), fmt(within_se), fmt(within_pct_decade, 1))])
writeLines(md, out_path("ants_results.md"))

# ---------------------------------------------------------------------------
# 7. Figures
# ---------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11))
top <- contrib[order(-n_events)][1:min(9, .N), fl.nm]
ev[, contrib_plot := fifelse(fl.nm %in% top, contributor, "other contributors")]
ev[, contrib_plot := factor(contrib_plot, levels = c(short_name(top), "other contributors"))]
p1 <- ggplot(ev[, .N, by = .(year, contrib_plot, Method)], aes(year, N, fill = contrib_plot)) +
  geom_col() + facet_wrap(~Method, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(scales::hue_pal()(length(top)), "grey70"), name = "Contributor file") +
  labs(x = "Sampling year", y = "Standardised campaign events",
       title = "ATLANTIC ANTS: standardised sampling events per year",
       subtitle = sprintf("%d campaign events, %d contributor files; Brazilian Atlantic Forest, 1990-2019", nrow(ev), res$events_campaign$n_contributors)) +
  theme(legend.position = "right", legend.text = element_text(size = 7))
ggsave(fig_path("ants_events_by_year.png"), p1, width = 10, height = 6.5, dpi = 200, bg = "white")

ye <- copy(year_effects)
ye[, label := fifelse(grepl("Mundlak", model), paste0(model, ": ", sub("year_", "", term)), model)]
ye[, label := paste0(fifelse(event_level == "plot", "[plot-date level] ", ""), label)]
ye[, label := factor(label, levels = rev(unique(label)))]
ye[, kind := fifelse(grepl("^abundance", model), "abundance", fifelse(event_level == "plot", "richness (plot-date sensitivity; unreliable unit)", "richness"))]
xlim_hi <- min(max(ye$pct_hi_decade, na.rm = TRUE), 150)
p2 <- ggplot(ye, aes(x = pct_decade, y = label, colour = kind)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo_decade, xmax = pct_hi_decade), width = 0.25, orientation = "y") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(richness = "black", `richness (plot-date sensitivity; unreliable unit)` = "grey55", abundance = "firebrick"), name = NULL) +
  coord_cartesian(xlim = c(-100, xlim_hi)) +
  labs(x = "% change per decade (95% Wald CI; axis truncated at 150 %)", y = NULL,
       title = "Year effect on litter-ant richness per standardised sampling event") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8))

pred_lines <- rbindlist(lapply(long_files, function(f) {
  dd <- evL[fl.nm == f]
  yrs <- seq(min(dd$year), max(dd$year), length.out = 20)
  meth <- names(which.max(table(dd$Method)))
  nd <- data.table(year = yrs, scaled_year = (yrs - yr_center) / yr_scale,
                   Method = factor(meth, levels = levels(evL$Method)),
                   log_effort = median(dd[Method == meth, log_effort]),
                   habitat3 = factor(names(which.max(table(dd$habitat3))), levels = levels(evL$habitat3)),
                   fl.nm = f, locality = "new", site = "new")
  nd[, year_between := mean(dd$scaled_year)][, year_within := scaled_year - year_between]
  nd[, fit := predict(m_cb, newdata = nd, type = "response", allow.new.levels = TRUE)]
  nd[, contributor := dd$contributor[1]]
  nd
}))
p3 <- ggplot(evL, aes(year, richness)) +
  geom_point(aes(colour = Method, size = effort), alpha = 0.45) +
  geom_line(data = pred_lines, aes(year, fit), colour = "black", linewidth = 0.9) +
  facet_wrap(~contributor, scales = "free", labeller = label_wrap_gen(30)) +
  scale_y_log10() + scale_size(range = c(0.6, 3.5), name = "Effort (samples/traps)") +
  labs(x = "Sampling year", y = "Species per campaign event (log scale)",
       title = "Monitoring programmes with >= 6 sampling years: within-contributor fit",
       subtitle = "Line: Mundlak model prediction at the programme's modal method and median effort") +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))
p_all <- p2 / p3 + plot_layout(heights = c(1.15, 1))
ggsave(fig_path("ants_year_effect.png"), p_all, width = 11, height = 12, dpi = 200, bg = "white")

message("[ants] year effects (per SD of year = ", fmt(yr_scale, 3), " yr; % per decade):")
print(year_effects[, .(event_level, model, term, n_events, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3),
                       pct_decade = round(pct_decade, 1), pct_lo = round(pct_lo_decade, 1), pct_hi = round(pct_hi_decade, 1))])
message("[ants] wrote ", derived_path("ants_events.rds"), ", ", out_path("ants_results.rds"),
        ", ", out_path("ants_results.md"), ", ", out_path("ants_annual_index.rds"), " and two figures.")
