# audit_provenance.R
# ---------------------------------------------------------------------------
# WHAT: Phase 0 of REVISION_PLAN.md — the provenance and design audit of the
# live-only analytical sample. No model fitting. Produces the tables and
# figures that (i) quantify the contributor / site / protocol / season /
# individual structure the referee flagged (§2.1, §2.2, §3.5, §3.6), (ii) put
# the published effect sizes on interpretable scales (§2.4), and (iii) check
# the reference numbers already quoted in REVISION_PLAN.md §0–§1, flagging any
# that differ (see `plan_checks` in every output object and the console).
#
# WHY: every number the revised Results / Supplement quote about sampling
# structure must come from one reproducible place, readable by index.qmd.
#
# INPUTS (from rebuild_passer90_live.R):
#   data/derived/passer90.rda         known-sex sample (12,571 rec, 73 spp)
#   data/derived/passer90_allsex.rda  same species + unknown-sex records
#   output/descriptive_summary.rds    published wing / bill estimates (optional)
#   output/mass_results.rds           published mass estimate (optional)
#   ../Manuscript/data/south_america.geojson  basemap (optional)
#
# OUTPUTS (Analysis/output/, all lists; scalars + data.frames):
#   audit_sources.rds      contributors per year, per-contributor table, spanning
#                          flags, wing-column mix, long-running contributors,
#                          contributors per species x period cell
#   audit_sites.rds        localities / municipalities / coordinate sites per
#                          period and shared, species x site shared, lon/lat/alt
#                          by period, map data
#   audit_individuals.rds  ring completeness, duplicated rings, recapture share
#   audit_season.rds       month / season / moult / reproductive stage by period,
#                          Date parsing, sex ratio by period
#   audit_traits.rds       completeness and provenance for the 6 traits + Hour
#   effect_scale.rds       SD(year), mean wing, mm/decade and % conversions,
#                          isometric expectation, comparator placeholders
# FIGURES (Analysis/figures/):
#   audit_records_per_year_by_source.png, audit_sites_map_shared.png,
#   audit_wingcol_by_year.png
#
# CONVENTIONS: early = Year <= 2006, late = Year >= 2013 (the quartile
# comparison of update_descriptive_stats.R); "wing records" = non-NA
# conc.wing.length. Counts labelled `_wing` are on wing records, the sample of
# the primary model; unlabelled counts are on all 12,571 known-sex records.
#
# RUN:  Rscript Analysis/scripts/audit_provenance.R      (< 1 min)
# Session: R 4.6.0, dplyr 1.1.x, tidyr 1.3.2, ggplot2 4.0.3, sf 1.1.2,
#          patchwork 1.3.2 (see R_session_info.txt)
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
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
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
dir.create(out_path(), showWarnings = FALSE, recursive = TRUE)
dir.create(fig_path(), showWarnings = FALSE, recursive = TRUE)

load(derived_path("passer90.rda"))          # passer90 (known sex)
load(derived_path("passer90_allsex.rda"))   # passer90_allsex
stopifnot(nrow(passer90) == 12571L, length(unique(passer90$Binomial)) == 73L)

EARLY_MAX <- 2006L; LATE_MIN <- 2013L
period_of <- function(y) ifelse(y <= EARLY_MAX, "early", ifelse(y >= LATE_MIN, "late", "mid"))
d  <- passer90 %>% mutate(period = period_of(Year), has_wing = !is.na(conc.wing.length))
da <- passer90_allsex %>% mutate(period = period_of(Year), has_wing = !is.na(conc.wing.length))
w  <- d %>% filter(has_wing)                      # primary-model sample (8,478)
wq <- w %>% filter(period != "mid")               # wing records in the two quartile periods

# --- plan-check bookkeeping ------------------------------------------------
plan_checks <- data.frame(section = character(), item = character(), plan = numeric(),
                          reproduced = numeric(), match = logical(), definition = character(),
                          stringsAsFactors = FALSE)
chk <- function(section, item, plan, reproduced, tol = 0, definition = "") {
  ok <- is.finite(plan) && is.finite(reproduced) && abs(plan - reproduced) <= tol
  plan_checks <<- rbind(plan_checks, data.frame(section = section, item = item, plan = plan,
                                                reproduced = reproduced, match = ok,
                                                definition = definition))
  invisible(ok)
}
pct <- function(x, y) 100 * x / y

# ===========================================================================
# A. SOURCES (contributors) — REVISION_PLAN §0 (1), §1 §2.1
# ===========================================================================
n_contrib_all  <- n_distinct(d$Main_researcher)
n_contrib_wing <- n_distinct(w$Main_researcher)
chk("sources", "contributors (wing records)", 42, n_contrib_wing,
    definition = "distinct Main_researcher among the 8,478 wing records (the plan's 'our sample'); all 12,571 records have 47")
chk("sources", "Main_researcher NA count", 0, sum(is.na(d$Main_researcher)))

contrib_year <- d %>%
  group_by(Year) %>%
  summarise(n_records = n(), n_wing = sum(has_wing),
            n_contributors = n_distinct(Main_researcher),
            n_contributors_wing = n_distinct(Main_researcher[has_wing]),
            n_municipalities = n_distinct(Municipality),
            n_species = n_distinct(Binomial), .groups = "drop")

contrib_tab <- d %>%
  group_by(Main_researcher) %>%
  summarise(
    n_records = n(), n_wing = sum(has_wing), n_mass = sum(!is.na(Body_mass.g.)),
    yr_min = min(Year), yr_max = max(Year), span = yr_max - yr_min, n_years = n_distinct(Year),
    in_early = any(Year <= EARLY_MAX), in_late = any(Year >= LATE_MIN),
    in_early_wing = any(Year <= EARLY_MAX & has_wing), in_late_wing = any(Year >= LATE_MIN & has_wing),
    n_species = n_distinct(Binomial), n_municipalities = n_distinct(Municipality),
    n_localities = n_distinct(Locality[!is.na(Locality)]),
    prop_right   = mean(wing_col[has_wing] == "right"),
    prop_left    = mean(wing_col[has_wing] == "left"),
    prop_generic = mean(wing_col[has_wing] == "generic"),
    mean_wing = mean(conc.wing.length, na.rm = TRUE),
    prop_recapture_yes = mean(Recapture == "Yes", na.rm = TRUE),
    prop_ring_present  = mean(!is.na(Ring)),
    .groups = "drop") %>%
  mutate(spans_both = in_early & in_late, spans_both_wing = in_early_wing & in_late_wing) %>%
  arrange(desc(n_records))

wing_contribs <- contrib_tab %>% filter(n_wing > 0)
n_span_wing   <- sum(wing_contribs$spans_both_wing)
rec_span_wing <- sum(wing_contribs$n_wing[wing_contribs$spans_both_wing])
chk("sources", "contributors spanning early & late (wing)", 6, n_span_wing)
chk("sources", "wing records from spanning contributors", 2828, rec_span_wing)
chk("sources", "wing records total", 8478, nrow(w))

# Wing-column (protocol proxy) mix by period, on wing records
wingcol_period <- wq %>%
  count(period, wing_col) %>% group_by(period) %>%
  mutate(prop = n / sum(n), n_period = sum(n)) %>% ungroup()
p_right_early   <- with(wingcol_period, prop[period == "early" & wing_col == "right"])
p_generic_late  <- with(wingcol_period, prop[period == "late"  & wing_col == "generic"])
chk("sources", "% right-wing column, early", 87, 100 * p_right_early, tol = 0.6)
chk("sources", "% generic wing column, late", 67, 100 * p_generic_late, tol = 0.6)
wingcol_year <- w %>% count(Year, wing_col) %>% group_by(Year) %>%
  mutate(prop = n / sum(n), n_year = sum(n)) %>% ungroup()

# Records per contributor per year (wing records) for the stacked figure
contrib_year_wing <- w %>% count(Year, Main_researcher, name = "n_wing")

# Long-running contributors (>= 8 distinct years, >= 200 wing records)
long_running <- wing_contribs %>% filter(n_years >= 8, n_wing >= 200)
chk("sources", "long-running contributors (>=8 yr, >=200 rec)", 6, nrow(long_running))
chk("sources", "wing records of long-running contributors", 3017, sum(long_running$n_wing))

# Contributors per species x period cell (mechanical heterogeneity, §0 variance).
# Computed under four definitions because the plan's 4.5 -> 6.6 was not
# reproduced exactly by any; closest is all records, species present in both
# quartile periods (the lnCVR cell set): 4.66 -> 6.65.
cells_fun <- function(x) x %>% group_by(Binomial, period) %>%
  summarise(n_contrib = n_distinct(Main_researcher), n = n(), .groups = "drop")
spp_both_fun <- function(x) x %>% group_by(Binomial) %>%
  summarise(e = any(period == "early"), l = any(period == "late"), .groups = "drop") %>%
  filter(e & l) %>% pull(Binomial)
dq <- d %>% filter(period != "mid")
cell_variants <- list(
  wing_all_species      = cells_fun(wq),
  wing_species_in_both  = cells_fun(wq %>% filter(Binomial %in% spp_both_fun(wq))),
  allrec_all_species    = cells_fun(dq),
  allrec_species_in_both = cells_fun(dq %>% filter(Binomial %in% spp_both_fun(dq))))
mean_contrib_cell <- bind_rows(lapply(names(cell_variants), function(nm)
  cell_variants[[nm]] %>% group_by(period) %>%
    summarise(mean_contributors = mean(n_contrib), median_contributors = median(n_contrib),
              n_cells = n(), .groups = "drop") %>% mutate(definition = nm)))
cell_contrib <- cell_variants$wing_all_species
mcc <- mean_contrib_cell %>% filter(definition == "allrec_species_in_both")
chk("sources", "mean contributors per species-period cell, early", 4.5,
    mcc$mean_contributors[mcc$period == "early"], tol = 0.2,
    definition = "all known-sex records, species present in both quartile periods (closest of 4 definitions; wing-only gives 3.7)")
chk("sources", "mean contributors per species-period cell, late", 6.6,
    mcc$mean_contributors[mcc$period == "late"], tol = 0.2,
    definition = "all known-sex records, species present in both quartile periods (wing-only gives 5.5)")

# Species attrition 73 -> in both quartile periods -> with >= 2 records per period
spp_period <- wq %>% group_by(Binomial) %>%
  summarise(n_early = sum(period == "early"), n_late = sum(period == "late"), .groups = "drop")
n_spp_both       <- sum(spp_period$n_early > 0 & spp_period$n_late > 0)
n_spp_both_ge2   <- sum(spp_period$n_early >= 2 & spp_period$n_late >= 2)
chk("sources", "species with wing records in both periods", 61, n_spp_both)

# Unknown-sex replication sample (M6): wing records, same 73 species
unk <- da %>% filter(!known_sex, has_wing)
unk_labelled <- unk %>% filter(!is.na(Sex))          # Sex == "Unknown" label only (NA-Sex excluded)
unk_complete <- unk %>% filter(!is.na(Main_researcher), !is.na(Municipality), !is.na(scaled_lat))
chk("sources", "unknown-sex wing records (M6 sample)", 3436, nrow(unk_labelled),
    definition = "wing records with Sex == 'Unknown'; adding the 96 wing records with Sex NA gives 3,532")

audit_sources <- list(
  generated = Sys.time(),
  definitions = c(early = sprintf("Year <= %d", EARLY_MAX), late = sprintf("Year >= %d", LATE_MIN),
                  wing_records = "non-NA conc.wing.length (coalesced right > left > generic)"),
  n_contributors_all = n_contrib_all, n_contributors_wing = n_contrib_wing,
  n_contributors_mass = n_distinct(d$Main_researcher[!is.na(d$Body_mass.g.)]),
  main_researcher_na = sum(is.na(d$Main_researcher)),
  n_span_both_wing = n_span_wing, n_wing_contributors = nrow(wing_contribs),
  wing_records_from_spanners = rec_span_wing, wing_records_total = nrow(w),
  spanning_contributors = wing_contribs$Main_researcher[wing_contribs$spans_both_wing],
  contributors_per_year = contrib_year,
  contributor_table = contrib_tab,
  records_per_contributor_year_wing = contrib_year_wing,
  wingcol_by_period = wingcol_period, wingcol_by_year = wingcol_year,
  prop_right_early = p_right_early, prop_generic_late = p_generic_late,
  long_running = long_running, long_running_n_wing = sum(long_running$n_wing),
  contributors_per_species_period_cell = cell_contrib,
  mean_contributors_per_cell = mean_contrib_cell,
  species_by_period = spp_period, n_spp_total = 73L, n_spp_both_periods = n_spp_both,
  n_spp_both_periods_ge2 = n_spp_both_ge2,
  unknown_sex = list(n_records_allsex = nrow(da), n_unknown_sex_records = sum(!da$known_sex),
                     n_unknown_sex_wing = nrow(unk), n_unknown_sex_wing_labelled = nrow(unk_labelled),
                     n_unknown_sex_wing_complete_covariates = nrow(unk_complete),
                     n_sex_na_raw = sum(is.na(da$Sex)),
                     n_sex_unknown_label = sum(da$Sex == "Unknown", na.rm = TRUE),
                     contributors = n_distinct(unk$Main_researcher),
                     species = n_distinct(unk$Binomial))
)

# ===========================================================================
# B. SITES — REVISION_PLAN §0 (2), §1 §2.2
# ===========================================================================
shared_count <- function(df, key) {
  key <- rlang::ensym(key)
  s <- df %>% filter(!is.na(!!key)) %>%
    group_by(!!key) %>%
    summarise(n = n(), early = any(period == "early"), late = any(period == "late"), .groups = "drop")
  list(n_units = nrow(s), n_shared = sum(s$early & s$late),
       n_records = sum(s$n), records_shared = sum(s$n[s$early & s$late]),
       n_units_early = sum(s$early), n_units_late = sum(s$late),
       n_na_records = sum(is.na(df[[rlang::as_string(key)]])))
}
loc_w <- shared_count(wq, Locality)
mun_w <- shared_count(wq, Municipality)
wq_sm <- wq %>% mutate(spp_mun = paste(Binomial, Municipality, sep = "|"))
sm_w  <- shared_count(wq_sm, spp_mun)
wq_xy <- wq %>% mutate(site_xy = paste(round(Longitude_decimal_degrees, 4),
                                       round(Latitude_decimal_degrees, 4), sep = "|"))
xy_w  <- shared_count(wq_xy, site_xy)
# The plan's 140 / 7 / 795 treated NA Locality (111 records) as one locality
# present in both periods. Named localities only: 139 / 6 / 684. Recomputed
# both ways so the source of the discrepancy is on record; the plan values
# are NOT reproduced and must be corrected in the manuscript.
loc_w_na_as_level <- {
  s <- wq %>% mutate(Locality = if_else(is.na(Locality), "<NA>", Locality)) %>% group_by(Locality) %>%
    summarise(n = n(), e = any(period == "early"), l = any(period == "late"), .groups = "drop")
  list(n_units = nrow(s), n_shared = sum(s$e & s$l), records_shared = sum(s$n[s$e & s$l]))
}
chk("sites", "wing records in quartile periods", 5240, nrow(wq))
chk("sites", "named localities (quartile wing records)", 140, loc_w$n_units,
    definition = "non-NA Locality; plan's 140 counted the NA level (139 + 1)")
chk("sites", "localities in both periods", 7, loc_w$n_shared,
    definition = "non-NA Locality; plan's 7 included the NA level (6 + 1)")
chk("sites", "records at shared localities", 795, loc_w$records_shared,
    definition = "non-NA Locality; plan's 795 included the 111 NA-locality records (684 + 111)")
chk("sites", "municipalities (quartile wing records)", 99, mun_w$n_units)
chk("sites", "municipalities in both periods", 15, mun_w$n_shared)
chk("sites", "records in shared municipalities", 1728, mun_w$records_shared)
chk("sites", "species x municipality combinations", 821, sm_w$n_units)
chk("sites", "species x municipality in both periods", 61, sm_w$n_shared)

# Same counts on ALL known-sex records (not only wing) for completeness
dq <- d %>% filter(period != "mid")
loc_all <- shared_count(dq, Locality); mun_all <- shared_count(dq, Municipality)

# Geography by period (wing records)
geo_period <- wq %>% group_by(period) %>%
  summarise(n = n(), n_sites_xy = n_distinct(paste(Longitude_decimal_degrees, Latitude_decimal_degrees)),
            lon_mean = mean(Longitude_decimal_degrees, na.rm = TRUE), lon_sd = sd(Longitude_decimal_degrees, na.rm = TRUE),
            lat_mean = mean(Latitude_decimal_degrees, na.rm = TRUE),  lat_sd = sd(Latitude_decimal_degrees, na.rm = TRUE),
            alt_mean = mean(Altitude, na.rm = TRUE), alt_median = median(Altitude, na.rm = TRUE),
            alt_sd = sd(Altitude, na.rm = TRUE), alt_na = sum(is.na(Altitude)),
            .groups = "drop")
# Correlation of year with geography (is "later" also "elsewhere"?)
geo_cor <- with(w, c(year_lon = cor(Year, Longitude_decimal_degrees, use = "complete.obs"),
                     year_lat = cor(Year, Latitude_decimal_degrees,  use = "complete.obs"),
                     year_alt = cor(Year, Altitude, use = "complete.obs")))
# Site x period data for the map (coordinate sites, wing records)
map_sites <- wq %>%
  group_by(Longitude_decimal_degrees, Latitude_decimal_degrees) %>%
  summarise(n_early = sum(period == "early"), n_late = sum(period == "late"),
            n_contrib = n_distinct(Main_researcher), .groups = "drop") %>%
  mutate(shared = n_early > 0 & n_late > 0)

audit_sites <- list(
  generated = Sys.time(), sample = "wing records in quartile periods unless suffixed _all",
  n_wing_quartile = nrow(wq),
  locality_wing = loc_w, locality_wing_NA_as_level = loc_w_na_as_level,
  municipality_wing = mun_w, species_municipality_wing = sm_w,
  coordinate_site_wing = xy_w,
  locality_all = loc_all, municipality_all = mun_all,
  n_localities_full_sample = n_distinct(d$Locality[!is.na(d$Locality)]),
  locality_na_full_sample = sum(is.na(d$Locality)),
  n_municipalities_full_sample = n_distinct(d$Municipality),
  n_coordinate_sites_full_sample = n_distinct(paste(d$Longitude_decimal_degrees, d$Latitude_decimal_degrees)),
  geography_by_period_wing = geo_period, cor_year_geography_wing = geo_cor,
  map_sites = map_sites
)

# ===========================================================================
# C. INDIVIDUALS — §2.1 recaptures
# ===========================================================================
ring_present <- mean(!is.na(d$Ring))
# An individual = Ring x Binomial: 36 ring strings occur on more than one
# species (ring-series re-use across contributors or data entry error), so the
# raw Ring string over-counts repeats (955 / 1,857 vs 927 / 1,813).
ring_tab     <- d %>% filter(!is.na(Ring)) %>% count(Ring, name = "n_rec")
ring_spp_tab <- d %>% filter(!is.na(Ring)) %>% count(Ring, Binomial, name = "n_rec")
rings_multi     <- sum(ring_spp_tab$n_rec > 1)
extra_records   <- sum(ring_spp_tab$n_rec[ring_spp_tab$n_rec > 1] - 1)
rings_multi_raw <- sum(ring_tab$n_rec > 1)
extra_raw       <- sum(ring_tab$n_rec[ring_tab$n_rec > 1] - 1)
rings_multi_spp <- sum((d %>% filter(!is.na(Ring)) %>% distinct(Ring, Binomial) %>% count(Ring))$n > 1)
chk("individuals", "% records with Ring", 91, 100 * ring_present, tol = 0.6)
chk("individuals", "rings with > 1 record", 927, rings_multi,
    definition = "Ring x Binomial (individual); raw Ring string gives 955")
chk("individuals", "extra (repeat) records", 1813, extra_records,
    definition = "Ring x Binomial; raw Ring string gives 1,857")
recap_fun <- function(x) x %>% filter(period != "mid") %>% group_by(period) %>%
  summarise(n = n(), recapture_yes = sum(Recapture == "Yes", na.rm = TRUE),
            recapture_na = sum(is.na(Recapture)),
            share_yes_of_all = recapture_yes / n,
            share_yes_of_recorded = mean(Recapture == "Yes", na.rm = TRUE), .groups = "drop")
recap_period      <- recap_fun(d)
recap_period_wing <- recap_fun(w)
chk("individuals", "% Recapture == Yes, early", 5, 100 * recap_period_wing$share_yes_of_all[recap_period_wing$period == "early"], tol = 0.6,
    definition = "wing records, Yes / all records (NA counted as not recapture); of recorded values 6.8 %; all records 7.7 % / 9.4 %")
chk("individuals", "% Recapture == Yes, late", 12, 100 * recap_period_wing$share_yes_of_all[recap_period_wing$period == "late"], tol = 0.6,
    definition = "wing records, Yes / all records; of recorded values 14.8 %; all records 20.9 % / 26.0 %")
# first-capture flag: earliest record per individual (Ring x Binomial) in the sample
d_first <- d %>% arrange(Year, month, day) %>% group_by(Ring, Binomial) %>%
  mutate(ring_seq = if_else(is.na(Ring), 1L, row_number())) %>% ungroup()
first_capture_n      <- sum(d_first$ring_seq == 1)
first_capture_n_wing <- sum(d_first$ring_seq == 1 & d_first$has_wing)
# Ring-based repeats across the sample vs the Recapture field (they measure different things:
# Recapture = previously ringed at capture, whether or not the earlier record is in ABT)
audit_individuals <- list(
  generated = Sys.time(), prop_ring_present = ring_present,
  individual_definition = "Ring x Binomial",
  n_rings_raw = nrow(ring_tab), n_individuals = nrow(ring_spp_tab),
  rings_on_more_than_one_species = rings_multi_spp,
  individuals_with_multiple_records = rings_multi, extra_records_from_repeats = extra_records,
  rings_with_multiple_records_raw = rings_multi_raw, extra_records_raw = extra_raw,
  records_per_individual_distribution = as.data.frame(table(n_rec = ring_spp_tab$n_rec)),
  recapture_by_period_all = recap_period, recapture_by_period_wing = recap_period_wing,
  recapture_overall = mean(d$Recapture == "Yes", na.rm = TRUE), recapture_na = sum(is.na(d$Recapture)),
  first_capture_only_n = first_capture_n, first_capture_only_n_wing = first_capture_n_wing,
  note = paste("first_capture_only keeps the earliest record per Ring x Binomial (records without Ring kept);",
               "Recapture field is NA for", sum(is.na(d$Recapture)), "records")
)

# ===========================================================================
# D. SEASON / DATE / MOULT / SEX RATIO — §2.1, §3.6
# ===========================================================================
date_na   <- sum(is.na(d$Date)); month_na <- sum(is.na(d$month))
date_nonstd <- unique(d$Date[!is.na(d$Date) & !grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$", d$Date)])
date_year_mismatch <- sum(!is.na(d$date_year) & d$date_year != d$Year)
chk("season", "records with Date missing", 6, date_na)
month_period <- d %>% filter(period != "mid", !is.na(month)) %>% count(period, month) %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
season_period <- d %>% filter(period != "mid") %>% count(period, season) %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
season_wing_period <- wq %>% count(period, season) %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
molt_period <- d %>% filter(period != "mid") %>% count(period, Molt) %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
molt_complete <- mean(!is.na(d$Molt))
chk("season", "% Molt recorded", 73, 100 * molt_complete, tol = 0.6)
repro_period <- d %>% filter(period != "mid") %>% count(period, Reproductive_stage) %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
molt_month <- d %>% filter(!is.na(month), !is.na(Molt)) %>% count(month, Molt) %>%
  group_by(month) %>% mutate(prop = n / sum(n)) %>% ungroup()
sex_period <- d %>% filter(period != "mid") %>% count(period, Sex) %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
sex_period_wing <- wq %>% count(period, Sex) %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
chk("season", "% Female, early", 46, 100 * sex_period_wing$prop[sex_period_wing$period == "early" & sex_period_wing$Sex == "Female"], tol = 0.6,
    definition = "wing records (all records: 47.8 %)")
chk("season", "% Female, late",  49, 100 * sex_period_wing$prop[sex_period_wing$period == "late"  & sex_period_wing$Sex == "Female"], tol = 0.6,
    definition = "wing records (all records: 48.5 %)")
# Crude season effect on wing (descriptive only: species-centred means)
season_wing_desc <- w %>% filter(!is.na(season)) %>% group_by(Binomial, Sex) %>%
  mutate(cwl_c = conc.wing.length - mean(conc.wing.length)) %>% ungroup() %>%
  group_by(season) %>% summarise(n = n(), mean_centred_wing = mean(cwl_c), se = sd(cwl_c) / sqrt(n()), .groups = "drop")
molt_wing_desc <- w %>% filter(!is.na(Molt)) %>% group_by(Binomial, Sex) %>%
  mutate(cwl_c = conc.wing.length - mean(conc.wing.length)) %>% ungroup() %>%
  group_by(Molt) %>% summarise(n = n(), mean_centred_wing = mean(cwl_c), se = sd(cwl_c) / sqrt(n()), .groups = "drop")
audit_season <- list(
  generated = Sys.time(),
  date_missing = date_na, month_missing = month_na, date_nonstandard_values = date_nonstd,
  date_year_vs_Year_mismatch = date_year_mismatch,
  hour_present = sum(!is.na(d$Hour)), hour_parsed = sum(!is.na(d$hour_num)),
  hour_unparsed_values = unique(d$Hour[!is.na(d$Hour) & is.na(d$hour_num)]),
  month_by_period = month_period, season_by_period = season_period,
  season_by_period_wing = season_wing_period,
  molt_by_period = molt_period, molt_completeness = molt_complete, molt_by_month = molt_month,
  reproductive_stage_by_period = repro_period,
  sex_by_period = sex_period, sex_by_period_wing = sex_period_wing,
  season_wing_descriptive = season_wing_desc, molt_wing_descriptive = molt_wing_desc,
  note = "descriptive species x sex-centred means only; the modelled season / moult effects belong to Phase 1"
)

# ===========================================================================
# E. TRAITS — completeness and provenance for 6 traits + Hour
# ===========================================================================
trait_cols <- c(wing = "conc.wing.length", body_mass = "Body_mass.g.", bill_width = "Bill_width.mm.",
                bill_length = "Bill_length.mm.", tail_length = "Tail_length.mm.",
                tarsus_length = "Tarsus_length.mm.", hour = "hour_num")
trait_tab <- bind_rows(lapply(names(trait_cols), function(nm) {
  v <- trait_cols[[nm]]; x <- d %>% filter(!is.na(.data[[v]]))
  ct <- x %>% group_by(Main_researcher) %>%
    summarise(e = any(Year <= EARLY_MAX), l = any(Year >= LATE_MIN), .groups = "drop")
  tibble(trait = nm, column = v, n = nrow(x), pct_of_sample = pct(nrow(x), nrow(d)),
         n_species = n_distinct(x$Binomial),
         n_contributors = n_distinct(x$Main_researcher), n_municipalities = n_distinct(x$Municipality),
         contributors_spanning_both = sum(ct$e & ct$l),
         n_early = sum(x$Year <= EARLY_MAX), n_late = sum(x$Year >= LATE_MIN),
         mean = mean(x[[v]]), sd = sd(x[[v]]),
         yr_min = min(x$Year), yr_max = max(x$Year))
}))
chk("traits", "body mass n", 11256, trait_tab$n[trait_tab$trait == "body_mass"])
chk("traits", "body mass contributors", 47, trait_tab$n_contributors[trait_tab$trait == "body_mass"])
chk("traits", "body mass municipalities", 150, trait_tab$n_municipalities[trait_tab$trait == "body_mass"])
chk("traits", "bill width n", 3209, trait_tab$n[trait_tab$trait == "bill_width"])
chk("traits", "bill width contributors", 22, trait_tab$n_contributors[trait_tab$trait == "bill_width"])
chk("traits", "bill width municipalities", 84, trait_tab$n_municipalities[trait_tab$trait == "bill_width"])
chk("traits", "bill width contributors spanning both", 2, trait_tab$contributors_spanning_both[trait_tab$trait == "bill_width"])
chk("traits", "bill length n", 7697, trait_tab$n[trait_tab$trait == "bill_length"])
chk("traits", "tail length n", 8872, trait_tab$n[trait_tab$trait == "tail_length"])
chk("traits", "tarsus length n", 4361, trait_tab$n[trait_tab$trait == "tarsus_length"])
chk("traits", "% records with Hour", 76.5, pct(sum(!is.na(d$Hour)), nrow(d)), tol = 0.6)
# Bill-width column provenance (98 % in the unspecified column, per plan)
bw_cols <- d %>% summarise(unspecified = sum(!is.na(Bill_width.mm.)))
raw_bw <- tryCatch({
  rb <- readr::read_csv(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"),
                        guess_max = 70000, show_col_types = FALSE,
                        col_select = c("ID_ABT", "Bill_width.mm.", "Bill_width_base.mm.", "Bill_width_nostril.mm."))
  rb <- rb[rb$ID_ABT %in% d$ID_ABT, ]
  c(unspecified = sum(!is.na(rb$Bill_width.mm.)), base = sum(!is.na(rb$Bill_width_base.mm.)),
    nostril = sum(!is.na(rb$Bill_width_nostril.mm.)))
}, error = function(e) c(unspecified = NA, base = NA, nostril = NA))
# Shared wing + mass records (isometry test sample)
n_shared_wing_mass <- sum(!is.na(d$conc.wing.length) & !is.na(d$Body_mass.g.))
chk("traits", "records with both wing and mass", 7577, n_shared_wing_mass)
# Trait x period x contributor table (for Supplementary)
trait_contrib_period <- bind_rows(lapply(names(trait_cols)[1:6], function(nm) {
  v <- trait_cols[[nm]]
  d %>% filter(!is.na(.data[[v]]), period != "mid") %>%
    group_by(Main_researcher, period) %>% summarise(n = n(), .groups = "drop") %>%
    mutate(trait = nm)
}))
# Contributor-level mean wing vs period (the between-contributor contrast)
contrib_wing_means <- w %>% group_by(Binomial, Sex) %>%
  mutate(cwl_c = conc.wing.length - mean(conc.wing.length)) %>% ungroup() %>%
  group_by(Main_researcher) %>%
  summarise(n = n(), mean_year = mean(Year), mean_centred_wing = mean(cwl_c), .groups = "drop")
audit_traits <- list(
  generated = Sys.time(), n_sample = nrow(d), trait_table = trait_tab,
  bill_width_columns_raw = raw_bw, n_shared_wing_mass = n_shared_wing_mass,
  n_shared_wing_mass_contributors = n_distinct(d$Main_researcher[!is.na(d$conc.wing.length) & !is.na(d$Body_mass.g.)]),
  trait_contributor_period = trait_contrib_period,
  contributor_centred_wing_means = contrib_wing_means,
  wing_columns = as.data.frame(table(wing_col = d$wing_col, useNA = "ifany")),
  note = "n counts are non-NA records in the 12,571-record known-sex sample; provenance = distinct Main_researcher / Municipality among those records"
)

# ===========================================================================
# F. EFFECT SCALE — §2.4
# ===========================================================================
sc_attr <- attributes(passer90$scaled_yr)
sd_year_scaling <- as.numeric(sc_attr[["scaled:scale"]])     # 5.0199, the SD the coefficient is actually per
mu_year_scaling <- as.numeric(sc_attr[["scaled:center"]])
sd_year_sample  <- sd(d$Year); sd_year_wing <- sd(w$Year)
mean_wing <- mean(w$conc.wing.length); mean_mass <- mean(d$Body_mass.g., na.rm = TRUE)
yr_range  <- range(d$Year); duration <- diff(yr_range)                  # 23 yr difference, 24 calendar years
chk("effect", "SD(year), analytical sample", 5.05, sd_year_sample, tol = 0.006)
chk("effect", "mean wing (mm)", 71.1, mean_wing, tol = 0.06)

# Published estimates: read from the pooled-summary files when present
ds <- if (file.exists(out_path("descriptive_summary.rds"))) readRDS(out_path("descriptive_summary.rds")) else NULL
mr <- if (file.exists(out_path("mass_results.rds"))) readRDS(out_path("mass_results.rds")) else NULL
wing_pub <- c(est = -0.92, lo = -1.40, hi = -0.44)         # descriptive_summary.rds (rounded)
if (!is.null(ds)) wing_pub <- c(est = ds$wing_yr, lo = ds$wing_yr_lo, hi = ds$wing_yr_hi)
mass_pub <- c(est = -0.0073, lo = -0.0168, hi = 0.0022)    # mass_results.rds rubin b_scaled_yr
if (!is.null(mr)) { r <- mr$rubin[mr$rubin$par == "b_scaled_yr", ]; mass_pub <- c(est = r$estimate, lo = r$lower, hi = r$upper) }
bill_pub <- c(est = 0.24, lo = 0.06, hi = 0.43)
if (!is.null(ds)) bill_pub <- c(est = ds$bill_yr, lo = ds$bill_yr_lo, hi = ds$bill_yr_hi)

convert_wing <- function(b, sd_yr, ref_mean = mean_wing, yrs = duration) {
  per_yr <- b / sd_yr
  c(per_sd_year = b, mm_per_year = per_yr, mm_per_decade = 10 * per_yr,
    pct_per_decade = 100 * 10 * per_yr / ref_mean,
    mm_over_record = yrs * per_yr, pct_over_record = 100 * yrs * per_yr / ref_mean)
}
convert_logmass <- function(b, sd_yr, yrs = duration) {
  per_yr <- b / sd_yr
  c(per_sd_year = b, log_per_year = per_yr,
    pct_per_decade = 100 * (exp(10 * per_yr) - 1),
    pct_over_record = 100 * (exp(yrs * per_yr) - 1))
}
wing_scale <- list(
  using_scaling_sd = t(sapply(wing_pub, convert_wing, sd_yr = sd_year_scaling)),
  using_sample_sd  = t(sapply(wing_pub, convert_wing, sd_yr = sd_year_sample))
)
mass_scale <- list(
  using_scaling_sd = t(sapply(mass_pub, convert_logmass, sd_yr = sd_year_scaling)),
  using_sample_sd  = t(sapply(mass_pub, convert_logmass, sd_yr = sd_year_sample))
)
bill_scale <- list(using_scaling_sd = t(sapply(bill_pub, convert_wing, sd_yr = sd_year_scaling,
                                                ref_mean = mean(d$Bill_width.mm., na.rm = TRUE))))
# Isometry: mass proportional to L^3 -> expected % mass change for a % wing change
wing_pct_record <- wing_scale$using_scaling_sd["est", "pct_over_record"]
iso_expect_pct  <- 100 * ((1 + wing_pct_record / 100)^3 - 1)
chk("effect", "-0.92 -> mm per decade", -1.8, wing_scale$using_sample_sd["est", "mm_per_decade"], tol = 0.06)
chk("effect", "-0.92 -> % per decade", -2.6, wing_scale$using_sample_sd["est", "pct_per_decade"], tol = 0.06)
chk("effect", "-0.92 -> mm over record", -4.2, wing_scale$using_sample_sd["est", "mm_over_record"], tol = 0.06)
chk("effect", "-0.92 -> % over record", -5.9, wing_scale$using_sample_sd["est", "pct_over_record"], tol = 0.06)
chk("effect", "mass -0.0073 -> % per decade", -1.4, mass_scale$using_sample_sd["est", "pct_per_decade"], tol = 0.06)
chk("effect", "mass -0.0073 -> % over record", -3.3, mass_scale$using_sample_sd["est", "pct_over_record"], tol = 0.06)
chk("effect", "mass lower CI -> % over record", -7.6, mass_scale$using_sample_sd["lo", "pct_over_record"], tol = 0.06,
    definition = "exact 100*(exp(b*23/SD)-1); the plan's -7.6 is the linear approximation 100*b*23/SD = -7.67; quote -7.4 %")
chk("effect", "isometric mass expectation (%)", -17, iso_expect_pct, tol = 0.6)

# Comparator table: values NOT verified here -> NA placeholders. Fill from the
# papers before use; `verified` must be set TRUE by whoever checks them.
comparators <- data.frame(
  study = c("Jirinec et al. 2021 (Sci. Adv.)", "Weeks et al. 2020 (Ecol. Lett.)",
            "Ryding et al. 2024", "This study, published baseline", "This study, controlled (Phase 1, to fill)"),
  system = c("Amazonian understorey birds, 1979-2019, live captures",
             "North American migrants, Chicago window kills, 1978-2016, museum",
             "TO VERIFY (appendage / shape-shifting synthesis)",
             "Atlantic Forest passerines, live captures, 1995-2018",
             "Atlantic Forest passerines, live captures, 1995-2018"),
  wing_pct_per_decade = c(NA, NA, NA, wing_scale$using_scaling_sd["est", "pct_per_decade"], NA),
  mass_pct_per_decade = c(NA, NA, NA, mass_scale$using_scaling_sd["est", "pct_per_decade"], NA),
  verified = c(FALSE, FALSE, FALSE, TRUE, FALSE),
  note = c("PLACEHOLDER: fill mass and wing % per decade from the paper (not verifiable offline here)",
           "PLACEHOLDER: fill tarsus/mass/wing % change from the paper (not verifiable offline here)",
           "PLACEHOLDER: confirm citation year, taxon set and metric before quoting",
           "conversions use the scaling SD actually used for scaled_yr (5.020)",
           "fill from controlled_wing_results.rds (P1) once available"),
  stringsAsFactors = FALSE)

effect_scale <- list(
  generated = Sys.time(),
  year_range = yr_range, duration_years_difference = duration, duration_calendar_years = duration + 1,
  sd_year_scaling = sd_year_scaling, mean_year_scaling = mu_year_scaling,
  sd_year_sample = sd_year_sample, sd_year_wing_sample = sd_year_wing,
  scaling_reference = attr(passer90, "scaling"),
  mean_wing_mm = mean_wing, sd_wing_mm = sd(w$conc.wing.length), n_wing = nrow(w),
  mean_mass_g = mean_mass, n_mass = sum(!is.na(d$Body_mass.g.)),
  mean_bill_width_mm = mean(d$Bill_width.mm., na.rm = TRUE),
  published = list(wing = wing_pub, mass_log = mass_pub, bill_width = bill_pub,
                   source = c(wing = "descriptive_summary.rds (brm0_multiphylo Rubin pooled)",
                              mass = "mass_results.rds (brm_mass_multiphylo Rubin pooled)",
                              bill = "descriptive_summary.rds (brm_bill_multiphylo Rubin pooled)")),
  wing = wing_scale, mass = mass_scale, bill_width = bill_scale,
  isometry = list(wing_pct_over_record = wing_pct_record, expected_mass_pct_isometric = iso_expect_pct,
                  observed_mass_pct_over_record = mass_scale$using_scaling_sd["est", "pct_over_record"],
                  observed_mass_pct_lower = mass_scale$using_scaling_sd["lo", "pct_over_record"],
                  note = "isometry: mass ~ L^3, so expected mass change = (1 + wing change)^3 - 1"),
  comparators = comparators,
  note = paste("scaled_yr was standardised on the PRE-threshold known-sex sample (SD 5.020), not the",
               "final 12,571-record sample (SD 5.047); REVISION_PLAN quotes 5.05. Conversions with both are",
               "given; the difference is 0.5 % of the effect. Use `using_scaling_sd` for the coefficient.")
)

# ===========================================================================
# G. Attach plan checks, save
# ===========================================================================
audit_sources$plan_checks     <- plan_checks
audit_sites$plan_checks       <- plan_checks[plan_checks$section == "sites", ]
audit_individuals$plan_checks <- plan_checks[plan_checks$section == "individuals", ]
audit_season$plan_checks      <- plan_checks[plan_checks$section == "season", ]
audit_traits$plan_checks      <- plan_checks[plan_checks$section == "traits", ]
effect_scale$plan_checks      <- plan_checks[plan_checks$section == "effect", ]

saveRDS(audit_sources,     out_path("audit_sources.rds"))
saveRDS(audit_sites,       out_path("audit_sites.rds"))
saveRDS(audit_individuals, out_path("audit_individuals.rds"))
saveRDS(audit_season,      out_path("audit_season.rds"))
saveRDS(audit_traits,      out_path("audit_traits.rds"))
saveRDS(effect_scale,      out_path("effect_scale.rds"))

# ===========================================================================
# H. Figures
# ===========================================================================
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank()))

# H1. Records per year by contributor (wing records), top contributors named
top_k <- 12
top_names <- contrib_tab %>% filter(n_wing > 0) %>% arrange(desc(n_wing)) %>% slice_head(n = top_k) %>% pull(Main_researcher)
f1_dat <- contrib_year_wing %>%
  mutate(contributor = if_else(Main_researcher %in% top_names, Main_researcher,
                               sprintf("other (%d contributors)", n_contrib_wing - top_k))) %>%
  group_by(Year, contributor) %>% summarise(n_wing = sum(n_wing), .groups = "drop") %>%
  mutate(contributor = factor(contributor, levels = c(top_names, setdiff(unique(contributor), top_names))))
pal <- c(scales::hue_pal()(top_k), "grey70")
p1 <- ggplot(f1_dat, aes(Year, n_wing, fill = contributor)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = pal, name = "Contributor (Main_researcher)") +
  scale_x_continuous(breaks = seq(1995, 2018, 2)) +
  labs(x = "Year", y = "Wing-length records",
       title = "Records per year by contributor (known-sex live sample, wing records)",
       subtitle = sprintf("%d records, %d contributors; %d contributors span both %s and %s (%d records)",
                          nrow(w), n_contrib_wing, n_span_wing,
                          sprintf("<= %d", EARLY_MAX), sprintf(">= %d", LATE_MIN), rec_span_wing)) +
  theme(legend.position = "right", legend.text = element_text(size = 7), legend.key.size = unit(0.4, "cm"))
p1b <- ggplot(contrib_year, aes(Year, n_contributors_wing)) +
  geom_col(fill = "grey50", width = 0.85) +
  scale_x_continuous(breaks = seq(1995, 2018, 2)) +
  labs(x = NULL, y = "Contributors")
ggsave(fig_path("audit_records_per_year_by_source.png"), p1 / p1b + plot_layout(heights = c(3, 1)),
       width = 10, height = 6.5, dpi = 200, bg = "white")

# H2. Sites map: early vs late, shared sites distinct from period-unique
basemap <- NULL
gj <- file.path(dirname(ANALYSIS_DIR), "Manuscript", "data", "south_america.geojson")
if (file.exists(gj) && requireNamespace("sf", quietly = TRUE)) {
  basemap <- tryCatch(sf::st_read(gj, quiet = TRUE), error = function(e) NULL)
}
map_long <- bind_rows(
  map_sites %>% filter(n_early > 0) %>% transmute(Longitude_decimal_degrees, Latitude_decimal_degrees,
                                                    n = n_early, shared, period = sprintf("Early (1995–%d)", EARLY_MAX)),
  map_sites %>% filter(n_late > 0)  %>% transmute(Longitude_decimal_degrees, Latitude_decimal_degrees,
                                                    n = n_late, shared, period = sprintf("Late (%d–2018)", LATE_MIN))) %>%
  mutate(site_type = if_else(shared, "Sampled in both periods", "Period-unique site"))
panel_lab <- map_long %>% group_by(period) %>%
  summarise(lab = sprintf("N = %d records, %d sites (%d shared)", sum(n), n(), sum(shared)), .groups = "drop")
p2 <- ggplot()
if (!is.null(basemap)) p2 <- p2 + geom_sf(data = basemap, fill = "grey93", colour = "grey60", linewidth = 0.2)
p2 <- p2 +
  geom_point(data = map_long %>% filter(!shared),
             aes(Longitude_decimal_degrees, Latitude_decimal_degrees, size = n, colour = site_type, shape = site_type), alpha = 0.7) +
  geom_point(data = map_long %>% filter(shared),
             aes(Longitude_decimal_degrees, Latitude_decimal_degrees, size = n, colour = site_type, shape = site_type), alpha = 0.95, stroke = 0.8) +
  geom_text(data = panel_lab, aes(x = -59, y = -33, label = lab), hjust = 0, vjust = 0, size = 3) +
  scale_colour_manual(values = c("Sampled in both periods" = "#B2182B", "Period-unique site" = "#4393C3"), name = NULL) +
  scale_shape_manual(values = c("Sampled in both periods" = 17, "Period-unique site" = 16), name = NULL) +
  scale_size_area(max_size = 7, name = "Wing records") +
  coord_sf(xlim = c(-60, -32), ylim = c(-34, -3), expand = FALSE) +
  facet_wrap(~ period) +
  labs(x = NULL, y = NULL,
       title = "Sampling sites by period: shared vs period-unique (wing records, known-sex live sample)",
       subtitle = sprintf("Coordinate sites: %d early, %d late, %d in both; named localities in both periods: %d of %d; municipalities: %d of %d",
                          sum(map_sites$n_early > 0), sum(map_sites$n_late > 0), sum(map_sites$shared),
                          loc_w$n_shared, loc_w$n_units, mun_w$n_shared, mun_w$n_units)) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "#EAF2F8", colour = NA))
ggsave(fig_path("audit_sites_map_shared.png"), p2, width = 10, height = 6.5, dpi = 200, bg = "white")

# H3. Wing-column proxy by year
wc_lab <- wingcol_year %>% group_by(Year) %>% summarise(n_year = first(n_year), .groups = "drop")
p3 <- ggplot(wingcol_year, aes(Year, prop, fill = factor(wing_col, levels = c("right", "left", "generic")))) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  geom_text(data = wc_lab, aes(Year, 1.02, label = n_year), inherit.aes = FALSE, size = 2.4, angle = 90, hjust = 0) +
  scale_fill_manual(values = c(right = "#2C7BB6", left = "#ABD9E9", generic = "#FDAE61"),
                    name = "Column populated", labels = c("Wing_length_right", "Wing_length_left", "Wing_length (unspecified)")) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  scale_x_continuous(breaks = seq(1995, 2018, 2)) +
  labs(x = "Year", y = "Share of wing records",
       title = "Which wing column supplied the coalesced wing length, by year",
       subtitle = sprintf("Protocol proxy: %.0f %% right-wing column in %s vs %.0f %% unspecified column in %s (numbers above bars = records)",
                          100 * p_right_early, sprintf("<= %d", EARLY_MAX), 100 * p_generic_late, sprintf(">= %d", LATE_MIN))) +
  theme(legend.position = "bottom")
ggsave(fig_path("audit_wingcol_by_year.png"), p3, width = 9, height = 5, dpi = 200, bg = "white")

# ===========================================================================
# I. Console report
# ===========================================================================
cat("\n================ PLAN CHECKS (REVISION_PLAN.md reference vs reproduced) ================\n")
print(plan_checks %>% mutate(reproduced = round(reproduced, 3)) %>% select(-definition), row.names = FALSE)
cat(sprintf("\n%d of %d reference numbers reproduced within tolerance.\n", sum(plan_checks$match), nrow(plan_checks)))
if (any(!plan_checks$match)) {
  cat("\nNOT reproduced (see `definition`):\n")
  print(plan_checks[!plan_checks$match, c("item", "plan", "reproduced", "definition")], row.names = FALSE)
}
cat("\nContributors per species x period cell (4 definitions):\n"); print(as.data.frame(mean_contrib_cell))
cat("\nTraits:\n"); print(as.data.frame(trait_tab %>% select(trait, n, pct_of_sample, n_contributors, n_municipalities, contributors_spanning_both, n_early, n_late)), digits = 4)
cat("\nBill width raw columns:\n"); print(raw_bw)
cat("\nEffect scale (wing, scaling SD):\n"); print(round(wing_scale$using_scaling_sd, 3))
cat("\nEffect scale (mass, scaling SD):\n"); print(round(mass_scale$using_scaling_sd, 3))
cat(sprintf("\nIsometric expectation for %.2f %% wing change: %.1f %% mass; observed %.1f %% [lower %.1f %%]\n",
            wing_pct_record, iso_expect_pct, effect_scale$isometry$observed_mass_pct_over_record,
            effect_scale$isometry$observed_mass_pct_lower))
cat(sprintf("\nUnknown-sex wing records (M6): %d (complete covariates %d); NA-Sex raw: %d; 'Unknown' label: %d\n",
            nrow(unk), nrow(unk_complete), audit_sources$unknown_sex$n_sex_na_raw, audit_sources$unknown_sex$n_sex_unknown_label))
cat(sprintf("Date: %d NA; non-standard values: %s; date_year != Year: %d\n", date_na,
            paste(date_nonstd, collapse = ", "), date_year_mismatch))
cat(sprintf("Geography x year correlations (wing records): lon %.2f, lat %.2f, alt %.2f\n", geo_cor[1], geo_cor[2], geo_cor[3]))
cat("\nWrote audit_*.rds, effect_scale.rds to", out_path(), "and 3 figures to", fig_path(), "\n")
