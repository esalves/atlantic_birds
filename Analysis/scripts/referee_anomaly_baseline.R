# referee_anomaly_baseline.R
# ---------------------------------------------------------------------------
# WHY (internal review before the Oikos submission, 2026-09-30): the local
# temperature anomaly dT is the record-year locality temperature minus the
# locality's mean over its own records (referee_reruns.R), i.e. over the years
# in which it was sampled, weighted by captures. The review asked whether that
# baseline depends on when birds were caught. This script
#   (1) compares it with a baseline that does not: the locality's unweighted mean
#       over every year 1995-2018 of the climate series (locality_year_tmean.rds),
#       and reports how each anomaly correlates with calendar year;
#   (2) refits the thermal models of block E3 (full M3 adjustment, glmmTMB REML,
#       phylogenetic species intercept, Rubin-pooled over the trees) with each
#       baseline, alone and with calendar year;
#   (3) describes the shared records that have no climate value (WorldClim land
#       mask), which the review suspected are not a random subset.
#
# RUN:  Rscript Analysis/scripts/referee_anomaly_baseline.R        (REFEREE_N_TREES, default 50)
# OUT:  Analysis/output/referee_reruns/anomaly_baseline.rds
# ---------------------------------------------------------------------------
local({
  f <- c("Analysis/scripts/referee_reruns.R", "scripts/referee_reruns.R", "referee_reruns.R")
  f <- f[file.exists(f)][1]
  if (is.na(f)) stop("Run from the repo root, Analysis/ or Analysis/scripts/.")
  ex <- parse(f)
  stop_at <- which(vapply(ex, function(e) grepl('note("E0', paste(deparse(e), collapse = " "), fixed = TRUE), TRUE))[1]
  for (e in ex[seq_len(stop_at - 1)]) eval(e, envir = globalenv())
})

# --- (1) the two baselines -----------------------------------------------------------
ly <- readRDS(derived_path("locality_year_tmean.rds"))
base_full <- ly %>% filter(year >= 1995, year <= 2018) %>% group_by(loc_id) %>%
  summarise(t_base_full = mean(tmean, na.rm = TRUE), n_yr_full = sum(!is.na(tmean)), .groups = "drop")
d3 <- shared %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season)) %>%
  left_join(base_full, by = "loc_id") %>%
  mutate(dT_full = rec_tmean - t_base_full)
stopifnot(!anyNA(d3$dT_full))
baselines <- data.frame(
  anomaly = c("dT (mean over the locality's records)", "dT_full (locality mean over 1995-2018)"),
  sd = c(sd(d3$dT), sd(d3$dT_full)),
  cor_year = c(cor(d3$dT, d3$Year), cor(d3$dT_full, d3$Year)),
  cor_with_dT = c(1, cor(d3$dT, d3$dT_full)))
note(sprintf("M3 thermal sample: %d records; full-series baseline over %s years", nrow(d3),
             paste(unique(range(base_full$n_yr_full)), collapse = "-")))
print(baselines, digits = 3)

# --- (2) E3 models with each baseline -------------------------------------------------
A3 <- phylo_A_list(sort(unique(d3$species_name)), n_trees = N_TREES)
rhs_adj <- "Sex + %s + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + %s || spp) + (1 | src) + (1 | site) + (1 | ind)"
specs <- list(
  list(label = "dT",                 x = "dT",                  slope = "dT",        terms = "dT"),
  list(label = "dT + year",          x = "dT + scaled_yr",      slope = "scaled_yr", terms = c("dT", "scaled_yr")),
  list(label = "dT_full",            x = "dT_full",             slope = "dT_full",   terms = "dT_full"),
  list(label = "dT_full + year",     x = "dT_full + scaled_yr", slope = "scaled_yr", terms = c("dT_full", "scaled_yr")))
responses <- list(c("wing", "wing length (mm)"), c("lnmass", "log body mass"), c("iso", "isometry contrast lnM-3lnL"))
rows <- list()
for (x in responses) for (s in specs) {
  f <- as.formula(sprintf(paste("%s ~", rhs_adj), x[1], s$x, s$slope))
  r <- safely(paste(x[1], s$label), run_phylo_trees(f, d3, A3, verbose = FALSE))
  if (is.null(r)) next
  rows[[length(rows) + 1]] <- pool_clean(r, s$terms) %>%
    transmute(response = x[2], model = s$label, term = par, estimate, se, lower, upper, p, trees_used, n = nrow(d3))
  note("  done: ", x[1], " ~ ", s$label)
}
fits <- do.call(rbind, rows)
print(fits, digits = 3)

# --- (3) shared records without a climate value ---------------------------------------
all_shared <- dat %>% filter(!is.na(wing), !is.na(lnmass))
missing_profile <- all_shared %>% mutate(no_climate = is.na(rec_tmean)) %>% group_by(no_climate) %>%
  summarise(records = n(), localities = n_distinct(loc_id), contributors = n_distinct(src), species = n_distinct(Binomial),
            mean_lat = mean(Latitude_decimal_degrees), mean_lon = mean(Longitude_decimal_degrees),
            median_alt_m = median(Altitude, na.rm = TRUE), mean_year = mean(Year), .groups = "drop")
missing_by_src <- all_shared %>% filter(is.na(rec_tmean)) %>% count(src, sort = TRUE)
print(missing_profile, width = 200); print(missing_by_src)

saveRDS(list(generated = Sys.time(), n_trees = N_TREES, n = nrow(d3), baselines = baselines, fits = fits,
             missing_profile = as.data.frame(missing_profile), missing_by_contributor = as.data.frame(missing_by_src),
             note = paste("dT: rec_tmean minus the locality's mean over its own records (referee_reruns.R);",
                          "dT_full: rec_tmean minus the locality's unweighted mean over 1995-2018 of the series.")),
        out_path("referee_reruns", "anomaly_baseline.rds"))
note("wrote ", out_path("referee_reruns", "anomaly_baseline.rds"))
