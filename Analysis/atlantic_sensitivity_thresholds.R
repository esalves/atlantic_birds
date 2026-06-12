# atlantic_sensitivity_thresholds.R
# ---------------------------------------------------------------------------
# Sensitivity analysis: does the year → wing-length result hold when the
# species-inclusion thresholds are relaxed?
#
# CURRENT thresholds (main analysis):
#   - Adults only
#   - Sex != "Unknown"
#   - n  >= 30 records per species
#   - temporal range >= 5 years per species
#
# This script fits the same model formula across a grid of alternative
# thresholds and compares the scaled_yr slope.  Because different thresholds
# yield different species sets, tree reconciliation is re-run per scenario.
# To keep runtime manageable, each scenario uses N_TREES trees (not 50);
# increase for the final run.
#
# Two-tier approach:
#   Tier 1 — all threshold combinations, non-phylogenetic brms (fast sweep).
#   Tier 2 — full phylogenetic model for the key comparison scenario
#             (sex filter dropped, n >= 30, range >= 5) vs. main analysis.
#
# Saves:
#   sensitivity_thresholds_tier1.rds  — Tier-1 pooled summaries
#   sensitivity_thresholds_tier2.rds  — Tier-2 full phylogenetic results
# ---------------------------------------------------------------------------

library(brms)
library(ape)
library(MCMCglmm)
library(prepR4pcm)
library(phytools)
library(future.apply)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

N_TREES <- 10L        # trees per scenario for Tier 2; Tier 1 ignores phylogeny
SEED    <- 20240202

set.seed(SEED)

# ── path helpers (mirrors atlantic_parallel.R) ─────────────────────────────
.find_file <- function(fname, subdir = "Analysis") {
  dirs <- getwd(); d <- getwd()
  for (i in 1:8) { d <- dirname(d); dirs <- c(dirs, d) }
  cand <- unique(c(file.path(dirs, fname), file.path(dirs, subdir, fname)))
  hit  <- cand[file.exists(cand)]
  if (!length(hit)) stop("Could not locate '", fname, "'. Run from inside the repo.")
  normalizePath(hit[1], winslash = "/")
}
ANALYSIS_DIR <- dirname(.find_file("passer90.rda"))
apath <- function(...) file.path(ANALYSIS_DIR, ...)

# ── raw data ───────────────────────────────────────────────────────────────
birds_raw <- read.csv(apath("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"))

base_filter <- function(df) {
  df %>%
    filter(Year >= 1990,
           Age  == "Adult",
           Order == "Passeriformes",
           AtlanticForests_20km_Buffer == "inside the 20 km polygon") %>%
    mutate(
      conc.wing.length  = coalesce(Wing_length_right.mm., Wing_length_left.mm., Wing_length.mm.),
      scaled_yr  = as.numeric(scale(Year)),
      scaled_lat = as.numeric(scale(Latitude_decimal_degrees * -1)),
      Binomial   = gsub(" ", "_", Binomial),
      spp        = Binomial,
      species_name = gsub("_", " ", Binomial)
    )
}

build_dataset <- function(df, min_n, min_range, drop_unknown_sex) {
  if (drop_unknown_sex) df <- filter(df, Sex != "Unknown")
  df %>%
    group_by(Binomial) %>%
    mutate(n_rec = n(), yr_min = min(Year), yr_max = max(Year),
           yr_range = yr_max - yr_min) %>%
    filter(n_rec >= min_n, yr_range >= min_range) %>%
    ungroup()
}

# ── eBird synonym map (same as atlantic_parallel.R) ────────────────────────
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
apply_synonyms <- function(df) {
  hit <- df$species_name %in% names(ebird_synonyms)
  df$species_name[hit] <- ebird_synonyms[df$species_name[hit]]
  df
}

# ── priors (same as main analysis) ─────────────────────────────────────────
priors_phylo <- c(prior(normal(71, 15), class = Intercept),
                  prior(normal(0, 10),  class = b),
                  prior(cauchy(0, 1),   class = sd),
                  prior(cauchy(0, 1),   class = sigma))

priors_nophylo <- priors_phylo   # same, no cov= term in this model

# ── threshold grid ─────────────────────────────────────────────────────────
threshold_grid <- expand.grid(
  min_n              = c(30L, 40L, 50L),
  min_range          = c(5L,  8L, 10L),
  drop_unknown_sex   = c(TRUE, FALSE),
  stringsAsFactors   = FALSE
)
# label for output
threshold_grid$label <- with(threshold_grid, paste0(
  "n>=", min_n, "_range>=", min_range,
  ifelse(drop_unknown_sex, "_sexfilter", "_nosexfilter")
))
message("Threshold grid: ", nrow(threshold_grid), " combinations")

# ============================================================================
# TIER 1: non-phylogenetic brms (species random intercept only)
# ============================================================================
message("\n=== TIER 1: non-phylogenetic sweep ===")

fit_nophylo <- function(dat_sc) {
  brm(conc.wing.length ~ 1 + Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp),
      data    = dat_sc,
      family  = gaussian(),
      prior   = priors_nophylo,
      iter    = 1500, warmup = 750, chains = 2, cores = 2,
      control = list(adapt_delta = 0.92),
      seed    = SEED,
      silent  = 2, refresh = 0)
}

base_dat <- base_filter(birds_raw)

tier1_results <- lapply(seq_len(nrow(threshold_grid)), function(i) {
  cfg <- threshold_grid[i, ]
  message("  [", i, "/", nrow(threshold_grid), "] ", cfg$label)

  dat_sc <- build_dataset(base_dat, cfg$min_n, cfg$min_range, cfg$drop_unknown_sex) %>%
    filter(!is.na(conc.wing.length), !is.na(Sex), Sex != "")

  n_spp <- n_distinct(dat_sc$Binomial)
  n_rec <- nrow(dat_sc)
  message("    N species = ", n_spp, ", N records = ", n_rec)

  if (n_spp < 10) {
    message("    Skipping — too few species")
    return(NULL)
  }

  fit  <- tryCatch(fit_nophylo(dat_sc), error = function(e) { message("    ERROR: ", e$message); NULL })
  if (is.null(fit)) return(NULL)

  fe <- fixef(fit)["scaled_yr", c("Estimate", "Q2.5", "Q97.5")]
  list(label    = cfg$label,
       min_n    = cfg$min_n,
       min_range= cfg$min_range,
       sex_filter = cfg$drop_unknown_sex,
       n_spp    = n_spp,
       n_rec    = n_rec,
       b_yr     = fe["Estimate"],
       b_yr_lo  = fe["Q2.5"],
       b_yr_hi  = fe["Q97.5"])
})

tier1_results <- Filter(Negate(is.null), tier1_results)
tier1_df <- bind_rows(lapply(tier1_results, as.data.frame))

message("\nTier 1 summary:")
print(tier1_df %>% arrange(b_yr))

saveRDS(tier1_df, apath("sensitivity_thresholds_tier1.rds"))
message("Saved: sensitivity_thresholds_tier1.rds")

# ── Tier 1 plot ─────────────────────────────────────────────────────────────
p_tier1 <- ggplot(tier1_df,
                  aes(x = b_yr, xmin = b_yr_lo, xmax = b_yr_hi,
                      y = reorder(label, b_yr),
                      colour = factor(min_n),
                      shape  = factor(sex_filter))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(height = 0.3) +
  geom_point(size = 2.5) +
  scale_colour_brewer(palette = "Dark2", name = "min records") +
  scale_shape_manual(values = c(16, 1), name = "Sex filter",
                     labels = c("FALSE" = "No", "TRUE" = "Yes")) +
  labs(x = "scaled_yr slope (mm per SD-year)",
       y = NULL,
       title = "Sensitivity: year → wing-length across inclusion thresholds",
       subtitle = "Non-phylogenetic models; dots = posterior mean, bars = 95% CrI") +
  theme_classic()

ggsave(apath("sensitivity_thresholds_tier1.png"), p_tier1,
       width = 9, height = 7, dpi = 150)
message("Saved: sensitivity_thresholds_tier1.png")

# ============================================================================
# TIER 2: full phylogenetic model for the most-relaxed scenario
#   n >= 30, range >= 5, sex filter dropped  (maximises species count)
# ============================================================================
message("\n=== TIER 2: full phylogenetic model (no sex filter, n>=30, range>=5) ===")

tryCatch({

dat_relaxed <- build_dataset(base_dat, min_n = 30L, min_range = 5L,
                             drop_unknown_sex = FALSE) %>%
  filter(!is.na(conc.wing.length)) %>%
  apply_synonyms()

message("Relaxed dataset: ", n_distinct(dat_relaxed$Binomial), " species, ",
        nrow(dat_relaxed), " records")

# retrieve and reconcile trees
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH")))
  clootl::get_avesdata_repo(path = ANALYSIS_DIR)

spp_relaxed <- unique(dat_relaxed$species_name)
# pr_get_tree only accepts n_tree = 100; subsample afterwards
got_relaxed <- pr_get_tree(spp_relaxed, source = "clootl",
                           n_tree = 100L, cache = TRUE)
rec_relaxed <- reconcile_tree(dat_relaxed, got_relaxed$tree[[1]],
                               x_species = "species_name",
                               fuzzy = TRUE, resolve = "flag")
print(reconcile_summary(rec_relaxed))

aligned_relaxed <- reconcile_apply(rec_relaxed, data = dat_relaxed,
                                   tree = got_relaxed$tree[[1]],
                                   species_col = "species_name",
                                   drop_unresolved = TRUE)
dat2  <- aligned_relaxed$data
dat2$species_name <- gsub(" ", "_", dat2$species_name)
dat2$spp          <- dat2$species_name

norm_us    <- function(x) gsub(" ", "_", x)
keep_tips  <- norm_us(aligned_relaxed$tree$tip.label)
trees_all  <- lapply(got_relaxed$tree, function(t) {
  t$tip.label <- norm_us(t$tip.label)
  ape::keep.tip(t, intersect(keep_tips, t$tip.label))
})
class(trees_all) <- "multiPhylo"
tree_samp2 <- sample(trees_all, N_TREES)

make_A <- function(tree) {
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A   <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

fit_one_phylo <- function(tree) {
  A <- make_A(tree)
  brm(conc.wing.length ~ 1 + Sex + scaled_yr + scaled_lat +
        (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)),
      data  = dat2, data2 = list(A = A),
      family  = gaussian(), prior = priors_phylo,
      iter    = 1500, warmup = 1000, chains = 1, cores = 1,
      control = list(max_treedepth = 12, adapt_delta = 0.95),
      seed    = SEED, silent = 2, refresh = 0)
}

plan(multisession)
fits_relaxed <- future_lapply(tree_samp2, fit_one_phylo, future.seed = TRUE)

pool_rubin <- function(fits,
                       pars = c("b_Intercept","b_SexMale","b_scaled_yr","b_scaled_lat")) {
  m     <- length(fits)
  draws <- lapply(fits, function(f) as_draws_df(f)[, pars, drop = FALSE])
  means <- t(sapply(draws, function(d) sapply(d, mean)))
  vars  <- t(sapply(draws, function(d) sapply(d, var)))
  qbar  <- colMeans(means)
  ubar  <- colMeans(vars)
  b     <- apply(means, 2, var)
  tot   <- ubar + (1 + 1/m) * b
  se    <- sqrt(tot)
  data.frame(par = pars, estimate = qbar, se = se,
             lower = qbar - 1.96*se, upper = qbar + 1.96*se)
}

rubin_relaxed <- pool_rubin(fits_relaxed)
message("\nTier 2 pooled results (no sex filter, n>=30, range>=5):")
print(rubin_relaxed)

# compare with main analysis pooled result
load(apath("brm0_multiphylo.rda"))   # → rubin_summary (main, n>=30, range>=5, sex filter)
message("\nMain (n>=30, range>=5, sex filter) vs No-sex-filter (n>=30, range>=5):")
print(rubin_summary)

saveRDS(
  list(
    relaxed_summary  = rubin_relaxed,
    relaxed_n_spp    = n_distinct(dat2$spp),
    relaxed_n_rec    = nrow(dat2),
    main_summary     = rubin_summary,
    tier1_df         = tier1_df
  ),
  apath("sensitivity_thresholds_tier2.rds")
)
message("Saved: sensitivity_thresholds_tier2.rds")

}, error = function(e) {
  message("ERROR in Tier 2: ", conditionMessage(e))
})
