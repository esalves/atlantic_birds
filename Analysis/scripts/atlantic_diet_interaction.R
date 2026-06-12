# atlantic_diet_interaction.R
# ---------------------------------------------------------------------------
# Tests whether the temporal wing-length decline interacts with diet:
# specifically, do species more reliant on invertebrates show a stronger
# body-size reduction over time?
#
# RATIONALE: if food limitation (declining arthropod prey) drives morphological
# change, the effect of year on wing length should be larger — more negative —
# in species with a higher proportion of invertebrates in their diet.
#
# TWO INTERACTION MODELS:
#   Model A (continuous): scaled_yr × diet_inv_std
#     diet_inv_std = z-scored proportion of diet from invertebrates (Diet-Inv
#     column from EltonTraits, 0–100 scale).
#   Model B (categorical): scaled_yr × diet_cat
#     diet_cat = Diet-5Cat from EltonTraits collapsed to Invertebrate vs. Other
#     (FruiNect, Omnivore, PlantSeed grouped as "Other").
#
# Both models include the same fixed effects as the main analysis
# (Sex, scaled_lat) and the same random structure (species intercept +
# slope, phylogenetic covariance across N_TREES trees).
#
# Saves:
#   diet_interaction_results.rds   — pooled summaries + draws for both models
#   diet_interaction_plot.png      — interaction plot (predicted year trend by
#                                    diet group / diet gradient)
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

N_TREES <- 10L
SEED    <- 20240303

set.seed(SEED)

# ── path helpers (repo layout: Analysis/{scripts,data,output,figures}) ───────
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts")))
      return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/ (need Analysis/data/derived and Analysis/scripts).")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)

# ── load data ────────────────────────────────────────────────────────────────
load(derived_path("passer90.rda"))   # → passer90 (already filtered; rerun Rmd to update thresholds)

elton_raw <- readr::read_tsv(raw_path("BirdFuncDat.txt"), show_col_types = FALSE)
elton_raw$Scientific <- gsub(" ", "_", elton_raw$Scientific)

# Join continuous invertebrate diet proportion (0–100)
diet_df <- elton_raw %>%
  select(Scientific, diet_inv = `Diet-Inv`, diet_cat = `Diet-5Cat`) %>%
  filter(Scientific %in% passer90$Binomial)

passer90 <- passer90 %>%
  left_join(diet_df, by = c("Binomial" = "Scientific")) %>%
  mutate(
    diet_cat2 = ifelse(diet_cat == "Invertebrate", "Invertebrate", "Other"),
    diet_cat2 = factor(diet_cat2, levels = c("Other", "Invertebrate"))
  )

# check coverage before dropping anything
n_missing_diet <- sum(is.na(passer90$diet_inv))
message("Records with missing Diet-Inv: ", n_missing_diet,
        " (", round(100 * n_missing_diet / nrow(passer90), 1), "%)")
missing_spp <- passer90 %>% filter(is.na(diet_inv)) %>% distinct(Binomial) %>% pull()
if (length(missing_spp)) {
  message("Species with no EltonTraits diet data (dropped from diet models):")
  message(paste(" -", missing_spp, collapse = "\n"))
}
message("Diet category breakdown (species with diet data):")
print(passer90 %>% filter(!is.na(diet_cat2)) %>%
        distinct(Binomial, diet_cat2) %>% count(diet_cat2))
message("Diet-Inv range (species with diet data):")
print(passer90 %>% filter(!is.na(diet_inv)) %>%
        distinct(Binomial, diet_inv) %>% pull(diet_inv) %>% summary())

# diet_inv_std computed after dropping NAs so scale() uses only present values
dat <- passer90 %>%
  filter(!is.na(diet_inv), !is.na(conc.wing.length),
         !is.na(Sex), Sex != "") %>%
  mutate(
    diet_inv_std = as.numeric(scale(diet_inv)),
    species_name = gsub("_", " ", Binomial)
  )

# ── eBird synonym map ────────────────────────────────────────────────────────
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
hit <- dat$species_name %in% names(ebird_synonyms)
dat$species_name[hit] <- ebird_synonyms[dat$species_name[hit]]

# ── phylogeny ────────────────────────────────────────────────────────────────
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH")))
  clootl::get_avesdata_repo(path = raw_path())

spp_data <- unique(dat$species_name)
.check <- pr_get_tree(spp_data, source = "clootl", n_tree = 1)
if (length(.check$unmatched) > 0) {
  message("Removing ", length(.check$unmatched), " species absent from eBird taxonomy: ",
          paste(.check$unmatched, collapse = ", "))
  spp_data <- setdiff(spp_data, .check$unmatched)
}
# pr_get_tree only accepts n_tree = 100; subsample afterwards
got  <- pr_get_tree(spp_data, source = "clootl", n_tree = 100L, cache = TRUE)
rec  <- reconcile_tree(dat, got$tree[[1]], x_species = "species_name",
                       fuzzy = TRUE, resolve = "flag")
print(reconcile_summary(rec))

aligned  <- reconcile_apply(rec, data = dat, tree = got$tree[[1]],
                             species_col = "species_name", drop_unresolved = TRUE)
dat      <- aligned$data
dat$species_name <- gsub(" ", "_", dat$species_name)
dat$spp          <- dat$species_name

norm_us   <- function(x) gsub(" ", "_", x)
keep_tips <- norm_us(aligned$tree$tip.label)
trees_all <- lapply(got$tree, function(t) {
  t$tip.label <- norm_us(t$tip.label)
  ape::keep.tip(t, intersect(keep_tips, t$tip.label))
})
class(trees_all) <- "multiPhylo"
tree_samp <- sample(trees_all, N_TREES)

make_A <- function(tree) {
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A   <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

# ── priors ───────────────────────────────────────────────────────────────────
priors <- c(prior(normal(71, 15), class = Intercept),
            prior(normal(0, 10),  class = b),
            prior(cauchy(0, 1),   class = sd),
            prior(cauchy(0, 1),   class = sigma))

# ── Rubin pooling helper ─────────────────────────────────────────────────────
pool_rubin <- function(fits, pars) {
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

# ============================================================================
# MODEL A: continuous interaction  scaled_yr × diet_inv_std
# ============================================================================
message("\n=== Model A: scaled_yr × diet_inv_std (continuous) ===")

fit_A <- function(tree) {
  A <- make_A(tree)
  brm(
    conc.wing.length ~ 1 + Sex + scaled_yr * diet_inv_std + scaled_lat +
      (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)),
    data = dat, data2 = list(A = A),
    family  = gaussian(), prior = priors,
    iter    = 2000, warmup = 1000, chains = 1, cores = 1,
    control = list(max_treedepth = 12, adapt_delta = 0.95),
    seed    = SEED, silent = 2, refresh = 0
  )
}

plan(multisession)
fits_A <- future_lapply(tree_samp, fit_A, future.seed = TRUE)

pars_A <- c("b_Intercept", "b_SexMale", "b_scaled_yr",
            "b_diet_inv_std", "b_scaled_lat", "b_scaled_yr:diet_inv_std")
rubin_A <- pool_rubin(fits_A, pars_A)
message("\nModel A pooled results:")
print(rubin_A)

# ============================================================================
# MODEL B: categorical interaction  scaled_yr × diet_cat2
# ============================================================================
message("\n=== Model B: scaled_yr × diet_cat2 (Invertebrate vs. Other) ===")

fit_B <- function(tree) {
  A <- make_A(tree)
  brm(
    conc.wing.length ~ 1 + Sex + scaled_yr * diet_cat2 + scaled_lat +
      (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)),
    data = dat, data2 = list(A = A),
    family  = gaussian(), prior = priors,
    iter    = 2000, warmup = 1000, chains = 1, cores = 1,
    control = list(max_treedepth = 12, adapt_delta = 0.95),
    seed    = SEED, silent = 2, refresh = 0
  )
}

fits_B <- future_lapply(tree_samp, fit_B, future.seed = TRUE)

pars_B <- c("b_Intercept", "b_SexMale", "b_scaled_yr",
            "b_diet_cat2Invertebrate", "b_scaled_lat",
            "b_scaled_yr:diet_cat2Invertebrate")
rubin_B <- pool_rubin(fits_B, pars_B)
message("\nModel B pooled results:")
print(rubin_B)

# ============================================================================
# Plots
# ============================================================================

# -- Model A: predicted year trend at low / mean / high diet_inv_std ----------
yr_seq   <- seq(-2, 2, length.out = 50)
diet_seq <- c(-1, 0, 1)   # -1 SD, mean, +1 SD invertebrate diet

# use one representative fit for visualisation
fit_A_rep <- fits_A[[1]]

pred_A <- expand.grid(
  scaled_yr     = yr_seq,
  diet_inv_std  = diet_seq,
  Sex           = "Female",
  scaled_lat    = 0,
  spp           = NA_character_,
  species_name  = NA_character_
)
pred_A$fitted <- predict(fit_A_rep, newdata = pred_A,
                          allow_new_levels = TRUE,
                          re_formula = NA)[, "Estimate"]
pred_A$diet_label <- factor(pred_A$diet_inv_std,
                             labels = c("−1 SD (low)", "Mean", "+1 SD (high)"))

# back-convert scaled_yr to calendar year for the x-axis
yr_mu <- mean(passer90$Year, na.rm = TRUE)
yr_sd <- sd(passer90$Year,  na.rm = TRUE)
pred_A$year <- pred_A$scaled_yr * yr_sd + yr_mu

p_A <- ggplot(pred_A, aes(x = year, y = fitted, colour = diet_label)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(
    values = c("−1 SD (low)" = "#1F77B4",
               "Mean"            = "#FF7F0E",
               "+1 SD (high)"    = "#D62728"),
    name = "Invertebrate\ndiet proportion"
  ) +
  labs(x = "Year",
       y = "Predicted wing length (mm)",
       title = "Interaction: year × invertebrate diet proportion",
       subtitle = paste0("Model A (continuous); interaction β = ",
                         round(rubin_A$estimate[rubin_A$par == "b_scaled_yr:diet_inv_std"], 3),
                         " [",
                         round(rubin_A$lower[rubin_A$par == "b_scaled_yr:diet_inv_std"], 3),
                         ", ",
                         round(rubin_A$upper[rubin_A$par == "b_scaled_yr:diet_inv_std"], 3),
                         "]")) +
  theme_classic()

# -- Model B: predicted year trend by diet_cat2 --------------------------------
pred_B <- expand.grid(
  scaled_yr    = yr_seq,
  diet_cat2    = levels(dat$diet_cat2),
  Sex          = "Female",
  scaled_lat   = 0,
  spp          = NA_character_,
  species_name = NA_character_
)
fit_B_rep <- fits_B[[1]]
pred_B$fitted <- predict(fit_B_rep, newdata = pred_B,
                          allow_new_levels = TRUE,
                          re_formula = NA)[, "Estimate"]
pred_B$year <- pred_B$scaled_yr * yr_sd + yr_mu

p_B <- ggplot(pred_B, aes(x = year, y = fitted, colour = diet_cat2)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Other" = "#1F77B4", "Invertebrate" = "#D62728"),
                      name = "Diet category") +
  labs(x = "Year",
       y = "Predicted wing length (mm)",
       title = "Interaction: year × diet category",
       subtitle = paste0("Model B (categorical); interaction β = ",
                         round(rubin_B$estimate[rubin_B$par == "b_scaled_yr:diet_cat2Invertebrate"], 3),
                         " [",
                         round(rubin_B$lower[rubin_B$par == "b_scaled_yr:diet_cat2Invertebrate"], 3),
                         ", ",
                         round(rubin_B$upper[rubin_B$par == "b_scaled_yr:diet_cat2Invertebrate"], 3),
                         "]")) +
  theme_classic()

# combine and save
library(patchwork)
p_combined <- p_A / p_B
img_dir <- file.path(dirname(ANALYSIS_DIR), "Manuscript", "images")
ggsave(file.path(img_dir, "diet_interaction_plot.png"), p_combined,
       width = 7, height = 9, dpi = 150)
message("Saved: Manuscript/images/diet_interaction_plot.png")

# ============================================================================
# Save all results
# ============================================================================
saveRDS(
  list(
    model_A_pooled = rubin_A,
    model_B_pooled = rubin_B,
    n_spp          = n_distinct(dat$spp),
    n_rec          = nrow(dat),
    diet_summary   = dat %>% distinct(spp, diet_inv, diet_cat2) %>%
                       group_by(diet_cat2) %>%
                       summarise(n = n(), mean_inv = mean(diet_inv), .groups = "drop")
  ),
  out_path("diet_interaction_results.rds")
)
message("Saved: diet_interaction_results.rds")

message("\n=== Key results ===")
message("Model A — scaled_yr:diet_inv_std interaction:")
message("  β = ", round(rubin_A$estimate[rubin_A$par == "b_scaled_yr:diet_inv_std"], 3),
        "  95% CI [",
        round(rubin_A$lower[rubin_A$par == "b_scaled_yr:diet_inv_std"], 3), ", ",
        round(rubin_A$upper[rubin_A$par == "b_scaled_yr:diet_inv_std"], 3), "]")
message("Model B — scaled_yr:Invertebrate interaction:")
message("  β = ", round(rubin_B$estimate[rubin_B$par == "b_scaled_yr:diet_cat2Invertebrate"], 3),
        "  95% CI [",
        round(rubin_B$lower[rubin_B$par == "b_scaled_yr:diet_cat2Invertebrate"], 3), ", ",
        round(rubin_B$upper[rubin_B$par == "b_scaled_yr:diet_cat2Invertebrate"], 3), "]")
message("\nInterpretation: a negative interaction term means insectivorous species")
message("show a stronger (more negative) year → wing-length trend, consistent")
message("with food-limitation driving body-size decline.")
