# atlantic_diet_interaction.R
# ---------------------------------------------------------------------------
# WHAT: Tests whether the temporal wing-length trend interacts with diet:
# do species more reliant on invertebrates show a stronger (more negative)
# year -> wing-length slope? Two interaction models (continuous Diet-Inv and a
# binary Invertebrate / Other category) with the fixed and random structure of
# the primary wing model, plus (2026-09 revision, REVISION_PLAN.md §1 "§3.3" and
# §3 "Phase 5") the descriptive material the referee asked for:
#
#   (a) a species-level distribution of Diet-Inv with the number of species and
#       records per 10-point bin (the covariate is bimodal: most species sit at
#       90-100 or <= 10, so a "+/-1 SD" projection lands on a handful of clades);
#   (b) predictions of the year slope at the OBSERVED 10th / 50th / 90th
#       percentiles of Diet-Inv instead of at +/-1 SD, with the species and
#       record support behind each percentile;
#   (c) a family-level random-effect sensitivity (Diet-Inv is strongly clustered
#       by family: Thamnophilidae / Conopophagidae at the top, Thraupidae /
#       Pipridae / Emberizidae at the bottom), so the interaction is re-fitted
#       with (1 | Family) and with (1 + scaled_yr || Family);
#   (d) an lme4 (REML, non-phylogenetic) analogue of every model so the new
#       pieces run in seconds without Stan or the tree cloud;
#   (e) (2026-09, phylogenetic tier) every model re-fitted WITH the phylogenetic
#       species term across the 50 trees of the published analysis using glmmTMB's
#       `propto` covariance structure (Williams, McGillycuddy, Drobniak, Bolker,
#       Warton & Nakagawa 2025, bioRxiv 10.64898/2025.12.20.695312) through the
#       shared engine scripts/_phylo_engine.R, Rubin-pooled over trees. glmmTMB
#       reproduces the published 50-tree brms wing model to 2-3 decimals
#       (glmmtmb_validation_wing.R) in seconds, so the interaction can now carry
#       phylogeny AND the Phase 1 provenance controls (contributor + municipality
#       intercepts) in the same fit, which brms could not afford locally.
#
# RATIONALE (unchanged): if food limitation (declining arthropod prey) drives
# morphological change, the year effect on wing length should be more negative
# in species with a higher proportion of invertebrates in the diet.
#
# MODELS (fixed effects Sex + scaled_yr * diet + scaled_lat; species random
# intercept + uncorrelated species year slope (1 + scaled_yr || spp) throughout):
#   Model A (continuous):  diet = diet_inv_std                 + phylogenetic term
#   Model B (categorical): diet = diet_cat2 (Invertebrate vs Other) + phylogenetic term
#   Model A_fam:           Model A + (1 | Family)               (clade sensitivity)
#   --- glmmTMB engine only (phylogenetic, 50 trees) ---
#   Model A_src_site:      Model A + (1 | Main_researcher) + (1 | Municipality)
#   Model A_fam_src_site:  Model A_fam + (1 | Main_researcher) + (1 | Municipality)
#   Model B_src_site:      Model B + (1 | Main_researcher) + (1 | Municipality)
#   --- lme4 only (non-phylogenetic, seconds) ---
#   Model A + (1 + scaled_yr || Family); Model A + contributor + municipality
#   Phylogenetic term: brms  (1 | gr(species_name, cov = A));
#                      glmmTMB propto(0 + species_name | g, A)  (appended by the engine)
#   diet_inv_std = z-scored Diet-Inv (EltonTraits, 0-100), standardised on the
#   RECORD-level analytical sample exactly as before (record-weighted mean/SD),
#   so the coefficients stay comparable with the published fit.
#
# INPUTS:  data/derived/passer90.rda            (live-only, known sex, 73 spp)
#          data/raw/BirdFuncDat.txt             (EltonTraits 1.0)
#          scripts/_phylo_engine.R              (glmmTMB engine; default path)
#          data/derived/phylo_A_50trees.rds     (the 50 published-tree correlation
#                                               matrices; built by the engine from
#                                               output/models/brm0_multiphylo.rda)
#          data/raw/AvesDataLite-main/          (clootl tree cloud; brms path only)
#          scripts/_sampling_config.R           (SAMPLING, SAMPLING_CONTROL; brms only)
# OUTPUTS: output/diet_distribution.rds         species table, bin counts,
#                                               quantiles and their support
#          figures/diet_distribution.png        species and records per bin
#          output/diet_lme4_results.rds         lme4 analogues + quantile slopes
#          figures/diet_quantile_predictions.png lme4 predicted trajectories at
#                                               the 10/50/90th percentiles
#          --- glmmTMB engine (default for --trees) ---
#          output/diet_interaction_phylo.rds    pooled fixed effects, per-tree
#                                               tables, variance components,
#                                               quantile slopes, comparison with
#                                               the published brms fit and lme4
#          figures/diet_interaction_plot.png    Model A at the observed percentiles
#                                               + Model B by category (50 trees)
#          figures/diet_quantile_predictions_phylo.png  A / A_fam / A_src_site
#          figures/diet_interaction_specifications.png  interaction coefficient
#                                               across engines and controls
#          figures/diet_species_slopes.png      species year slopes vs Diet-Inv
#          --- brms engine only (Bayesian cross-check; Totoro) ---
#          output/diet_interaction_results.rds  model_A_pooled, model_B_pooled,
#                                               n_spp, n_rec, diet_summary (as
#                                               before) + model_A_family_pooled,
#                                               quantile_slopes, quantile_support
#          figures/diet_interaction_plot_brms.png
#
# RUN:  Rscript Analysis/scripts/atlantic_diet_interaction.R                  # glmmTMB, 50 trees (~5 min)
#       Rscript Analysis/scripts/atlantic_diet_interaction.R --trees 50       # same, explicit
#       Rscript Analysis/scripts/atlantic_diet_interaction.R --no-brms        # descriptive + lme4 only
#       Rscript Analysis/scripts/atlantic_diet_interaction.R --engine brms --trees 10  # published
#                                                                             # brms path (Totoro)
#       Rscript Analysis/scripts/atlantic_diet_interaction.R --engine brms --smoke    # 1 tree,
#                                                                             # 2 chains, 400 iter:
#                                                                             # a SMOKE TEST
#       --no-brms (alias --no-phylo) and the glmmTMB engine never touch
#       diet_interaction_results.rds (the published 10-tree brms object).
#
# Session (2026-09-09): R 4.6.0, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0,
# lme4 2.0.1, glmmTMB 1.1.14, ggplot2 4.0.3, patchwork 1.3.2 (glmmTMB path and
# --no-brms path executed locally in full). brms path last run with brms 2.23.0 /
# prepR4pcm 0.5.0.9000 / clootl 0.1.4 (10 trees, server). Versions are recorded
# in the output rds files.
# ---------------------------------------------------------------------------

# ── command-line flags ───────────────────────────────────────────────────────
.args     <- commandArgs(trailingOnly = TRUE)
RUN_PHYLO <- !any(c("--no-brms", "--no-phylo") %in% .args)
SMOKE     <- "--smoke" %in% .args
ENGINE    <- "glmmTMB"                       # default for --trees mode
if (any(grepl("^--engine(=|$)", .args))) {
  i <- grep("^--engine(=|$)", .args)[1]
  v <- sub("^--engine=?", "", .args[i])
  if (!nzchar(v) && length(.args) > i) v <- .args[i + 1]
  ENGINE <- match.arg(tolower(v), c("glmmtmb", "brms"))
  ENGINE <- if (ENGINE == "glmmtmb") "glmmTMB" else "brms"
}
RUN_BRMS <- RUN_PHYLO && ENGINE == "brms"
N_TREES  <- if (ENGINE == "brms") 10L else 50L   # 10 = published brms fit; 50 = full tree sample
if (any(grepl("^--trees(=|$)", .args))) {
  i <- grep("^--trees(=|$)", .args)[1]
  v <- sub("^--trees=?", "", .args[i])
  if (!nzchar(v) && length(.args) > i) v <- .args[i + 1]
  N_TREES <- as.integer(v)
  if (is.na(N_TREES) || N_TREES < 1L) stop("--trees needs a positive integer")
}
if (SMOKE) N_TREES <- 1L
SEED <- 20240303
set.seed(SEED)

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
})
if (RUN_BRMS) suppressMessages({
  library(brms)
  library(ape)
  library(MCMCglmm)
  library(prepR4pcm)
  library(phytools)
  library(future.apply)
  library(posterior)
})

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
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
if (RUN_BRMS) source(script_path("_sampling_config.R"))  # SAMPLING settings
if (!dir.exists(fig_path())) dir.create(fig_path(), recursive = TRUE)

.pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
session_info <- list(
  date = as.character(Sys.Date()), R = R.version.string,
  packages = sapply(c("dplyr", "tidyr", "lme4", "glmmTMB", "ggplot2", "patchwork",
                      "brms", "prepR4pcm", "clootl", "posterior"), .pkg_ver),
  flags = list(run_phylo = RUN_PHYLO, engine = ENGINE, run_brms = RUN_BRMS,
               smoke = SMOKE, n_trees = N_TREES)
)

# ── load data ────────────────────────────────────────────────────────────────
load(derived_path("passer90.rda"))   # → passer90 (live-only, known sex, 73 species)
# scaled_yr / scaled_lat are 1-column scale() matrices in passer90.rda; lme4 and
# expand.grid want plain numerics.
passer90$scaled_yr  <- as.numeric(passer90$scaled_yr)
passer90$scaled_lat <- as.numeric(passer90$scaled_lat)
# Year centre / SD actually used for scaled_yr (pre-threshold reference sample,
# 2009.487 / 5.020; see REVISION_NOTES_P0.md). Read from the scale() attributes
# first, then the "scaling" attribute, then fall back to the sample itself.
.sc <- attributes(passer90$scaled_yr)
yr_mu <- .sc[["scaled:center"]]; yr_sd <- .sc[["scaled:scale"]]
if (is.null(yr_mu) || is.null(yr_sd)) {
  yr_mu <- attr(passer90, "scaling")$year[["center"]]
  yr_sd <- attr(passer90, "scaling")$year[["scale"]]
}
if (is.null(yr_mu) || is.null(yr_sd)) {          # older builds: recover from the data
  yr_mu <- mean(passer90$Year, na.rm = TRUE); yr_sd <- sd(passer90$Year, na.rm = TRUE)
}
yr_mu <- as.numeric(yr_mu); yr_sd <- as.numeric(yr_sd)

elton_raw <- readr::read_tsv(raw_path("BirdFuncDat.txt"), show_col_types = FALSE)
elton_raw$Scientific <- gsub(" ", "_", elton_raw$Scientific)

# Join continuous invertebrate diet proportion (0–100) + EltonTraits certainty
diet_df <- elton_raw %>%
  select(Scientific, diet_inv = `Diet-Inv`, diet_cat = `Diet-5Cat`,
         diet_certainty = `Diet-Certainty`) %>%
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

# diet_inv_std computed after dropping NAs so scale() uses only present values.
# NOTE: scale() here is RECORD-weighted (one row per record), as in the published
# fit; the species-weighted mean/SD are reported in diet_distribution.rds.
dat <- passer90 %>%
  filter(!is.na(diet_inv), !is.na(conc.wing.length),
         !is.na(Sex), Sex != "") %>%
  mutate(
    diet_inv_std = as.numeric(scale(diet_inv)),
    species_name = gsub("_", " ", Binomial)
  )
diet_center <- mean(dat$diet_inv); diet_scale <- sd(dat$diet_inv)
z_of <- function(x) (x - diet_center) / diet_scale     # Diet-Inv (0-100) -> diet_inv_std

# ============================================================================
# (a) Descriptive: species-level Diet-Inv distribution with record support
# ============================================================================
message("\n=== Diet-Inv distribution (species level) ===")

spp_tab <- dat %>%
  group_by(Binomial, Family, diet_inv, diet_cat, diet_cat2, diet_certainty) %>%
  summarise(n_wing_records = n(), n_years = n_distinct(Year), .groups = "drop") %>%
  left_join(passer90 %>% count(Binomial, name = "n_all_records"), by = "Binomial") %>%
  arrange(diet_inv, Binomial)

bin_breaks <- seq(0, 100, by = 10)
spp_tab$diet_bin <- cut(spp_tab$diet_inv, breaks = bin_breaks, include.lowest = TRUE,
                        right = TRUE, labels = paste0(head(bin_breaks, -1), "-", bin_breaks[-1]))
bin_tab <- spp_tab %>%
  group_by(diet_bin, .drop = FALSE) %>%
  summarise(n_spp = n(), n_records = sum(n_wing_records),
            families = paste(sort(unique(Family)), collapse = ", "), .groups = "drop") %>%
  mutate(pct_spp = 100 * n_spp / sum(n_spp), pct_records = 100 * n_records / sum(n_records))
print(as.data.frame(bin_tab %>% select(diet_bin, n_spp, n_records, families)))

# Observed quantiles: species-weighted (each species once) is the primary set
# because Diet-Inv is a species trait; record-weighted is reported for
# completeness (it is what scale() used).
probs <- c(0.10, 0.50, 0.90)
q_spp <- quantile(spp_tab$diet_inv, probs, type = 7, names = FALSE)
q_rec <- quantile(dat$diet_inv,     probs, type = 7, names = FALSE)
# Support: species / records within +/- 5 Diet-Inv points of each quantile
support_at <- function(q, half_width = 5) {
  s <- spp_tab %>% filter(abs(diet_inv - q) <= half_width)
  data.frame(diet_inv = q, n_spp = nrow(s), n_records = sum(s$n_wing_records),
             families = paste(names(sort(table(s$Family), decreasing = TRUE)), collapse = ", "))
}
quantile_support <- bind_rows(lapply(seq_along(probs), function(i)
  cbind(quantile = paste0("p", probs[i] * 100), support_at(q_spp[i]))))
old_sd_points <- c(minus1sd = diet_center - diet_scale, mean = diet_center,
                   plus1sd = diet_center + diet_scale)
old_sd_support <- bind_rows(lapply(names(old_sd_points), function(nm)
  cbind(point = nm, support_at(old_sd_points[[nm]]))))
message("Species-weighted Diet-Inv percentiles (10/50/90): ",
        paste(round(q_spp, 1), collapse = " / "),
        "; record-weighted: ", paste(round(q_rec, 1), collapse = " / "))
message("Old +/-1 SD projection points (record-weighted): ",
        paste(round(old_sd_points, 1), collapse = " / "))
print(quantile_support)

# Family composition of the Diet-Inv scale
fam_tab <- spp_tab %>%
  group_by(Family) %>%
  summarise(n_spp = n(), n_records = sum(n_wing_records),
            diet_inv_median = median(diet_inv), diet_inv_min = min(diet_inv),
            diet_inv_max = max(diet_inv), .groups = "drop") %>%
  arrange(desc(n_spp))
print(as.data.frame(fam_tab))

# --- figure: species per bin (stacked by family) with record counts ---------
top_fam <- head(fam_tab$Family, 6)
plot_spp <- spp_tab %>%
  mutate(Family_plot = ifelse(Family %in% top_fam, Family, "Other families"),
         Family_plot = factor(Family_plot, levels = c(top_fam, "Other families")))
fam_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#999999")
names(fam_cols) <- levels(plot_spp$Family_plot)
bin_lab <- bin_tab %>% mutate(lab = paste0(n_spp, " spp\n", format(n_records, big.mark = ","), " rec"))

p_dist <- ggplot(plot_spp, aes(x = diet_bin, fill = Family_plot)) +
  geom_bar(colour = "grey30", linewidth = 0.2) +
  geom_text(data = bin_lab, aes(x = diet_bin, y = n_spp, label = lab),
            inherit.aes = FALSE, vjust = -0.15, size = 2.6, lineheight = 0.85) +
  scale_fill_manual(values = fam_cols, name = "Family", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Diet-Inv (EltonTraits, % of diet from invertebrates)",
       y = "Number of species",
       title = "Species-level distribution of the diet covariate",
       subtitle = sprintf("%d species, %s wing records; labels give species and records per bin",
                          nrow(spp_tab), format(nrow(dat), big.mark = ","))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# Second panel: where the projection points fall. Rug of species (jittered),
# observed percentiles (solid) vs the old +/-1 SD points (dotted).
marks <- bind_rows(
  data.frame(x = q_spp, lab = paste0("p", probs * 100), type = "observed percentile"),
  data.frame(x = old_sd_points, lab = c("-1 SD", "mean", "+1 SD"), type = "old +/-1 SD projection"))
p_marks <- ggplot(spp_tab, aes(x = diet_inv, y = n_wing_records)) +
  geom_vline(data = marks, aes(xintercept = x, linetype = type), colour = "grey40") +
  geom_point(aes(colour = Family %in% top_fam), alpha = 0.7, size = 2, show.legend = FALSE) +
  geom_text(data = marks, aes(x = x, y = max(spp_tab$n_wing_records) * 1.05, label = lab),
            inherit.aes = FALSE, size = 2.6, angle = 90, hjust = 1, vjust = -0.3, colour = "grey30") +
  scale_colour_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey60")) +
  scale_linetype_manual(values = c("observed percentile" = "solid",
                                   "old +/-1 SD projection" = "dotted"), name = NULL) +
  scale_y_log10() +
  labs(x = "Diet-Inv (%)", y = "Wing records per species (log scale)",
       title = "Where the projections sit on the observed covariate") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")

ggsave(fig_path("diet_distribution.png"), p_dist / p_marks + plot_layout(heights = c(1.2, 1)),
       width = 8, height = 8.5, dpi = 150)
message("Saved: ", fig_path("diet_distribution.png"))

diet_distribution <- list(
  species_table   = spp_tab,
  bins            = bin_tab,
  families        = fam_tab,
  quantiles       = data.frame(quantile = paste0("p", probs * 100), prob = probs,
                               diet_inv_species_weighted = q_spp,
                               diet_inv_record_weighted  = q_rec,
                               diet_inv_std_species_weighted = z_of(q_spp)),
  quantile_support = quantile_support,
  standardisation = list(center_record_weighted = diet_center, sd_record_weighted = diet_scale,
                         center_species_weighted = mean(spp_tab$diet_inv),
                         sd_species_weighted = sd(spp_tab$diet_inv),
                         old_sd_points = old_sd_points, old_sd_support = old_sd_support),
  n_spp = nrow(spp_tab), n_rec = nrow(dat),
  n_spp_total = n_distinct(passer90$Binomial), missing_spp = missing_spp,
  certainty = as.data.frame(table(certainty = spp_tab$diet_certainty)),
  session = session_info
)
saveRDS(diet_distribution, out_path("diet_distribution.rds"))
message("Saved: ", out_path("diet_distribution.rds"))

# ============================================================================
# (d) lme4 analogues (fast, non-phylogenetic; run on every path)
# ============================================================================
message("\n=== lme4 analogues (REML; species intercept + uncorrelated year slope) ===")
dat_l <- dat %>% mutate(spp = Binomial, src = Main_researcher, site = Municipality)

fml <- list(
  A_baseline   = conc.wing.length ~ Sex + scaled_yr * diet_inv_std + scaled_lat +
                   (1 + scaled_yr || spp),
  A_family_int = conc.wing.length ~ Sex + scaled_yr * diet_inv_std + scaled_lat +
                   (1 + scaled_yr || spp) + (1 | Family),
  A_family_slp = conc.wing.length ~ Sex + scaled_yr * diet_inv_std + scaled_lat +
                   (1 + scaled_yr || spp) + (1 + scaled_yr || Family),
  A_src_site   = conc.wing.length ~ Sex + scaled_yr * diet_inv_std + scaled_lat +
                   (1 + scaled_yr || spp) + (1 | src) + (1 | site),
  A_src_site_family_slp = conc.wing.length ~ Sex + scaled_yr * diet_inv_std + scaled_lat +
                   (1 + scaled_yr || spp) + (1 + scaled_yr || Family) + (1 | src) + (1 | site),
  B_baseline   = conc.wing.length ~ Sex + scaled_yr * diet_cat2 + scaled_lat +
                   (1 + scaled_yr || spp),
  B_family_slp = conc.wing.length ~ Sex + scaled_yr * diet_cat2 + scaled_lat +
                   (1 + scaled_yr || spp) + (1 + scaled_yr || Family)
)
lme4_fits <- lapply(fml, function(f)
  lmer(f, data = dat_l, REML = TRUE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))))

tidy_fixef <- function(fit, model) {
  cf <- summary(fit)$coefficients
  data.frame(model = model, par = rownames(cf), estimate = cf[, "Estimate"],
             se = cf[, "Std. Error"], t = cf[, "t value"],
             lower = cf[, "Estimate"] - 1.96 * cf[, "Std. Error"],
             upper = cf[, "Estimate"] + 1.96 * cf[, "Std. Error"], row.names = NULL)
}
lme4_fixef <- bind_rows(lapply(names(lme4_fits), function(m) tidy_fixef(lme4_fits[[m]], m)))
lme4_convergence <- sapply(lme4_fits, function(f) {
  msg <- f@optinfo$conv$lme4$messages; if (is.null(msg)) "ok" else paste(msg, collapse = "; ")
})
lme4_varcomp <- bind_rows(lapply(names(lme4_fits), function(m) {
  vc <- as.data.frame(VarCorr(lme4_fits[[m]]))
  data.frame(model = m, grp = vc$grp, var1 = vc$var1, sd = vc$sdcor)
}))
int_rows <- lme4_fixef %>% filter(grepl("^scaled_yr:", par) | par == "scaled_yr") %>%
  select(model, par, estimate, se, t, lower, upper)
print(int_rows, digits = 3)
message("Convergence: ", paste(names(lme4_convergence), lme4_convergence, sep = " = ", collapse = "; "))

# Year slope at the observed Diet-Inv percentiles: slope(q) = b_yr + b_int * z(q),
# Var = v_yy + z^2 v_ii + 2 z v_yi (fixed-effect vcov; Wald 95% interval).
slope_at_quantiles <- function(fit, model, z, labels, diet_vals) {
  b <- fixef(fit); V <- as.matrix(vcov(fit))
  est <- b["scaled_yr"] + b["scaled_yr:diet_inv_std"] * z
  se  <- sqrt(V["scaled_yr", "scaled_yr"] + z^2 * V["scaled_yr:diet_inv_std", "scaled_yr:diet_inv_std"] +
              2 * z * V["scaled_yr", "scaled_yr:diet_inv_std"])
  data.frame(model = model, quantile = labels, diet_inv = diet_vals, diet_inv_std = z,
             slope_per_sd_year = unname(est), se = unname(se),
             lower = unname(est - 1.96 * se), upper = unname(est + 1.96 * se),
             mm_per_decade = unname(est) * 10 / yr_sd, row.names = NULL)
}
cont_models <- grep("^A_", names(lme4_fits), value = TRUE)
lme4_quantile_slopes <- bind_rows(lapply(cont_models, function(m)
  slope_at_quantiles(lme4_fits[[m]], m, z_of(q_spp), paste0("p", probs * 100), q_spp)))
message("\nYear slope (mm per SD-year) at species-weighted Diet-Inv percentiles, lme4:")
print(lme4_quantile_slopes %>% select(model, quantile, diet_inv, slope_per_sd_year, lower, upper, mm_per_decade),
      digits = 3)

# Predicted trajectories (fixed effects only; females at mean latitude) for the
# quantile figure, baseline and family-slope models side by side.
yr_seq <- seq(min(dat$scaled_yr), max(dat$scaled_yr), length.out = 50)
pred_lme4 <- function(fit, model) {
  nd <- expand.grid(scaled_yr = yr_seq, diet_inv = q_spp, Sex = "Female", scaled_lat = 0,
                    stringsAsFactors = FALSE)
  nd$Sex <- factor(nd$Sex, levels = c("Female", "Male"))   # full levels for model.matrix
  nd$diet_inv_std <- z_of(nd$diet_inv)
  X  <- model.matrix(~ Sex + scaled_yr * diet_inv_std + scaled_lat, nd)
  X  <- X[, names(fixef(fit)), drop = FALSE]
  b  <- fixef(fit); V <- as.matrix(vcov(fit))
  nd$fitted <- as.numeric(X %*% b)
  nd$se     <- sqrt(rowSums((X %*% V) * X))
  # Change relative to the record midpoint (scaled_yr = 0, i.e. 2009.5): the
  # ribbon then shows the uncertainty of the SLOPE at each percentile rather than
  # the (much larger) uncertainty of the grand intercept across species.
  z   <- z_of(nd$diet_inv)
  se_slope <- sqrt(V["scaled_yr", "scaled_yr"] + z^2 * V["scaled_yr:diet_inv_std", "scaled_yr:diet_inv_std"] +
                   2 * z * V["scaled_yr", "scaled_yr:diet_inv_std"])
  nd$delta    <- (b["scaled_yr"] + b["scaled_yr:diet_inv_std"] * z) * nd$scaled_yr
  nd$delta_se <- abs(nd$scaled_yr) * se_slope
  nd$year   <- nd$scaled_yr * yr_sd + yr_mu
  nd$quantile <- factor(paste0("p", probs * 100)[match(nd$diet_inv, q_spp)],
                        levels = paste0("p", probs * 100))
  nd$model <- model; nd
}
pred_q <- bind_rows(pred_lme4(lme4_fits$A_baseline, "A_baseline"),
                    pred_lme4(lme4_fits$A_family_slp, "A_family_slp"),
                    pred_lme4(lme4_fits$A_src_site, "A_src_site"))
q_labels <- sprintf("%s: Diet-Inv %.0f (%d spp, %s rec)", paste0("p", probs * 100), q_spp,
                    quantile_support$n_spp, format(quantile_support$n_records, big.mark = ","))
names(q_labels) <- paste0("p", probs * 100)
p_q <- ggplot(pred_q, aes(x = year, y = delta, colour = quantile, fill = quantile)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_ribbon(aes(ymin = delta - 1.96 * delta_se, ymax = delta + 1.96 * delta_se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ model, labeller = as_labeller(c(
    A_baseline = "lme4 analogue of Model A",
    A_family_slp = "+ (1 + year || Family)",
    A_src_site = "+ (1 | contributor) + (1 | municipality)"))) +
  scale_colour_manual(values = c(p10 = "#1F77B4", p50 = "#FF7F0E", p90 = "#D62728"),
                      labels = q_labels, name = "Diet-Inv percentile\n(species-weighted)") +
  scale_fill_manual(values = c(p10 = "#1F77B4", p50 = "#FF7F0E", p90 = "#D62728"),
                    labels = q_labels, name = "Diet-Inv percentile\n(species-weighted)") +
  labs(x = "Year", y = sprintf("Predicted change in wing length (mm) relative to %.1f", yr_mu),
       title = "Year x diet interaction evaluated at observed Diet-Inv percentiles",
       subtitle = "Fixed-effect slope predictions +/- 1.96 SE (lme4 REML, non-phylogenetic); support = species / records within +/-5 Diet-Inv points") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(fig_path("diet_quantile_predictions.png"), p_q, width = 11, height = 5.5, dpi = 150)
message("Saved: ", fig_path("diet_quantile_predictions.png"))

saveRDS(list(
  fixef = lme4_fixef, varcomp = lme4_varcomp, convergence = lme4_convergence,
  quantile_slopes = lme4_quantile_slopes, quantiles = diet_distribution$quantiles,
  quantile_support = quantile_support, predictions = pred_q,
  n_spp = n_distinct(dat_l$spp), n_rec = nrow(dat_l), n_family = n_distinct(dat_l$Family),
  n_src = n_distinct(dat_l$src), n_site = n_distinct(dat_l$site),
  year_scaling = c(center = yr_mu, sd = yr_sd),
  formulas = sapply(fml, function(f) paste(deparse(f), collapse = " ")),
  note = paste("lme4 REML analogues of the brms diet models, fitted on all species with",
               "EltonTraits data (no tree matching, so N species may exceed the brms fit by 1).",
               "A_family_int / A_family_slp add a family intercept / uncorrelated family year slope;",
               "A_src_site adds contributor (Main_researcher) and municipality intercepts (Phase 1 controls)."),
  session = session_info), out_path("diet_lme4_results.rds"))
message("Saved: ", out_path("diet_lme4_results.rds"))

if (!RUN_PHYLO) {
  message("\n--no-brms / --no-phylo: stopping before the phylogenetic fits. ",
          "diet_interaction_results.rds and diet_interaction_phylo.rds were NOT touched.")
  quit(save = "no", status = 0)
}

# ============================================================================
# glmmTMB engine (default): phylogenetic models across the 50 published trees
# ============================================================================
# Uses scripts/_phylo_engine.R (sourced AFTER the path helpers; it needs
# derived_path / out_path / raw_path). Nothing phylogenetic is re-implemented
# here: fit_phylo_glmmtmb() appends propto(0 + species_name | g, A) to each
# formula, run_phylo_trees() loops over the tree sample and Rubin-pools the fixed
# effects, tidy_phylo_fit() gives the variance components, species_slopes_phylo()
# the species-specific slopes. Derived quantities (year slope at the observed
# Diet-Inv percentiles and by diet category) are computed per tree from the
# fixed-effect vcov and pooled with the SAME Rubin rule (pool_rubin_df).
if (ENGINE == "glmmTMB") {
  source(script_path("_phylo_engine.R"))
  message("\n=== glmmTMB engine (propto phylogenetic term), ", N_TREES, " tree(s) ===")
  message("glmmTMB ", as.character(packageVersion("glmmTMB")))

  # -- data: tree tip names (eBird synonyms), factors, plain numerics ----------
  dat_g <- dat %>%
    mutate(species_name = phylo_species_name(Binomial),
           spp  = species_name,                        # non-phylogenetic species RE
           Sex  = factor(Sex, levels = c("Female", "Male")),
           src  = Main_researcher, site = Municipality,
           scaled_yr = as.numeric(scaled_yr), scaled_lat = as.numeric(scaled_lat)) %>%
    filter(species_name != "Herpsilochmus_sellowi")    # absent from the tree (engine note)
  stopifnot(!anyNA(dat_g$src), !anyNA(dat_g$site), !anyNA(dat_g$Family))
  spp_g  <- sort(unique(dat_g$species_name))
  A_list <- phylo_A_list(spp_g, n_trees = N_TREES)
  stopifnot(length(A_list) == N_TREES, all(spp_g %in% rownames(A_list[[1]])))
  message("Sample: ", nrow(dat_g), " wing records, ", length(spp_g), " species, ",
          n_distinct(dat_g$src), " contributors, ", n_distinct(dat_g$site), " municipalities, ",
          n_distinct(dat_g$Family), " families; ", length(A_list), " trees from the published cache")

  # -- formulas (the engine appends the phylogenetic term) ---------------------
  # Two species random-effect specifications are fitted for every model:
  #   "published": (1 + scaled_yr || spp)  iid species intercept + slope, as in the
  #                brms fits. The iid intercept is REDUNDANT with the phylogenetic
  #                intercept here: in the validated wing model its SD is 0.005 mm
  #                (glmmtmb_validation_wing.rds) and the phylogenetic term carries
  #                95 % of the variance. Without provenance terms the optimizer parks
  #                it at ~0 and converges; with contributor + municipality intercepts
  #                it lands on the mirror-image degenerate solution on most trees
  #                (iid SD 16.1, phylogenetic SD 0.09, singular Hessian, logLik NA),
  #                with fixed effects unchanged to 3 decimals.
  #   "reduced":   (0 + scaled_yr | spp)  species slope only; the species intercept
  #                variance is carried entirely by the phylogenetic term. Converges
  #                on every tree, identical log-likelihood and fixed effects where
  #                both converge. THIS IS THE PRIMARY TIER; the published structure
  #                is reported as a check, pooled over its converged trees only.
  fixed_part <- c(A = "conc.wing.length ~ Sex + scaled_yr * diet_inv_std + scaled_lat",
                  B = "conc.wing.length ~ Sex + scaled_yr * diet_cat2 + scaled_lat")
  extra_re <- c(A = "", A_fam = "+ (1 | Family)", A_src_site = "+ (1 | src) + (1 | site)",
                A_fam_src_site = "+ (1 | Family) + (1 | src) + (1 | site)",
                B = "", B_src_site = "+ (1 | src) + (1 | site)")
  spp_term <- c(reduced = "+ (0 + scaled_yr | spp)", published = "+ (1 + scaled_yr || spp)")
  make_fml <- function(mdl, spec) as.formula(paste(fixed_part[[substr(mdl, 1, 1)]], spp_term[[spec]], extra_re[[mdl]]))
  phylo_fml <- setNames(lapply(names(extra_re), make_fml, spec = "reduced"), names(extra_re))
  model_labels <- c(
    A = "Model A (phylogenetic)", A_fam = "A + (1 | Family)",
    A_src_site = "A + (1 | contributor) + (1 | municipality)",
    A_fam_src_site = "A + Family + contributor + municipality",
    B = "Model B (phylogenetic)", B_src_site = "B + (1 | contributor) + (1 | municipality)")
  int_term <- c(A = "scaled_yr:diet_inv_std", A_fam = "scaled_yr:diet_inv_std",
                A_src_site = "scaled_yr:diet_inv_std", A_fam_src_site = "scaled_yr:diet_inv_std",
                B = "scaled_yr:diet_cat2Invertebrate", B_src_site = "scaled_yr:diet_cat2Invertebrate")

  # -- derived quantity from one fit: slope of year at diet value z ------------
  # slope(z) = b_yr + b_int * z ; Var = v_yy + z^2 v_ii + 2 z v_yi (Wald; z = 0/1
  # for Model B gives the Other / Invertebrate slopes).
  slope_at <- function(fit, z, int_par) {
    b <- fixef(fit)$cond; V <- as.matrix(vcov(fit)$cond)
    est <- b["scaled_yr"] + b[int_par] * z
    se  <- sqrt(V["scaled_yr", "scaled_yr"] + z^2 * V[int_par, int_par] + 2 * z * V["scaled_yr", int_par])
    data.frame(z = z, estimate = unname(est), se = unname(se))
  }
  z_q <- z_of(q_spp)                                    # species-weighted p10 / p50 / p90

  # -- fit every model x species specification across the trees ---------------
  # Pooling uses ONLY trees with a positive-definite Hessian and finite SEs: the
  # non-converged fits are of two kinds (recorded in $per_tree_interaction), a
  # benign boundary case (a redundant SD at 0, fixed effects unchanged) and a
  # genuine bad optimum (fixed effects off, SEs tiny), and Rubin's rule must not
  # see the second. Per model the PRIMARY specification is the one with more
  # converged trees (ties -> reduced); both are saved.
  # NB: pool_rubin_df() returns a column named `m` (number of trees), so the loop
  # variable is `mdl` and is referenced as .env$mdl inside dplyr verbs.
  phylo_runs <- list()
  spp_re_lookup <- dat_g %>% distinct(spp, Binomial, Family, diet_inv, diet_cat2) %>%
    left_join(dat_g %>% count(spp, name = "n_records"), by = "spp")
  for (mdl in names(extra_re)) for (spec in names(spp_term)) {
    f <- make_fml(mdl, spec)
    message("\n--- ", mdl, ": ", model_labels[[mdl]], " [", spec, " spp spec] ---")
    run <- suppressWarnings(run_phylo_trees(f, dat_g, A_list, keep_fits = TRUE, verbose = TRUE))
    fin <- run$per_tree_fixed %>% group_by(tree) %>% summarise(finite = all(is.finite(se)), .groups = "drop")
    ok  <- run$varcomp$tree[run$varcomp$converged & run$varcomp$tree %in% fin$tree[fin$finite]]
    run$converged_trees <- ok; run$n_converged <- length(ok)
    run$naive_pool_all_trees <- run$pooled                     # naive pooling over every tree (diagnostic only)
    # exact [[ ]] indexing and a stored NULL: `$` would partial-match another field
    run["pooled"] <- list(if (length(ok)) pool_rubin_df(run$per_tree_fixed %>% filter(tree %in% ok)) else NULL)
    run$per_tree_interaction <- run$per_tree_fixed %>% filter(par == int_term[[mdl]]) %>%
      left_join(run$varcomp %>% select(tree, converged), by = "tree") %>%
      group_by(converged) %>% summarise(n = n(), est_min = min(estimate), est_max = max(estimate),
                                        se_min = min(se), se_max = max(se), .groups = "drop")
    # derived slopes (year slope at diet values) per converged tree -> Rubin-pooled
    is_A <- grepl("^A", mdl)
    zs   <- if (is_A) z_q else c(0, 1)
    labs <- if (is_A) paste0("p", probs * 100) else c("Other", "Invertebrate")
    run$derived <- if (length(ok)) {
      per_tree <- bind_rows(lapply(ok, function(i) {
        d <- slope_at(run$fits[[i]], zs, int_term[[mdl]])
        data.frame(tree = i, component = "derived", par = labs, estimate = d$estimate, se = d$se)
      }))
      pool_rubin_df(per_tree) %>%
        mutate(model = .env$mdl, spec = .env$spec, quantile = par, z = zs[match(par, labs)],
               diet_inv = if (is_A) q_spp[match(par, labs)] else NA_real_,
               mm_per_decade = estimate * 10 / yr_sd) %>%
        select(model, spec, quantile, diet_inv, diet_inv_std = z, slope_per_sd_year = estimate, se, lower, upper, m, mm_per_decade)
    } else NULL
    # species-specific year slopes (fixed + BLUP) from the first converged tree
    run$species_slopes <- if (mdl %in% c("A", "A_src_site") && length(ok)) {
      fit1 <- run$fits[[ok[1]]]
      sl <- species_slopes_phylo(fit1, "scaled_yr", "spp") %>%
        mutate(model = .env$mdl, spec = .env$spec, tree = ok[1]) %>% left_join(spp_re_lookup, by = "spp")
      # The engine's `slope` is the year MAIN effect + BLUP and ignores the
      # interaction, so it is orthogonal to diet by construction. Add each
      # species' own diet contribution b_int * z(Diet-Inv) (and its Wald
      # variance) so the point is the species-specific total year slope.
      fx <- slope_at(fit1, z_of(sl$diet_inv), int_term[[mdl]])
      sl %>% mutate(slope_main_plus_blup = slope, se_main_plus_blup = se_total,
                    slope = fx$estimate + ranef, se_total = sqrt(fx$se^2 + se_ranef^2),
                    lo = slope - 1.96 * se_total, hi = slope + 1.96 * se_total)
    } else NULL
    run$fits <- NULL; gc(verbose = FALSE)                # 50 fits x ~1.4 MB; drop before the next fit
    run$formula <- paste(deparse(f), collapse = " "); run$spec <- spec; run$label <- model_labels[[mdl]]
    phylo_runs[[mdl]][[spec]] <- run
    message(sprintf("%s [%s]: %d/%d trees converged (pdHess & finite SE), %.0f s; interaction over converged trees %s",
                    mdl, spec, run$n_converged, run$n_trees, run$secs,
                    if (!is.null(run[["pooled"]])) { r <- run[["pooled"]][run[["pooled"]]$par == int_term[[mdl]], ]
                      sprintf("%+.3f [%+.3f, %+.3f]", r$estimate, r$lower, r$upper) } else "NA"))
    print(as.data.frame(run$per_tree_interaction), digits = 4)
  }
  primary_spec <- sapply(names(extra_re), function(mdl) {
    nc <- sapply(phylo_runs[[mdl]], `[[`, "n_converged")
    if (nc[["published"]] > nc[["reduced"]]) "published" else "reduced" })
  primary <- setNames(lapply(names(extra_re), function(mdl) phylo_runs[[mdl]][[primary_spec[[mdl]]]]), names(extra_re))
  message("\nPrimary species specification per model (more converged trees): ",
          paste(names(primary_spec), primary_spec, sep = " = ", collapse = "; "))

  quantile_slopes_phylo <- bind_rows(lapply(primary[grepl("^A", names(primary))], `[[`, "derived"))
  category_slopes_phylo <- bind_rows(lapply(primary[grepl("^B", names(primary))], `[[`, "derived"))
  quantile_slopes_all_specs <- bind_rows(lapply(unlist(phylo_runs, recursive = FALSE), `[[`, "derived"))
  message("\nYear slope (mm per SD-year) at species-weighted Diet-Inv percentiles, glmmTMB phylogenetic, Rubin-pooled (primary spec):")
  print(as.data.frame(quantile_slopes_phylo %>% select(model, spec, quantile, diet_inv, slope_per_sd_year, lower, upper, m, mm_per_decade)), digits = 3)
  message("Year slope by diet category (Model B):")
  print(as.data.frame(category_slopes_phylo %>% select(model, spec, quantile, slope_per_sd_year, lower, upper, m, mm_per_decade)), digits = 3)

  # -- pooled fixed effects ----------------------------------------------------
  runs_flat <- unlist(phylo_runs, recursive = FALSE)     # names "A.reduced", ...
  fixed_pooled_all <- bind_rows(lapply(names(extra_re), function(mdl) bind_rows(lapply(names(spp_term), function(spec) {
    r <- phylo_runs[[mdl]][[spec]]; if (is.null(r[["pooled"]])) return(NULL)
    r[["pooled"]] %>% mutate(model = .env$mdl, label = model_labels[[.env$mdl]], spec = .env$spec,
                        primary = .env$spec == primary_spec[[.env$mdl]], n_trees = r$n_trees,
                        n_converged = r$n_converged, secs = r$secs) %>% relocate(model, label, spec, primary) }))))
  fixed_pooled <- fixed_pooled_all %>% filter(primary)
  make_int_table <- function(fp) fp %>%
    filter(component == "cond", par %in% c("scaled_yr", unname(int_term))) %>%
    mutate(term = ifelse(par == "scaled_yr", "year", "interaction")) %>%
    select(model, label, spec, primary, term, par, estimate, se, lower, upper, z, between_tree_var, n_trees, n_converged, secs)
  interaction_table <- make_int_table(fixed_pooled)
  interaction_table_all_specs <- make_int_table(fixed_pooled_all)
  message("\nInteraction and year main effect, glmmTMB phylogenetic, Rubin-pooled over converged trees (primary spec per model):")
  print(as.data.frame(interaction_table %>% select(model, spec, term, estimate, se, lower, upper, z, n_converged, n_trees, secs)), digits = 3)
  message("Same, every specification:")
  print(as.data.frame(interaction_table_all_specs %>% select(model, spec, primary, term, estimate, se, lower, upper, n_converged, n_trees)), digits = 3)

  # -- variance components (mean over converged trees) -------------------------
  varcomp_summary <- bind_rows(lapply(names(extra_re), function(mdl) bind_rows(lapply(names(spp_term), function(spec) {
    r <- phylo_runs[[mdl]][[spec]]; vc <- r$varcomp[r$varcomp$tree %in% r$converged_trees, , drop = FALSE]
    if (!nrow(vc)) return(NULL)
    vc %>% summarise(across(where(is.numeric) & !tree, mean), n_trees_used = n()) %>%
      mutate(model = .env$mdl, spec = .env$spec, primary = .env$spec == primary_spec[[.env$mdl]]) %>% relocate(model, spec, primary) }))))
  message("\nVariance components (SD; mean over converged trees):")
  print(as.data.frame(varcomp_summary), digits = 3)

  # -- comparison: published brms (10 trees) vs lme4 (no phylogeny) vs glmmTMB ---
  brms_pub <- if (file.exists(out_path("diet_interaction_results.rds")))
    readRDS(out_path("diet_interaction_results.rds")) else NULL   # read only; never rewritten here
  pick <- function(tab, par, model = NA, engine, label) {
    r <- tab[tab$par == par, ]; if (!is.na(model)) r <- r[r$model == model, ]
    if (!nrow(r)) return(NULL)
    data.frame(engine = engine, model = label, par = par, estimate = r$estimate[1], se = r$se[1],
               lower = r$lower[1], upper = r$upper[1])
  }
  comparison <- bind_rows(
    if (!is.null(brms_pub)) list(
      pick(brms_pub$model_A_pooled, "b_scaled_yr:diet_inv_std", engine = "brms (published, 10 trees)", label = "A"),
      pick(brms_pub$model_A_pooled, "b_scaled_yr", engine = "brms (published, 10 trees)", label = "A"),
      pick(brms_pub$model_B_pooled, "b_scaled_yr:diet_cat2Invertebrate", engine = "brms (published, 10 trees)", label = "B"),
      pick(brms_pub$model_B_pooled, "b_scaled_yr", engine = "brms (published, 10 trees)", label = "B")),
    # lme4 analogues: A_fam_src_site is matched to A_src_site_family_slp (family
    # SLOPE, the only lme4 spec with family + provenance terms)
    bind_rows(lapply(list(c("A", "A_baseline"), c("A_fam", "A_family_int"), c("A_src_site", "A_src_site"),
                          c("A_fam_src_site", "A_src_site_family_slp"), c("B", "B_baseline")), function(p)
      bind_rows(pick(lme4_fixef, "scaled_yr", p[2], "lme4 (no phylogeny)", p[1]),
                pick(lme4_fixef, int_term[[p[1]]], p[2], "lme4 (no phylogeny)", p[1])))),
    interaction_table_all_specs %>%
      transmute(engine = ifelse(primary, "glmmTMB phylogenetic (primary spp spec)", "glmmTMB phylogenetic (alternative spp spec)"),
                model, par, estimate, se, lower, upper, spec, n_converged, n_trees)
  ) %>% mutate(par = sub("^b_", "", par),
               term = ifelse(par == "scaled_yr", "year", "interaction")) %>%
    relocate(engine, model, term)
  # phylogeny vs no phylogeny on the same model: glmmTMB (primary spec) minus lme4
  phylo_vs_lme4 <- comparison %>% filter(grepl("^lme4|primary", engine)) %>%
    mutate(engine = ifelse(grepl("^lme4", engine), "lme4", "glmmTMB")) %>%
    select(engine, model, term, estimate, se) %>%
    pivot_wider(names_from = engine, values_from = c(estimate, se)) %>%
    mutate(diff_estimate = estimate_glmmTMB - estimate_lme4, se_ratio = se_glmmTMB / se_lme4) %>%
    filter(!is.na(estimate_lme4))
  message("\nComparison across engines (interaction and year main effect):")
  print(as.data.frame(comparison), digits = 3)
  message("\nPhylogeny (glmmTMB) minus no phylogeny (lme4), same model:")
  print(as.data.frame(phylo_vs_lme4), digits = 3)

  # ── figures (Analysis/figures only) ───────────────────────────────────────────
  engine_tag <- sprintf("glmmTMB propto phylogenetic term, %d tree%s, Rubin-pooled over converged trees; %d species / %s records",
                        N_TREES, ifelse(N_TREES > 1, "s", ""), length(spp_g), format(nrow(dat_g), big.mark = ","))
  conv_lab <- function(mdl, base) sprintf("%s\n%s spp spec, %d/%d trees", base, primary_spec[[mdl]],
                                          primary[[mdl]]$n_converged, primary[[mdl]]$n_trees)
  .fmt_ci <- function(tab, mdl, par) { r <- tab[tab$model == mdl & tab$par == par, ]
    sprintf("%+.3f [%+.3f, %+.3f]", r$estimate, r$lower, r$upper) }
  # trajectories relative to the record midpoint: delta = slope * scaled_yr,
  # SE = |scaled_yr| * SE(slope), from the POOLED slopes (slope is linear in year)
  traj <- function(sl, mdl) {
    s <- sl %>% filter(model == mdl)
    expand.grid(scaled_yr = yr_seq, quantile = s$quantile, stringsAsFactors = FALSE) %>%
      left_join(s, by = "quantile") %>%
      mutate(delta = slope_per_sd_year * scaled_yr, delta_se = abs(scaled_yr) * se,
             year = scaled_yr * yr_sd + yr_mu, model = mdl)
  }
  pred_A_phylo <- bind_rows(lapply(c("A", "A_fam", "A_src_site"), function(m) traj(quantile_slopes_phylo, m))) %>%
    mutate(quantile = factor(quantile, levels = paste0("p", probs * 100)))
  pred_B_phylo <- bind_rows(lapply(c("B", "B_src_site"), function(m) traj(category_slopes_phylo, m))) %>%
    mutate(quantile = factor(quantile, levels = c("Other", "Invertebrate")))
  q_cols <- c(p10 = "#1F77B4", p50 = "#FF7F0E", p90 = "#D62728")

  # (1) manuscript-style figure: Model A at percentiles / Model B by category
  p_A_g <- ggplot(pred_A_phylo %>% filter(model == "A"), aes(x = year, y = delta, colour = quantile, fill = quantile)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_ribbon(aes(ymin = delta - 1.96 * delta_se, ymax = delta + 1.96 * delta_se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = q_cols, labels = q_labels, name = "Diet-Inv percentile (species-weighted)") +
    scale_fill_manual(values = q_cols, labels = q_labels, name = "Diet-Inv percentile (species-weighted)") +
    labs(x = "Year", y = sprintf("Change in wing length (mm) relative to %.1f", yr_mu),
         title = "Year x invertebrate diet proportion (Model A)",
         subtitle = paste0("interaction = ", .fmt_ci(interaction_table, "A", "scaled_yr:diet_inv_std"),
                           " mm per SD-year per SD Diet-Inv\nwith (1 | Family): ",
                           .fmt_ci(interaction_table, "A_fam", "scaled_yr:diet_inv_std"),
                           "\nwith contributor + municipality intercepts: ",
                           .fmt_ci(interaction_table, "A_src_site", "scaled_yr:diet_inv_std"))) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")
  p_B_g <- ggplot(pred_B_phylo, aes(x = year, y = delta, colour = quantile, fill = quantile, linetype = model)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_ribbon(data = pred_B_phylo %>% filter(model == "B"),
                aes(ymin = delta - 1.96 * delta_se, ymax = delta + 1.96 * delta_se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = c(Other = "#1F77B4", Invertebrate = "#D62728"), name = "Diet category (EltonTraits Diet-5Cat)") +
    scale_fill_manual(values = c(Other = "#1F77B4", Invertebrate = "#D62728"), name = "Diet category (EltonTraits Diet-5Cat)") +
    scale_linetype_manual(values = c(B = "solid", B_src_site = "dashed"),
                          labels = c(B = "Model B", B_src_site = "B + contributor + municipality"), name = NULL) +
    labs(x = "Year", y = sprintf("Change in wing length (mm) relative to %.1f", yr_mu),
         title = "Year x diet category (Model B)",
         subtitle = paste0("interaction = ", .fmt_ci(interaction_table, "B", "scaled_yr:diet_cat2Invertebrate"),
                           "; with contributor + municipality: ",
                           .fmt_ci(interaction_table, "B_src_site", "scaled_yr:diet_cat2Invertebrate"),
                           "\nribbons: Model B only")) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")
  p_comb_g <- (p_A_g / p_B_g) + plot_annotation(caption = engine_tag)
  fig_main <- fig_path(if (SMOKE) "diet_interaction_plot_SMOKE.png" else "diet_interaction_plot.png")
  ggsave(fig_main, p_comb_g, width = 7.5, height = 11, dpi = 150)
  message("Saved: ", fig_main)

  # (2) percentile trajectories across the three continuous specifications
  p_q_g <- ggplot(pred_A_phylo, aes(x = year, y = delta, colour = quantile, fill = quantile)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_ribbon(aes(ymin = delta - 1.96 * delta_se, ymax = delta + 1.96 * delta_se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~ model, labeller = as_labeller(c(A = conv_lab("A", "Model A (phylogenetic)"), A_fam = conv_lab("A_fam", "+ (1 | Family)"),
                                                 A_src_site = conv_lab("A_src_site", "+ (1 | contributor) + (1 | municipality)")))) +
    scale_colour_manual(values = q_cols, labels = q_labels, name = "Diet-Inv percentile\n(species-weighted)") +
    scale_fill_manual(values = q_cols, labels = q_labels, name = "Diet-Inv percentile\n(species-weighted)") +
    labs(x = "Year", y = sprintf("Predicted change in wing length (mm) relative to %.1f", yr_mu),
         title = "Year x diet interaction at observed Diet-Inv percentiles, phylogenetic models",
         subtitle = paste0(engine_tag, "\nribbons = Rubin-pooled slope uncertainty (+/- 1.96 SE); facet labels give the species-RE spec and converged trees")) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave(fig_path(if (SMOKE) "diet_quantile_predictions_phylo_SMOKE.png" else "diet_quantile_predictions_phylo.png"),
         p_q_g, width = 11, height = 5.5, dpi = 150)
  message("Saved: ", fig_path("diet_quantile_predictions_phylo.png"))

  # (3) the interaction coefficient across engines and control sets
  spec_tab <- comparison %>% filter(term == "interaction") %>%
    mutate(diet_form = ifelse(grepl("^A", model), "continuous (per SD Diet-Inv)", "categorical (Invertebrate - Other)"),
           spec = ifelse(grepl("^glmmTMB", engine),
                         sprintf("%s  |  %s [%s spp RE, %d/%d trees]", model, engine, spec, n_converged, n_trees),
                         paste0(model, "  |  ", engine))) %>%
    arrange(diet_form, model, engine)
  spec_tab$spec <- factor(spec_tab$spec, levels = rev(unique(spec_tab$spec)))
  p_spec <- ggplot(spec_tab, aes(x = estimate, y = spec, colour = engine)) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, orientation = "y") +
    geom_point(size = 2.4) +
    facet_wrap(~ diet_form, ncol = 1, scales = "free") +
    scale_colour_manual(values = c("#7570B3", "#1B9E77", "#D95F02", "#E7298A"), name = NULL) +
    labs(x = "Year x diet interaction (mm per SD-year per unit diet), 95% interval", y = NULL,
         title = "Diet x year interaction across engines and provenance controls",
         subtitle = "brms = published 10-tree fit; lme4 = REML, no phylogeny;\nglmmTMB = phylogenetic (propto), Rubin-pooled over converged trees [species-RE spec, converged/total trees]") +
    theme_classic(base_size = 11) + theme(legend.position = "bottom", plot.title.position = "plot")
  ggsave(fig_path(if (SMOKE) "diet_interaction_specifications_SMOKE.png" else "diet_interaction_specifications.png"),
         p_spec, width = 14, height = 7, dpi = 150)
  message("Saved: ", fig_path("diet_interaction_specifications.png"))

  # (4) species-specific year slopes (fixed + BLUP, tree 1) against Diet-Inv
  ss <- bind_rows(lapply(primary[c("A", "A_src_site")], `[[`, "species_slopes"))
  p_ss <- ggplot(ss, aes(x = diet_inv, y = slope * 10 / yr_sd, colour = Family)) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_errorbar(aes(ymin = lo * 10 / yr_sd, ymax = hi * 10 / yr_sd), width = 0, alpha = 0.35) +
    geom_point(aes(size = n_records), alpha = 0.85) +
    geom_abline(data = quantile_slopes_phylo %>% filter(model %in% c("A", "A_src_site"), quantile == "p50") %>%
                  left_join(fixed_pooled %>% filter(component == "cond", par == "scaled_yr:diet_inv_std") %>% select(model, b_int = estimate), by = "model") %>%
                  mutate(slope = b_int * 10 / yr_sd / diet_scale,
                         intercept = (slope_per_sd_year * 10 / yr_sd) - slope * diet_inv),
                aes(slope = slope, intercept = intercept), colour = "grey30") +
    facet_wrap(~ model, labeller = as_labeller(c(A = conv_lab("A", "Model A (phylogenetic)"), A_src_site = conv_lab("A_src_site", "A + contributor + municipality")))) +
    scale_size_area(max_size = 6, name = "Wing records") +
    labs(x = "Diet-Inv (EltonTraits, % invertebrates)", y = "Species year slope (mm per decade), fixed + BLUP",
         title = "Species-specific wing-length trends against the diet covariate",
         subtitle = sprintf("Points: species slopes from one tree (tree %s; +/- 1.96 SE incl. BLUP conditional SD); line: pooled fixed-effect interaction",
                            paste(unique(ss$tree), collapse = "/"))) +
    theme_classic(base_size = 10) + theme(legend.position = "right")
  ggsave(fig_path(if (SMOKE) "diet_species_slopes_SMOKE.png" else "diet_species_slopes.png"), p_ss, width = 12, height = 6, dpi = 150)
  message("Saved: ", fig_path("diet_species_slopes.png"))

  # ── save ──────────────────────────────────────────────────────────────────────
  phylo_results <- list(
    generated = Sys.time(), engine = "glmmTMB", glmmTMB_version = as.character(packageVersion("glmmTMB")),
    n_trees = N_TREES, tree_source = "data/derived/phylo_A_50trees.rds (the 50 trees of brm0_multiphylo.rda)",
    smoke_test = SMOKE,
    n_spp = length(spp_g), n_rec = nrow(dat_g), n_src = n_distinct(dat_g$src), n_site = n_distinct(dat_g$site),
    n_family = n_distinct(dat_g$Family), species = spp_g,
    species_re_spec = list(
      reduced   = "(0 + scaled_yr | spp) + propto phylogenetic species intercept (species intercept variance carried by the phylogenetic term)",
      published = "(1 + scaled_yr || spp) + propto phylogenetic species intercept (iid species intercept redundant with the phylogenetic one; SD ~ 0 when it converges)",
      rule      = "primary spec per model = the one with more trees converged (pdHess TRUE and finite SEs); ties -> reduced. Pooling always over converged trees only.",
      primary_spec = primary_spec),
    models = lapply(phylo_runs, function(by_spec) lapply(by_spec, function(r)
      r[c("label", "spec", "formula", "pooled", "naive_pool_all_trees", "per_tree_fixed", "per_tree_interaction", "varcomp",
          "converged_trees", "n_trees", "n_converged", "secs", "derived")])),
    fixed_pooled = fixed_pooled,                       # primary spec per model
    fixed_pooled_all_specs = fixed_pooled_all,
    interaction_table = interaction_table,             # primary spec per model
    interaction_table_all_specs = interaction_table_all_specs,
    varcomp_summary = varcomp_summary,
    quantile_slopes = quantile_slopes_phylo,           # primary spec per model
    quantile_slopes_all_specs = quantile_slopes_all_specs,
    category_slopes = category_slopes_phylo,
    quantiles = diet_distribution$quantiles, quantile_support = quantile_support,
    comparison = comparison, phylo_vs_lme4 = phylo_vs_lme4,
    brms_published = if (!is.null(brms_pub)) brms_pub[c("model_A_pooled", "model_B_pooled", "n_spp", "n_rec")] else NULL,
    predictions_A = pred_A_phylo, predictions_B = pred_B_phylo,
    species_slopes_first_tree = ss,
    year_scaling = c(center = yr_mu, sd = yr_sd),
    diet_standardisation = c(center = diet_center, sd = diet_scale),
    total_secs = sum(sapply(runs_flat, `[[`, "secs")),
    note = paste("Phylogenetic diet-interaction models fitted with glmmTMB propto (scripts/_phylo_engine.R)",
                 "across the tree sample of the published brms analysis and Rubin-pooled over the trees whose",
                 "Hessian is positive definite with finite SEs (species_re_spec$rule).",
                 "Intervals are Wald +/- 1.96 SE with the Rubin between-tree variance added.",
                 "quantile_slopes / category_slopes are derived per tree from the fixed-effect vcov and pooled",
                 "with the same rule. species_slopes_first_tree conditions on one tree (BLUPs).",
                 "diet_interaction_results.rds (published 10-tree brms) is read for the comparison and never rewritten."),
    session = session_info
  )
  out_phylo <- out_path(if (SMOKE) "diet_interaction_phylo_SMOKE.rds" else "diet_interaction_phylo.rds")
  saveRDS(phylo_results, out_phylo)
  message("Saved: ", out_phylo)

  message("\n=== Key results (glmmTMB phylogenetic, ", N_TREES, " trees; primary spp spec per model) ===")
  for (mdl in names(primary)) message(sprintf("%-16s [%-9s] interaction %s | year %s | %d/%d trees converged | %.0f s  [other spec: %d/%d]",
    mdl, primary_spec[[mdl]], .fmt_ci(interaction_table, mdl, int_term[[mdl]]), .fmt_ci(interaction_table, mdl, "scaled_yr"),
    primary[[mdl]]$n_converged, primary[[mdl]]$n_trees, primary[[mdl]]$secs,
    phylo_runs[[mdl]][[setdiff(names(spp_term), primary_spec[[mdl]])]]$n_converged, N_TREES))
  message("Total glmmTMB wall time (both specs): ", round(phylo_results$total_secs), " s")
  quit(save = "no", status = 0)
}

# ============================================================================
# brms engine (--engine brms; Totoro): phylogeny, Models A, B, A_fam, quantile
# slopes, figure. Kept as the Bayesian cross-check of the glmmTMB engine.
# ============================================================================
if (SMOKE) {
  message("\n*** SMOKE TEST: 1 tree, 2 chains, 400 iterations. Nothing below is a result. ***")
  SAMPLING$chains <- 2L; SAMPLING$iter <- 400L; SAMPLING$warmup <- 200L
  SAMPLING$workers <- 1L
}

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

# ── Rubin pooling helpers ────────────────────────────────────────────────────
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
# Same rule applied to a derived quantity computed from the draws of each tree
# (list of numeric matrices, draws x quantities, one per tree).
pool_rubin_derived <- function(mats, names) {
  m     <- length(mats)
  means <- t(sapply(mats, colMeans))
  vars  <- t(sapply(mats, function(x) apply(x, 2, var)))
  if (m == 1L) { means <- matrix(means, nrow = 1); vars <- matrix(vars, nrow = 1) }
  qbar <- colMeans(means); ubar <- colMeans(vars)
  b    <- if (m > 1L) apply(means, 2, var) else 0
  se   <- sqrt(ubar + (1 + 1/m) * b)
  data.frame(par = names, estimate = qbar, se = se,
             lower = qbar - 1.96*se, upper = qbar + 1.96*se, row.names = NULL)
}

fit_formula <- function(formula, tree) {
  A <- make_A(tree)
  brm(
    formula, data = dat, data2 = list(A = A),
    family  = gaussian(), prior = priors,
    iter = SAMPLING$iter, warmup = SAMPLING$warmup, chains = SAMPLING$chains,
    cores = SAMPLING$cores, backend = SAMPLING$backend,
    control = SAMPLING_CONTROL,
    seed    = SEED, silent = 2, refresh = 0
  )
}

# ============================================================================
# MODEL A: continuous interaction  scaled_yr × diet_inv_std   (unchanged)
# ============================================================================
message("\n=== Model A: scaled_yr × diet_inv_std (continuous) ===")
fit_A <- function(tree) fit_formula(
  conc.wing.length ~ 1 + Sex + scaled_yr * diet_inv_std + scaled_lat +
    (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)), tree)

plan(multisession, workers = SAMPLING$workers)
fits_A <- future_lapply(tree_samp, fit_A, future.seed = TRUE)

pars_A <- c("b_Intercept", "b_SexMale", "b_scaled_yr",
            "b_diet_inv_std", "b_scaled_lat", "b_scaled_yr:diet_inv_std")
rubin_A <- pool_rubin(fits_A, pars_A)
message("\nModel A pooled results:")
print(rubin_A)

# ============================================================================
# MODEL B: categorical interaction  scaled_yr × diet_cat2     (unchanged)
# ============================================================================
message("\n=== Model B: scaled_yr × diet_cat2 (Invertebrate vs. Other) ===")
fit_B <- function(tree) fit_formula(
  conc.wing.length ~ 1 + Sex + scaled_yr * diet_cat2 + scaled_lat +
    (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)), tree)

fits_B <- future_lapply(tree_samp, fit_B, future.seed = TRUE)

pars_B <- c("b_Intercept", "b_SexMale", "b_scaled_yr",
            "b_diet_cat2Invertebrate", "b_scaled_lat",
            "b_scaled_yr:diet_cat2Invertebrate")
rubin_B <- pool_rubin(fits_B, pars_B)
message("\nModel B pooled results:")
print(rubin_B)

# ============================================================================
# MODEL A_fam: Model A + (1 | Family)   (clade-confounding sensitivity, new)
# ============================================================================
message("\n=== Model A_fam: Model A + (1 | Family) ===")
fit_A_fam <- function(tree) fit_formula(
  conc.wing.length ~ 1 + Sex + scaled_yr * diet_inv_std + scaled_lat +
    (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)) + (1 | Family), tree)

fits_A_fam <- future_lapply(tree_samp, fit_A_fam, future.seed = TRUE)
rubin_A_fam <- pool_rubin(fits_A_fam, pars_A)
message("\nModel A_fam pooled results:")
print(rubin_A_fam)

# ============================================================================
# Year slope at the observed Diet-Inv percentiles (posterior, Rubin-pooled)
# ============================================================================
z_q <- z_of(q_spp)
slope_draws <- function(fits) lapply(fits, function(f) {
  d <- as_draws_df(f)
  sapply(z_q, function(z) d[["b_scaled_yr"]] + d[["b_scaled_yr:diet_inv_std"]] * z)
})
quantile_slopes <- bind_rows(
  cbind(model = "A", pool_rubin_derived(slope_draws(fits_A), paste0("p", probs * 100))),
  cbind(model = "A_fam", pool_rubin_derived(slope_draws(fits_A_fam), paste0("p", probs * 100)))) %>%
  rename(quantile = par, slope_per_sd_year = estimate) %>%
  mutate(diet_inv = rep(q_spp, 2), diet_inv_std = rep(z_q, 2),
         mm_per_decade = slope_per_sd_year * 10 / yr_sd)
message("\nYear slope at species-weighted Diet-Inv percentiles (brms, Rubin-pooled):")
print(quantile_slopes)

# ============================================================================
# Plots
# ============================================================================
# -- Model A: predicted year trend at the observed 10/50/90th percentiles ------
# (replaces the +/-1 SD projection, which fell on 11 species / 1 clade group)
fit_A_rep <- fits_A[[1]]      # one representative tree for the visual only
pred_A <- expand.grid(
  scaled_yr     = yr_seq,
  diet_inv      = q_spp,
  Sex           = "Female",
  scaled_lat    = 0,
  spp           = NA_character_,
  species_name  = NA_character_,
  Family        = NA_character_
)
pred_A$diet_inv_std <- z_of(pred_A$diet_inv)
pe <- fitted(fit_A_rep, newdata = pred_A, allow_new_levels = TRUE, re_formula = NA)
pred_A$fitted <- pe[, "Estimate"]; pred_A$lower <- pe[, "Q2.5"]; pred_A$upper <- pe[, "Q97.5"]
pred_A$quantile <- factor(paste0("p", probs * 100)[match(pred_A$diet_inv, q_spp)],
                          levels = paste0("p", probs * 100))
pred_A$year <- pred_A$scaled_yr * yr_sd + yr_mu

.fmt_int <- function(tab, par) sprintf("%.3f [%.3f, %.3f]", tab$estimate[tab$par == par],
                                       tab$lower[tab$par == par], tab$upper[tab$par == par])
p_A <- ggplot(pred_A, aes(x = year, y = fitted, colour = quantile, fill = quantile)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(p10 = "#1F77B4", p50 = "#FF7F0E", p90 = "#D62728"),
                      labels = q_labels, name = "Diet-Inv percentile") +
  scale_fill_manual(values = c(p10 = "#1F77B4", p50 = "#FF7F0E", p90 = "#D62728"),
                    labels = q_labels, name = "Diet-Inv percentile") +
  labs(x = "Year", y = "Predicted wing length (mm)",
       title = "Interaction: year × invertebrate diet proportion",
       subtitle = paste0("Model A (continuous); interaction β = ",
                         .fmt_int(rubin_A, "b_scaled_yr:diet_inv_std"),
                         "; with (1 | Family): ", .fmt_int(rubin_A_fam, "b_scaled_yr:diet_inv_std"))) +
  theme_classic() + theme(legend.position = "bottom", legend.direction = "vertical")

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
peB <- fitted(fit_B_rep, newdata = pred_B, allow_new_levels = TRUE, re_formula = NA)
pred_B$fitted <- peB[, "Estimate"]; pred_B$lower <- peB[, "Q2.5"]; pred_B$upper <- peB[, "Q97.5"]
pred_B$year <- pred_B$scaled_yr * yr_sd + yr_mu

p_B <- ggplot(pred_B, aes(x = year, y = fitted, colour = diet_cat2, fill = diet_cat2)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Other" = "#1F77B4", "Invertebrate" = "#D62728"),
                      name = "Diet category") +
  scale_fill_manual(values = c("Other" = "#1F77B4", "Invertebrate" = "#D62728"),
                    name = "Diet category") +
  labs(x = "Year",
       y = "Predicted wing length (mm)",
       title = "Interaction: year × diet category",
       subtitle = paste0("Model B (categorical); interaction β = ",
                         .fmt_int(rubin_B, "b_scaled_yr:diet_cat2Invertebrate"))) +
  theme_classic()

# combine and save. Written to Analysis/figures only (2026-09): the manuscript
# figure is produced by P6's make_figures pipeline from the files on disk, and
# the glmmTMB engine owns figures/diet_interaction_plot.png; the brms version
# is the cross-check.
p_combined <- p_A / p_B
ggsave(fig_path(if (SMOKE) "diet_interaction_plot_brms_SMOKE.png" else "diet_interaction_plot_brms.png"),
       p_combined, width = 7, height = 10, dpi = 150)
message("Saved: ", fig_path("diet_interaction_plot_brms.png"))

# ============================================================================
# Save all results (existing fields kept; new fields appended)
# ============================================================================
results <- list(
  model_A_pooled = rubin_A,
  model_B_pooled = rubin_B,
  n_spp          = n_distinct(dat$spp),
  n_rec          = nrow(dat),
  diet_summary   = dat %>% distinct(spp, diet_inv, diet_cat2) %>%
                     group_by(diet_cat2) %>%
                     summarise(n = n(), mean_inv = mean(diet_inv), .groups = "drop"),
  # --- new (2026-09) ---
  model_A_family_pooled = rubin_A_fam,
  quantile_slopes       = quantile_slopes,
  quantile_support      = quantile_support,
  quantiles             = diet_distribution$quantiles,
  n_trees               = N_TREES,
  smoke_test            = SMOKE,
  session               = session_info
)
out_file <- out_path(if (SMOKE) "diet_interaction_results_SMOKE.rds" else "diet_interaction_results.rds")
saveRDS(results, out_file)
message("Saved: ", out_file)

message("\n=== Key results ===")
message("Model A — scaled_yr:diet_inv_std interaction: ", .fmt_int(rubin_A, "b_scaled_yr:diet_inv_std"))
message("Model A + (1|Family) — interaction:           ", .fmt_int(rubin_A_fam, "b_scaled_yr:diet_inv_std"))
message("Model B — scaled_yr:Invertebrate interaction:  ", .fmt_int(rubin_B, "b_scaled_yr:diet_cat2Invertebrate"))
message("\nInterpretation: a negative interaction term means insectivorous species")
message("show a stronger (more negative) year → wing-length trend, consistent")
message("with food-limitation driving body-size decline. Read the slopes at the")
message("observed percentiles (quantile_slopes) rather than at +/-1 SD.")
