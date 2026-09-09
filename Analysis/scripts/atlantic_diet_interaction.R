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
#       pieces run in seconds without Stan or the tree cloud.
#
# RATIONALE (unchanged): if food limitation (declining arthropod prey) drives
# morphological change, the year effect on wing length should be more negative
# in species with a higher proportion of invertebrates in the diet.
#
# MODELS (brms, N_TREES trees, Rubin-pooled; lme4 analogues with the same
# fixed effects and species random intercept + uncorrelated slope):
#   Model A (continuous):  wing ~ Sex + scaled_yr * diet_inv_std + scaled_lat
#                          + (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A))
#   Model B (categorical): as A with diet_cat2 (Invertebrate vs Other) in
#                          place of diet_inv_std
#   Model A_fam:           Model A + (1 | Family)          (clade sensitivity)
#   lme4 only:             Model A + (1 + scaled_yr || Family)
#                          Model A + (1 | Main_researcher) + (1 | Municipality)
#                          (provenance controls of Phase 1; not phylogenetic)
#   diet_inv_std = z-scored Diet-Inv (EltonTraits, 0-100), standardised on the
#   RECORD-level analytical sample exactly as before (record-weighted mean/SD),
#   so the brms coefficients stay comparable with the published fit.
#
# INPUTS:  data/derived/passer90.rda            (live-only, known sex, 73 spp)
#          data/raw/BirdFuncDat.txt             (EltonTraits 1.0)
#          data/raw/AvesDataLite-main/          (clootl tree cloud; brms path only)
#          scripts/_sampling_config.R           (SAMPLING, SAMPLING_CONTROL)
# OUTPUTS: output/diet_distribution.rds         species table, bin counts,
#                                               quantiles and their support
#          figures/diet_distribution.png        species and records per bin
#          output/diet_lme4_results.rds         lme4 analogues + quantile slopes
#          figures/diet_quantile_predictions.png lme4 predicted trajectories at
#                                               the 10/50/90th percentiles
#          --- brms path only (unchanged files, extended contents) ---
#          output/diet_interaction_results.rds  model_A_pooled, model_B_pooled,
#                                               n_spp, n_rec, diet_summary (as
#                                               before) + model_A_family_pooled,
#                                               quantile_slopes, quantile_support
#          ../Manuscript/images/diet_interaction_plot.png (quantile version)
#
# RUN:  Rscript Analysis/scripts/atlantic_diet_interaction.R              # full brms (Totoro)
#       Rscript Analysis/scripts/atlantic_diet_interaction.R --no-brms    # descriptive + lme4 only
#       Rscript Analysis/scripts/atlantic_diet_interaction.R --trees 10   # brms with 10 trees
#       Rscript Analysis/scripts/atlantic_diet_interaction.R --smoke      # 1 tree, 2 chains,
#                                                                         # 400 iter: a SMOKE
#                                                                         # TEST, not a result
#       The --no-brms path never touches diet_interaction_results.rds or the
#       manuscript figure, so the published brms numbers survive a local run.
#
# Session (2026-09-09, --no-brms path executed locally): R 4.6.0, dplyr 1.2.1,
# tidyr 1.3.2, readr 2.2.0, lme4 2.0.1, ggplot2 4.0.3, patchwork 1.3.2.
# brms path last run with brms 2.23.0 / prepR4pcm 0.5.0.9000 / clootl 0.1.4
# (10 trees, server). Versions are recorded in the output rds files.
# ---------------------------------------------------------------------------

# ── command-line flags ───────────────────────────────────────────────────────
.args    <- commandArgs(trailingOnly = TRUE)
RUN_BRMS <- !("--no-brms" %in% .args)
SMOKE    <- "--smoke" %in% .args
N_TREES  <- 10L
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
  packages = sapply(c("dplyr", "tidyr", "lme4", "ggplot2", "patchwork",
                      "brms", "prepR4pcm", "clootl", "posterior"), .pkg_ver),
  flags = list(run_brms = RUN_BRMS, smoke = SMOKE, n_trees = N_TREES)
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

if (!RUN_BRMS) {
  message("\n--no-brms: stopping before the phylogenetic brms fits. ",
          "diet_interaction_results.rds and Manuscript/images/diet_interaction_plot.png were NOT touched.")
  quit(save = "no", status = 0)
}

# ============================================================================
# brms path (server): phylogeny, Models A, B, A_fam, quantile slopes, figure
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

# combine and save (the manuscript figure; Analysis/figures gets a copy)
p_combined <- p_A / p_B
img_dir <- file.path(dirname(ANALYSIS_DIR), "Manuscript", "images")
if (!SMOKE) {
  ggsave(file.path(img_dir, "diet_interaction_plot.png"), p_combined,
         width = 7, height = 10, dpi = 150)
  message("Saved: Manuscript/images/diet_interaction_plot.png")
}
ggsave(fig_path(if (SMOKE) "diet_interaction_plot_SMOKE.png" else "diet_interaction_plot.png"),
       p_combined, width = 7, height = 10, dpi = 150)

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
