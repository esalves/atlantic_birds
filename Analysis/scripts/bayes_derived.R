# bayes_derived.R
# ---------------------------------------------------------------------------
# Quantities the manuscript quotes that are functions of several parameters,
# computed draw by draw from the Stan tier (mixture posterior over trees):
#   - species-specific wing year slopes under "Wing: fully adjusted"
#     (fixed year slope + species deviation), with counts by direction
#   - wing year slope at the p10/p50/p90 invertebrate-diet percentiles under
#     "Wing x diet: + contributor & municipality"
# Needs the per-tree fits in output/bayes/fits/ (not tracked in git).
#
# RUN:  Rscript Analysis/scripts/bayes_derived.R      (after bayes_aggregate.R)
# OUT:  output/bayes/bayes_derived.rds
# ---------------------------------------------------------------------------
suppressMessages({ library(dplyr) })
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/")
  repeat {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived"))) return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/.")
}
ANALYSIS_DIR <- .find_analysis_dir()
out_path <- function(...) file.path(ANALYSIS_DIR, "output", ...)

# per-decade factor: scaled_yr is calendar year / its SD
SD_YR <- as.numeric(attr(readRDS(file.path(ANALYSIS_DIR, "data", "derived", "passer90_climate.rds"))$scaled_yr, "scaled:scale"))
DEC <- 10 / SD_YR

pooled_draws <- function(fit_id) {
  files <- list.files(out_path("bayes", "fits", fit_id), pattern = "^tree_[0-9]+\\.rds$", full.names = TRUE)
  if (!length(files)) return(NULL)
  do.call(rbind, lapply(files, function(f) readRDS(f)$draws))
}
summ <- function(x) c(estimate = mean(x), lower = unname(quantile(x, 0.025)), upper = unname(quantile(x, 0.975)),
                      p_neg = mean(x < 0))

out <- list(generated = Sys.time(), sd_yr = SD_YR, dec = DEC)

# species slopes, fully adjusted wing model (mm per decade)
dr <- pooled_draws("atlantic_parallel_controlled__011__conc.wing.length__M3")
if (!is.null(dr)) {
  dev <- grep("^slope_dev_", colnames(dr), value = TRUE)
  sl <- (dr[, "b_scaled_yr"] + dr[, dev, drop = FALSE]) * DEC
  sp <- as.data.frame(t(apply(sl, 2, summ)))
  sp$species <- sub("^slope_dev_", "", dev); rownames(sp) <- NULL
  out$species_slopes <- sp
  out$species_counts <- with(sp, c(n = nrow(sp), neg_point = sum(estimate < 0),
                                   cri_spans_zero = sum(lower < 0 & upper > 0),
                                   cri_neg = sum(upper < 0), cri_pos = sum(lower > 0),
                                   median_mm_decade = median(estimate)))
  out$community_mm_decade <- summ(dr[, "b_scaled_yr"] * DEC)
}

# diet percentiles (diet_inv_std values of the p10/p50/p90 species-weighted percentiles)
dr <- pooled_draws("atlantic_diet_interaction__006__conc.wing.length")
if (!is.null(dr)) {
  q <- readRDS(out_path("diet_interaction_phylo.rds"))$quantile_slopes_all_specs
  q <- q[q$model == "A_src_site" & q$spec == "reduced", c("quantile", "diet_inv", "diet_inv_std")]
  out$diet_quantile_slopes <- bind_rows(lapply(seq_len(nrow(q)), function(i) {
    s <- (dr[, "b_scaled_yr"] + dr[, "b_scaled_yr:diet_inv_std"] * q$diet_inv_std[i]) * DEC
    data.frame(quantile = q$quantile[i], diet_inv = q$diet_inv[i], t(summ(s)))
  }))
}

# species slopes, baseline wing model (mm per decade): the hollow comparison
# markers of the species-slope figure
dr <- pooled_draws("atlantic_parallel_controlled__001__conc.wing.length__M0")
if (!is.null(dr)) {
  dev <- grep("^slope_dev_", colnames(dr), value = TRUE)
  sl <- (dr[, "b_scaled_yr"] + dr[, dev, drop = FALSE]) * DEC
  sp <- as.data.frame(t(apply(sl, 2, summ)))
  sp$species <- sub("^slope_dev_", "", dev); rownames(sp) <- NULL
  out$species_slopes_M0 <- sp
  out$community_mm_decade_M0 <- summ(dr[, "b_scaled_yr"] * DEC)
}

# wing year slope by diet category (Other vs Invertebrate) under
# "Wing x diet category: + contributor & municipality"
dr <- pooled_draws("atlantic_diet_interaction__012__conc.wing.length")
if (!is.null(dr)) {
  out$diet_category_slopes <- rbind(
    data.frame(category = "Other", t(summ(dr[, "b_scaled_yr"] * DEC))),
    data.frame(category = "Invertebrate",
               t(summ((dr[, "b_scaled_yr"] + dr[, "b_scaled_yr:diet_cat2Invertebrate"]) * DEC))))
}

# sampler diagnostics per tree for the coefficients the manuscript quotes
# (year, anomaly and interaction terms, and the phylogenetic share and
# correlation of the species slopes) and, separately, for the species
# intercept SD, which mixes most slowly. rhat / ess are per tree (4 chains);
# ess_pooled sums the per-tree ESS over trees, since the trees are independent runs.
focal <- c("b_scaled_yr", "b_yr_within_src", "b_yr_src_mean", "b_dT", "b_dT_detr",
           "b_scaled_yr:diet_inv_std", "b_scaled_yr:diet_cat2Invertebrate",
           "sigma_scaled_yr", "sigma_scaled_tmean", "h2_slope", "cor_phylo_int_slope")
diag_rows <- list()
for (fid in list.files(out_path("bayes", "fits"))) {
  for (f in list.files(out_path("bayes", "fits", fid), pattern = "^tree_[0-9]+\\.rds$", full.names = TRUE)) {
    x <- readRDS(f); pars <- intersect(c(focal, "sd_spp_int"), colnames(x$draws))
    # the slope h2 and correlation are estimated only in the correlated (__pcor) fits
    if (!grepl("__pcor$", fid)) pars <- setdiff(pars, c("h2_slope", "cor_phylo_int_slope"))
    if (!length(pars)) next
    ch <- split(seq_along(x$chain), x$chain); nd <- min(lengths(ch))
    arr <- array(NA_real_, c(nd, length(ch), length(pars)), dimnames = list(NULL, NULL, pars))
    for (k in seq_along(ch)) arr[, k, ] <- x$draws[ch[[k]][seq_len(nd)], pars, drop = FALSE]
    s <- posterior::summarise_draws(posterior::as_draws_array(arr), "rhat", "ess_bulk", "ess_tail")
    diag_rows[[length(diag_rows) + 1]] <- data.frame(fit = fid, tree = basename(f), s)
  }
}
if (length(diag_rows)) {
  d <- bind_rows(diag_rows)
  pooled <- d %>% group_by(fit, variable) %>%
    summarise(rhat_max = max(rhat), ess_bulk_pooled = sum(ess_bulk), ess_tail_pooled = sum(ess_tail),
              n_trees = n(), .groups = "drop")
  out$diag_terms <- d
  out$diag_focal <- c(rhat_max = max(d$rhat[d$variable %in% focal]),
                      n_tree_terms = sum(d$variable %in% focal),
                      n_rhat_gt_101 = sum(d$rhat[d$variable %in% focal] > 1.01),
                      ess_bulk_pooled_min = min(pooled$ess_bulk_pooled[pooled$variable %in% focal]),
                      ess_tail_pooled_min = min(pooled$ess_tail_pooled[pooled$variable %in% focal]))
  out$diag_sd_spp_int <- c(rhat_max = max(d$rhat[d$variable == "sd_spp_int"]),
                           n_trees = sum(d$variable == "sd_spp_int"),
                           n_rhat_gt_101 = sum(d$rhat[d$variable == "sd_spp_int"] > 1.01),
                           ess_bulk_min = min(d$ess_bulk[d$variable == "sd_spp_int"]))
}

saveRDS(out, out_path("bayes", "bayes_derived.rds"))
print(out$species_counts); print(out$community_mm_decade); print(out$diet_quantile_slopes)
print(out$diet_category_slopes); print(out$diag_focal); print(out$diag_sd_spp_int)
