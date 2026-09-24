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

saveRDS(out, out_path("bayes", "bayes_derived.rds"))
print(out$species_counts); print(out$community_mm_decade); print(out$diet_quantile_slopes)
