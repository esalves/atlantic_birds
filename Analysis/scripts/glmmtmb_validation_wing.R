# glmmtmb_validation_wing.R
# ---------------------------------------------------------------------------
# WHAT: refit the published 50-tree brms wing-length model (atlantic_parallel.R,
#       brm0_multiphylo.rda) in glmmTMB using the phylogenetic `propto` covariance
#       structure (Williams, McGillycuddy, Drobniak, Bolker, Warton & Nakagawa 2025,
#       bioRxiv 10.64898/2025.12.20.695312), on the IDENTICAL data and per-tree
#       covariance matrices stored in each brmsfit ($data / $data2$A), pool across
#       trees with Rubin's rules, and compare with the brms Rubin summary.
# WHY:  establishes whether glmmTMB (seconds per fit) can replace brms (hours per
#       50-tree run) as the phylogenetic engine for the revision. Result 2026-09-09:
#       year -0.922 [-1.400, -0.444] (glmmTMB) vs -0.922 [-1.403, -0.440] (brms);
#       Sex, latitude, species-slope SD, phylogenetic SD and sigma agree to 2-3
#       decimals; 50 fits in ~24 s. See GLMMTMB_ENGINE.md.
# INPUTS:  Analysis/output/models/brm0_multiphylo.rda (git-ignored; on Totoro too)
# OUTPUTS: Analysis/output/glmmtmb_validation_wing.rds
# RUN:     Rscript Analysis/scripts/glmmtmb_validation_wing.R
# NOTE:  glmmTMB >= 1.1.14 (CRAN) already carries `propto`; the 1.1.15 GitHub build
#        fails against TMB 1.9.21 on macOS and is not needed.
# ---------------------------------------------------------------------------
suppressMessages({library(glmmTMB); library(dplyr)})
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived"))) return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/.")
}
ANALYSIS_DIR <- .find_analysis_dir()
out_path <- function(...) file.path(ANALYSIS_DIR, "output", ...)

t0 <- Sys.time()
load(out_path("models", "brm0_multiphylo.rda"))          # fits (50 brmsfit), rubin_summary
dat <- fits[[1]]$data; dat$g <- factor(1)
stopifnot(nrow(dat) == 8478)
cat("brms formula:", deparse(fits[[1]]$formula$formula), "\n")

fit_one <- function(i) {
  A <- fits[[i]]$data2$A
  lv <- levels(factor(dat$species_name)); A <- A[lv, lv]
  t1 <- Sys.time()
  m <- glmmTMB(conc.wing.length ~ 1 + Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) +
                 propto(0 + species_name | g, A), data = dat, REML = TRUE)
  secs <- as.numeric(Sys.time() - t1, units = "secs")
  fe <- summary(m)$coefficients$cond
  vc <- VarCorr(m)$cond
  spp_sd <- attr(vc$spp, "stddev")                        # diag block: (Intercept), scaled_yr
  data.frame(tree = i, par = rownames(fe), estimate = fe[, "Estimate"], se = fe[, "Std. Error"],
             sd_phylo = attr(vc$g, "stddev")[1], sd_spp_int = spp_sd[["(Intercept)"]],
             sd_spp_yr = spp_sd[["scaled_yr"]], sigma = sigma(m), secs = secs, row.names = NULL)
}
res <- do.call(rbind, lapply(seq_along(fits), fit_one))
pool <- res %>% group_by(par) %>%
  summarise(m = n(), qbar = mean(estimate), ubar = mean(se^2), b = var(estimate), .groups = "drop") %>%
  mutate(se = sqrt(ubar + (1 + 1/m) * b), estimate = qbar, lower = qbar - 1.96 * se, upper = qbar + 1.96 * se) %>%
  select(par, m, estimate, se, lower, upper, ubar, b)
vcomp <- res %>% filter(par == "(Intercept)") %>%
  mutate(phylo_prop = sd_phylo^2 / (sd_phylo^2 + sd_spp_int^2 + sd_spp_yr^2 + sigma^2))
brms_tree1 <- list(fixef = fixef(fits[[1]])[, c("Estimate", "Est.Error")],
                   sd = lapply(VarCorr(fits[[1]]), function(x) if (!is.null(x$sd)) x$sd[, "Estimate"]))
cat("\n== glmmTMB Rubin-pooled over", length(fits), "trees ==\n"); print(as.data.frame(pool), digits = 4)
cat("\n== brms Rubin-pooled (published) ==\n"); print(rubin_summary, digits = 4)
cat("\nvariance components, mean over trees: phylo sd", round(mean(vcomp$sd_phylo), 3),
    "| spp intercept sd", round(mean(vcomp$sd_spp_int), 3), "| spp year-slope sd", round(mean(vcomp$sd_spp_yr), 3),
    "| sigma", round(mean(vcomp$sigma), 3), "| phylo proportion", round(mean(vcomp$phylo_prop), 3),
    "[", round(min(vcomp$phylo_prop), 3), ",", round(max(vcomp$phylo_prop), 3), "] (brms 0.95 [0.94, 0.96])\n")
cat("brms tree 1: phylo sd", round(brms_tree1$sd$species_name, 3), "| spp sds", round(brms_tree1$sd$spp, 3), "| sigma", round(brms_tree1$sd$residual__, 3), "\n")
cat("median seconds per glmmTMB fit:", round(median(res$secs), 2), "| total", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
saveRDS(list(generated = Sys.time(), glmmTMB_version = as.character(packageVersion("glmmTMB")),
             per_tree = res, pooled = pool, varcomp = vcomp, brms_rubin = rubin_summary, brms_tree1 = brms_tree1),
        out_path("glmmtmb_validation_wing.rds"))
