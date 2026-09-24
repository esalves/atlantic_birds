# stan_validation_wing.R
# ---------------------------------------------------------------------------
# Checks the custom Stan engine (_stan_engine.R + stan/phylo_lmm_marginal.stan)
# against glmmTMB on the same data and the same tree before any Bayesian run
# is trusted. Model: the M3-level thermal-sample wing model of referee_reruns.R
# (contributor, municipality and ring intercepts; species intercept + year
# slope; phylogenetic intercept), plus the phylogenetic-slope variant.
#
# Pass criteria: fixed effects within ~0.1 posterior SD of the REML estimates,
# posterior SDs within ~15 % of the REML SEs, variance components close, and
# no divergences.
#
# RUN: Rscript Analysis/scripts/stan_validation_wing.R [--trees 1] [--iter 500] [--subset 0]
# OUTPUT: output/stan_validation_wing.rds
# ---------------------------------------------------------------------------
suppressMessages({ library(dplyr); library(glmmTMB) })
# Totoro is shared: cap BLAS/OpenMP threads (_phylo_engine.R also caps TMB at 1 thread).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
set.seed(20240101)

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
N_TREES <- as.integer(get_flag("--trees", "1"))
ITER    <- as.integer(get_flag("--iter", "500"))
SUBSET  <- as.integer(get_flag("--subset", "0"))   # >0: random subset of species, for a quick laptop check

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
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
source(file.path(ANALYSIS_DIR, "scripts", "_phylo_engine.R"))
source(file.path(ANALYSIS_DIR, "scripts", "_stan_engine.R"))

pc <- readRDS(derived_path("passer90_climate.rds"))
d <- pc %>%
  transmute(wing = conc.wing.length, Sex = factor(Sex, levels = c("Female", "Male")),
            scaled_yr = as.numeric(scaled_yr), scaled_lat = as.numeric(scaled_lat),
            scaled_lon = as.numeric(scaled_lon), scaled_alt = as.numeric(scaled_alt),
            season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")),
            molt = factor(ifelse(is.na(Molt), "unknown", as.character(Molt))),
            src = factor(Main_researcher), site = factor(Municipality),
            ind = factor(paste(Ring, Binomial, sep = "_")),
            species_name = phylo_species_name(Binomial), spp = factor(species_name)) %>%
  filter(complete.cases(.))
A_list <- phylo_A_list(sort(unique(d$species_name)), n_trees = N_TREES)
d <- d[d$species_name %in% rownames(A_list[[1]]), ]
if (SUBSET > 0) {
  keep <- sample(unique(d$species_name), SUBSET)
  d <- d[d$species_name %in% keep, ]; d$spp <- droplevels(d$spp)
  A_list <- lapply(A_list, function(A) A[sort(keep), sort(keep)])
}
message(sprintf("[validation] %d records, %d species, %d rings (%d with repeats)",
                nrow(d), n_distinct(d$species_name), n_distinct(d$ind), sum(table(d$ind) > 1)))

f <- wing ~ Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season + molt +
  (1 + scaled_yr || spp) + (1 | src) + (1 | site) + (1 | ind)

out <- list(generated = Sys.time(), n = nrow(d), n_trees = N_TREES, iter = ITER, subset = SUBSET, trees = list())
for (k in seq_len(N_TREES)) {
  A <- A_list[[k]]
  tmb <- fit_phylo_glmmtmb(f, d, A)
  cf <- summary(tmb)$coefficients$cond
  vc <- VarCorr(tmb)$cond
  tmb_sd <- c(sd_spp_int_phylo = attr(vc$g, "stddev")[1], sd_spp_int_iid = attr(vc$spp, "stddev")[1],
              sd_spp_slope = attr(vc$spp, "stddev")[2], sd_src = attr(vc$src, "stddev"),
              sd_site = attr(vc$site, "stddev"), sd_ind = attr(vc$ind, "stddev"), sigma = sigma(tmb))

  st  <- fit_stan_tree(f, d, A, iter_warmup = ITER, iter_sampling = ITER)
  st2 <- fit_stan_tree(f, d, A, iter_warmup = ITER, iter_sampling = ITER, phylo_slope = TRUE)

  D <- st$draws
  cmp <- data.frame(par = rownames(cf), reml = cf[, 1], reml_se = cf[, 2],
                    stan = colMeans(D[, paste0("b_", rownames(cf))]),
                    stan_sd = apply(D[, paste0("b_", rownames(cf))], 2, sd), row.names = NULL) %>%
    mutate(diff_in_sd = (stan - reml) / stan_sd, sd_ratio = stan_sd / reml_se)
  vc_cmp <- data.frame(par = names(tmb_sd), reml = unname(tmb_sd),
                       stan = c(mean(D[, "sd_spp_int"] * sqrt(D[, "h2_int"])), mean(D[, "sd_spp_int"] * sqrt(1 - D[, "h2_int"])),
                                mean(D[, "sd_spp_slope"]), mean(D[, "sd_src"]), mean(D[, "sd_site"]),
                                mean(D[, "sd_ind"]), mean(exp(D[, "sigma_(Intercept)"]))))
  D2 <- st2$draws
  out$trees[[k]] <- list(fixed = cmp, varcomp = vc_cmp, diag = st$diag, secs = st$secs,
                         phylo_slope = list(yr = quantile(D2[, "b_scaled_yr"], c(.025, .5, .975)),
                                            h2_slope = quantile(D2[, "h2_slope"], c(.025, .5, .975)),
                                            sd_slope = quantile(D2[, "sd_spp_slope"], c(.025, .5, .975)),
                                            diag = st2$diag, secs = st2$secs))
  message(sprintf("[validation] tree %d: stan %.0f s (phylo-slope %.0f s)", k, st$secs, st2$secs))
  print(cmp, digits = 3); print(vc_cmp, digits = 3); str(st$diag); str(out$trees[[k]]$phylo_slope)
}
saveRDS(out, out_path(if (SUBSET > 0) "stan_validation_wing_subset.rds" else "stan_validation_wing.rds"))
