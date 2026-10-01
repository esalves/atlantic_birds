# wing_year_effects.R
# ---------------------------------------------------------------------------
# WHY (Shinichi Nakagawa's comment on the wing-trend figure, 2026-09-30): the
# figure showed raw annual means +/- 1.96 SE, which treat every record as
# independent. Records from the same contributor, municipality and species are
# not, so those intervals are too narrow, most visibly in the well-sampled years
# around 2010. This script estimates one effect per capture year from the model
# instead: the M0 and M3 wing models with calendar year as a factor in place of
# the linear year term (species random slopes on scaled_yr are kept, as in M3).
#
# Year effects are centred on their record-weighted mean, the same reference as
# the fitted linear trends in the figure (zero at the record-weighted mean year),
# with standard errors from the model covariance of the centred contrasts.
# Engine: glmmTMB REML with the phylogenetic species intercept on each tree,
# Rubin-pooled per year (as the other REML checks).
#
# RUN:  Rscript Analysis/scripts/wing_year_effects.R --trees 50 --cores 8
# OUT:  output/figure_data/fig_wingtrend_year_effects.csv (read by make_figures.py)
#       output/wing_year_effects.rds
# ---------------------------------------------------------------------------
local({
  f <- c("Analysis/scripts/atlantic_parallel_controlled.R", "scripts/atlantic_parallel_controlled.R",
         "atlantic_parallel_controlled.R")
  f <- f[file.exists(f)][1]
  if (is.na(f)) stop("Run from the repo root, Analysis/ or Analysis/scripts/.")
  ex <- parse(f)
  stop_at <- which(vapply(ex, function(e) startsWith(deparse(e)[1], "YEAR_TERMS <-"), TRUE))[1]
  for (e in ex[seq_len(stop_at - 1)]) eval(e, envir = globalenv())
})
suppressMessages(library(parallel))
source(script_path("_phylo_engine.R"))

opt <- function(flag, default) { a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i)) default else if (i == length(a) || startsWith(a[i + 1], "--")) TRUE else a[i + 1] }
N_TREES_Y <- as.integer(opt("--trees", 50))
CORES <- as.integer(opt("--cores", max(1L, detectCores() - 1L)))

# primary wing frame on the tree names, as in section 4 of atlantic_parallel_controlled.R
tree_species <- readRDS(derived_path("phylo_A_50trees.rds"))$species
d <- w_all %>% mutate(species_name = phylo_species_name(Binomial)) %>%
  filter(species_name %in% setdiff(tree_species, "Herpsilochmus_sellowi")) %>%
  mutate(spp = factor(species_name), yr_f = factor(Year))
stopifnot(nrow(d) == nrow(w_all))
years <- levels(d$yr_f); w <- as.numeric(table(d$yr_f)[years])
A_list <- phylo_A_list(sort(unique(as.character(d$species_name))), n_trees = N_TREES_Y)

# year factor in place of the linear year term
f_yr <- function(s) as.formula(paste("conc.wing.length ~ 1 +", sub("scaled_yr", "yr_f", s$fx), "+", s$re))
FORMULAS <- list(M0 = f_yr(SPECS$M0), M3 = f_yr(SPECS$M3))

# centred year effects of one fit: b = K beta (reference year = 0), centred on the
# record-weighted mean, var = C V C'
K <- rbind(0, diag(length(years) - 1))
C <- (diag(length(years)) - matrix(w / sum(w), length(years), length(years), byrow = TRUE)) %*% K
centred <- function(fit) {
  co <- fixef(fit)$cond; V <- vcov(fit)$cond
  nm <- paste0("yr_f", years[-1]); stopifnot(all(nm %in% names(co)))
  b <- drop(C %*% co[nm]); se <- sqrt(diag(C %*% V[nm, nm] %*% t(C)))
  data.frame(year = as.integer(years), estimate = b, se = se)
}
jobs <- expand.grid(model = names(FORMULAS), tree = seq_len(N_TREES_Y), stringsAsFactors = FALSE)
t0 <- Sys.time()
rows <- mclapply(seq_len(nrow(jobs)), function(k) {
  fit <- tryCatch(suppressWarnings(fit_phylo_glmmtmb(FORMULAS[[jobs$model[k]]], d, A_list[[jobs$tree[k]]])),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$sdr$pdHess)) return(NULL)
  cbind(model = jobs$model[k], tree = jobs$tree[k], centred(fit))
}, mc.cores = CORES)
per_tree <- do.call(rbind, rows[!vapply(rows, function(x) is.null(x) || inherits(x, "try-error"), TRUE)])
pooled <- per_tree %>% group_by(model, year) %>%
  summarise(trees = n(), qbar = mean(estimate), ubar = mean(se^2), b = if (n() > 1) var(estimate) else 0, .groups = "drop") %>%
  mutate(se = sqrt(ubar + (1 + 1 / trees) * b), lo = qbar - 1.96 * se, hi = qbar + 1.96 * se) %>%
  left_join(data.frame(year = as.integer(years), n_records = w), by = "year") %>%
  transmute(model, year, estimate = qbar, se, lo, hi, trees, n_records)
print(as.data.frame(pooled), digits = 3)

dir.create(out_path("figure_data"), showWarnings = FALSE)
write.csv(pooled, out_path("figure_data", "fig_wingtrend_year_effects.csv"), row.names = FALSE)
saveRDS(list(generated = Sys.time(), n_trees = N_TREES_Y, formulas = lapply(FORMULAS, deparse1),
             centring = "record-weighted mean of the year effects = 0", pooled = as.data.frame(pooled),
             per_tree = per_tree, secs = as.numeric(Sys.time() - t0, units = "secs")),
        out_path("wing_year_effects.rds"))
message(sprintf("[done] %.1f min", as.numeric(Sys.time() - t0, units = "mins")))
