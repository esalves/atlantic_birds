# bayes_fit_job.R
# ---------------------------------------------------------------------------
# Bayesian (Stan) refit of ONE exported glmmTMB model on ONE tree.
#
# Jobs are the model specifications the glmmTMB scripts write when they run with
# PHYLO_EXPORT_DIR set (see run_phylo_trees() in _phylo_engine.R): formula,
# dispformula, the exact model frame and the pooled glmmTMB result. This keeps the
# Bayesian tier on the identical data and formulas of the published REML tier.
#
# RUN:  Rscript Analysis/scripts/bayes_fit_job.R --job output/bayes/jobs/<id>.rds --tree 7
#         [--phylo-slope] [--phylo-cor] [--chains 4] [--warmup 1000] [--sampling 500] [--thin 1]
# OUT:  output/bayes/fits/<id>[__pslope|__pcor]/tree_07.rds
#       (draws on the response scale, sampler diagnostics, timing)
#       Skipped when tree_07.rds exists or a tree_07.rds.claimed marker does.
# Driver: run_bayes_totoro.sh runs job x tree pairs in parallel and skips done ones.
# ---------------------------------------------------------------------------
suppressMessages({ library(dplyr); library(posterior) })
# Totoro is shared: cap BLAS/OpenMP threads (_phylo_engine.R also caps TMB at 1 thread).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
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

job_file <- get_flag("--job"); tree <- as.integer(get_flag("--tree", "1"))
if (is.null(job_file)) stop("--job is required")
if (!file.exists(job_file)) job_file <- file.path(ANALYSIS_DIR, job_file)
phylo_slope <- isTRUE(get_flag("--phylo-slope", FALSE))
phylo_cor   <- isTRUE(get_flag("--phylo-cor", FALSE))
if (phylo_cor) phylo_slope <- TRUE
job <- readRDS(job_file)
variant <- if (phylo_cor) "__pcor" else if (phylo_slope) "__pslope" else ""
dest_dir <- out_path("bayes", "fits", paste0(job$id, variant))
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
dest <- file.path(dest_dir, sprintf("tree_%02d.rds", tree))
if (file.exists(dest)) { message("done already: ", dest); quit(save = "no") }
# A tree_NN.rds.claimed marker means another host fits this tree (work split
# between Totoro and Kohaku); its tree_NN.rds is rsynced in afterwards.
if (file.exists(paste0(dest, ".claimed"))) { message("claimed by another host: ", dest); quit(save = "no") }

A <- phylo_A_list(sort(unique(as.character(job$data[[job$species_col]]))), n_trees = max(tree, 50))[[tree]]
r <- fit_stan_tree(job$formula, job$data, A, dispformula = job$dispformula,
                   species_col = job$species_col, phylo_slope = phylo_slope, phylo_cor = phylo_cor,
                   chains = as.integer(get_flag("--chains", "4")),
                   iter_warmup = as.integer(get_flag("--warmup", "1000")),
                   iter_sampling = as.integer(get_flag("--sampling", "500")),
                   seed = 20240101 + tree)
r <- thin_tree_draws(r, as.integer(get_flag("--thin", "1")))
r$tree <- tree; r$job <- job$id; r$variant <- variant
saveRDS(r, dest)
message(sprintf("[bayes] %s%s tree %d: %.0f s, rhat_max %.3f, ess_min %.0f, divergent %d",
                job$id, variant, tree, r$secs, r$diag$rhat_max, r$diag$ess_bulk_min, r$diag$divergent))
