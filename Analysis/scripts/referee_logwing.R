# referee_logwing.R
# ---------------------------------------------------------------------------
# WHY (internal review before the Oikos submission, 2026-09-30): wing length is
# modelled in millimetres while body mass and the isometry contrast are modelled
# on the log scale. An assemblage year slope in mm weights large-winged species
# more, and the manuscript's "% per decade" for wing is a conversion. This script
# refits the wing ladder with log(wing length) as the response, side by side
# with the mm response, so the proportional trend is estimated directly.
#
# Models (same formulas as atlantic_parallel_controlled.R, response swapped):
#   M0, M1, M3, M5src (year within contributor) and M3_noanom.
# Engine: glmmTMB REML with the phylogenetic species intercept (propto) on each
# of --trees trees (default 50), Rubin-pooled; the same engine as the REML checks
# quoted in the manuscript. --trees 0 gives a fast lme4 (no phylogeny) version.
#
# RUN (Totoro): Rscript Analysis/scripts/referee_logwing.R --trees 50 --cores 20
#      smoke:   Rscript Analysis/scripts/referee_logwing.R --trees 2 --models M0,M3
# OUT: Analysis/output/referee_reruns/logwing.rds (+ logwing_summary.csv)
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

opt <- function(flag, default) { a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i)) default else if (i == length(a) || startsWith(a[i + 1], "--")) TRUE else a[i + 1] }
N_TREES_LW <- as.integer(opt("--trees", 50))
CORES  <- as.integer(opt("--cores", max(1L, detectCores() - 1L)))
MODELS <- strsplit(opt("--models", "M0,M1,M3,M5src,M3_noanom"), ",")[[1]]
stopifnot(all(MODELS %in% names(SPECS)))
OUT_DIR <- out_path("referee_reruns"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
tag <- if (N_TREES_LW < 50) sprintf("_%dtrees", N_TREES_LW) else ""
MMDEC <- 10 / SD_YEAR
YEAR_TERM <- c(M5src = "yr_within_src")                         # default: scaled_yr
year_term <- function(m) if (m %in% names(YEAR_TERM)) YEAR_TERM[[m]] else "scaled_yr"
RESPONSES <- c(mm = "conc.wing.length", log = "lnwing")
formula_for <- function(s, resp) as.formula(paste(resp, "~ 1 +", s$fx, "+", s$re))

frames_lw <- lapply(frames, function(x) { x$lnwing <- log(x$conc.wing.length); x })
t0 <- Sys.time()

if (N_TREES_LW > 0) {
  source(script_path("_phylo_engine.R"))
  # species on the tree, as in section 4 of atlantic_parallel_controlled.R
  tree_species <- readRDS(derived_path("phylo_A_50trees.rds"))$species
  NOT_ON_TREE <- "Herpsilochmus_sellowi"
  to_tree <- function(x) {
    x$species_name <- phylo_species_name(x$Binomial)
    x %>% filter(species_name %in% setdiff(tree_species, NOT_ON_TREE)) %>%
      mutate(spp = factor(species_name)) %>%
      select(-yr_site_mean, -yr_src_mean, -lat_spp_mean, -yr_within_site, -yr_within_src, -lat_within_spp) %>%
      finish()
  }
  frames_lw <- lapply(frames_lw, to_tree)
}

fit_one <- function(m, resp, tree) {
  s <- guard_spec(SPECS[[m]], frames_lw[[SPECS[[m]]$data]]); dat <- frames_lw[[s$data]]
  f <- formula_for(s, RESPONSES[[resp]]); term <- year_term(m)
  if (N_TREES_LW == 0) {
    fit <- lmer(f, data = dat, REML = TRUE, control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
    co <- summary(fit)$coefficients
    return(data.frame(model = m, response = resp, tree = NA, estimate = co[term, 1], se = co[term, 2], pd = TRUE))
  }
  A <- phylo_A_list(sort(unique(as.character(dat$species_name))), n_trees = N_TREES_LW)[[tree]]
  fit <- tryCatch(suppressWarnings(fit_phylo_glmmtmb(f, dat, A)), error = function(e) NULL)
  if (is.null(fit)) return(data.frame(model = m, response = resp, tree = tree, estimate = NA, se = NA, pd = FALSE))
  co <- summary(fit)$coefficients$cond
  data.frame(model = m, response = resp, tree = tree, estimate = co[term, 1], se = co[term, 2], pd = isTRUE(fit$sdr$pdHess))
}

jobs <- expand.grid(model = MODELS, response = names(RESPONSES), tree = if (N_TREES_LW > 0) seq_len(N_TREES_LW) else NA,
                    stringsAsFactors = FALSE)
message(sprintf("[logwing] %d fits (%s; %s) on %d cores", nrow(jobs), paste(MODELS, collapse = "/"),
                if (N_TREES_LW > 0) sprintf("glmmTMB propto, %d trees", N_TREES_LW) else "lme4", CORES))
rows <- mclapply(seq_len(nrow(jobs)), function(k) fit_one(jobs$model[k], jobs$response[k], jobs$tree[k]), mc.cores = CORES)
bad <- vapply(rows, inherits, TRUE, "try-error"); if (any(bad)) warning(sum(bad), " fits failed")
per_tree <- do.call(rbind, rows[!bad])

# Rubin pooling over trees; non-PD fits are dropped (as pool_rubin_df does)
pooled <- per_tree %>% filter(pd, is.finite(estimate), is.finite(se)) %>%
  group_by(model, response) %>%
  summarise(m_used = n(), qbar = mean(estimate), ubar = mean(se^2), b = if (n() > 1) var(estimate) else 0, .groups = "drop") %>%
  mutate(se = sqrt(ubar + (1 + 1 / m_used) * b), lo = qbar - 1.96 * se, hi = qbar + 1.96 * se)
mean_wing <- mean(frames_lw$w_all$conc.wing.length)
summary_tbl <- pooled %>% transmute(
  model, response, m_used,
  per_decade = ifelse(response == "mm", qbar * MMDEC, 100 * (exp(qbar * MMDEC) - 1)),
  lo = ifelse(response == "mm", lo * MMDEC, 100 * (exp(lo * MMDEC) - 1)),
  hi = ifelse(response == "mm", hi * MMDEC, 100 * (exp(hi * MMDEC) - 1)),
  unit = ifelse(response == "mm", "mm/decade", "%/decade"),
  # mm slope as a share of the mean wing, for a rough side-by-side with the log fit
  pct_of_mean_wing = ifelse(response == "mm", 100 * qbar * MMDEC / mean_wing, NA_real_)) %>%
  arrange(match(model, MODELS), response)
print(as.data.frame(summary_tbl), digits = 3)

res <- list(generated = Sys.time(), n_trees = N_TREES_LW, models = MODELS,
            engine = if (N_TREES_LW > 0) "glmmTMB REML, propto phylogenetic species intercept, Rubin-pooled" else "lme4 REML",
            mean_wing_mm = mean_wing, per_tree = per_tree, pooled = pooled, summary = as.data.frame(summary_tbl),
            secs = as.numeric(Sys.time() - t0, units = "secs"))
saveRDS(res, file.path(OUT_DIR, paste0("logwing", tag, ".rds")))
write.csv(summary_tbl, file.path(OUT_DIR, paste0("logwing_summary", tag, ".csv")), row.names = FALSE)
message(sprintf("[done] %.1f min", res$secs / 60))
