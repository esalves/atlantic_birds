# bayes_aggregate.R
# ---------------------------------------------------------------------------
# Pools the per-tree Stan fits written by bayes_fit_job.R (mixture posterior over
# trees) and sets each against the pooled glmmTMB (REML, Rubin) result stored in
# the job file.
#
# RUN: Rscript Analysis/scripts/bayes_aggregate.R
# OUT: output/bayes/bayes_summary.rds  ($fixed, $varcomp, $diag, $comparison, $comparison_sigma)
#      output/bayes/BAYES_SUMMARY.md
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
source(file.path(ANALYSIS_DIR, "scripts", "_stan_engine.R"))

# Readable model names. Job ids (<script>__<nnn>__<response>__<code>) stay as the
# on-disk keys; `name` is what the summary shows, `code` the label used in the
# source script. Variant suffixes on a fit id (__site, __pcor) are appended.
MODEL_NAMES <- tribble(
  ~job_id,                                                   ~name,                                          ~code,
  "atlantic_parallel_controlled__001__conc.wing.length__M0",    "Wing: baseline",                               "M0",
  "atlantic_parallel_controlled__011__conc.wing.length__M3",    "Wing: fully adjusted",                         "M3",
  "atlantic_parallel_controlled__012__conc.wing.length__M4",    "Wing: fully adjusted + per-contributor trends", "M4",
  "atlantic_parallel_controlled__015__conc.wing.length__M5src", "Wing: year within vs between contributors",    "M5src",
  "atlantic_parallel_controlled__038__conc.wing.length__M6",    "Wing: fully adjusted, unknown-sex birds",      "M6",
  "atlantic_parallel_controlled__017__conc.wing.length__M7",    "Wing: fully adjusted, first captures only",    "M7",
  "atlantic_multitrait__003__y",                                "Body mass (ln): fully adjusted",               "multitrait M3_ring",
  "atlantic_multitrait__007__y",                                "Bill width: fully adjusted",                   "multitrait M3_ring",
  "atlantic_diet_interaction__002__conc.wing.length",           "Wing x diet: baseline",                        "diet A_baseline",
  "atlantic_diet_interaction__006__conc.wing.length",           "Wing x diet: + contributor & municipality",    "diet A_src_site",
  "atlantic_variance_sigma__002__conc.wing.length__S1",         "Wing (mean + SD model): baseline",             "S1",
  "atlantic_variance_sigma__004__conc.wing.length__S3",         "Wing (mean + SD model): + contributor & municipality", "S3",
  "referee_reruns__005__iso",                                   "Isometry vs temperature anomaly",              "E3",
  "referee_reruns__010__iso",                                   "Isometry vs year: fully adjusted",                          "E4b",
  # main-text models added 2026-09-24 (names agreed with ESAS)
  "atlantic_parallel_controlled__023__conc.wing.length__M3_cc",    "Wing: fully adjusted, complete cases",          "M3_cc",
  "atlantic_parallel_controlled__027__conc.wing.length__M5src_cc", "Wing: year within vs between contributors, complete cases", "M5src_cc",
  "atlantic_multitrait__001__y",                                "Body mass (ln): baseline",                     "multitrait M0_baseline",
  "atlantic_multitrait__004__y",                                "Body mass (ln): year within vs between contributors", "multitrait MWsrc",
  "atlantic_multitrait__021__y",                                "Body mass (ln): capture-time sample",          "multitrait HOUR0",
  "atlantic_multitrait__022__y",                                "Body mass (ln): capture-time sample + capture hour", "multitrait HOUR1",
  "atlantic_multitrait__005__y",                                "Bill width: baseline",                         "multitrait M0_baseline",
  "atlantic_multitrait__006__y",                                "Bill width: + contributor & municipality",     "multitrait M1_src_site",
  "atlantic_multitrait__008__y",                                "Bill width: year within vs between contributors", "multitrait MWsrc",
  "atlantic_variance_sigma__009__conc.wing.length__S7",         "Wing (mean + SD model): + local temperature",  "S7",
  "atlantic_diet_interaction__012__conc.wing.length",           "Wing x diet category: + contributor & municipality", "diet B_src_site",
  "referee_reruns__001__wing",                                  "Wing vs temperature anomaly",                  "E3",
  "referee_reruns__003__lnmass",                                "Body mass (ln) vs temperature anomaly",        "E3",
  "referee_reruns__006__iso",                                   "Isometry vs temperature anomaly + year",       "E3b",
  "referee_reruns__009__iso",                                   "Isometry vs year: parsimonious adjustment",    "E4a",
  "referee_reruns__013__iso",                                   "Isometry vs detrended temperature anomaly",    "E1b (phylo, M3)",
  "referee_reruns__014__iso",                                   "Isometry vs temperature anomaly + locality-year", "E2a (phylo, M3)",
  "referee_reruns__011__lnwing",                                "Wing (ln): fully adjusted, shared records",    "E4c",
  "referee_reruns__012__lnmass",                                "Body mass (ln): fully adjusted, shared records", "E4c")
VARIANT_NAMES <- c(site = "+ municipality", pcor = "correlated phylo intercept-slope", pslope = "phylo slope")
TERM_NAMES <- c(b_scaled_yr = "year", b_yr_within_src = "year within contributor", b_yr_src_mean = "contributor mean year",
                b_dT = "temperature anomaly", `b_scaled_yr:diet_inv_std` = "year x diet (invertebrate share)",
                b_yr_mean_src = "contributor mean year", b_dT_detr = "detrended temperature anomaly",
                b_hour_num = "capture hour", sigma_scaled_yr = "year", sigma_SexMale = "male (vs female)",
                sigma_scaled_tmean = "local temperature")
fit_name <- function(id) vapply(id, function(x) {
  parts <- strsplit(x, "__")[[1]]
  sfx <- intersect(parts, names(VARIANT_NAMES)); base <- paste(setdiff(parts, names(VARIANT_NAMES)), collapse = "__")
  nm <- MODEL_NAMES$name[match(base, MODEL_NAMES$job_id)]
  if (is.na(nm)) nm <- x
  if (length(sfx)) nm <- paste0(nm, " [", paste(VARIANT_NAMES[sfx], collapse = ", "), "]")
  nm
}, "", USE.NAMES = FALSE)
fit_code <- function(id) MODEL_NAMES$code[match(gsub("__(site|pcor|pslope)", "", id), MODEL_NAMES$job_id)]
term_name <- function(par) ifelse(par %in% names(TERM_NAMES), TERM_NAMES[par], par)

fit_dirs <- list.dirs(out_path("bayes", "fits"), recursive = FALSE)
res <- lapply(fit_dirs, function(dd) {
  files <- list.files(dd, pattern = "^tree_[0-9]+\\.rds$", full.names = TRUE)
  if (!length(files)) return(NULL)
  per_tree <- lapply(files, readRDS)
  id <- basename(dd); job_id <- sub("__p(slope|cor)$", "", id)
  job <- readRDS(out_path("bayes", "jobs", paste0(job_id, ".rds")))
  pooled <- pool_stan_draws(per_tree) %>% mutate(fit = id, .before = 1)
  diag <- do.call(rbind, lapply(per_tree, function(r) data.frame(fit = id, tree = r$tree, secs = r$secs, r$diag)))
  cmp <- NULL
  if (!is.null(job$glmmtmb_pooled)) {
    g <- job$glmmtmb_pooled %>% filter(component == "cond") %>%
      transmute(par = paste0("b_", par), reml = estimate, reml_se = se, reml_lo = lower, reml_hi = upper)
    cmp <- inner_join(pooled %>% select(fit, par, estimate, sd, lower, upper, p_neg), g, by = "par") %>%
      mutate(diff_in_sd = (estimate - reml) / sd, sd_ratio = sd / reml_se)
  }
  # residual-SD (dispersion) terms: glmmTMB `disp` and Stan `sigma_` are both on the log-sigma scale
  cmp_sigma <- NULL
  if (!is.null(job$glmmtmb_pooled) && any(job$glmmtmb_pooled$component == "disp" & job$glmmtmb_pooled$par != "(Intercept)")) {
    g <- job$glmmtmb_pooled %>% filter(component == "disp", par != "(Intercept)") %>%
      transmute(par = paste0("sigma_", par), reml = estimate, reml_se = se, reml_lo = lower, reml_hi = upper)
    cmp_sigma <- inner_join(pooled %>% select(fit, par, estimate, sd, lower, upper, p_neg), g, by = "par") %>%
      mutate(diff_in_sd = (estimate - reml) / sd, sd_ratio = sd / reml_se)
  }
  list(pooled = pooled, diag = diag, comparison = cmp, comparison_sigma = cmp_sigma,
       meta = data.frame(fit = id, name = fit_name(id), code = fit_code(id), script = job$script, label = job$label, formula = paste(deparse(job$formula), collapse = " "),
                         n = nrow(job$data), n_trees = length(per_tree)))
})
res <- Filter(Negate(is.null), res)
out <- list(generated = Sys.time(),
            meta = bind_rows(lapply(res, `[[`, "meta")),
            fixed = bind_rows(lapply(res, function(r) r$pooled %>% filter(startsWith(par, "b_") | startsWith(par, "sigma_")))),
            varcomp = bind_rows(lapply(res, function(r) r$pooled %>% filter(!startsWith(par, "b_"), !startsWith(par, "sigma_"), !startsWith(par, "slope_dev_")))),
            species_slope_dev = bind_rows(lapply(res, function(r) r$pooled %>% filter(startsWith(par, "slope_dev_")))),
            diag = bind_rows(lapply(res, `[[`, "diag")),
            comparison = bind_rows(lapply(res, `[[`, "comparison")),
            comparison_sigma = bind_rows(lapply(res, `[[`, "comparison_sigma")))
saveRDS(out, out_path("bayes", "bayes_summary.rds"))

fmt <- function(x) formatC(x, digits = 3, format = "f")
key <- out$comparison %>% filter(grepl("scaled_yr|dT|yr_within|yr_src|yr_mean|hour_num|:", par)) %>% arrange(fit)
diag_tab <- out$diag %>% group_by(fit) %>% arrange(fit) %>%
  summarise(n = n(), r = max(rhat_max), e = min(ess_bulk_min), dv = sum(divergent), s = median(secs))
legend <- out$meta %>% arrange(fit)
lines <- c("# Bayesian (Stan) tier: summary", "",
           sprintf("Generated %s. %d fits. Posterior pooled over trees (mixture); REML = glmmTMB Rubin-pooled over the same trees.",
                   format(out$generated), nrow(out$meta)), "",
           "## How to read the tables", "",
           "- **P(<0)** is the posterior probability that the coefficient is negative: the share of posterior draws,",
           "  pooled over all trees, that fall below zero. P(<0) = 0.96 means that, given the model, the priors and the",
           "  data, there is a 96 % probability that the effect is negative (e.g. that wing length declined over time).",
           "  The probability of a positive effect is 1 - P(<0). Values near 0.5 mean no evidence for either direction;",
           "  values near 1 are strong evidence for a negative effect, values near 0 strong evidence for a positive one.",
           "  It is not a p-value: it is a direct statement about the effect, not about data under a null hypothesis.",
           "  For a roughly symmetric posterior, the 95 % CrI excludes zero exactly when P(<0) > 0.975 or < 0.025.",
           "- **shift (post. SD)**: (posterior mean - REML estimate) / posterior SD. |shift| < 0.1 means the two engines agree.",
           "- **SD ratio**: posterior SD / REML SE. Slightly above 1 is expected, because the posterior also carries the",
           "  uncertainty in the variance components that REML plugs in as known.",
           "- Year coefficients are per SD of calendar year (`scaled_yr`), on the response scale (mm for wing and bill,",
           "  ln units for mass and the isometry contrast ln M - 3 ln L).", "",
           "## Models", "",
           "| model | code in source script | fit id | n | trees |", "|---|---|---|---|---|",
           with(legend, sprintf("| %s | %s | `%s` | %d | %d |", name, ifelse(is.na(code), "", code), fit, n, n_trees)),
           "", "## Sampler diagnostics (worst tree per fit)", "",
           "| model | trees | max R-hat | min bulk-ESS | divergences | median s/tree |", "|---|---|---|---|---|---|",
           with(diag_tab, sprintf("| %s | %d | %.3f | %.0f | %d | %.0f |", fit_name(fit), n, r, e, dv, s)),
           "", "## Time terms: posterior vs REML", "",
           "| model | term | posterior mean [95% CrI] | P(<0) | REML [95% CI] | shift (post. SD) | SD ratio |", "|---|---|---|---|---|---|---|",
           with(key, sprintf("| %s | %s | %s [%s, %s] | %.3f | %s [%s, %s] | %.2f | %.2f |", fit_name(fit), term_name(par), fmt(estimate), fmt(lower), fmt(upper),
                             p_neg, fmt(reml), fmt(reml_lo), fmt(reml_hi), diff_in_sd, sd_ratio)),
           "", "## Residual-SD terms (mean + SD models): posterior vs REML", "",
           "Coefficients are on the log(residual SD) scale; `% SD change` = 100 (exp(b) - 1), the percentage change",
           "in within-species residual SD per SD of year (or males vs females). P(<0) here is the probability that",
           "the residual SD shrinks (for year: that wing length became less variable over time).", "",
           "| model | term | posterior mean [95% CrI] | % SD change [95% CrI] | P(<0) | REML [95% CI] | shift (post. SD) | SD ratio |",
           "|---|---|---|---|---|---|---|---|",
           with(out$comparison_sigma %>% arrange(fit),
                sprintf("| %s | %s | %s [%s, %s] | %+.1f [%+.1f, %+.1f] | %.3f | %s [%s, %s] | %.2f | %.2f |", fit_name(fit), term_name(par),
                        fmt(estimate), fmt(lower), fmt(upper), 100 * (exp(estimate) - 1), 100 * (exp(lower) - 1), 100 * (exp(upper) - 1),
                        p_neg, fmt(reml), fmt(reml_lo), fmt(reml_hi), diff_in_sd, sd_ratio)))
writeLines(lines, out_path("bayes", "BAYES_SUMMARY.md"))
message("wrote ", out_path("bayes", "BAYES_SUMMARY.md"))
