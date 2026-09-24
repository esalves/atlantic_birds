# bayes_export_thermal.R
# ---------------------------------------------------------------------------
# Writes Stan-tier job files for the thermal-anomaly and isometry models of
# referee_reruns.R that the manuscript quotes, all with the municipality
# intercept (1 | site). The E1/E2 robustness checks were plain glmmTMB fits
# (no phylogeny, parsimonious adjustment); here they become phylogenetic models
# at the M3 adjustment set, so every quoted thermal number comes from one engine.
#
# The data are rebuilt by evaluating referee_reruns.R up to its E0 block, so the
# model frame is exactly the one the REML tier used. Where referee_reruns.rds
# holds a matching 50-tree glmmTMB fit, its pooled row is attached for the
# posterior-vs-REML comparison; the two new specifications have none.
#
# RUN:  Rscript Analysis/scripts/bayes_export_thermal.R
# OUT:  output/bayes/jobs/referee_reruns__0NN__<resp>__site.rds
#       output/bayes/jobs_thermal.txt (job list for run_bayes_totoro.sh)
# ---------------------------------------------------------------------------
Sys.setenv(REFEREE_N_TREES = "1")          # the A matrices built during data prep are not used here
local({
  here <- normalizePath(getwd(), winslash = "/")
  f <- if (file.exists("Analysis/scripts/referee_reruns.R")) "Analysis/scripts/referee_reruns.R" else "scripts/referee_reruns.R"
  ex <- parse(f)
  stop_at <- which(vapply(ex, function(e) grepl('note("E0', paste(deparse(e), collapse = " "), fixed = TRUE), TRUE))[1]
  for (e in ex[seq_len(stop_at - 1)]) eval(e, envir = globalenv())
})

d3 <- shared %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season))
rr <- readRDS(out_path("referee_reruns", "referee_reruns.rds"))

m3_rhs <- "Sex + %s + scaled_lat + scaled_lon + scaled_alt + season + molt + (1 + %s || spp) + (1 | src) + (1 | site) + (1 | ind)"
specs <- list(
  list(id = "referee_reruns__001__wing__site",   f = sprintf(paste("wing ~", m3_rhs), "dT", "dT"),
       rr = list(tbl = "E3", model = "E3_M3_adjusted_phylo", response = "wing length (mm)")),
  list(id = "referee_reruns__003__lnmass__site", f = sprintf(paste("lnmass ~", m3_rhs), "dT", "dT"),
       rr = list(tbl = "E3", model = "E3_M3_adjusted_phylo", response = "log body mass")),
  list(id = "referee_reruns__006__iso__site",    f = sprintf(paste("iso ~", m3_rhs), "dT + scaled_yr", "scaled_yr"),
       rr = list(tbl = "E3", model = "E3b_M3_adjusted_phylo_plus_year", response = "isometry contrast lnM-3lnL")),
  list(id = "referee_reruns__009__iso__site",
       f = "iso ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp) + (1 | src) + (1 | site)",
       rr = list(tbl = "E4", model = "E4a_iso_parsimonious", response = "isometry contrast lnM-3lnL")),
  list(id = "referee_reruns__013__iso__site",    f = sprintf(paste("iso ~", m3_rhs), "dT_detr", "dT_detr"), rr = NULL),
  list(id = "referee_reruns__014__iso__site",    f = paste(sprintf(paste("iso ~", m3_rhs), "dT", "dT"), "+ (1 | loc_year)"), rr = NULL),
  # separate wing and mass trends on the shared records (the two sides of the contrast)
  list(id = "referee_reruns__011__lnwing__site", f = sprintf(paste("lnwing ~", m3_rhs), "scaled_yr", "scaled_yr"),
       rr = list(tbl = "E4", model = "E4c_separate_M3_slopes", response = "lnwing")),
  list(id = "referee_reruns__012__lnmass__site", f = sprintf(paste("lnmass ~", m3_rhs), "scaled_yr", "scaled_yr"),
       rr = list(tbl = "E4", model = "E4c_separate_M3_slopes", response = "lnmass"))
)

rr_pooled <- function(s) {
  if (is.null(s)) return(NULL)
  t <- rr[[s$tbl]]; t <- t[t$model == s$model & t$response == s$response, , drop = FALSE]
  if (!nrow(t)) return(NULL)
  tibble::tibble(component = "cond", par = t$term, estimate = t$estimate, se = t$se, lower = t$lower, upper = t$upper)
}

dir <- out_path("bayes", "jobs")
for (s in specs) {
  f <- as.formula(s$f)
  vars <- unique(c(all.vars(f), "species_name"))
  job <- list(id = s$id, script = "referee_reruns", label = NA_character_, formula = f, dispformula = ~1,
              species_col = "species_name", data = as.data.frame(d3)[, vars, drop = FALSE],
              glmmtmb_pooled = rr_pooled(s$rr), exported = Sys.time())
  saveRDS(job, file.path(dir, paste0(s$id, ".rds")))
  message(sprintf("%s: n = %d, REML comparison: %s", s$id, nrow(d3), if (is.null(job$glmmtmb_pooled)) "none" else "yes"))
}
writeLines(file.path("Analysis/output/bayes/jobs", paste0(vapply(specs, `[[`, "", "id"), ".rds")),
           out_path("bayes", "jobs_thermal.txt"))
