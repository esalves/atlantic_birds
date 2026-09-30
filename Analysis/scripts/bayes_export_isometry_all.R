# bayes_export_isometry_all.R
# ---------------------------------------------------------------------------
# Writes Stan-tier job files for the year-only isometry models on ALL shared
# wing + mass records, without requiring a temperature anomaly (co-author
# comment C16: the year effect does not need the climate data, so the temporal
# isometry test should not lose the records outside the WorldClim land mask).
#
# The formulas are those of the climate-restricted jobs quoted in the main text
# (referee_reruns__0{09,10,11,12}__*__site, exported by bayes_export_thermal.R);
# only the data change. The data are rebuilt by evaluating referee_reruns.R up
# to its E0 block, as bayes_export_thermal.R does, and the same M3 guards
# (longitude, altitude and season recorded) are applied.
#
# RUN:  Rscript Analysis/scripts/bayes_export_isometry_all.R
# OUT:  output/bayes/jobs/referee_reruns__0{15,16,17,18}__*__all.rds
#       output/bayes/jobs_isometry_all.txt (job list for run_bayes_totoro.sh;
#       run 016 a second time with --phylo-cor for the correlated-slope check)
# ---------------------------------------------------------------------------
Sys.setenv(REFEREE_N_TREES = "1")          # the A matrices built during data prep are not used here
local({
  f <- if (file.exists("Analysis/scripts/referee_reruns.R")) "Analysis/scripts/referee_reruns.R" else "scripts/referee_reruns.R"
  ex <- parse(f)
  stop_at <- which(vapply(ex, function(e) grepl('note("E0', paste(deparse(e), collapse = " "), fixed = TRUE), TRUE))[1]
  for (e in ex[seq_len(stop_at - 1)]) eval(e, envir = globalenv())
})

shared_all <- dat %>% filter(!is.na(wing), !is.na(lnmass)) %>%
  mutate(iso = lnmass - 3 * lnwing)
d_all <- shared_all %>% filter(!is.na(scaled_lon), !is.na(scaled_alt), !is.na(season))

dir <- out_path("bayes", "jobs")
pairs <- c(referee_reruns__015__iso__all    = "referee_reruns__009__iso__site",
           referee_reruns__016__iso__all    = "referee_reruns__010__iso__site",
           referee_reruns__017__lnwing__all = "referee_reruns__011__lnwing__site",
           referee_reruns__018__lnmass__all = "referee_reruns__012__lnmass__site")
for (id in names(pairs)) {
  src <- readRDS(file.path(dir, paste0(pairs[[id]], ".rds")))
  f <- src$formula
  stopifnot(!any(c("dT", "dT_detr") %in% all.vars(f)))
  vars <- unique(c(all.vars(f), src$species_col))
  job <- list(id = id, script = "referee_reruns", label = NA_character_, formula = f,
              dispformula = src$dispformula, species_col = src$species_col,
              data = as.data.frame(d_all)[, vars, drop = FALSE],
              glmmtmb_pooled = NULL, source_job = pairs[[id]], exported = Sys.time())
  saveRDS(job, file.path(dir, paste0(id, ".rds")))
  message(sprintf("%s (formula of %s): n = %d (climate-restricted job: n = %d), %d species",
                  id, pairs[[id]], nrow(job$data), nrow(src$data),
                  length(unique(job$data[[job$species_col]]))))
}
writeLines(file.path("Analysis/output/bayes/jobs", paste0(names(pairs), ".rds")),
           out_path("bayes", "jobs_isometry_all.txt"))
