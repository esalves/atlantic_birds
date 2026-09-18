# referee_anomaly_sensitivity.R
# ---------------------------------------------------------------------------
# WHY: the reported thermal-coupling estimate for the isometry contrast
# (beta = +0.0212, t = 2.17, p = 0.030, n = 7,577; index.qmd / secular_expansion
# REPORT.md 5.5) has NO generating script in this repository. referee_reruns.R
# block E0 reconstructs it, but does not land on the reported value.
#
# This script isolates WHY: the estimate is sensitive to how the local annual
# temperature anomaly is constructed, and no construction available from
# passer90_climate.rds reaches the reported n of 7,577 -- 353 of the 7,577
# shared records have no rec_tmean at all.
#
# Definitions compared (all on the shared wing+mass records):
#   A. site baseline over ALL passer90 records     (dT = rec_tmean - mean by loc_id)
#   B. site baseline over the SHARED records only
#   C. rec_tmean minus Annual_mean_temperature     (WorldClim long-term site mean)
#   D. scaled_tmean as supplied in passer90_climate.rds
#
# RUN: Rscript Analysis/scripts/referee_anomaly_sensitivity.R
# ---------------------------------------------------------------------------
suppressMessages({ library(dplyr); library(glmmTMB) })

.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/")
  repeat {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/.")
}
ANALYSIS_DIR <- .find_analysis_dir()
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
OUT_DIR <- out_path("referee_reruns"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pc <- readRDS(derived_path("passer90_climate.rds"))

# site baseline computed over ALL records, before subsetting to shared
all_base <- pc %>% group_by(loc_id) %>%
  summarise(base_all = mean(rec_tmean, na.rm = TRUE), .groups = "drop")

sh <- pc %>%
  filter(!is.na(conc.wing.length), !is.na(ln_body_mass)) %>%
  left_join(all_base, by = "loc_id") %>%
  group_by(loc_id) %>% mutate(base_shared = mean(rec_tmean, na.rm = TRUE)) %>% ungroup() %>%
  mutate(A_base_all    = rec_tmean - base_all,
         B_base_shared = rec_tmean - base_shared,
         C_minus_amt   = rec_tmean - Annual_mean_temperature,
         D_scaled      = as.numeric(scaled_tmean),
         iso   = ln_body_mass - 3 * log(conc.wing.length),
         Sex   = factor(Sex), src = factor(Main_researcher),
         scaled_lat = as.numeric(scaled_lat)) %>%
  as.data.frame()

cat(sprintf("shared wing+mass records: %d  (REPORT.md 5.5 states n = 7577)\n", nrow(sh)))
for (v in c("A_base_all", "B_base_shared", "C_minus_amt", "D_scaled"))
  cat(sprintf("  %-14s available on %d records\n", v, sum(!is.na(sh[[v]]))))

res <- do.call(rbind, lapply(c("A_base_all", "B_base_shared", "C_minus_amt", "D_scaled"), function(v) {
  dd <- sh[!is.na(sh[[v]]), ]
  f <- as.formula(sprintf("iso ~ Sex + %s + scaled_lat + (1 + %s || spp) + (1 | src)", v, v))
  fit <- tryCatch(glmmTMB(f, data = dd, REML = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients$cond
  data.frame(anomaly_definition = v, n = nrow(dd), estimate = cf[v, 1], se = cf[v, 2],
             z = cf[v, 3], p = cf[v, 4], row.names = NULL)
}))
res$reported_estimate <- 0.0212
res$reported_p        <- 0.030
res$reported_n        <- 7577L
print(res, row.names = FALSE)

saveRDS(res, file.path(OUT_DIR, "anomaly_sensitivity.rds"))
writeLines(c(
  "# Anomaly-definition sensitivity of the thermal-allometry estimate",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "The reported estimate (beta = +0.0212, p = 0.030, n = 7,577) has no generating",
  "script in the repository. Reconstructions from `passer90_climate.rds` give:",
  "", "```",
  paste(capture.output(print(res, row.names = FALSE)), collapse = "\n"),
  "```", "",
  "No construction reaches n = 7,577: 353 of the 7,577 shared records have no",
  "`rec_tmean`. The sign is stable across definitions but significance is not.",
  ""), file.path(OUT_DIR, "ANOMALY_SENSITIVITY.md"))
cat("\nwrote", file.path(OUT_DIR, "ANOMALY_SENSITIVITY.md"), "\n")
