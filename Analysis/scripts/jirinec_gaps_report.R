# jirinec_gaps_report.R -- render output/jirinec_gaps/jirinec_gaps_results.rds
# to a plain-markdown numbers dump (output/jirinec_gaps/jirinec_gaps_tables.md).
# Companion to jirinec_gap_analysis.R. Not part of the manuscript pipeline.
suppressMessages({ library(dplyr) })
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    p <- dirname(d); if (identical(p, d)) break; d <- p
  }
  stop("Could not locate Analysis/.")
}
AD <- .find_analysis_dir()
gp <- function(...) file.path(AD, "output", "jirinec_gaps", ...)
res <- readRDS(gp("jirinec_gaps_results.rds"))

con <- file(gp("jirinec_gaps_tables.md"), open = "wt")
w  <- function(...) cat(..., "\n", sep = "", file = con)
tb <- function(df, digits = 4) {
  if (is.null(df) || !nrow(df)) { w("_(empty)_\n"); return(invisible()) }
  df <- as.data.frame(df)
  num <- sapply(df, is.numeric); df[num] <- lapply(df[num], function(x) round(x, digits))
  w("| ", paste(names(df), collapse = " | "), " |")
  w("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  for (i in seq_len(nrow(df))) w("| ", paste(format(unlist(df[i, ]), trim = TRUE), collapse = " | "), " |")
  w("")
}
w("# Jirinec-gap analysis - numbers\n")
w("Generated: ", format(res$generated), "  | trees: ", res$n_trees, "  | runtime: ", res$elapsed_min, " min\n")

w("## Climate coverage\n");        tb(res$climate_coverage)
w("## G1 VIF (their max was 2.01)\n")
tb(data.frame(set = rep(names(res$G1_vif), sapply(res$G1_vif, length)),
              term = unlist(lapply(res$G1_vif, names)),
              VIF = unlist(res$G1_vif), row.names = NULL), 3)
w("## G1 year-slope correlation with climate anomalies\n")
tb(data.frame(term = names(res$G1_year_climate_cor), r = as.numeric(res$G1_year_climate_cor)), 3)
w("## G1 year slope, with and without lagged climate\n")
tb(res$G1_year_shift %>% select(trait, model, estimate, lower, upper,
                                pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n), 3)
w("## G1 AIC across climate variants (ML refits)\n"); tb(res$G1_aic, 2)
w("## G1 climate coefficients (all tiers)\n")
tb(res$G1_table %>% filter(term != "scaled_yr") %>%
     select(trait, model, term, estimate, se, lower, upper, t), 5)
w("## G2 species-slope tallies\n");  tb(res$G2_tally, 3)
w("## G3 mass:wing ratio, pooled slope\n")
tb(res$G3_table %>% select(trait, model, estimate, lower, upper,
                           pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n), 4)
w("## G3 mass:wing species tally\n"); tb(res$G3_tally, 3)
w("## G4 phylogenetic signal in slopes\n")
tb(res$G4_slope_signal %>% select(metric, n_spp, moran_I, moran_expected, moran_p_median,
                                  moran_p_max, lambda, lambda_lo, lambda_hi, lambda_p_median), 4)
w("## G5 foraging-stratum interactions\n")
tb(res$G5_table %>% select(trait, model, term, estimate, se, lower, upper, t, n), 5)
w("## G5 species per dominant stratum\n"); tb(res$G5_strata_counts)
w("## G6 design-matched restricted tier\n")
w("Criteria: span >= ", res$G6_criteria$span_min, " yr, n >= ", res$G6_criteria$n_min,
  "; kept ", res$G6_criteria$n_contributors_kept, " contributors\n")
tb(as.data.frame(t(res$G6_n)))
tb(res$G6_table %>% select(trait, model, estimate, lower, upper,
                           pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n), 4)
w("## G6 contributor coverage (top 15 by span)\n"); tb(head(res$G6_contributor_coverage, 15))
w("## G7 symmetric-coverage filter\n"); tb(res$G7_summary)
tb(res$G7_table %>% select(trait, model, estimate, lower, upper,
                           pct_per_decade, pct_per_decade_lower, pct_per_decade_upper, n), 4)
close(con)
message("[write] ", gp("jirinec_gaps_tables.md"))
