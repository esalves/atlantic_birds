# update_descriptive_stats.R
# ---------------------------------------------------------------------------
# Recomputes every descriptive / quartile / lnCVR number quoted in the
# manuscript from the canonical analytical sample (passer90.rda, the 73-species
# live-only build with thresholds n>=30 records and span>=5 years), mirroring the exact
# methodology in atlantic_birds_ms.Rmd. Reads the pooled model objects
# (brm0_multiphylo.rda wing, brm_bill_multiphylo.rda bill) for fixed effects,
# phylogenetic signal, and analysis N.
#
# Writes:
#   descriptive_summary.rds  - all quoted numbers, for index.qmd to read inline
#   passer90_export.csv      - the analytical sample, for make_figures.py
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr)
  library(data.table)
  library(metafor)
  library(posterior)
  library(brms)
})

# --- robust path resolution (repo layout: Analysis/{scripts,data,output,figures}) ---
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts")))
      return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/ (need Analysis/data/derived and Analysis/scripts).")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)

source(script_path("functions.R"))   # s2.lnCVR()
load(derived_path("passer90.rda"))    # canonical analytical sample

out <- list()

# --- raw, pre-threshold counts (Methods sentence) --------------------------
raw <- fread(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"))
out$raw_total       <- nrow(raw)
out$raw_passerine   <- sum(raw$Order == "Passeriformes", na.rm = TRUE)
out$raw_pass_pct    <- round(100 * out$raw_passerine / out$raw_total)
out$n_adult         <- sum(raw$Age == "Adult", na.rm = TRUE)
# adult + passerine + 1990-2018 (pre species/sex/buffer thresholds)
out$n_adult_pass_yr <- raw[Age == "Adult" & Order == "Passeriformes" &
                           Year >= 1990 & Year <= 2018, .N]
rm(raw); invisible(gc())

# --- analytical sample sizes ------------------------------------------------
out$n_spp <- length(unique(passer90$Binomial))
out$n_rec <- nrow(passer90)
out$yr_min <- min(passer90$Year, na.rm = TRUE)
out$yr_max <- max(passer90$Year, na.rm = TRUE)

# --- quartile means table (1st vs 3rd quartile of the record) --------------
passer.quart.years <- passer90 %>% filter(Year <= 2006 | Year >= 2013)
passer.quart.years$period <- ifelse(passer.quart.years$Year <= 2006,
                                     "1990-2006", "2013-2018")

wing.by.spp <- passer.quart.years %>%
  group_by(spp = Binomial, period) %>%
  summarise(mean.wing       = mean(conc.wing.length, na.rm = TRUE),
            sd.wing         = sd(conc.wing.length,   na.rm = TRUE),
            sample.n.wing   = sum(!is.na(conc.wing.length)),
            mean.bill.width = mean(Bill_width.mm.,   na.rm = TRUE),
            sd.bill         = sd(Bill_width.mm.,     na.rm = TRUE),
            sample.n.bill   = sum(!is.na(Bill_width.mm.)),
            .groups = "drop_last") %>%
  mutate(size.trend      = mean.wing - lag(mean.wing),
         size.trend.bill = mean.bill.width - lag(mean.bill.width)) %>%
  ungroup()

# Direction counts: a species "decreased" if late-period mean < early-period
# mean (size.trend computed as late - early within species; NA when a species
# is present in only one quartile -> excluded, as in the Rmd/figures).
out$wing_dec <- sum(wing.by.spp$size.trend      < 0, na.rm = TRUE)
out$wing_inc <- sum(wing.by.spp$size.trend      > 0, na.rm = TRUE)
out$bill_dec <- sum(wing.by.spp$size.trend.bill < 0, na.rm = TRUE)
out$bill_inc <- sum(wing.by.spp$size.trend.bill > 0, na.rm = TRUE)
out$wing_quartile_n <- out$wing_dec + out$wing_inc
out$bill_quartile_n <- out$bill_dec + out$bill_inc

# --- lnCVR meta-analyses (wing & bill) -------------------------------------
dat.lnCVR <- dcast(setDT(wing.by.spp), spp ~ period,
                   value.var = c("mean.wing", "sd.wing", "mean.bill.width",
                                 "sd.bill", "sample.n.wing", "sample.n.bill"))

cor.wing.e <- cor(log(wing.by.spp$mean.wing[wing.by.spp$period == "1990-2006"]),
                  log(wing.by.spp$sd.wing[wing.by.spp$period == "1990-2006"]),
                  use = "complete.obs")
cor.wing.l <- cor(log(wing.by.spp$mean.wing[wing.by.spp$period == "2013-2018"]),
                  log(wing.by.spp$sd.wing[wing.by.spp$period == "2013-2018"]),
                  use = "complete.obs")
cor.bill.e <- cor(log(wing.by.spp$mean.bill.width[wing.by.spp$period == "1990-2006"]),
                  log(wing.by.spp$sd.bill[wing.by.spp$period == "1990-2006"]),
                  use = "complete.obs")
cor.bill.l <- cor(log(wing.by.spp$mean.bill.width[wing.by.spp$period == "2013-2018"]),
                  log(wing.by.spp$sd.bill[wing.by.spp$period == "2013-2018"]),
                  use = "complete.obs")

dat.lnCVR$lnCVR.wing <- log((dat.lnCVR$`sd.wing_2013-2018`/dat.lnCVR$`mean.wing_2013-2018`)/
                            (dat.lnCVR$`sd.wing_1990-2006`/dat.lnCVR$`mean.wing_1990-2006`)) +
  (1/(2*(dat.lnCVR$`sample.n.wing_2013-2018`-1))) - (1/(2*(dat.lnCVR$`sample.n.wing_1990-2006`-1)))
dat.lnCVR$s2.lnCVR.wing <- s2.lnCVR(dat.lnCVR$`mean.wing_1990-2006`, dat.lnCVR$`sd.wing_1990-2006`,
                                    dat.lnCVR$`sample.n.wing_1990-2006`, cor.wing.e,
                                    dat.lnCVR$`mean.wing_2013-2018`, dat.lnCVR$`sd.wing_2013-2018`,
                                    dat.lnCVR$`sample.n.wing_2013-2018`, cor.wing.l)

dat.lnCVR$lnCVR.bill <- log((dat.lnCVR$`sd.bill_2013-2018`/dat.lnCVR$`mean.bill.width_2013-2018`)/
                            (dat.lnCVR$`sd.bill_1990-2006`/dat.lnCVR$`mean.bill.width_1990-2006`)) +
  (1/(2*(dat.lnCVR$`sample.n.bill_2013-2018`-1))) - (1/(2*(dat.lnCVR$`sample.n.bill_1990-2006`-1)))
dat.lnCVR$s2.lnCVR.bill <- s2.lnCVR(dat.lnCVR$`mean.bill.width_1990-2006`, dat.lnCVR$`sd.bill_1990-2006`,
                                    dat.lnCVR$`sample.n.bill_1990-2006`, cor.bill.e,
                                    dat.lnCVR$`mean.bill.width_2013-2018`, dat.lnCVR$`sd.bill_2013-2018`,
                                    dat.lnCVR$`sample.n.bill_2013-2018`, cor.bill.l)

m.wing <- rma(yi = lnCVR.wing, vi = s2.lnCVR.wing, method = "REML", data = dat.lnCVR)
m.bill <- rma(yi = lnCVR.bill, vi = s2.lnCVR.bill, method = "REML", data = dat.lnCVR)

out$lncvr_wing  <- round(as.numeric(m.wing$b), 3)
out$lncvr_wing_lo <- round(m.wing$ci.lb, 3)
out$lncvr_wing_hi <- round(m.wing$ci.ub, 3)
out$lncvr_wing_k  <- m.wing$k
out$lncvr_wing_I2 <- round(m.wing$I2, 1)
out$lncvr_bill  <- round(as.numeric(m.bill$b), 3)
out$lncvr_bill_lo <- round(m.bill$ci.lb, 3)
out$lncvr_bill_hi <- round(m.bill$ci.ub, 3)
out$lncvr_bill_k  <- m.bill$k
out$lncvr_bill_I2 <- round(m.bill$I2, 1)

# Attach variance results slots if available (REVISION_NOTES_P3.md §6)
var_file <- out_path("variance_results.rds")
if (file.exists(var_file)) {
  var_res <- readRDS(var_file)
  if (!is.null(var_res$lncvr)) {
    # within-contributor lnCVR
    r_src <- var_res$lncvr[var_res$lncvr$analysis == "within species x sex x contributor, n >= 5" & var_res$lncvr$model == "rma", ]
    if (nrow(r_src) > 0) {
      out$lncvr_wing_src <- round(r_src$estimate[1], 3)
      out$lncvr_wing_src_lo <- round(r_src$lower[1], 3)
      out$lncvr_wing_src_hi <- round(r_src$upper[1], 3)
      out$lncvr_wing_src_k <- r_src$k[1]
    }
  }
  if (!is.null(var_res$contributors_per_cell_summary)) {
    c_sum <- as.data.frame(var_res$contributors_per_cell_summary)
    out$contrib_per_cell_early <- c_sum$mean_contrib[c_sum$period == "early"]
    out$contrib_per_cell_late  <- c_sum$mean_contrib[c_sum$period == "late"]
  }
  if (!is.null(var_res$variance_decomposition_summary)) {
    v_sum <- as.data.frame(var_res$variance_decomposition_summary)
    out$var_share_early <- v_sum$median_share_between[v_sum$period == "early"]
    out$var_share_late  <- v_sum$median_share_between[v_sum$period == "late"]
  }
}

# --- WING model: pooled fixed effects + phylogenetic signal ----------------
load(out_path("models", "brm0_multiphylo.rda"))          # fits, rubin_summary, clootl_version
pe <- function(df, par, col) round(df[[col]][df$par == par], 2)
out$wing_yr      <- pe(rubin_summary, "b_scaled_yr",  "estimate")
out$wing_yr_lo   <- pe(rubin_summary, "b_scaled_yr",  "lower")
out$wing_yr_hi   <- pe(rubin_summary, "b_scaled_yr",  "upper")
out$wing_sex     <- pe(rubin_summary, "b_SexMale",    "estimate")
out$wing_sex_lo  <- pe(rubin_summary, "b_SexMale",    "lower")
out$wing_sex_hi  <- pe(rubin_summary, "b_SexMale",    "upper")
out$wing_lat     <- pe(rubin_summary, "b_scaled_lat", "estimate")
out$wing_lat_lo  <- pe(rubin_summary, "b_scaled_lat", "lower")
out$wing_lat_hi  <- pe(rubin_summary, "b_scaled_lat", "upper")
out$n_rec_wing   <- nobs(fits[[1]])

wing_sig <- sapply(fits, function(f)
  hypothesis(f,
    "sd_species_name__Intercept^2 / (sd_species_name__Intercept^2 + sd_spp__Intercept^2 + sd_spp__scaled_yr^2 + sigma^2) = 0",
    class = NULL)$hypothesis$Estimate)
out$wing_phylo    <- round(mean(wing_sig), 2)
out$wing_phylo_lo <- round(quantile(wing_sig, .025), 2)
out$wing_phylo_hi <- round(quantile(wing_sig, .975), 2)

# --- BILL model: pooled fixed effects + phylogenetic signal ----------------
load(out_path("models", "brm_bill_multiphylo.rda"))      # bill_fits, bill_rubin_summary, bill_sig
out$bill_yr      <- pe(bill_rubin_summary, "b_scaled_yr",  "estimate")
out$bill_yr_lo   <- pe(bill_rubin_summary, "b_scaled_yr",  "lower")
out$bill_yr_hi   <- pe(bill_rubin_summary, "b_scaled_yr",  "upper")
out$bill_lat     <- pe(bill_rubin_summary, "b_scaled_lat", "estimate")
out$bill_lat_lo  <- pe(bill_rubin_summary, "b_scaled_lat", "lower")
out$bill_lat_hi  <- pe(bill_rubin_summary, "b_scaled_lat", "upper")
out$n_rec_bill   <- nobs(bill_fits[[1]])
out$bill_phylo    <- round(mean(bill_sig), 2)
out$bill_phylo_lo <- round(quantile(bill_sig, .025), 2)
out$bill_phylo_hi <- round(quantile(bill_sig, .975), 2)
out$bill_yr_excl_zero <- (out$bill_yr_lo > 0) || (out$bill_yr_hi < 0)

saveRDS(out, out_path("descriptive_summary.rds"))

# --- export analytical sample for make_figures.py --------------------------
exp_cols <- passer90 %>%
  transmute(Binomial, Year,
            cwl = conc.wing.length,
            Bill_width.mm.,
            Longitude_decimal_degrees, Latitude_decimal_degrees)
data.table::fwrite(exp_cols, derived_path("passer90_export.csv"))

# --- print everything ------------------------------------------------------
cat("\n================ DESCRIPTIVE SUMMARY ================\n")
str(out)
cat("\nWING: beta_yr", out$wing_yr, "[", out$wing_yr_lo, out$wing_yr_hi, "] N", out$n_rec_wing,
    " phylo", out$wing_phylo, "[", out$wing_phylo_lo, out$wing_phylo_hi, "]\n")
cat("BILL: beta_yr", out$bill_yr, "[", out$bill_yr_lo, out$bill_yr_hi, "] N", out$n_rec_bill,
    " excl0:", out$bill_yr_excl_zero, " phylo", out$bill_phylo, "\n")
cat("Quartile wing: dec", out$wing_dec, "inc", out$wing_inc, " (of", out$wing_quartile_n, ")\n")
cat("Quartile bill: dec", out$bill_dec, "inc", out$bill_inc, " (of", out$bill_quartile_n, ")\n")
cat("lnCVR wing:", out$lncvr_wing, "[", out$lncvr_wing_lo, out$lncvr_wing_hi, "] k", out$lncvr_wing_k, "I2", out$lncvr_wing_I2, "\n")
cat("lnCVR bill:", out$lncvr_bill, "[", out$lncvr_bill_lo, out$lncvr_bill_hi, "] k", out$lncvr_bill_k, "I2", out$lncvr_bill_I2, "\n")
