# secular_threshold_sensitivity.R
# ---------------------------------------------------------------------------
# Follow-up to SCRIPT_CONSTANTS_AUDIT.md section B.
#
# secular_museum_expansion.R applies several undocumented thresholds that decide
# which data enter the secular models. This script tests whether the reported
# secular results depend on them. Nothing here changes the main analysis; it
# only quantifies the exposure.
#
#   S1. SPECIES-INCLUSION RULE. Tier 4 (museum only) uses n >= 10 & span >= 10;
#       Tier 5 (integrated live + museum) uses n >= 30 & span >= 15; the main
#       analysis uses n >= 30 & span >= 5. Refit both tiers under all three
#       rules.
#   S2. 4-SD TRIMS. The secular annual-mean series for wing, mass and the
#       isometry contrast are trimmed at |deviation| <= 4 SD before the
#       breakpoint regressions. Refit those regressions untrimmed, and at 3 and
#       5 SD, with the breakpoint fixed at 1980 as in the source script.
#   S3. 1960 START. The mass and allometry secular series begin at 1960 while
#       the wing series begins at 1880. Refit with alternative start years.
#
# OUTPUT: output/referee_reruns/secular_threshold_sensitivity.{rds,md}
# RUN:    Rscript Analysis/scripts/secular_threshold_sensitivity.R
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(readr); library(lme4) })
options(width = 150)

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
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
OUT <- out_path("referee_reruns"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

local({
  e <- new.env(); load(derived_path("passer90.rda"), envir = e)
  SD_YR  <<- as.numeric(attr(e$passer90$scaled_yr, "scaled:scale"))
  CTR_YR <<- as.numeric(attr(e$passer90$scaled_yr, "scaled:center"))
})
DEC <- 10 / SD_YR
note <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

note("loading raw ATLANTIC BIRD TRAITS ...")
birds_raw <- read_csv(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"),
                      guess_max = 70000, show_col_types = FALSE)
d_base <- birds_raw %>%
  filter(AtlanticForests_20km_Buffer == "inside the 20 km polygon", !is.na(Year)) %>%
  mutate(conc_wing = coalesce(Wing_length_right.mm., Wing_length_left.mm., Wing_length.mm.),
         ln_mass   = log(Body_mass.g.),
         Binomial  = factor(gsub(" ", "_", Binomial)),
         Status    = factor(ifelse(Status %in% c("live", "museum"), Status, NA_character_)),
         site      = factor(Municipality),
         Collection = factor(Collection),
         scaled_yr = (Year - CTR_YR) / SD_YR)

CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
ff <- function(f, d) tryCatch(lmer(f, data = d, control = CTRL), error = function(e) NULL)
grab <- function(fit, id, rule, n_spp, mw) {
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients
  if (!"scaled_yr" %in% rownames(cf)) return(NULL)
  e <- cf["scaled_yr", 1]; s <- cf["scaled_yr", 2]
  data.frame(analysis = id, rule = rule, n = nobs(fit), n_species = n_spp,
             estimate = e, se = s, t = cf["scaled_yr", 3],
             mm_decade = e * DEC, pct_decade = 100 * e * DEC / mw,
             pct_lo = 100 * (e - 1.96 * s) * DEC / mw,
             pct_hi = 100 * (e + 1.96 * s) * DEC / mw, row.names = NULL)
}

# ===========================================================================
# S1. species-inclusion rule
# ===========================================================================
note("S1: species-inclusion thresholds ...")
RULES <- list("n>=10, span>=10 (as run, Tier 4)" = c(10, 10),
              "n>=30, span>=15 (as run, Tier 5)" = c(30, 15),
              "n>=30, span>=5  (main analysis)"  = c(30, 5))

s1 <- list()
for (tier in c("T4_museum_only", "T5_integrated")) {
  pool <- if (tier == "T4_museum_only")
    d_base %>% filter(Status == "museum", Year >= 1880, !is.na(conc_wing))
  else
    d_base %>% filter(Status %in% c("live", "museum"), Year >= 1880, !is.na(conc_wing))
  for (rn in names(RULES)) {
    r <- RULES[[rn]]
    keep <- pool %>% group_by(Binomial) %>%
      summarise(n = n(), span = max(Year) - min(Year), .groups = "drop") %>%
      filter(n >= r[1], span >= r[2]) %>% pull(Binomial)
    d <- pool %>% filter(Binomial %in% keep)
    if (!nrow(d)) next
    f <- if (tier == "T4_museum_only")
      conc_wing ~ scaled_yr + (1 + scaled_yr || Binomial) + (1 | Collection)
    else
      conc_wing ~ scaled_yr + Status + (1 + scaled_yr || Binomial) + (1 | site)
    s1[[length(s1) + 1]] <- grab(ff(f, d), tier, rn, length(keep), mean(d$conc_wing))
  }
  note("  done ", tier)
}
S1 <- bind_rows(s1); print(S1, row.names = FALSE)

# ===========================================================================
# S2. 4-SD trims on the secular annual-mean series
# ===========================================================================
note("S2: SD trims on the annual-mean breakpoint regressions ...")
bp_val <- 1980   # fixed a priori, as in the source script

trim_series <- function(d, value_col, k) {
  d <- d %>% group_by(Binomial) %>%
    mutate(centered = .data[[value_col]] - mean(.data[[value_col]], na.rm = TRUE)) %>% ungroup()
  if (is.finite(k)) d <- d %>% filter(abs(centered) <= k * sd(centered, na.rm = TRUE))
  d %>% group_by(Year) %>% summarise(n = n(), mean_c = mean(centered), .groups = "drop")
}
bp_fit <- function(ann) {
  # a series that starts at (or after) the breakpoint has no pre-period, so the
  # pre-slope term is dropped by lm(); return NA for it rather than failing.
  m <- lm(mean_c ~ pmin(Year - bp_val, 0) + pmax(Year - bp_val, 0), data = ann, weights = ann$n)
  cf <- summary(m)$coefficients
  gv <- function(i, j) if (nrow(cf) >= i) cf[i, j] else NA_real_
  has_pre <- sum(ann$Year < bp_val) > 1
  data.frame(pre_slope_decade  = if (has_pre) gv(2, 1) * 10 else NA_real_,
             pre_t             = if (has_pre) gv(2, 3) else NA_real_,
             pre_p             = if (has_pre) gv(2, 4) else NA_real_,
             post_slope_decade = if (has_pre) gv(3, 1) * 10 else gv(2, 1) * 10,
             post_t            = if (has_pre) gv(3, 3) else gv(2, 3),
             post_p            = if (has_pre) gv(3, 4) else gv(2, 4),
             n_years = nrow(ann), row.names = NULL)
}
series <- list(
  wing = d_base %>% filter(Status %in% c("live", "museum"), !is.na(conc_wing), Year >= 1880) %>%
    group_by(Binomial) %>% filter(n() >= 10) %>% ungroup(),
  mass = d_base %>% filter(Order == "Passeriformes", !is.na(ln_mass), Year >= 1960, Year <= 2018) %>%
    group_by(Binomial) %>% filter(n() >= 10) %>% ungroup(),
  iso  = d_base %>% filter(Order == "Passeriformes", !is.na(conc_wing), !is.na(ln_mass),
                           Year >= 1960, Year <= 2018) %>%
    mutate(iso_contrast = (ln_mass - 3 * log(conc_wing)) * 100) %>%
    group_by(Binomial) %>% filter(n() >= 10) %>% ungroup())
cols <- c(wing = "conc_wing", mass = "ln_mass", iso = "iso_contrast")

s2 <- list()
for (nm in names(series)) for (k in c(3, 4, 5, Inf)) {
  ann <- trim_series(series[[nm]], cols[[nm]], k)
  r <- bp_fit(ann)
  r$series <- nm; r$trim <- if (is.finite(k)) paste0(k, " SD") else "none"
  r$n_records <- if (is.finite(k))
    nrow(series[[nm]] %>% group_by(Binomial) %>%
           mutate(c = .data[[cols[[nm]]]] - mean(.data[[cols[[nm]]]], na.rm = TRUE)) %>% ungroup() %>%
           filter(abs(c) <= k * sd(c, na.rm = TRUE))) else nrow(series[[nm]])
  s2[[length(s2) + 1]] <- r
}
S2 <- bind_rows(s2) %>% select(series, trim, n_records, n_years, everything())
print(S2, row.names = FALSE)

# ===========================================================================
# S3. start year for the mass / allometry secular series
# ===========================================================================
note("S3: start year for the mass and allometry series ...")
s3 <- list()
for (nm in c("mass", "iso")) for (y0 in c(1900, 1940, 1960, 1980)) {
  base <- d_base %>% filter(Order == "Passeriformes", Year >= y0, Year <= 2018)
  base <- if (nm == "mass") base %>% filter(!is.na(ln_mass))
          else base %>% filter(!is.na(conc_wing), !is.na(ln_mass)) %>%
                 mutate(iso_contrast = (ln_mass - 3 * log(conc_wing)) * 100)
  base <- base %>% group_by(Binomial) %>% filter(n() >= 10) %>% ungroup()
  if (nrow(base) < 50) next
  ann <- trim_series(base, cols[[nm]], 4)
  r <- bp_fit(ann); r$series <- nm; r$start_year <- y0; r$n_records <- nrow(base)
  s3[[length(s3) + 1]] <- r
}
S3 <- bind_rows(s3) %>% select(series, start_year, n_records, n_years, everything())
print(S3, row.names = FALSE)

# ===========================================================================
res <- list(S1_species_rule = S1, S2_sd_trim = S2, S3_start_year = S3,
            meta = list(generated = Sys.time(), SD_YR = SD_YR, CTR_YR = CTR_YR,
                        bp_val = bp_val,
                        note = "Sensitivity of the secular results to undocumented thresholds in secular_museum_expansion.R (SCRIPT_CONSTANTS_AUDIT.md section B)."))
saveRDS(res, file.path(OUT, "secular_threshold_sensitivity.rds"))

fmt <- function(d) { d <- as.data.frame(d)
  d[] <- lapply(d, function(x) if (is.numeric(x)) signif(x, 3) else x)
  paste(capture.output(print(d, row.names = FALSE)), collapse = "\n") }
writeLines(c(
  "# Secular threshold sensitivity", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
  "Follow-up to SCRIPT_CONSTANTS_AUDIT.md section B. Tests whether the secular",
  "results depend on thresholds that were set without a stated rationale.", "",
  "## S1. Species-inclusion rule", "", "```", fmt(S1), "```", "",
  "## S2. SD trim on the annual-mean series (breakpoint fixed at 1980)", "", "```", fmt(S2), "```", "",
  "## S3. Start year for the mass / allometry series", "", "```", fmt(S3), "```", ""),
  file.path(OUT, "SECULAR_THRESHOLD_SENSITIVITY.md"))
note("wrote ", file.path(OUT, "SECULAR_THRESHOLD_SENSITIVITY.md"))
note("DONE")
