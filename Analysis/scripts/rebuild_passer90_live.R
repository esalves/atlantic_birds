# rebuild_passer90_live.R
# ---------------------------------------------------------------------------
# WHAT: Regenerates the LIVE-ONLY analytical sample of adult Atlantic Forest
# passerines from the raw ATLANTIC BIRD TRAITS csv and writes TWO derived files:
#
#   data/derived/passer90.rda         -> `passer90`        known-sex records only
#                                        (12,571 records, 73 species; the sample
#                                        every published model was fitted on)
#   data/derived/passer90_allsex.rda  -> `passer90_allsex` the SAME 73 species,
#                                        unknown-sex records included, with the
#                                        SAME scaled_yr / scaled_lat centring
#                                        (for the unknown-sex replication, M6)
#
# WHY (2026-06, live only): museum specimens shrink/shift on preservation and are
# concentrated in the early years (see body_mass_descriptive.qmd), so all main
# analyses use Status == "live". The filter is applied BEFORE the n >= 30 /
# span >= 5 species thresholds so every retained species has >= 30 live records
# (89 -> 73 species, 15,332 -> 12,571 records).
#
# WHY (2026-09, revision Phase 0; REVISION_PLAN.md §3 "Data pipeline update"):
# the referee asked for provenance (contributor), site, season, protocol and
# individual (recapture) controls, plus an unknown-sex replication and a
# multi-trait check. The old build dropped Sex == "Unknown" at the first filter
# and kept only 15 columns, so none of that could be run from the derived data.
# This build keeps all adult live passerines inside the buffer, flags
# known_sex, and carries every metadata field and morphological trait the
# revision needs. The KNOWN-SEX sample is unchanged: same rows, same species,
# same values in every column the downstream scripts already read.
#
# INVARIANTS (asserted below):
#   * species thresholds (n >= 30, year span >= 5) are evaluated on the KNOWN-SEX
#     sample only, exactly as before, so passer90.rda keeps 12,571 rows / 73 spp;
#   * scaled_yr / scaled_lat are scale()d on the PRE-THRESHOLD known-sex sample
#     (14,048 records; centre 2009.487 / SD 5.020 for Year), exactly as before,
#     and the identical centre/SD is applied to the unknown-sex records so year
#     coefficients from the two files are on the same scale;
#   * scaled_lon / scaled_alt (new) use the same reference sample.
#
# INPUTS:  data/raw/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv
# OUTPUTS: data/derived/passer90.rda, data/derived/passer90_allsex.rda
#
# Downstream readers of passer90.rda (all just load() it): atlantic_parallel.R,
# atlantic_parallel_bill.R, atlantic_diet_interaction.R, climate_extraction.R
# (-> passer90_climate.rds -> atlantic_drmsem.R), update_descriptive_stats.R,
# audit_provenance.R, atlantic_parallel_controlled.R (Phase 1). The mass and
# threshold-sensitivity scripts re-read the raw csv with the same filters.
#
# NON-OBVIOUS DECISIONS (also in REVISION_NOTES_P0.md):
#   * Sex is NA (not "Unknown") for 406 pre-threshold records. The old filter
#     `Sex != "Unknown"` silently dropped NA as well, so known_sex is defined as
#     !is.na(Sex) & Sex != "Unknown". Sex is left as recorded (NA stays NA).
#   * month is parsed from the FIRST numeric field of Date (m/d/yyyy in this
#     file), so "1/1/2006-2007" yields month 1 instead of NA; date_year is the
#     4-digit year in Date and is compared with Year in audit_provenance.R.
#   * hour_num is decimal hours from "H:MM" or "H:MM:SS"; the free-text values
#     "morning" / "afternoon" become NA (they are kept verbatim in Hour).
#   * wing_col records WHICH column supplied conc.wing.length (coalesce order
#     right > left > generic), the only measurement-protocol proxy in ABT.
#   * scaled_yr / scaled_lat stay 1-column matrices (scale() output) in
#     passer90.rda for byte-compatibility with the old file; in
#     passer90_allsex.rda and for the new scaled_lon / scaled_alt they are plain
#     numeric vectors.
#
# RUN:  Rscript Analysis/scripts/rebuild_passer90_live.R
#       (works from the repo root, Analysis/, or Analysis/scripts/)
# Session: R 4.6.0, dplyr 1.1.x, readr 2.1.x (see R_session_info.txt).
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(readr) })
set.seed(20260909)   # nothing stochastic below; recorded for the convention

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
  stop("Could not locate Analysis/ (need Analysis/data/derived and Analysis/scripts). ",
       "Run from inside the atlantic_birds repo.")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)

LIVE_ONLY <- TRUE   # the entire purpose of this script; FALSE reproduces the old 89-spp sample

# Expected counts for the live-only known-sex build (LIVE_ONLY_MIGRATION.md)
EXPECT_N_KNOWN   <- if (LIVE_ONLY) 12571L else 15332L
EXPECT_SPP_KNOWN <- if (LIVE_ONLY) 73L    else 89L

# ---------------------------------------------------------------------------
# 1. Raw data and the population-level filters (sex NOT filtered here)
# ---------------------------------------------------------------------------
birds <- read_csv(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"),
                  guess_max = 70000, show_col_types = FALSE)

base <- birds %>%
  filter(Year >= 1990, Age == "Adult",
         Order == "Passeriformes",
         AtlanticForests_20km_Buffer == "inside the 20 km polygon")
if (LIVE_ONLY) base <- base %>% filter(Status == "live")

# ---------------------------------------------------------------------------
# 2. Helpers: date, hour, season
# ---------------------------------------------------------------------------
# Date is "m/d/yyyy" (US order, verified against Year). A handful of values are
# non-standard, e.g. "1/1/2006-2007"; take month = first field, day = second,
# year = first 4-digit run, so month is not lost when it is parseable.
parse_date_parts <- function(x) {
  x   <- trimws(as.character(x))
  m   <- regmatches(x, regexec("^([0-9]{1,2})/([0-9]{1,2})/([0-9]{4})", x))
  mon <- vapply(m, function(v) if (length(v) == 4) as.integer(v[2]) else NA_integer_, integer(1))
  day <- vapply(m, function(v) if (length(v) == 4) as.integer(v[3]) else NA_integer_, integer(1))
  yr  <- vapply(m, function(v) if (length(v) == 4) as.integer(v[4]) else NA_integer_, integer(1))
  mon[!is.na(mon) & (mon < 1 | mon > 12)] <- NA_integer_
  data.frame(month = mon, day = day, date_year = yr)
}
# Hour is "H:MM" or "H:MM:SS"; free text ("morning", "afternoon") -> NA.
parse_hour <- function(x) {
  x  <- trimws(as.character(x))
  ok <- grepl("^[0-9]{1,2}:[0-9]{2}(:[0-9]{2})?$", x)
  h  <- rep(NA_real_, length(x))
  if (any(ok)) {
    p <- strsplit(x[ok], ":", fixed = TRUE)
    h[ok] <- vapply(p, function(v) as.numeric(v[1]) + as.numeric(v[2]) / 60, numeric(1))
  }
  h[!is.na(h) & (h < 0 | h >= 24)] <- NA_real_
  h
}
season_of <- function(m) {
  factor(dplyr::case_when(m %in% c(12, 1, 2) ~ "DJF", m %in% 3:5 ~ "MAM",
                          m %in% 6:8 ~ "JJA",        m %in% 9:11 ~ "SON",
                          TRUE ~ NA_character_),
         levels = c("DJF", "MAM", "JJA", "SON"))
}

# ---------------------------------------------------------------------------
# 3. Record-level derived variables (identical for both output files)
# ---------------------------------------------------------------------------
dp <- parse_date_parts(base$Date)

base <- base %>%
  mutate(
    known_sex           = !is.na(Sex) & Sex != "Unknown",
    conc.wing.length    = coalesce(Wing_length_right.mm., Wing_length_left.mm., Wing_length.mm.),
    wing_col            = case_when(!is.na(Wing_length_right.mm.) ~ "right",
                                    !is.na(Wing_length_left.mm.)  ~ "left",
                                    !is.na(Wing_length.mm.)       ~ "generic",
                                    TRUE ~ NA_character_),
    ln_conc_wing_length = log(conc.wing.length),
    ln_wing_length      = log(Wing_length_right.mm.),
    ln_body_mass        = log(Body_mass.g.),
    month               = dp$month,
    day                 = dp$day,
    date_year           = dp$date_year,
    season              = season_of(month),
    hour_num            = parse_hour(Hour),
    Binomial            = gsub(" ", "_", Binomial),
    spp                 = Binomial
  )

# ---------------------------------------------------------------------------
# 4. Standardisation constants: PRE-THRESHOLD KNOWN-SEX sample, as before.
#    (The original pipeline called scale() before the species thresholds.)
# ---------------------------------------------------------------------------
ref <- base %>% filter(known_sex)
sc  <- function(v) c(center = mean(v, na.rm = TRUE), scale = sd(v, na.rm = TRUE))
scaling <- list(
  year = sc(ref$Year),
  lat  = sc(ref$Latitude_decimal_degrees * -1),   # southern latitude as positive
  lon  = sc(ref$Longitude_decimal_degrees),
  alt  = sc(ref$Altitude),
  reference_sample = "pre-threshold known-sex live adult passerines inside buffer",
  n_reference      = nrow(ref)
)
z <- function(v, s) (v - s[["center"]]) / s[["scale"]]

# ---------------------------------------------------------------------------
# 5. Species thresholds on the KNOWN-SEX sample (unchanged rule)
# ---------------------------------------------------------------------------
keep_spp <- ref %>%
  group_by(Binomial) %>%
  summarise(n = n(), range = max(Year) - min(Year), .groups = "drop") %>%
  filter(n >= 30, range >= 5) %>%
  pull(Binomial)

# ---------------------------------------------------------------------------
# 6. Column selection shared by both files
# ---------------------------------------------------------------------------
keep_cols <- c(
  "ID_ABT", "spp", "Binomial", "Family", "Sex", "known_sex", "Status", "Year", "scaled_yr",
  "Date", "date_year", "month", "day", "season", "Hour", "hour_num",
  "Main_researcher", "State", "Municipality", "Locality", "Ring", "Recapture",
  "Altitude", "scaled_alt", "Longitude_decimal_degrees", "scaled_lon",
  "Latitude_decimal_degrees", "scaled_lat",
  "wing_col", "Wing_length_right.mm.", "ln_wing_length", "conc.wing.length", "ln_conc_wing_length",
  "Body_mass.g.", "ln_body_mass", "Bill_width.mm.", "Bill_length.mm.",
  "Tail_length.mm.", "Tarsus_length.mm.", "Molt", "Reproductive_stage",
  "Annual_mean_temperature"
)

# ---------------------------------------------------------------------------
# 7a. passer90 (known sex): scale() applied on the pre-threshold known-sex
#     sample -> 1-column matrices with scaled:center / scaled:scale attributes,
#     byte-compatible with the previous build.
# ---------------------------------------------------------------------------
passer90 <- ref %>%
  mutate(scaled_lat = scale(Latitude_decimal_degrees * -1),
         scaled_yr  = scale(Year),
         scaled_lon = z(Longitude_decimal_degrees, scaling$lon),
         scaled_alt = z(Altitude, scaling$alt)) %>%
  filter(Binomial %in% keep_spp) %>%
  dplyr::select(all_of(keep_cols)) %>%
  as.data.frame()

# ---------------------------------------------------------------------------
# 7b. passer90_allsex: same species, unknown-sex records included, same centring
# ---------------------------------------------------------------------------
passer90_allsex <- base %>%
  filter(Binomial %in% keep_spp) %>%
  mutate(scaled_yr  = z(Year, scaling$year),
         scaled_lat = z(Latitude_decimal_degrees * -1, scaling$lat),
         scaled_lon = z(Longitude_decimal_degrees, scaling$lon),
         scaled_alt = z(Altitude, scaling$alt)) %>%
  dplyr::select(all_of(keep_cols)) %>%
  as.data.frame()
attr(passer90_allsex, "scaling") <- scaling
attr(passer90, "scaling")        <- scaling

# ---------------------------------------------------------------------------
# 8. Assertions
# ---------------------------------------------------------------------------
stopifnot(!LIVE_ONLY || all(passer90$Status == "live"))
stopifnot(all(passer90$known_sex))
stopifnot(all(passer90$Sex %in% c("Female", "Male")))
stopifnot(nrow(passer90) == EXPECT_N_KNOWN)
stopifnot(length(unique(passer90$Binomial)) == EXPECT_SPP_KNOWN)
stopifnot(setequal(unique(passer90_allsex$Binomial), unique(passer90$Binomial)))
stopifnot(nrow(passer90_allsex) >= nrow(passer90))
stopifnot(sum(passer90_allsex$known_sex) == nrow(passer90))
# Same centring in both files: scaled_yr of known-sex rows must agree
stopifnot(isTRUE(all.equal(
  as.numeric(passer90$scaled_yr),
  passer90_allsex$scaled_yr[match(passer90$ID_ABT, passer90_allsex$ID_ABT)])))
stopifnot(!anyDuplicated(passer90$ID_ABT), !anyDuplicated(passer90_allsex$ID_ABT))
# month must be parseable wherever Date carries m/d/yyyy at the start
stopifnot(sum(is.na(passer90$month)) == sum(is.na(passer90$Date) |
                                            !grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}", passer90$Date)))

# Regression check against the previous build on disk (columns in common).
prev_file <- derived_path("passer90.rda")
if (file.exists(prev_file)) {
  e <- new.env(); load(prev_file, envir = e)
  old <- e$passer90
  common <- intersect(names(old), names(passer90))
  same <- isTRUE(all.equal(
    old[order(old$Year, old$Binomial, old$conc.wing.length, old$Latitude_decimal_degrees), common],
    passer90[order(passer90$Year, passer90$Binomial, passer90$conc.wing.length,
                   passer90$Latitude_decimal_degrees), common],
    check.attributes = FALSE))
  cat(sprintf("Regression check vs previous passer90.rda (%d rows, %d shared cols): %s\n",
              nrow(old), length(common), if (same) "IDENTICAL" else "DIFFERS (inspect!)"))
  if (!same && nrow(old) == nrow(passer90)) {
    for (cc in common) {
      a <- old[[cc]]; b <- passer90[[cc]]
      if (!isTRUE(all.equal(sort(as.character(a)), sort(as.character(b)))))
        cat("  column differs:", cc, "\n")
    }
  }
}

cat(sprintf("passer90 (%s, known sex): %d records, %d species; wing non-NA %d, bill-width non-NA %d, mass non-NA %d\n",
            if (LIVE_ONLY) "LIVE-only" else "all-status",
            nrow(passer90), length(unique(passer90$Binomial)),
            sum(!is.na(passer90$conc.wing.length)), sum(!is.na(passer90$Bill_width.mm.)),
            sum(!is.na(passer90$Body_mass.g.))))
cat(sprintf("passer90_allsex: %d records (%d known sex, %d unknown/NA sex), %d species; wing non-NA %d\n",
            nrow(passer90_allsex), sum(passer90_allsex$known_sex), sum(!passer90_allsex$known_sex),
            length(unique(passer90_allsex$Binomial)), sum(!is.na(passer90_allsex$conc.wing.length))))
cat(sprintf("Scaling (reference n = %d): Year centre %.3f SD %.4f | -lat centre %.3f SD %.3f | lon centre %.3f SD %.3f | alt centre %.1f SD %.1f\n",
            scaling$n_reference, scaling$year[1], scaling$year[2], scaling$lat[1], scaling$lat[2],
            scaling$lon[1], scaling$lon[2], scaling$alt[1], scaling$alt[2]))
cat(sprintf("Date: %d NA, %d unparsed month among non-NA; Hour: %d non-NA, %d parsed to hour_num\n",
            sum(is.na(passer90_allsex$Date)),
            sum(is.na(passer90_allsex$month) & !is.na(passer90_allsex$Date)),
            sum(!is.na(passer90_allsex$Hour)), sum(!is.na(passer90_allsex$hour_num))))
print(table(wing_col = passer90$wing_col, useNA = "ifany"))

save(passer90,        file = derived_path("passer90.rda"))
save(passer90_allsex, file = derived_path("passer90_allsex.rda"))
cat("Wrote", derived_path("passer90.rda"), "\n")
cat("Wrote", derived_path("passer90_allsex.rda"), "\n")
