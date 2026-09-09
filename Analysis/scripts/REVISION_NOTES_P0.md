# Revision notes — Phase 0 (provenance and design audit), 2026-09-09

Companion to `REVISION_PLAN.md` §3 "Phase 0". Records what changed in the data
pipeline, the counts, the reference numbers that did or did not reproduce, and
what downstream agents (P1–P7) need to know. Style follows
`LIVE_ONLY_MIGRATION.md`.

Executed locally with R 4.6.0 (dplyr, readr, tidyr, ggplot2 4.0.3, sf 1.1.2,
patchwork 1.3.2). Both scripts run in under a minute; nothing here is a smoke
test — every number below is from the full run.

---

## 1. `rebuild_passer90_live.R` — what changed

The old build dropped `Sex == "Unknown"` at the first filter and kept 15
columns, which made the unknown-sex replication (M6), the provenance / site /
season / individual controls and the multi-trait check impossible from the
derived data. The script now:

1. Filters `Year >= 1990, Age == "Adult", Order == "Passeriformes",
   Status == "live", AtlanticForests_20km_Buffer == "inside the 20 km polygon"`
   — **no sex filter**.
2. Adds `known_sex = !is.na(Sex) & Sex != "Unknown"`. **Note the NA case**:
   406 pre-threshold records have `Sex` NA (not the string "Unknown"); the old
   `filter(Sex != "Unknown")` silently dropped those as well, so `known_sex`
   must treat NA as unknown to reproduce the 12,571-row sample. `Sex` itself is
   left as recorded (NA stays NA; 130 NA-sex records survive the species filter).
3. Computes species thresholds (n ≥ 30, span ≥ 5 yr) on the **known-sex**
   sample only, as before → the same 73 species.
4. Standardises `scaled_yr` / `scaled_lat` with `scale()` on the
   **pre-threshold known-sex sample (n = 14,048)**, exactly as the original
   `newdataframes` chunk did (centre 2009.487, SD 5.0199 for Year; 19.991 /
   7.902 for −latitude). The same centre/SD is applied to the unknown-sex
   records, so year coefficients from the two files are directly comparable.
   New `scaled_lon` (centre −44.584, SD 5.840) and `scaled_alt` (centre 412.4 m,
   SD 389.2 m) use the same reference sample. Constants are stored as
   `attr(passer90, "scaling")` and in `effect_scale.rds$scaling_reference`.
5. Saves two objects:

| File | Object | Rows | Species | Sex |
|---|---|---:|---:|---|
| `data/derived/passer90.rda` | `passer90` | **12,571** | **73** | Female / Male only (`known_sex` all TRUE) |
| `data/derived/passer90_allsex.rda` | `passer90_allsex` | **18,708** | 73 | 12,571 known + 6,007 "Unknown" + 130 NA |

**Regression check (built into the script):** all 15 columns of the previous
`passer90.rda` are reproduced **identically** (values and column types;
`scaled_yr` / `scaled_lat` remain 1-column `scale()` matrices with the
`scaled:center` / `scaled:scale` attributes). Wing non-NA 8,478; bill-width
3,209; mass 11,256 — unchanged. Downstream scripts that `load(passer90.rda)`
(`atlantic_parallel.R`, `atlantic_parallel_bill.R`, `atlantic_diet_interaction.R`,
`climate_extraction.R`, `update_descriptive_stats.R`, `atlantic_drmsem.R` via
`passer90_climate.rds`) need no edits.

### Columns now carried by both files (42)

| Group | Columns |
|---|---|
| Identity | `ID_ABT` (raw record id, unique), `spp`, `Binomial` (underscored), `Family`, `Sex`, `known_sex`, `Status` |
| Time | `Year`, `scaled_yr`, `Date` (raw string), `date_year`, `month`, `day`, `season` (factor DJF/MAM/JJA/SON), `Hour` (raw), `hour_num` (decimal hours) |
| Provenance / site | `Main_researcher`, `State`, `Municipality`, `Locality`, `Ring`, `Recapture` |
| Geography | `Altitude`, `scaled_alt`, `Longitude_decimal_degrees`, `scaled_lon`, `Latitude_decimal_degrees`, `scaled_lat` |
| Wing | `wing_col` ("right" / "left" / "generic" = which column supplied the coalesced value; NA if no wing), `Wing_length_right.mm.`, `ln_wing_length`, `conc.wing.length`, `ln_conc_wing_length` |
| Other traits | `Body_mass.g.`, `ln_body_mass`, `Bill_width.mm.`, `Bill_length.mm.`, `Tail_length.mm.`, `Tarsus_length.mm.` |
| Condition | `Molt`, `Reproductive_stage` |
| Legacy | `Annual_mean_temperature` (static WorldClim; drmSEM fallback) |

`ID_ABT`, `Family`, `State`, `day`, `date_year` are extras beyond the plan's
list (cheap, and `Family` is needed for the P5 clade check).

### Parsing decisions

- `Date` is `m/d/yyyy`. `month` is taken from the first numeric field by regex,
  so `"1/1/2006-2007"` (the one non-standard value) gives month 1 rather than
  NA. In the known-sex sample `Date` is NA for 6 records and `month` is NA for
  exactly those 6 (no parse losses). **62 records have a 4-digit year in `Date`
  that differs from `Year`**; `Year` is retained as the analysis variable, and
  the mismatch is recorded in `audit_season.rds$date_year_vs_Year_mismatch`
  for P1/P7 to mention if they wish.
- `Hour`: `"H:MM"` / `"H:MM:SS"` → `hour_num`; the free-text values
  `"morning"` and `"afternoon"` become NA (kept verbatim in `Hour`). In the
  known-sex sample 9,615 records have `Hour`, 9,580 parse.
- `wing_col` follows the coalesce order right > left > generic. 410 records
  have both right and left; none have right and generic.

---

## 2. `audit_provenance.R` — outputs and headline numbers

Conventions: early = `Year <= 2006`, late = `Year >= 2013` (the quartile
comparison in `update_descriptive_stats.R`); "wing records" = non-NA
`conc.wing.length` (8,478, the primary-model sample); 5,240 wing records fall
in the two quartile periods (2,370 early, 2,870 late).

Every `audit_*.rds` / `effect_scale.rds` is a list of scalars and data.frames
and carries a `plan_checks` table (plan value, reproduced value, match,
definition). `audit_sources.rds$plan_checks` holds the full 53-row table.
**49 of 53 reference numbers in REVISION_PLAN.md reproduce; the 4 that do not
are the locality counts and one mass conversion (see §3).**

### `audit_sources.rds`
- 47 contributors (`Main_researcher`) in the 12,571-record sample, **42 among
  wing records**, 0 NA.
- **6 of 42 wing contributors span both periods** (E. Carrano, A. Ross,
  M. Alves, C. Fontana, A. Piratelli, L. Bugoni), contributing **2,828 of 8,478**
  wing records.
- Long-running contributors (≥ 8 distinct years, ≥ 200 wing records): 6,
  **3,017 wing records** (E. Carrano 1,040; A. Ross 618; M. Alves 433;
  C. Fontana 388; A. Piratelli 295; A. Bispo 243). Five of the six used the
  right-wing column exclusively.
- Wing-column protocol proxy: **87.4 % right-wing column early → 67.5 %
  unspecified column late** (`wingcol_by_period`; per-year table
  `wingcol_by_year`).
- Contributors per species × period cell (four definitions in
  `mean_contributors_per_cell`): wing records 3.7 → 5.5; all records,
  species present in both periods 4.66 → 6.65 (the plan's 4.5 → 6.6 is
  closest to the latter).
- Species with wing records in both quartile periods: **61 of 73** (`species_by_period`).
- Unknown-sex replication sample (`unknown_sex`): 6,137 unknown/NA-sex records,
  **3,532 with wing** (3,436 labelled "Unknown" + 96 with `Sex` NA); 3,532 have
  complete `Main_researcher` / `Municipality` / latitude.
- Tables for figures: `contributors_per_year`, `contributor_table` (per
  contributor: n, span, years, period flags, wing-column mix, mean wing,
  recapture share, ring completeness), `records_per_contributor_year_wing`.

### `audit_sites.rds` (quartile wing records unless suffixed `_all`)
- Named localities: **139, 6 in both periods, 684 records at shared
  localities**; 111 quartile wing records have `Locality` NA.
- Municipalities: **99, 15 shared, 1,728 records**.
- Species × municipality: **821 combinations, 61 in both periods**.
- Coordinate sites (exact lon/lat pairs): 233; 93 early, 144 late, **4 in
  both** (543 records).
- Full sample: 192 named localities (1,384 NA), 153 municipalities, 455
  coordinate sites.
- Geography by period (`geography_by_period_wing`): mean latitude −24.8° early →
  −21.9° late (SD 3.7 → 6.7); mean longitude −47.5 → −46.2; median altitude
  175 m → 489 m (mean 391 → 445; `Altitude` NA for 122 early, 10 late records).
  Record-level correlations with year (wing records): lon 0.10, lat 0.17, alt
  0.07 — "later" is also "further north and higher".
- `map_sites`: one row per coordinate site with `n_early`, `n_late`, `shared`
  — the data behind the redesigned Fig. 1.

### `audit_individuals.rds`
- `Ring` present for **90.9 %** of records (9,571 distinct strings).
- **36 ring strings occur on more than one species**, so an individual is
  defined as `Ring × Binomial` (9,615 individuals). **927 individuals have > 1
  record (1,813 repeat records)**; on the raw ring string it is 955 / 1,857.
  P1 should build the ring random effect on `interaction(Ring, Binomial)`.
- `Recapture == "Yes"` share by period, wing records, NA counted as not a
  recapture: **4.8 % → 12.3 %** (of recorded values 6.8 % → 14.8 %; all
  records 7.7 % → 20.9 %). `Recapture` is NA for 3,206 records.
- First-capture-only sample (earliest record per individual; unringed kept):
  10,758 records, **7,734 with wing**.

### `audit_season.rds`
- `Date` NA for 6 records; `month` NA for the same 6; one non-standard value
  (`1/1/2006-2007`) parsed; 62 `date_year ≠ Year`.
- Season of wing captures shifts: DJF 12.6 % early → 21.6 % late; SON 28.9 %
  → 23.8 % (`season_by_period_wing`).
- Descriptive species × sex-centred wing means by season (no model): DJF −0.55
  mm, MAM +0.21, JJA −0.04, SON +0.25 (SE ≈ 0.1). Moulting birds −0.24 vs
  non-moulting +0.10 mm. **These are descriptive only; the modelled effects are
  Phase 1.**
- `Molt` recorded for **73.0 %** overall, but 43 % NA early vs 13 % NA late.
- Sex ratio (wing records): **45.6 % female early → 48.7 % late** (all records
  47.8 % → 48.5 %).
- `Hour` unparsed values: "morning", "afternoon".

### `audit_traits.rds` (`trait_table`)

| Trait | n | % of 12,571 | contributors | municipalities | contributors in both periods |
|---|---:|---:|---:|---:|---:|
| Wing (coalesced) | 8,478 | 67.4 | 42 | 130 | 6 |
| Body mass | 11,256 | 89.5 | 47 | 150 | 6 |
| Bill width | 3,209 | 25.5 | 22 | 84 | **2** |
| Bill length | 7,697 | 61.2 | 42 | 120 | 4 |
| Tail length | 8,872 | 70.6 | 42 | 129 | 6 |
| Tarsus length | 4,361 | 34.7 | 27 | 62 | 3 |
| Hour (parsed) | 9,580 | 76.2 | 35 | 112 | 4 |

Bill width: all 3,209 values come from the unspecified `Bill_width.mm.`
column; the raw file has 367 `Bill_width_base` and 45 `Bill_width_nostril`
values for these records that the pipeline does not use. Records with both wing
and mass: **7,577** (42 contributors).

### `effect_scale.rds`
- `scaled_yr` SD actually used = **5.020** (pre-threshold known-sex sample);
  SD(Year) of the final 12,571 sample = 5.047 (the plan's "5.05"). Both
  conversions are stored (`using_scaling_sd`, `using_sample_sd`); the
  difference is 0.5 % of the effect. **Use `using_scaling_sd`** — that is the
  unit the coefficient is per.
- Published wing β = −0.92 [−1.40, −0.44] per SD-year (read from
  `descriptive_summary.rds`) → **−0.18 mm/yr, −1.83 mm/decade, −2.58 %/decade,
  −4.2 mm (−5.9 %) over the 23-year difference 1995→2018** (interval −9.0 to
  −2.8 % over record). Mean wing 71.13 mm (n = 8,478).
- Published log-mass β = −0.0073 [−0.0168, +0.0022] (read from
  `mass_results.rds`) → **−1.44 %/decade, −3.3 % over record, interval −7.4 %
  to +1.0 %**.
- Isometric expectation for a −5.93 % wing change: **−16.7 % mass**
  (`isometry`).
- Bill width: +0.24 per SD-year → +0.48 mm/decade (mean 7.93 mm) — for the
  retraction text.
- `duration_years_difference` = 23, `duration_calendar_years` = 24: the
  over-record conversions use 23 (2018 − 1995). The plan's "24 years" is the
  inclusive count; say "1995–2018 (24 calendar years)" but multiply slopes by 23.
- `comparators`: Jirinec 2021, Weeks 2020, Ryding 2024 rows are **NA
  placeholders with `verified = FALSE`** — I could not check the papers
  offline and refused to type recollected numbers into a file the manuscript
  reads. P7 must fill them from the papers and flip `verified`.

### Figures (`Analysis/figures/`)
- `audit_records_per_year_by_source.png` — stacked wing records per year for
  the 12 largest contributors + "other (30)", with contributors per year below.
- `audit_sites_map_shared.png` — early vs late panels on the
  `Manuscript/data/south_america.geojson` basemap; shared coordinate sites as
  red triangles, period-unique as blue circles sized by records; N records /
  sites / shared per panel. Candidate replacement for Fig. 1 (P6 may restyle
  from `audit_sites.rds$map_sites`).
- `audit_wingcol_by_year.png` — share of wing records by column populated,
  per year, with n above bars.

---

## 3. Reference numbers in REVISION_PLAN.md that did NOT reproduce

| Plan | Reproduced | Why | Action |
|---|---|---|---|
| 140 named localities, 7 shared, 795 records | **139, 6, 684** | the plan counted NA `Locality` (111 records) as a locality present in both periods | use 139 / 6 / 684 in the manuscript; `locality_wing_NA_as_level` in `audit_sites.rds` documents the inflated version |
| mass lower CI → −7.6 % over record | **−7.4 %** | the plan used the linear approximation 100·b·23/SD = −7.67; the exact `exp()` conversion is −7.38 | quote −7.4 % |

Numbers that reproduce only under a specific definition (now recorded in the
`definition` column): 42 contributors = wing records (all records: 47); 927 /
1,813 repeat rings = Ring × Binomial (raw string 955 / 1,857); recapture
5 % → 12 % = wing records with NA as "no" (4.8 → 12.3 %); sex ratio 46/54 →
49/51 = wing records; 3,436 unknown-sex records = the "Unknown" label only
(3,532 including NA sex); 4.5 → 6.6 contributors per cell = all records,
species in both periods (4.66 → 6.65; wing-only 3.7 → 5.5).

---

## 4. What downstream agents must know

- **P1 (`atlantic_parallel_controlled.R`)**: `load(derived_path("passer90.rda"))`
  now gives all covariates: `Main_researcher` (src), `Municipality` (site;
  `Locality` is 11 % NA, 1,384 records), `wing_col`, `scaled_lon`, `scaled_alt`
  (`Altitude` NA for 205 known-sex records — check the complete-case N),
  `season`, `Molt`, `hour_num`, `Ring`. For M6 load `passer90_allsex.rda` and
  filter `!known_sex` (3,532 wing records) or `Sex == "Unknown"` (3,436);
  `scaled_yr` is on the same scale as the known-sex file. Ring random effect:
  `interaction(Ring, Binomial)`. First captures: earliest record per
  Ring × Binomial, unringed kept (7,734 wing records). `scaled_yr` / `scaled_lat`
  in `passer90` are 1-column matrices (as before); in `passer90_allsex` and for
  `scaled_lon` / `scaled_alt` they are plain numeric. Coerce with
  `as.numeric()` if a formula helper complains.
- **P2 (mass / multitrait)**: `Body_mass.g.`, `ln_body_mass` and the other
  traits are now in `passer90.rda`; `passer90_mass.rda` (from
  `atlantic_parallel_mass.R`) is no longer the only route to mass. 7,577
  shared wing + mass records.
- **P3 (variance)**: `audit_sources.rds$mean_contributors_per_cell` and
  `contributors_per_species_period_cell` give the cell heterogeneity numbers;
  sex ratio by period in `audit_season.rds`.
- **P5**: `Family` is in both files for the clade sensitivity.
- **P6 / P7**: every number above is in the rds files; the plan's
  locality counts and the −7.6 % must be corrected as in §3; comparators need
  filling. `effect_scale.rds$note` explains the 5.02 vs 5.05 SD.
- `passer90_climate.rds` remains **stale** (pre-live build); `climate_extraction.R`
  must be re-run (P5) — it will inherit the new columns automatically.
- Run order (Step 0–1 of REVISION_PLAN §4) is unchanged:
  `Rscript Analysis/scripts/rebuild_passer90_live.R` then
  `Rscript Analysis/scripts/audit_provenance.R`.
