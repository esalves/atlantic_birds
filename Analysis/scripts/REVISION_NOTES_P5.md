# Revision notes — Phase 5 (climate context and diet), 2026-09-09

Companion to `REVISION_PLAN.md` §1 "§3.3" / "§3.4" and §3 "Phase 5", and to
`Analysis/scripts/climate_extraction.R` and
`Analysis/scripts/atlantic_diet_interaction.R`. Records what was run, what was
not, every design decision that is not obvious from the code, and what P6/P7
need. Style follows `REVISION_NOTES_P0.md`.

Executed locally with R 4.6.0 (terra 1.9.27, dplyr 1.2.1, tidyr 1.3.2, lme4
2.0.1, ggplot2 4.0.3, patchwork 1.3.2, readr 2.2.0; SPEI 1.8.1 installed for
the cross-check only). **Everything in the climate script and the `--no-brms`
path of the diet script ran in full; every number below is from a full run.
The brms path of the diet script (Models A, B, A_fam with phylogeny) was NOT
run here — not even as a smoke test — because it needs the tree cloud and
~hours of Stan; it is for Totoro (§6).** Every number quoted here is in
`Analysis/output/climate_trends.rds`, `diet_distribution.rds` or
`diet_lme4_results.rds`.

Provenance of this phase: a first P5 run (interrupted by a spend limit) edited
both scripts, ran the diet `--no-brms` path, regenerated
`passer90_climate.rds` on the live-only sample and produced a first
`climate_trends.rds` / `.png`. This run verified that work, fixed two
statistical problems in the locality-level block (§2.3, §2.4), re-ran both
scripts, and wrote these notes. Nothing from the first run was discarded.

---

## 0. Bottom line

**Climate (Fig. S-climate, `figures/climate_trends.png`).** At the 357 sampled
localities with raster data, annual mean temperature rose **+0.25 °C per
decade** (95 % CI +0.08 to +0.41; +0.57 °C over 1995–2018, CI +0.18 to +0.95),
every one of the 357 per-locality OLS slopes being positive. Warm-quarter
maximum temperature rose by a similar amount but with an interval that touches
zero (**+0.24 °C/decade, CI −0.01 to +0.48**). There is **no supportable
drying trend**: SPEI-12 −0.26 per decade (CI −0.62 to +0.10), annual
precipitation −45 mm/decade (CI −123 to +34). These intervals come from a model
that lets the localities share the regional weather (a year random effect); the
naive locality-replicate model gives t ≈ 47 for Tmean and would have been an
overclaim by a factor of ~17 in the SE.

**"Later" is also "warmer" in the bird sample, mostly because of where later
sampling happened, not when.** Record-level r(Year, rec_tmean) = 0.30
(n = 11,799); split by locality this is 0.13 *within* (the warming a fixed
site experienced) and 0.25–0.32 *between* (later records come from warmer
places — the geographic drift P0 documented, mean latitude −24.8 → −21.9).
Any record-level temperature effect in the drmSEM is therefore mostly a
spatial contrast.

**Diet (`figures/diet_distribution.png`, `diet_quantile_predictions.png`).**
The covariate is bimodal: of 65 species, 20 sit at Diet-Inv 90–100 (1,871
records; Thamnophilidae, Conopophagidae, Tyrannidae, Parulidae) and 11 at 0–10
(1,520 records; Thraupidae, Pipridae, Fringillidae). The old ±1 SD projection
points (23.6 / 57.9 / 92.2) are supported within ±5 Diet-Inv points by 6 / 3 /
2 species (374 / 1,064 / 607 records; the +1 SD end is Thamnophilidae only).
The lme4 analogue reproduces the published brms interaction exactly (+0.665 vs
+0.660 mm per SD-year per SD-diet) and it does **not** move when a family
intercept or family year-slope is added (+0.664; the family slope variance is
estimated at zero). It **halves under the Phase 1 provenance controls**
(contributor + municipality intercepts: +0.321, 95 % CI −0.015 to +0.658,
t = 1.87). At the observed percentiles (species-weighted Diet-Inv 10 / 50 /
100), the year slope is −3.6 / −2.0 / −0.1 mm per decade in the baseline and
−1.8 / −1.1 / −0.1 mm per decade with provenance controls; the p90 (obligate
insectivore) slope is indistinguishable from zero in every specification. The
manuscript's reading — the decline is *not* concentrated in insectivores — is
robust; the *strength* of the "frugivores/omnivores decline more" contrast is
not, once contributor turnover is controlled.

---

## 1. What ran, what did not

| Script / path | Ran here? | Time | Notes |
|---|---|---|---|
| `climate_extraction.R` (both blocks) | **Yes, in full** | ~4 min | All 3 × 348 WorldClim 2.5′ monthly layers cached under `data/raw/worldclim/` since June (13 GB); nothing was downloaded. |
| `atlantic_diet_interaction.R --no-brms` | **Yes, in full** | 10 s | Descriptive block + 7 lme4 models. |
| `atlantic_diet_interaction.R` (brms, 10 trees) | **No** | — | Totoro. Flags `--trees N`, `--smoke` added; see §6. |
| TerraClimate SPEI route | **Cannot run** | — | `data/raw/terraclimate/` holds only `TerraClimate_tmax_1990.nc` (97 MB) and a truncated `_1991.nc` (17 MB) plus two 0-byte files. SPEI needs `ppt` and `pet`; the code path (`spei_terraclimate()`) is written and self-activates when `TerraClimate_{ppt,pet}_1990…2018.nc` are present. |

---

## 2. `climate_extraction.R`

### 2.1 Per-record block (output schema unchanged)

`data/derived/passer90_climate.rds` is now the **live-only build**: 12,571
records / 73 species / 46 columns (the 42 columns of the 2026-09 `passer90.rda`
+ `species_name`, `loc_id`, `rec_tmean`, `scaled_tmean`; `spp` is overwritten
with `species_name` as before). `scaled_yr` / `scaled_lat` stay 1-column
`scale()` matrices; `ID_ABT` is unique, so P3 (`atlantic_variance_sigma.R`)
and P4ef (`atlantic_drmsem.R`), which refuse the stale file and join on
`ID_ABT`, can now use it. `rec_tmean` is NA for **772 records** (98
localities; §2.5). Content is identical to the first P5 run of 10:47 — only
the locality-level block changed.

The previous file (15,332 rows / 89 species, pre-live build, no `ID_ABT`) is
gone; the drmSEM results in `output/drmsem_results_*` that were produced on a
server copy must be regenerated by P4ef/Totoro from this file.

Within-locality interannual range of `rec_tmean`: median 0.50 °C, max 1.58 °C
(156 localities sampled in > 1 year) — the time-resolved series does vary
within site, unlike the static `Annual_mean_temperature`.

### 2.2 Locality-level series

455 unique coordinate localities; 357 have raster data. Per locality and year
1990–2018: annual mean temperature (mean of monthly (tmax+tmin)/2, identical
to `rec_tmean`), warm-quarter Tmax (maximum over the centred 3-month running
mean of monthly Tmax whose central month falls in the year — so a DJF austral
summer is assigned to the January year rather than split by the calendar),
annual precipitation, SPEI-12 (December value = calendar-year water balance;
and the annual mean of the 12 monthly SPEI-12 values).

Trend window 1995–2018 (the bird record), `decade = (year − 2006.5)/10`.

### 2.3 SPEI — what changed and why

*Source.* TerraClimate `ppt`/`pet` are absent (§1), so SPEI is computed from
the WorldClim monthly series: Hargreaves–Samani PET (FAO-56 eq. 52; tmin,
tmax, extraterrestrial radiation from latitude at mid-month), D = P − PET,
12-month sums, standardised per calendar month and locality over 1990–2018
(29 years; a 30-year normal is not available in the cached files).
`climate_trends.rds$spei_source` records this; `spei_terraclimate_available`
is FALSE.

*Distribution — the fix.* The first run fitted the strict three-parameter
log-logistic of Vicente-Serrano et al. (2010) by PWMs and returned NA whenever
the fit was invalid. That estimator requires β = 1/τ₃ > 1, i.e. **positive
sample L-skewness**, and a humid-region 12-month water balance is often
negatively skewed: **993 of 4,284 valid locality × calendar-month cells
(23 %) failed, and 71 of 357 localities dropped out of the SPEI trend** (286
remained). The script now fits the **generalized logistic** from sample
L-moments (Hosking 1990; Hosking & Wallis 1997 A.7), which is exactly what
`SPEI::spei(distribution = "log-Logistic", fit = "ub-pwm")` does internally
(`lmom::pelglo` / `cdfglo`, verified in the 1.8.1 source); the strict
log-logistic is its κ < 0 branch and negative skew is handled by reflection
(κ > 0). Result: **0 real failures** (the 1,176 "failed" cells are exactly 12 ×
the 98 no-raster localities), **1,520 of 4,284 cells (35 %) have κ > 0** — the
cells the old code could not fit — and all 357 localities carry SPEI.
Cross-check on the locality with most bird records (loc 241) against
`SPEI::spei`: r = 1.0000, max |difference| = 9 × 10⁻¹⁵. A synthetic test
(positive-skew, negative-skew, symmetric series) also matched to 10⁻¹⁵.

### 2.4 Trend inference — what changed and why

The first run's pooled model `y ~ decade + (1 + decade || loc_id)` treats 357
localities as independent replicates of the same 24 years. They are not —
they share the regional weather (the black line in the figure *is* the shared
signal) — so its SE is grossly optimistic (Tmean t = 47). Two additions:

1. `+ (1 | year)` in the pooled model (rows with `weighting` ending in
   `"+ year RE"`; `settings$pooled_headline`). The year term absorbs the shared
   anomaly (SD 0.29 °C for Tmean, 0.43 °C for warm-quarter Tmax, 0.62 for
   SPEI) and the trend SE grows 17-fold for Tmean.
2. OLS on the across-locality regional-mean series (`regional_trend`; 24
   points, t-based CI). It reproduces the year-RE estimates and SEs to three
   decimals (as it must for a balanced panel) and reports the lag-1 residual
   autocorrelation, which is ≈ 0 for all four variables (−0.007, −0.029,
   0.054, 0.126), so no autocorrelation correction is needed.

The naive rows are kept in `pooled` for comparison. **Quote the year-RE rows.**

| Variable (per decade, 1995–2018) | Year-RE pooled trend | 95 % CI | t | Change over window (CI) | Naive t |
|---|---:|---|---:|---|---:|
| Annual mean T (°C) | **+0.246** | +0.076, +0.415 | 2.84 | +0.57 (+0.18, +0.95) | 47.0 |
| Warm-quarter Tmax (°C) | +0.235 | −0.014, +0.484 | 1.85 | +0.54 (−0.03, +1.11) | 26.1 |
| SPEI-12 (December) | −0.256 | −0.617, +0.105 | −1.39 | −0.59 (−1.42, +0.24) | −15.8 |
| Annual precipitation (mm) | −44.7 | −123.1, +33.7 | −1.12 | −103 (−283, +78) | −12.6 |

Record-weighted variants (localities weighted by bird records) are within
0.01 °C for the temperatures and slightly larger for SPEI/precipitation
(−0.29, −53 mm), all with the same conclusion.

Per-locality OLS slopes (357 localities): Tmean mean +0.246, median +0.265,
2.5–97.5 % range +0.089 to +0.317, **100 % positive**; warm-quarter Tmax mean
+0.235, median +0.252, range +0.019 to +0.428, 98.3 % positive; SPEI-12 mean
−0.256, range −0.80 to +0.32, 17 % positive; precipitation mean −45 mm, range
−133 to +87, 19 % positive. The locality-slope heterogeneity is real (SD of
Tmean slopes 0.03–0.05 °C/decade) but small relative to the shared trend.

### 2.5 Localities outside the WorldClim land mask

98 of 455 coordinate localities (772 records) fall in cells the CRU-TS
downscaled grid leaves empty: Espírito Santo 91 localities / 376 records
(the Guarapari–Setiba coastal restinga cluster, −40.41, −20.57), Santa
Catarina 4 / 249, Paraíba 2 / 144, Rio de Janeiro 1 / 3. Median distance to
the nearest valid cell 10.8 km (4 within 5 km, 13 within 10 km). These points
are on land; the mask is coarse at the coast. **Decision: no nearest-cell
fill** — the per-record file keeps the same 772 NA as every previous build so
downstream consumers are not silently changed; the list with distances is in
`climate_trends.rds$no_raster_localities` (`loc_id`, coordinates, State,
Municipality, `n_records`, `km_to_nearest_valid_cell`). If P4ef/P7 want these
records in the drmSEM, a nearest-valid-cell fill is a five-line change in the
per-record block and should be recorded as a separate decision.

### 2.6 Record-level year–temperature correlation

`climate_trends.rds$record_level_cor`: r(Year, rec_tmean) = **0.302**
(n = 11,799 records with temperature; 357 localities); within-locality
(both centred on locality means) **0.131**; between-locality (locality means,
record-weighted) **0.317**; between-locality with each locality once 0.250.
The plan's sentence that year and temperature are "uncorrelated at record
level" is not right on this build; the correlation is moderate and mostly
spatial. For P7: a record-level temperature coefficient (drmSEM) is
identified mainly by *where*, not *when*, birds were measured.

### 2.7 Caveats to carry into Methods

- WorldClim historical monthly weather = CRU-TS 4.09 downscaled with WorldClim
  2.1 at 2.5′ (~4.6 km); the interannual signal is CRU's station interpolation,
  sparse in the interior Atlantic Forest, so the series are spatially smooth
  (hence the singular locality-slope variances in the naive models).
- Exposure-window mismatch: capture-year temperature cannot cause a fixed
  adult wing length; these series are regional CONTEXT, not per-individual
  exposure.
- SPEI reference period is 29 years (1990–2018), not a 30-year normal, and
  PET is Hargreaves–Samani (temperature-based), not Penman–Monteith; the
  TerraClimate route (Penman–Monteith PET) is coded for when the files exist.
- Duration in the figure and tables is 1995–2018 = 24 years.

---

## 3. `atlantic_diet_interaction.R`

### 3.1 Structure (first P5 run, verified here)

Flags: `--no-brms` (descriptive + lme4 only; **never touches**
`diet_interaction_results.rds` or `Manuscript/images/diet_interaction_plot.png`),
`--trees N`, `--smoke` (1 tree, 2 chains × 400 iterations; writes
`*_SMOKE.rds/.png` only). The brms block was refactored into
`fit_formula()` but formula, priors (`normal(71,15)` intercept, `normal(0,10)`
b, `cauchy(0,1)` sd/sigma), `SAMPLING`, `SAMPLING_CONTROL`, seed 20240303 and
the tree subsample call (`sample(trees_all, N_TREES)`, no random draws added
before it) are unchanged, so a Totoro re-run reproduces the published 10-tree
fit up to Stan. New brms pieces: Model A_fam (`+ (1 | Family)`), Rubin-pooled
year slopes at the observed percentiles (`pool_rubin_derived()`), the
manuscript figure redrawn at the percentiles instead of ±1 SD, and new fields
appended to `diet_interaction_results.rds` (existing fields kept).

Sample: 65 species / 8,086 wing records with Diet-Inv. Of the 73 species, 7
have no EltonTraits row (*Tangara sayaca, Ceratopipra rubrocapilla,
Myiothlypis leucoblephara, Microspingus cabanisi, Rhopias gularis,
Myrmoderus squamosus, Pyriglena pernambucensis*; 521 records, 4.1 %) and 1
(*Herpsilochmus pectoralis*) has diet data but no wing measurement (73
mass-only records). EltonTraits certainty: A 53, B 2, C 1, D1 8, D2 1
species — "D" means the diet was inferred from congeners; P7 should qualify
EltonTraits with this.

### 3.2 Distribution (`diet_distribution.rds`, `figures/diet_distribution.png`)

| Diet-Inv bin | Species | Wing records | Families |
|---|---:|---:|---|
| 0–10 | 11 | 1,520 | Thraupidae, Pipridae, Fringillidae |
| 10–20 | 6 | 374 | Thraupidae, Pipridae, Turdidae, Passeridae |
| 20–30 | 5 | 309 | Emberizidae, Thraupidae, Fringillidae |
| 30–40 | 5 | 320 | Thraupidae, Pipridae, Turdidae |
| 40–50 | 8 | 998 | Thraupidae, Pipridae, Turdidae, Tityridae |
| 50–60 | 3 | 1,064 | Thraupidae, Thamnophilidae |
| 60–70 | 3 | 599 | Thraupidae, Turdidae, Icteridae |
| 70–80 | 2 | 424 | Thamnophilidae, Cardinalidae |
| 80–90 | 2 | 607 | Thamnophilidae |
| 90–100 | 20 | 1,871 | Thamnophilidae, Conopophagidae, Tyrannidae, Parulidae, Polioptilidae, Thraupidae |

Standardisation used by the models (record-weighted, as in the published fit):
centre 57.9, SD 34.3 (species-weighted: 56.2 / 36.4). Old ±1 SD projection
points 23.6 / 57.9 / 92.2, supported (±5 points) by 6 / 3 / 2 species and
374 / 1,064 / 607 records; the +1 SD point is Thamnophilidae only. The plan's
"11 species / 1,520 records" for −1 SD is the 0–10 bin, not the ±5 window.

Observed percentiles (species-weighted, primary — Diet-Inv is a species
trait): p10 = 10, p50 = 50, p90 = 100 (record-weighted 10 / 60 / 100).
Support within ±5 points: p10 5 species / 735 records (Pipridae, Thraupidae,
Passeridae); p50 8 / 998 (Pipridae, Thraupidae, Tityridae, Turdidae); p90
20 / 1,871 (Thamnophilidae, Conopophagidae, Tyrannidae, Parulidae,
Polioptilidae, Thraupidae).

### 3.3 lme4 analogues (`diet_lme4_results.rds`; REML, `(1 + scaled_yr || spp)`, non-phylogenetic)

Interaction `scaled_yr:diet_inv_std` (mm per SD-year per SD Diet-Inv):

| Model | β_int | SE | t | 95 % CI | β_year | Convergence |
|---|---:|---:|---:|---|---:|---|
| A baseline (lme4 analogue of the published model) | **+0.665** | 0.227 | 2.93 | +0.219, +1.111 | −0.856 | ok |
| + `(1 | Family)` | +0.664 | 0.226 | 2.94 | +0.221, +1.107 | −0.858 | ok |
| + `(1 + scaled_yr || Family)` | +0.664 | 0.226 | 2.94 | +0.221, +1.107 | −0.858 | singular (family slope SD = 0) |
| + `(1 | contributor) + (1 | municipality)` | **+0.321** | 0.172 | 1.87 | −0.015, +0.658 | −0.465 | ok |
| + contributor + municipality + family slope | +0.321 | 0.171 | 1.88 | −0.014, +0.655 | −0.465 | singular (family slope SD = 0) |
| B baseline (`scaled_yr:Invertebrate`) | +1.046 | 0.496 | 2.11 | +0.074, +2.018 | −1.382 | ok |
| B + family slope | +1.006 | 0.518 | 1.94 | −0.009, +2.021 | −1.335 | ok |

Published brms (10 trees, Rubin): A +0.660 [+0.215, +1.106]; B +1.025
[+0.046, +2.004] — the lme4 baselines reproduce them (the phylogenetic term
adds nothing to the interaction). Family intercept SD 15.6 mm absorbs most of
the species-intercept variance (16.5 → 10.6 mm) but leaves the interaction
untouched; the family × year slope variance is zero in both continuous models
(0.28 in B_family_slp), so clade-level slope confounding is not what drives
the interaction. **Contributor (SD 3.1 mm) and municipality (SD 1.6 mm)
intercepts halve it and halve the main year effect** — the same attenuation
Phase 1 documents for the primary model; the interaction inherits the
provenance confound.

Year slope at the species-weighted percentiles (mm per SD-year; mm per decade
uses SD(year) = 5.020):

| Model | p10 (Diet-Inv 10) | p50 (50) | p90 (100) |
|---|---|---|---|
| A baseline | −1.785 [−2.534, −1.037]; −3.56 mm/dec | −1.009 [−1.481, −0.537]; −2.01 | −0.040 [−0.789, +0.710]; −0.08 |
| + family | −1.786 [−2.530, −1.043]; −3.56 | −1.011 [−1.480, −0.542]; −2.01 | −0.043 [−0.787, +0.702]; −0.08 |
| + contributor + municipality | −0.914 [−1.512, −0.316]; −1.82 | −0.539 [−0.940, −0.138]; −1.07 | −0.070 [−0.661, +0.520]; −0.14 |

`figures/diet_quantile_predictions.png` shows these three panels as change
relative to the record midpoint (2009.5) so the ribbon is the slope
uncertainty, not the grand-intercept uncertainty across species.

### 3.4 Decisions

- Percentiles are species-weighted because Diet-Inv is a species trait;
  record-weighted values are stored alongside. p90 = 100 coincides with the
  maximum (20 species sit there).
- `diet_inv_std` keeps the record-weighted `scale()` of the published fit so
  the brms coefficients stay comparable; species-weighted centre/SD stored.
- Support is counted within ±5 Diet-Inv points (EltonTraits is in 10-point
  steps, so this is "the same value").
- The lme4 provenance model uses `Main_researcher` and `Municipality`
  intercepts (Phase 1's M-series terms), not the wing-column proxy or season;
  P1's controlled script is the reference for the full control set.

---

## 4. For P6 (figures)

- `figures/climate_trends.png` (11 × 8 in, 150 dpi) is Fig. S-climate: four
  panels — Tmean anomaly, warm-quarter Tmax anomaly, SPEI-12 (December), and
  the histogram of per-locality slopes. Panel subtitles carry the year-RE
  pooled trend and CI. Replotting from `climate_trends.rds`: `annual`
  (locality × year), `regional_mean`, `pooled` (filter `weighting ==
  "unweighted + year RE"`), `locality_slopes`.
- `figures/diet_distribution.png` (8 × 8.5 in) and
  `figures/diet_quantile_predictions.png` (11 × 5.5 in) are supplementary.
  The manuscript's `images/diet_interaction_plot.png` is only redrawn by the
  brms path on Totoro.

## 5. For P7 (manuscript)

- Replace the ±1 SD text (`index.qmd` lines ~132–133, 202: `slope_low`,
  `slope_high`) with the percentile slopes. Until Totoro runs the brms path,
  `diet_interaction_results.rds$quantile_slopes` does not exist; the lme4
  values are in `diet_lme4_results.rds$quantile_slopes` and reproduce the
  brms baseline. After the Totoro run, read `quantile_slopes` (Rubin-pooled,
  models A and A_fam) and `model_A_family_pooled` from
  `diet_interaction_results.rds`.
- State that the interaction attenuates by half under contributor +
  municipality controls (t = 1.87) and that clade (family) controls do not
  change it; keep the substantive point (insectivores show no steeper decline;
  their slope is ≈ 0 in every specification).
- Qualify EltonTraits: 9 of 65 species have certainty D1/D2 (inferred).
- Climate section: quote the year-RE trends (§2.4 table), duration 1995–2018 =
  24 years, no drying signal, exposure-window caveat, and the within/between
  decomposition of r(Year, rec_tmean) (§2.6). The title's "climate change"
  is supported for temperature (+0.57 °C at the sampled sites) and **not** for
  drought/precipitation.
- `passer90_climate.rds` is live-only; the drmSEM must be re-run on it before
  any of its numbers are quoted (P4ef).

## 6. Totoro run

```bash
cd Analysis/scripts
Rscript climate_extraction.R                       # ~4 min with cached rasters; idempotent
Rscript atlantic_diet_interaction.R --trees 10     # brms A, B, A_fam + percentile slopes + manuscript figure
# optional, if TerraClimate_{ppt,pet}_1990..2018.nc are placed in data/raw/terraclimate/:
Rscript climate_extraction.R                       # SPEI switches to the TerraClimate route automatically
```

The WorldClim cache (13 GB, `data/raw/worldclim/wc2.1_cruts4.09_2.5m_{tmax,tmin,prec}_{1990-1999,2000-2009,2010-2019}/`)
is needed on Totoro or the script will download it (~3 × 3 zips from
geodata.ucdavis.edu).

## 7. Files

| File | Status |
|---|---|
| `Analysis/scripts/climate_extraction.R` | edited (locality block; §2.3–2.6) |
| `Analysis/scripts/atlantic_diet_interaction.R` | edited by the first P5 run; verified, unchanged here |
| `Analysis/data/derived/passer90_climate.rds` | **regenerated, live-only** (12,571 / 73 / 46 cols; 772 NA `rec_tmean`) |
| `Analysis/output/climate_trends.rds` | new (slots listed in §2; `$settings`, `$session`) |
| `Analysis/figures/climate_trends.png` | new |
| `Analysis/output/diet_distribution.rds`, `diet_lme4_results.rds` | new |
| `Analysis/figures/diet_distribution.png`, `diet_quantile_predictions.png` | new |
| `Analysis/output/diet_interaction_results.rds`, `Manuscript/images/diet_interaction_plot.png` | **untouched** (published 10-tree brms; Totoro) |
| `Analysis/scripts/REVISION_NOTES_P5.md` | this file |

Not touched: raw data, `passer90.rda`, any fitted model object,
`Manuscript/index.qmd`. No git operations.
