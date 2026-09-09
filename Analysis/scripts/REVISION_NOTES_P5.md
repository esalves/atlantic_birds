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

> **Update 2026-09-09 (agent E4, §8).** The phylogenetic tier has since been
> run in full locally with glmmTMB's `propto` covariance structure across the
> 50 published trees (`--engine glmmTMB`, now the default;
> `output/diet_interaction_phylo.rds`, 694 s). It reproduces the published brms
> interaction (+0.663 [+0.217, +1.108] vs +0.660 [+0.215, +1.106]), shows that
> phylogeny leaves every diet coefficient unchanged relative to the lme4 tier,
> and that the interaction halves and its interval includes zero under
> contributor + municipality controls (+0.319 [−0.019, +0.657]). The brms path
> remains selectable (`--engine brms`) as the Bayesian cross-check for Totoro.

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
| `atlantic_diet_interaction.R` (brms, 10 trees) | **No** | — | Totoro (`--engine brms --trees 10`). Flags `--trees N`, `--smoke` added; see §6. |
| `atlantic_diet_interaction.R --trees 50` (glmmTMB phylogenetic, default engine) | **Yes, in full** (E4, §8) | 694 s (12 runs × 50 trees) | 6 models × 2 species-RE specs; `output/diet_interaction_phylo.rds`, four `figures/diet_*.png`. |
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
  `slope_high`) with the percentile slopes. **Read them from
  `diet_interaction_phylo.rds$quantile_slopes` (phylogenetic, 50 trees,
  Rubin-pooled; §8.3)**; `diet_lme4_results.rds$quantile_slopes` is the
  non-phylogenetic fast tier and agrees to 2 decimals. The brms
  `diet_interaction_results.rds$quantile_slopes` will only exist after a
  Totoro `--engine brms` run and is a cross-check, not a prerequisite.
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
Rscript atlantic_diet_interaction.R --engine brms --trees 10   # brms A, B, A_fam + percentile slopes (Bayesian cross-check; figures/diet_interaction_plot_brms.png)
# the default engine (glmmTMB, 50 trees, ~13 min) has already been run locally: Rscript atlantic_diet_interaction.R
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
| `Analysis/scripts/atlantic_diet_interaction.R` | edited by the first P5 run; verified here; **E4 added the glmmTMB engine (§8.5)** |
| `Analysis/output/diet_interaction_phylo.rds`, `Analysis/figures/diet_interaction_plot.png`, `diet_quantile_predictions_phylo.png`, `diet_interaction_specifications.png`, `diet_species_slopes.png` | **new (E4, §8)** |
| `Analysis/data/derived/passer90_climate.rds` | **regenerated, live-only** (12,571 / 73 / 46 cols; 772 NA `rec_tmean`) |
| `Analysis/output/climate_trends.rds` | new (slots listed in §2; `$settings`, `$session`) |
| `Analysis/figures/climate_trends.png` | new |
| `Analysis/output/diet_distribution.rds`, `diet_lme4_results.rds` | new |
| `Analysis/figures/diet_distribution.png`, `diet_quantile_predictions.png` | new |
| `Analysis/output/diet_interaction_results.rds`, `Manuscript/images/diet_interaction_plot.png` | **untouched** (published 10-tree brms; the glmmTMB engine reads the former for the comparison and never writes it; the manuscript image is P6's) |
| `Analysis/scripts/REVISION_NOTES_P5.md` | this file |

Not touched: raw data, `passer90.rda`, any fitted model object,
`Manuscript/index.qmd`. No git operations.

---

## 8. Phylogenetic tier of the diet interaction (glmmTMB, 50 trees; added 2026-09-09, agent E4)

Companion to `atlantic_diet_interaction.R` (default engine now `glmmTMB`) and
`Analysis/output/diet_interaction_phylo.rds`. Every number below is in that rds
(`$interaction_table`, `$interaction_table_all_specs`, `$quantile_slopes`,
`$category_slopes`, `$varcomp_summary`, `$comparison`, `$phylo_vs_lme4`,
`$species_slopes_first_tree`, `$models[[model]][[spec]]`). Executed locally in
full: R 4.6.0, glmmTMB 1.1.14 (`propto` covariance structure, Williams et al.
2025, bioRxiv 10.64898/2025.12.20.695312), engine `scripts/_phylo_engine.R`,
the **50 species correlation matrices of the published brms analysis**
(`data/derived/phylo_A_50trees.rds`, cached from `brm0_multiphylo.rda`; the
published diet fit used 10 of these). Sample 65 species / 8,086 wing records /
42 contributors / 130 municipalities / 14 families (identical to the lme4 tier;
*Herpsilochmus sellowi* is not in the diet sample, so nothing was dropped for
the tree). Wall time 694 s for 12 model × specification runs of 50 trees (final run,
14:24); the six primary runs took 25–107 s each (`$models[[m]][[spec]]$secs`).

### 8.0 Bottom line

1. **glmmTMB reproduces the published brms interaction.** Model A, 50 trees:
   year × Diet-Inv **+0.663 [+0.217, +1.108]** mm per SD-year per SD Diet-Inv
   (published brms, 10 trees: +0.660 [+0.215, +1.106]); year main effect −0.858
   [−1.327, −0.389] (brms −0.855 [−1.324, −0.387]). Model B: +1.024 [+0.056,
   +1.991] (brms +1.025 [+0.046, +2.004]). All 50 trees converged.
2. **Phylogeny changes nothing.** glmmTMB-with-phylogeny minus lme4-without,
   same model: |Δ estimate| ≤ 0.022 for every year and interaction coefficient
   (largest: Model B interaction 1.024 vs 1.046), SE ratios 0.995–1.004
   (`$phylo_vs_lme4`). The phylogenetic term carries 96 % of the species-level
   variance (SD 24.1 mm vs residual 4.60 mm), yet the slopes are unchanged
   because the interaction is identified within species over time.
3. **Family adds nothing on top of phylogeny.** Model A + (1 | Family):
   +0.665 [+0.221, +1.108]; the family SD (13.6 mm) is carved out of the
   phylogenetic SD (24.1 → 21.0 mm; phylogenetic proportion 0.96 → 0.68) and
   the fixed effects do not move.
4. **The interaction does not survive contributor + municipality controls at
   the 95 % level.** With (1 | Main_researcher) + (1 | Municipality) and
   phylogeny: **+0.319 [−0.019, +0.657]**, z = 1.85 (49/50 trees); with family
   as well: +0.321 [−0.014, +0.656], z = 1.88 (48/50). Model B with controls:
   **+0.387 [−0.335, +1.109]**, z = 1.05 (50/50). The point estimate halves,
   exactly as in the lme4 tier (+0.321), and the year main effect halves with
   it (−0.858 → −0.466 [−0.865, −0.066]; Model B −1.377 → −0.665 [−1.187,
   −0.143]). The interaction inherits the provenance confound documented for
   the primary model in Phase 1; phylogeny does not rescue it.
5. **What is robust and what is not.** In every one of the 12 specifications
   the year slope of obligate insectivores (Diet-Inv 100, 20 species) is
   indistinguishable from zero (−0.09 to −0.15 mm per decade, intervals ±≈ 3
   mm/decade). The manuscript's substantive reading — *insectivores show no
   steeper decline* — stands. The *contrast* (frugivore/omnivore species
   declining faster) is the part that is not robust: it is +0.66 without and
   +0.32 with provenance controls, and only the former excludes zero. Note that
   the low-Diet-Inv slope itself stays negative with the controls
   (p10: −1.82 mm per decade [−3.01, −0.62]); it is the *difference* between
   diet groups whose interval includes zero.

### 8.1 Models and species random-effect specification (a finding about the fit, not the biology)

Fixed effects Sex + scaled_yr × diet + scaled_lat throughout; species year
slope (spp); phylogenetic species intercept `propto(0 + species_name | g, A)`
appended by the engine. Two species specifications were fitted for every model
because the published structure `(1 + scaled_yr || spp)` carries an **iid
species intercept that is redundant with the phylogenetic intercept**: in the
validated wing model its SD is 0.005 mm (`glmmtmb_validation_wing.rds`), and
here 0.005–0.028 mm whenever the fit converges. The likelihood is flat along
the trade-off between the two, so the Hessian is frequently singular in the
variance parameters:

| Model | spec `(1 + yr ‖ spp)` published | spec `(0 + yr | spp)` reduced | Primary (more converged trees) |
|---|---:|---:|---|
| A | 50/50 | 50/50 | reduced |
| A + Family | 34/50 | 28/50 | published |
| A + contributor + municipality | **7/50** | 49/50 | reduced |
| A + Family + contributor + municipality | 48/50 | **0/50** | published |
| B | 50/50 | 50/50 | reduced |
| B + contributor + municipality | **9/50** | 50/50 | reduced |

"Converged" = positive-definite Hessian (`pdHess`) and finite SEs. Two kinds of
failure occur (`$models[[m]][[spec]]$per_tree_interaction`, which tabulates the
interaction estimate and SE by convergence status):

- *Benign boundary*: a redundant SD sits at 0 (iid species intercept, or the
  family intercept in the reduced A + Family + provenance fit, Family SD 0.001),
  fixed effects unchanged to 3 decimals (e.g. A_fam_src_site reduced, 50
  non-pdHess trees: interaction 0.318–0.320, SE 0.172, vs 0.320–0.322 on the
  48 converged published-spec trees).
- *Genuine bad optimum*: the optimizer lands on the mirror-image solution in
  which the iid species SD takes the species variance (16.1 mm) and the
  phylogenetic SD collapses to 0.09 mm, logLik NA (A_src_site published, 43
  trees; fixed effects still 0.3214 but no usable likelihood), or on a
  degenerate family solution with interaction 0.467–0.663 and SE 0.057 (A_fam
  reduced, 22 trees).

**Rule applied**: Rubin pooling uses converged trees only; the primary
specification per model is the one with more converged trees (ties → reduced).
Both specifications are saved (`$interaction_table_all_specs`,
`$quantile_slopes_all_specs`); where both have ≥ 7 converged trees their
pooled interactions agree to ≤ 0.002. The naive all-tree pooling is kept as a
diagnostic (`$models[[m]][[spec]]$naive_pool_all_trees`) and should not be
quoted. Wherever the published structure converges, the iid species SD is
< 0.03 mm, i.e. the two specifications are the same model at the MLE.

### 8.2 Interaction and year main effect (primary spec; mm per SD-year; Rubin-pooled over converged trees)

| Model | Spec | Trees | Interaction | 95 % CI | z | Year | 95 % CI |
|---|---|---:|---:|---|---:|---:|---|
| A (continuous) | reduced | 50 | **+0.663** | +0.217, +1.108 | 2.91 | −0.858 | −1.327, −0.389 |
| A + (1 \| Family) | published | 34 | +0.665 | +0.221, +1.108 | 2.94 | −0.852 | −1.319, −0.385 |
| A + contributor + municipality | reduced | 49 | **+0.319** | −0.019, +0.657 | 1.85 | −0.466 | −0.865, −0.066 |
| A + Family + contributor + municipality | published | 48 | +0.321 | −0.014, +0.656 | 1.88 | −0.463 | −0.860, −0.067 |
| B (Invertebrate vs Other) | reduced | 50 | **+1.024** | +0.056, +1.991 | 2.07 | −1.377 | −2.025, −0.729 |
| B + contributor + municipality | reduced | 50 | **+0.387** | −0.335, +1.109 | 1.05 | −0.665 | −1.187, −0.143 |

Comparators (`$comparison`): published brms A +0.660 [+0.215, +1.106], year
−0.855 [−1.324, −0.387]; B +1.025 [+0.046, +2.004], year −1.376 [−2.030,
−0.722]. lme4 tier (§3.3): A +0.665, A_family_int +0.664, A_src_site +0.321
[−0.015, +0.658], B +1.046 [+0.074, +2.018]. The between-tree variance
component of Rubin's rule is negligible everywhere (interaction estimates
range 0.662–0.664 across the 50 trees in Model A).

Variance components, mean over converged trees (`$varcomp_summary`): Model A
phylogenetic SD 24.1 mm, species year-slope SD 1.74, residual 4.60, phylogenetic
proportion 0.960; with contributor + municipality: phylogenetic 23.4,
contributor 3.12, municipality 1.59, species slope SD 1.24, residual 4.18 (the
species-slope SD shrinks by 30 % once contributor and site intercepts absorb
between-source differences). **Caveat**: with municipality intercepts the
latitude coefficient is no longer identified (scaled_lat −0.65 [SE 0.11] → +0.24
[SE 0.42]); latitude is a municipality-level variable and the site intercepts
absorb it — the same issue Phase 1 handles with the Mundlak decomposition.

### 8.3 Year slope at the observed Diet-Inv percentiles (species-weighted p10 / p50 / p90 = Diet-Inv 10 / 50 / 100; Rubin-pooled; `$quantile_slopes`)

| Model | p10 (5 spp, 735 rec) | p50 (8 spp, 998 rec) | p90 (20 spp, 1,871 rec) |
|---|---|---|---|
| A | −1.784 [−2.533, −1.036]; **−3.55 mm/dec** | −1.011 [−1.483, −0.539]; −2.01 | −0.045 [−0.794, +0.705]; **−0.09** |
| A + Family | −1.781 [−2.526, −1.037]; −3.55 | −1.006 [−1.475, −0.536]; −2.00 | −0.036 [−0.781, +0.710]; −0.07 |
| A + contributor + municipality | −0.912 [−1.512, −0.312]; **−1.82** | −0.539 [−0.942, −0.137]; −1.07 | −0.074 [−0.667, +0.519]; **−0.15** |
| A + Family + contributor + municipality | −0.912 [−1.507, −0.316]; −1.82 | −0.537 [−0.937, −0.138]; −1.07 | −0.069 [−0.658, +0.519]; −0.14 |

Model B by category (`$category_slopes`): Other −1.377 [−2.025, −0.729]
(−2.74 mm/dec) vs Invertebrate −0.353 [−1.071, +0.365] (−0.70); with
contributor + municipality: Other −0.665 [−1.187, −0.143] (−1.33) vs
Invertebrate −0.278 [−0.840, +0.283] (−0.55). mm per decade = slope × 10 /
5.020. These reproduce the lme4 tier (§3.3) to 2 decimals.

Species-specific total year slopes (fixed main effect + the species' own diet
term b_int × z(Diet-Inv) + BLUP, tree 1; `$species_slopes_first_tree`,
`figures/diet_species_slopes.png`): Model A 53 of 65 species negative (13 with
intervals below zero, 2 above), median −1.03 mm/decade, r(slope, Diet-Inv) =
0.41; with contributor + municipality 49 of 65 negative (8 below, 1 above),
median −0.66 mm/decade, r = 0.31. Family means (≥ 3 species), baseline →
with controls, mm/decade: Thraupidae (20 spp, Diet-Inv ≈ 33) −3.26 → −1.27;
Pipridae (7) −1.28 → −0.62; Turdidae (5) −1.85 → −1.53; Emberizidae (3) −1.79 →
−2.05; Thamnophilidae (15, Diet-Inv ≈ 95) −0.50 → −0.28; Conopophagidae (3)
−0.74 → −0.85; Tyrannidae (3) +0.29 → +0.82. The steep-declining low-Diet-Inv
group is essentially the Thraupidae, and it is the group whose decline the
provenance controls shrink most. (The engine's `species_slopes_phylo()`
quantity, year main effect + BLUP only, is kept as `slope_main_plus_blup`; see
§8.5 iv.)

### 8.4 Figures (all in `Analysis/figures/`; nothing written to `Manuscript/images/`)

- `diet_interaction_plot.png` — Model A trajectories at p10/p50/p90 and Model
  B by category (solid: baseline; dashed: + contributor + municipality), change
  relative to 2009.5 with Rubin-pooled slope ribbons; subtitles carry the
  pooled interactions. Replaces the brms-drawn version for P6.
- `diet_quantile_predictions_phylo.png` — A / + Family / + contributor +
  municipality side by side (facet labels give spec and converged trees).
- `diet_interaction_specifications.png` — forest plot of the interaction across
  brms (published), lme4, and glmmTMB (both spp specs, with converged/total).
- `diet_species_slopes.png` — species slopes vs Diet-Inv, with and without the
  provenance terms.
- `diet_quantile_predictions.png` (lme4) and `diet_distribution.png` unchanged.
  The brms engine (`--engine brms`) now writes `diet_interaction_plot_brms.png`
  here instead of `Manuscript/images/diet_interaction_plot.png` (P6 owns the
  manuscript images).

### 8.5 Script changes (`atlantic_diet_interaction.R`)

- `--engine glmmTMB|brms` (default glmmTMB); default `--trees` 50 for glmmTMB,
  10 for brms (the published fit); `--no-brms` / `--no-phylo` unchanged
  (descriptive + lme4 only). The brms code path is byte-for-byte the previous
  one apart from the figure destination, so `--engine brms --trees 10` on
  Totoro remains the Bayesian cross-check and still writes
  `diet_interaction_results.rds`. The glmmTMB path never touches that file.
- The engine is sourced after the path helpers; nothing phylogenetic is
  re-implemented. Derived quantities (slope at a diet value) are computed per
  tree from the fixed-effect vcov and pooled with the engine's `pool_rubin_df`.
- Requests for `_phylo_engine.R` (not edited; owner: Synthesis): (i)
  `run_phylo_trees()` should carry a converged flag into pooling (or expose a
  `converged_only` argument) — `pool_rubin_df()` currently propagates NA SEs
  from singular-Hessian fits into the pooled SE, and a naive pool can include
  bad optima whose fixed effects differ (A_fam: 0.467 vs 0.664); (ii) document
  that the iid species intercept in `(1 + x || spp)` is redundant with the
  `propto` intercept and often fails `pdHess` once further intercepts are added
  — `(0 + x | spp)` is the safer default when the phylogenetic term is present;
  (iii) `tidy_phylo_fit()` column names such as `spp..Intercept.` are awkward
  to select downstream; (iv) `species_slopes_phylo()` returns the year MAIN
  effect + BLUP and ignores any fixed-effect interaction with the slope
  covariate, so in an interaction model its `slope` is orthogonal to the
  moderator by construction (cor with Diet-Inv = 0 to machine precision here);
  the diet script adds b_int × z(Diet-Inv) and its Wald variance per species
  (`slope`), keeping the engine's quantity as `slope_main_plus_blup`. A
  `moderator =` argument (or a `newdata` interface) in the engine would make
  this generic.

### 8.6 For P7

- Quote the phylogenetic tier (§8.2–8.3) in place of the lme4 numbers in §3.3;
  the published brms values may stay as the Bayesian comparator. State: 50
  trees, glmmTMB `propto`, Rubin-pooled over converged trees, and that
  phylogeny leaves every diet coefficient unchanged (|Δ| ≤ 0.02 vs lme4).
- Wording: "insectivorous species show no detectable wing-length trend in any
  specification (−0.1 mm per decade, CI ±3); species with lower invertebrate
  diet fractions decline faster in the baseline model (interaction +0.66 mm per
  SD-year per SD Diet-Inv, CI +0.22 to +1.11), but the contrast halves and its
  interval includes zero once contributor and municipality intercepts are
  added (+0.32, CI −0.02 to +0.66); a family random intercept does not change
  it." Do not describe the interaction as robust.
- Keep the EltonTraits certainty qualification (§3.1) and the bimodality of the
  covariate (§3.2).
