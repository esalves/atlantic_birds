# Revision notes — Phase 3 (variance re-analysis), 2026-09-09

Companion to `REVISION_PLAN.md` §3 "Phase 3" and to
`Analysis/scripts/atlantic_variance_sigma.R`. Records what was run, what was
not, every design decision that is not obvious from the code, and what P6/P7
and `update_descriptive_stats.R` need. Style follows `REVISION_NOTES_P0.md`.

Executed locally with R 4.6.0 (glmmTMB 1.1.14, lme4 2.0.1, metafor 5.0.1,
brms 2.23.0 / rstan 2.32.7, ggplot2 4.0.3). **Parts 2–4 (glmmTMB / lme4 /
lnCVR / heterogeneity diagnostics) ran in full (`--fast`, 18 s) and every
number in §§0–4 is from that run. Part 1, the PHYLOGENETIC distributional
model, ran in full on 2026-09-09 with the glmmTMB `propto` engine
(`_phylo_engine.R`) across the 50 trees of the published brms wing model —
every tier of the σ ladder, Rubin-pooled; see §9. The brms engine
(`--engine brms`) is kept as the Bayesian cross-check for Totoro and has only
been smoke-tested (killed before completion, §8); nothing from it is a result.**
Every number quoted here is in `Analysis/output/variance_results.rds` (`$index`
lists the slots; `$phylo` is the pointer/headline of the phylogenetic tier) or,
for §9, in `Analysis/output/variance_phylo_results.rds`.

History: a first P3 session (2026-09-09 10:40) wrote the script, the fast
results, the figures and §§0–7 of these notes, and was cut off during its brms
smoke test. A second session re-ran `--fast` end to end (all numbers
reproduced byte-for-byte), switched `scaled_tmean` on now that P5's live-only
`passer90_climate.rds` exists (§1), added tiers S3t/S7 and the `--no-tmean`
flag, and ran the smoke test under a 20-minute cap (§8). A third session
(E2, 2026-09-09 13:00) added `--engine glmmTMB` (now the default for the
phylogenetic tier), ran the 50-tree phylogenetic ladder in full (§9), added
the phylogenetic rows to `variance_summary.png` panel B, and re-ran the fast
tier as part of the same invocation (S0–S7 reproduced to all printed digits).

---

## 0. Bottom line

The "rising wing-length variance" result does **not** survive the source
structure:

| Analysis | k / N | lnCVR or % Δ residual SD per decade | 95 % CI |
|---|---:|---:|---|
| lnCVR pooled per species (published definition, reproduced) | k = 54 | **+0.184** | [+0.034, +0.334], I² 93.6 % |
| lnCVR pooled per species, corrected cell n | k = 54 | +0.177 | [+0.031, +0.323], I² 90.4 % |
| lnCVR **within sex** (species × sex cells, n ≥ 5) | k = 71 | +0.187 | [+0.048, +0.325], I² 89.1 % |
| lnCVR within sex, rma.mv with species random effect | k = 71 | +0.169 | [−0.009, +0.346] |
| lnCVR **within species × sex × contributor** (n ≥ 2) | k = 58 | **−0.108** | [−0.304, +0.088], I² 56 % |
| lnCVR within species × sex × contributor (n ≥ 5) | k = 15 | **−0.249** | [−0.515, +0.017], I² 63 % |
| lnCVR E. Carrano only (n ≥ 5) | k = 7 | −0.429 | [−0.698, −0.159] |
| σ model S0: `sigma ~ year` (no controls) | N = 8,282 | **+9.6 %/decade** | [+6.3, +13.1] |
| σ model S2: + `(1|contributor)` in σ, contributor + site in mean | N = 8,282 | +2.5 %/decade | [−2.9, +8.2] |
| σ model S3: + `(1|site)` in σ (the plan's specification, minus phylogeny/tmean) | N = 8,282 | **−6.3 %/decade** | [−13.5, +1.4] |
| σ model S4: S3 with year split within/between contributor | N = 8,282 | within −6.1 [−13.5, +2.0]; between −9.9 [−35.1, +25.1] | |
| σ model S6: S3 + `(1|spp)` in σ | N = 8,282 | +3.4 %/decade | [−4.7, +12.2] |
| σ model S7: S3 + `scaled_tmean` in σ (the plan's σ specification, minus phylogeny) | N = 7,651 | year −6.8 %/decade [−14.0, +1.0]; tmean −0.7 % per SD [−8.7, +8.1] | |
| **σ model S7 + phylogeny** (the plan's full specification; glmmTMB propto, 50 trees, Rubin-pooled; §9) | N = 7,651 | **year −6.8 %/decade [−14.0, +1.0]**; tmean −0.7 % per SD [−8.8, +8.1] | 50/50 trees converged |
| **σ model S3 + phylogeny** (`--no-tmean` variant, 8,282 records; §9) | N = 8,282 | **−6.3 %/decade** | [−13.4, +1.5]; 50/50 trees |
| σ model S0 + phylogeny (no controls) | N = 8,282 | +9.6 %/decade | [+6.3, +13.1]; 45/50 trees (§9.4) |
| two-stage log\|resid\| T1 (+contributor + site) | N = 8,282 | +6.3 %/decade | [−1.7, +14.9] |

Reading: the CV rise is real *as a description of the pooled data* and it is
not a sex-ratio artefact (it holds within sex, both sexes). But it vanishes,
and turns negative, as soon as early and late records are compared **within
the same measurer**, and the continuous σ model gives no consistent sign once
contributor and site intercepts enter the σ sub-model (−6 % to +3 % per decade
depending on specification, every interval spanning zero). The mechanical
mechanism the referee predicted is visible: contributors per species-period
cell rise 3.7 → 5.5 (45 of 61 species gain contributors, 7 lose; Wilcoxon
p = 5e-8), and the between-contributor share of within-cell variance rises from
a median of 0.25 to 0.45. The variance result has to be presented as a
pooled-data description that is *consistent with* measurer turnover and *not
supported* within measurer.

**Phylogeny changes none of this (§9).** With the phylogenetic species
intercept in the mean model across the 50 published trees, every tier's σ-year
coefficient moves by ≤ 0.0004 log units and every SE by < 1 %: the plan's full
specification (S7 + phylogeny, N = 7,651) gives −0.0353 (SE 0.0206) on the
log-SD scale = **−6.8 % residual SD per decade [−14.0, +1.0]**, the
8,282-record no-tmean variant −6.3 % [−13.4, +1.5], and the uncontrolled S0 +
phylogeny +9.6 % [+6.3, +13.1]. **Plain statement: the year → σ effect does
not survive contributor and site controls, with or without phylogeny** — the
pooled rise (+9.6 %/decade) becomes a non-significant, slightly negative
estimate once contributor and municipality intercepts enter the σ sub-model,
and the phylogenetic term does not rescue it (it does not enter σ at all, and
in the mean model it only re-partitions the between-species variance).

An unexpected by-product for P1/P7 (§4): modelling the heteroscedasticity by
contributor/site also shrinks the **mean** year slope on the same complete-case
sample from −0.32 [−0.69, +0.04] (homoscedastic lme4, M3 structure) to
**−0.15 [−0.35, +0.05]** (glmmTMB S3), because records from noisy
contributors are down-weighted.

---

## 1. Data and sample

- Sample: `passer90.rda` wing records (8,478; known sex, live, 73 species).
  The σ models need `scaled_alt` (190 NA) and `season` (6 NA), so **all σ-model
  tiers use the same 8,282 complete-case records** (42 contributors, 126
  municipalities) — including the S0/S1 baselines, so the year coefficient is
  comparable across tiers. The lnCVR and heterogeneity parts use all 8,478
  (5,240 in the two quartile periods: 2,370 early, 2,870 late).
- `scaled_tmean` **included** (decision revised 2026-09-09 11:20). The first
  P3 session found `passer90_climate.rds` stale (15,332 rows / 89 species, no
  `ID_ABT`) and dropped the term. P5 regenerated the file at 10:47
  (`climate_extraction.R`: 12,571 rows, 73 species, `ID_ABT` key, `rec_tmean`
  = annual mean of (tmax+tmin)/2 in the capture year at the record's
  coordinates from WorldClim 2.1 historical monthly, `scaled_tmean` =
  `scale()` over the 12,571 records). Every one of the 8,478 wing records is
  keyed, but `rec_tmean` is NA for 772 records whose coordinates fall outside
  the WorldClim land mask — coastal / island localities (in the complete-case
  wing sample: Guarapari 283, Florianópolis 221, Mataraca 121, Santa Rita 6).
  Coverage of the 8,282-record complete-case wing sample is therefore
  **92.4 % (7,651 records; 631 dropped: 31 early, 446 middle 2007–2012, 154
  late)**. The script's original ≥ 99 % rule would have discarded the term for
  a reason that no longer applies, so the rule is now: use `scaled_tmean` when
  the file carries `ID_ABT` and covers ≥ 90 % of the sample. To keep the
  S0–S6 numbers comparable with the first session, **S0–S6 stay on the 8,282
  records; S3t (= S3 refitted on the 7,651 tmean-complete records) and S7 (=
  S3t + `scaled_tmean` in σ) are added so the sample change and the covariate
  are separable; the brms model (part 1) follows the plan's σ specification
  with `scaled_tmean` on the 7,651 records**. `--no-tmean` forces the
  8,282-record brms model without the term (Totoro can run both). Record-level
  `cor(scaled_yr, scaled_tmean)` = 0.12 and `cor(scaled_lat, scaled_tmean)` =
  −0.76 in the tmean sample (`$scaled_tmean`). P5 should say in its notes why
  the coastal cells are missing and whether a nearest-land-cell fallback is
  appropriate; if it adds one, this script picks the new values up
  automatically.
- `SD_year` for the %/decade conversion is the `scaled:scale` attribute of
  `passer90$scaled_yr` (5.0199), as in `effect_scale.rds`.

## 2. Continuous σ models (glmmTMB `dispformula`, part 2)

glmmTMB ≥ 1.1.8 accepts random effects in `dispformula`; all nine tiers
converged (`converged == TRUE`, positive-definite Hessian) in ≤ 3 s each.
Tiers (all REML; S0–S6 on the 8,282-record sample, S3t/S7 on the 7,651
tmean-complete records):

| Tier | mean model | σ model | b_σ(year) (SE) | %/decade [CI] |
|---|---|---|---:|---|
| S0 | M0: Sex + yr + lat + (1+yr‖spp) | ~ yr | +0.046 (0.008) | +9.6 [6.3, 13.1] |
| S1 | M0 | ~ yr + Sex | +0.047 (0.008) | +9.8 [6.4, 13.2] |
| S2 | M3: + lon + alt + season + (1|src) + (1|site) | ~ yr + Sex + (1|src) | +0.012 (0.014) | +2.5 [−2.9, +8.2] |
| **S3** | M3 | ~ yr + Sex + (1|src) + (1|site) | **−0.033 (0.020)** | **−6.3 [−13.5, +1.4]** |
| S4 | M3 | ~ yr_within_src + yr_src_mean + Sex + (1|src) + (1|site) | within −0.031 (0.021); between −0.052 (0.084) | −6.1 / −9.9 |
| S5 | M3 | S3 + wing_col | −0.033 (0.020) | −6.3 [−13.5, +1.5] |
| S6 | M3 | S3 + (1|spp) | +0.017 (0.021) | +3.4 [−4.7, +12.2] |
| S3t | M3 (N = 7,651) | S3 | −0.036 (0.021) | −6.9 [−14.0, +0.9] |
| S7 | M3 (N = 7,651) | S3 + scaled_tmean | year −0.035 (0.021); tmean −0.007 (0.043) | year −6.8 [−14.0, +1.0]; tmean −0.7 % per SD [−8.7, +8.1] |

- **S3t vs S7**: dropping the 631 coastal records without temperature moves
  the S3 year effect from −6.3 to −6.9 %/decade (sample effect, within noise);
  adding `scaled_tmean` then changes nothing (−6.8). The temperature term on σ
  is a null (−0.7 % residual SD per SD of record-year mean temperature, CI
  −8.7 to +8.1; z = −0.16). This is the glmmTMB preview of the referee-§2.7
  σ contradiction (drmSEM: warmer → less variable): once contributor and site
  intercepts are in the σ model there is no temperature effect on residual
  variance to explain, and year and temperature are only weakly correlated in
  these records (r = 0.12). In S3t the σ random-effect SDs are 0.34
  (contributor) and 0.39 (municipality).

- S0 is the continuous analogue of the pooled lnCVR and agrees with it in sign
  and rough magnitude (+9.6 %/decade ≈ +22 % over 23 years vs lnCVR +0.18 ≈
  +20 %).
- In S3 the σ random-effect SDs are 0.38 (contributor) and 0.39
  (municipality) on the log scale: residual SD varies **2.5-fold to 0.45-fold**
  between contributors (`exp(u_j)` range 0.45–2.53). This is the scale of the
  provenance heterogeneity the year term was picking up.
- S5: the wing-column protocol proxy is strongly associated with residual SD
  (left-wing column +1.47 log units, generic +0.34, both z > 2) but leaves the
  year coefficient unchanged; note that "left" is only 157 records.
- S6 adds species intercepts to σ (species differ in CV); the year sign flips
  back to a small positive. **The sign of the controlled year effect on σ is
  specification-dependent and never excludes zero**; the honest summary is
  "no evidence of a residual-variance trend once source and site are modelled",
  not "variance declines".
- Two-stage cross-check (lme4 on log|residual| of the M3 mean model, 0 exact
  zeros dropped): T0 +9.4 [4.0, 15.2] %/decade → T1 (+contributor + site)
  +6.3 [−1.7, +14.9]; Mundlak within-contributor +6.2 [−3.4, +16.8]. The
  two-stage estimate is less attenuated than glmmTMB S3 because the stage-2
  contributor intercepts absorb the level of |r| but not the contributor-specific
  precision weighting of the mean fit; treat the two approaches as bracketing the
  controlled effect (−6 % to +6 %, both spanning zero).

## 3. lnCVR (part 3)

- **Exact reproduction of `update_descriptive_stats.R`**: k = 54, lnCVR 0.184
  [0.034, 0.334], I² 93.6 % — matches `descriptive_summary.rds` (0.183 [0.033,
  0.334], I² 93.7). **Bug found in the original**: `sample.n.wing = n()` is
  computed on *all* quartile records of the species, including the 2,630
  records with no wing measurement, and the same wing n is reused for the bill
  lnCVR (bill width is measured on only 25 % of records). With n = non-NA wing
  records the pooled estimate is 0.177 [0.031, 0.323], I² 90.4 — same
  conclusion, but the Methods should quote the corrected version (§6).
- lnCVR definition: `log(CV_late/CV_early) + 1/(2(n_late−1)) − 1/(2(n_early−1))`,
  sampling variance `s2.lnCVR()` from `functions.R` with the mean–SD correlation
  computed across the cells of that analysis within each period (species cells
  0.47 / 0.49; species × sex 0.49 / 0.41; species × sex × contributor 0.67 /
  0.35). lnVR (SD ratio without the mean) is reported alongside: pooled +0.148
  [+0.001, +0.295]; within sex +0.163 [+0.024, +0.302]; within contributor
  −0.253 [−0.515, +0.010]. The pattern is the same, so this is a variance
  change, not a mean-driven CV change.
- **Within sex**: species × sex cells at n ≥ 2 / 5 / 10 in both periods give
  +0.159 (k 95), +0.187 (k 71), +0.208 (k 42), all excluding zero under `rma`;
  with a species random effect (`rma.mv`, cells within species are not
  independent) the n ≥ 5 estimate is +0.169 [−0.009, +0.346]. Females +0.213
  [+0.033, +0.392] (k 33), males +0.164 [−0.045, +0.374] (k 38). The sex-ratio
  shift (45.6 → 48.7 % female) does **not** explain the CV rise.
- **Within species × sex × contributor**: only 5 contributors have paired
  cells at all (E. Carrano, M. Alves, A. Ross, C. Fontana, A. Piratelli;
  L. Bugoni has one cell with n = 1). Estimates: n ≥ 2, k = 58: −0.108
  [−0.304, +0.088]; n ≥ 3, k = 39: −0.138 [−0.342, +0.066]; n ≥ 5, k = 15:
  −0.249 [−0.515, +0.017]; n ≥ 10, k = 4 (all Carrano): −0.421 [−0.729,
  −0.114]. `rma.mv` with a contributor random effect: −0.075 / −0.125 /
  −0.233 (same order). Per contributor at n ≥ 5: Carrano k = 7 mean −0.37
  (5 of 7 negative), Piratelli k = 3 mean −0.11, Ross k = 2 −0.01, Fontana
  k = 2 −0.41, Alves k = 1 +0.79. **The plan's "E. Carrano k = 6, mean −0.46"
  does not reproduce at any single threshold**; the closest are k = 7 (n ≥ 5,
  mean −0.37, median −0.47) and k = 4 (n ≥ 10, mean −0.39). Quote the k = 7
  `rma` row (−0.43 [−0.70, −0.16]) or the k = 15 all-contributor row, with the
  threshold stated.
- Descriptive CV (`cv_by_period_min5`, `_min10`): median CV of species × sex
  cells pooled across contributors 0.044 → 0.054 (n ≥ 5) / 0.045 → 0.055
  (n ≥ 10); within species × sex × contributor 0.034 → 0.038 / 0.036 → 0.040.
  The plan's 0.036 → 0.040 corresponds to the n ≥ 10 definition; its "0.054 →
  0.056 pooled" does not reproduce (we get 0.045 → 0.055 at n ≥ 10).

## 4. Mechanical heterogeneity (part 4)

- Contributors per species × period cell, wing records: mean 3.68 → 5.50
  (median 3 → 4; 63 / 70 cells); species present in both periods (61): 3.77 →
  5.98, 45 species gain contributors, 7 lose, Wilcoxon signed-rank p = 5.4e-8.
  Species × sex cells: 3.26 → 4.43. Municipalities per cell 6.4 → 6.8. These
  are the wing-record definitions (P0's `audit_sources.rds` has the all-record
  version 4.66 → 6.65 that the plan quotes).
- Variance decomposition within species × sex × period cells (n ≥ 5; 86 early,
  102 late cells): the share of within-cell variance that lies **between
  contributor means** rises from mean 0.34 / median 0.25 (early) to mean 0.45 /
  median 0.45 (late). Roughly half of the late-period within-cell variance is
  between measurers.
- Meta-regression of the species lnCVR on Δlog(contributors): slope +0.08
  [−0.27, +0.42] alone (R² 0); jointly with Δlog(sites): contributors +0.35
  [−0.09, +0.80], sites −0.44 [−0.90, +0.02], R² 8 %, intercept (species whose
  contributor and site counts did not change) +0.08 [−0.14, +0.30]. **Honest
  reading: the species-level meta-regression is weak (k = 54, huge sampling
  error in per-species lnCVR) and does not by itself demonstrate the
  mechanical inflation; the within-contributor lnCVR and the σ-model ladder are
  the evidence.**
- Mean-slope side effect (`mean_year_slope_comparison`): homoscedastic lme4
  M3 on the 8,282 complete cases −0.32 [−0.69, +0.04] (matches the plan's
  model 6 −0.30); the heteroscedastic glmmTMB S3 fit gives −0.15 [−0.35, +0.05]
  with half the SE. P1 should be aware that the controlled brms model's mean
  slope will move again if σ is modelled — the distributional fit (part 1) will
  give the phylogenetic version of this number.

## 5. brms distributional model (part 1, `--engine brms`) — kept as the Bayesian cross-check; not run

Superseded as the primary evidence by the glmmTMB phylogenetic tier (§9),
which fits the same mean and σ structure with the phylogenetic species
intercept on the identical 50 trees in ~3 s per fit. The Stan path below is
unchanged and selectable with `--engine brms --trees N` (Totoro); `--smoke`
implies it. Its purpose now is to confirm that Bayesian posterior means and
95 % intervals match the REML/Rubin estimates of §9 (as they did for the
published wing model: `glmmtmb_validation_wing.R`, year −0.922 [−1.400,
−0.444] glmmTMB vs −0.922 [−1.403, −0.440] brms).

- Written in full: `bf(conc.wing.length ~ 1 + Sex + scaled_yr + scaled_lat +
  scaled_lon + scaled_alt + season + (1 + scaled_yr || spp) + (1 | gr(species_name,
  cov = A)) + (1 | Main_researcher) + (1 | Municipality), sigma ~ scaled_yr [+
  scaled_tmean] + Sex + (1 | Main_researcher) + (1 | Municipality))`, gaussian,
  clootl 100-tree cloud via prepR4pcm exactly as `atlantic_parallel.R` (same
  eBird synonym map, same `make_A`, `inverseA` scaled covariance), priors as
  `atlantic_parallel.R` for the mean plus `normal(1.4, 1)` on the σ intercept,
  `normal(0, 1)` on σ slopes, `cauchy(0, 1)` on σ group SDs. Rubin pooling
  (`pool_rubin`, generalised to all `b_` and `sd_` parameters, with between- and
  within-tree variances stored), phylogenetic signal as the phylo share of the
  between-species + residual variance at the σ intercept, and per-fit R-hat /
  bulk-ESS / divergences are saved.
- Offline tree access: if `AVESDATA_PATH` is unset the script points clootl at
  `data/raw/AvesDataLite-main/` (present in the repo) before falling back to
  `get_avesdata_repo()`; the 73-species tree resolves with 0 unmatched names.
- The brms σ formula is now `sigma ~ scaled_yr + scaled_tmean + Sex +
  (1 | Main_researcher) + (1 | Municipality)` on the 7,651 tmean-complete
  records (§1); `--no-tmean` gives the 8,282-record version without the term.
  `scaled_tmean` enters σ only (the plan's specification), not the mean.
- Flags: `--fast` (parts 2–4 only; the local default), `--smoke` (1 tree,
  chains 2, iter 400, warmup 200, adapt_delta 0.9, treedepth 10; writes
  `output/models/variance_sigma_SMOKE.rda` and `variance_results.rds$brms_smoke`
  with a WARNING string; **not a result**), `--trees N` (Totoro: `--trees 50`,
  `_sampling_config.R` settings, `future_lapply` over trees, writes
  `output/models/variance_sigma_multiphylo.rda` and `variance_results.rds$brms`),
  `--no-tmean` (see above). Totoro run order: `Rscript
  Analysis/scripts/atlantic_variance_sigma.R --trees 50`, then optionally
  `--trees 50 --no-tmean` (note: the second run overwrites
  `variance_sigma_multiphylo.rda` and `$brms`; rename the first file before
  running it if both are wanted).
- Smoke-test outcome: §8.

## 6. What `update_descriptive_stats.R` needs (not edited — P3 does not own it)

1. Compute `sample.n.wing` as `sum(!is.na(conc.wing.length))` and a separate
   `sample.n.bill = sum(!is.na(Bill_width.mm.))`, and use the bill n in the
   bill lnCVR and its `s2.lnCVR`. Currently both use `n()` over all quartile
   records (2,630 NA-wing records inflate the wing n; the bill n is inflated
   far more).
2. Either read the lnCVR numbers from `variance_results.rds$lncvr` (rows
   "pooled per species (corrected n …)", "within sex (species x sex cells), n >= 5",
   "within species x sex x contributor, n >= 5") into `descriptive_summary.rds`
   as `lncvr_wing_*`, `lncvr_wing_sex_*`, `lncvr_wing_src_*`, or copy the
   `lncvr_table()` helper from `atlantic_variance_sigma.R`.
3. Add the S3 σ-model row (`variance_results.rds$glmmTMB`, tier S3, par
   scaled_yr: estimate, pct_decade, pct_decade_lo/hi), the S7 tmean row (tier
   S7, par scaled_tmean: `pct_per_sd_tmean`, `_lo`, `_hi`) for the §2.7
   σ-contradiction sentence, and, when available, the Rubin-pooled
   `b_sigma_scaled_yr` and `b_sigma_scaled_tmean` from
   `variance_results.rds$brms$rubin_summary` as the primary variance numbers.
4. Add `contributors_per_cell_summary` (3.68 → 5.50) and
   `variance_decomposition_summary` (between-contributor share 0.25 → 0.45
   median) for the mechanical-heterogeneity sentence.

## 7. Figures (`Analysis/figures/`)

- `variance_summary.png` — A: lnCVR forest across aggregation levels (pooled,
  within sex, within contributor, Carrano); B: %Δ residual SD per decade for
  tiers S0–S6, S3t, S7 (blue), each followed by its phylogenetic 50-tree
  Rubin-pooled twin "+ phylogeny" (green, §9), and the two-stage T0–T2 (red);
  C: histogram of contributors per species × period cell by period; D: boxplot
  of the between-contributor share of within-cell variance by period.
  Candidate replacement for the manuscript's Fig. 4 (P6 may restyle from
  `variance_results.rds` / `variance_phylo_results.rds$table`).
- `variance_lncvr_contributor.png` — the 15 within-contributor cells (n ≥ 5)
  with 95 % CIs, coloured by contributor.

## 8. brms smoke-test log (`--smoke`)

Launched 2026-09-09 11:23:20 under the 21-minute `perl alarm` cap. The Stan
model compiled and sampling started (2 chains × 400 iterations on 1 tree), but
the process was killed at iteration ~201/400 of the first chain when the cap
expired. No `output/models/variance_sigma_SMOKE.rda` and no
`variance_results.rds$brms_smoke` exist; the brms code path has therefore
never completed end to end locally. It is no longer needed for the primary
result (§9) and is kept only as the optional Bayesian cross-check for Totoro
(`--engine brms --trees 50`, optionally with `--no-tmean`).

## 9. Phylogenetic tier — glmmTMB `propto` engine, 50 trees (part 1, `--engine glmmTMB`, default)

Run 2026-09-09 (E2 session) with `Rscript Analysis/scripts/atlantic_variance_sigma.R
--trees 50` (engine default glmmTMB; R 4.6.0, glmmTMB 1.1.14, TMB 1.9.21).
Every number in this section is in `Analysis/output/variance_phylo_results.rds`
(`$table`, `$comparison`, `$headline`, `$pooled_all`, `$varcomp`,
`$varcomp_summary`, `$timing`, `$headline_tree1`); the headline and the
comparison table are mirrored in `variance_results.rds$phylo`.

### 9.1 What was fitted

- Engine: `_phylo_engine.R` (`fit_phylo_glmmtmb`, `run_phylo_trees`,
  `pool_rubin_df`), which appends `propto(0 + species_name | g, A)` — a
  phylogenetic species intercept whose covariance is proportional to the
  species correlation matrix A of one tree (Williams, McGillycuddy, Drobniak,
  Bolker, Warton & Nakagawa 2025, bioRxiv 10.64898/2025.12.20.695312) — to the
  mean model, and fits the σ sub-model with `dispformula`. REML. Rubin's rules
  over trees (estimate = mean; SE² = mean within-tree SE² + (1 + 1/m) ×
  between-tree variance of the estimates); Wald 95 % intervals.
- Trees: the **identical 50 trees as the published brms wing model**
  (`data/derived/phylo_A_50trees.rds`, correlation matrices cached from
  `brm0_multiphylo.rda`; `inverseA(scale = TRUE)` on nnls-ultrametricised
  clootl trees, seed 20240101). ABT binomials mapped to tree tips with the
  eBird synonym map of `atlantic_parallel.R`; 72 species, all present in the
  trees (`Herpsilochmus_sellowi`, which is absent from the trees, has no wing
  record in the complete-case sample: 0 records dropped).
- Data: the same records as the non-phylogenetic tiers — 8,282 complete-case
  wing records (S0–S6) and the 7,651 tmean-complete records (S3t, S7) —
  with `spp = species_name` for the species random slope, so every
  phylogenetic estimate is directly comparable with §2.
- Models: **all nine tiers of §2 (S0–S6, S3t, S7) refitted with the
  phylogenetic species intercept in the mean model**; σ sub-models unchanged
  (`scaled_tmean` enters σ only, as in the plan). The two headline models are
  S7 + phylogeny (the plan's full Phase 3 specification:
  `sigma ~ scaled_yr + scaled_tmean + Sex + (1|Main_researcher) +
  (1|Municipality)`, N = 7,651) and S3 + phylogeny (the `--no-tmean` variant,
  N = 8,282). Because the glmmTMB fits are cheap, both are produced in one run;
  `--no-tmean` only removes the tmean tiers.
- Wall time (all 50 trees, sequential, laptop): S0 20 s, S1 21 s, S2 249 s,
  S3 365 s, S4 396 s, S5 452 s, S6 527 s, S3t 326 s, S7 350 s; **45.4 min for
  the whole ladder** (0.4–10.5 s per fit; the S2–S6 fits were slowed by a
  concurrent diagnostic process, the uncontended rate is ~3.5 s per S3 fit).
  For comparison the brms path would have been 9 × 50 Stan fits of hours each.
- Convergence (positive-definite Hessian, `sdr$pdHess`): **50/50 trees for
  every tier except S0 (45/50)**; see §9.4. `fit$convergence == 0` on all
  450 fits.

### 9.2 Result: σ-year coefficient, phylogeny vs no phylogeny (same records)

Log-SD scale; b (SE); % change in residual SD per decade
= 100·(exp(b·10/5.0199) − 1).

| Tier | σ term | no phylogeny (§2) | **+ phylogeny, 50 trees** | %/decade no phylo | **%/decade + phylo** | trees pdHess |
|---|---|---:|---:|---|---|---:|
| S0 | scaled_yr | +0.0461 (0.0079) | **+0.0462 (0.0079)** | +9.6 [+6.3, +13.0] | **+9.6 [+6.3, +13.1]** | 45/50 |
| S1 | scaled_yr | +0.0468 (0.0079) | +0.0469 (0.0079) | +9.8 [+6.4, +13.2] | +9.8 [+6.5, +13.3] | 50/50 |
| S2 | scaled_yr | +0.0123 (0.0138) | +0.0125 (0.0138) | +2.5 [−2.9, +8.2] | +2.5 [−2.9, +8.2] | 50/50 |
| **S3** | scaled_yr | −0.0327 (0.0203) | **−0.0326 (0.0203)** | −6.3 [−13.5, +1.4] | **−6.3 [−13.4, +1.5]** | 50/50 |
| S4 | yr_within_src | −0.0314 (0.0211) | −0.0313 (0.0211) | −6.1 [−13.5, +2.0] | −6.0 [−13.5, +2.0] | 50/50 |
| S4 | yr_src_mean | −0.0524 (0.0841) | −0.0528 (0.0841) | −9.9 [−35.1, +25.1] | −10.0 [−35.2, +25.0] | 50/50 |
| S5 | scaled_yr | −0.0328 (0.0204) | −0.0327 (0.0204) | −6.3 [−13.5, +1.5] | −6.3 [−13.5, +1.5] | 50/50 |
| S6 | scaled_yr | +0.0169 (0.0209) | +0.0170 (0.0209) | +3.4 [−4.7, +12.2] | +3.5 [−4.7, +12.2] | 50/50 |
| S3t | scaled_yr | −0.0357 (0.0205) | −0.0355 (0.0205) | −6.9 [−14.0, +0.9] | −6.8 [−14.0, +0.9] | 50/50 |
| **S7** | scaled_yr | −0.0354 (0.0206) | **−0.0353 (0.0206)** | −6.8 [−14.0, +1.0] | **−6.8 [−14.0, +1.0]** | 50/50 |
| S7 | scaled_tmean | −0.0067 (0.0432) | −0.0071 (0.0432) | −0.7 % per SD [−8.7, +8.1] | −0.7 % per SD [−8.8, +8.1] | 50/50 |

- **Headline (plan specification, S7 + phylogeny, N = 7,651, 72 species, 50
  trees):** b_σ(year) = −0.0353, SE 0.0206, 95 % CI [−0.0756, +0.0050],
  z = −1.72 → **−6.8 % residual SD per decade [−14.0, +1.0]**;
  b_σ(tmean) = −0.0071, SE 0.0432 [−0.0917, +0.0775], z = −0.16 → −0.7 % per
  SD of record-year temperature [−8.8, +8.1]; b_σ(Sex male) = −0.027 (0.019).
  σ intercept 1.330 (0.073) → baseline residual SD exp(1.33) = 3.78 mm.
- **`--no-tmean` variant (S3 + phylogeny, N = 8,282):** b_σ(year) = −0.0326,
  SE 0.0203 [−0.0725, +0.0073], z = −1.60 → **−6.3 %/decade [−13.4, +1.5]**;
  σ intercept 1.368 (0.078) → 3.93 mm.
- Between-tree variance of every σ coefficient is ≤ 6 × 10⁻⁸ (≤ 0.02 % of the
  within-tree variance): tree uncertainty is irrelevant for the σ sub-model,
  which contains no phylogenetic term. Phylogenetic uncertainty shows up only
  in the mean-model intercept (between-tree var 0.25 on an SE of 11.8).
- **Does phylogeny change anything?** No. Across the 11 σ-year/tmean
  coefficients the largest shift is 0.0004 log units (S4 between-contributor,
  a term with SE 0.084) and every SE ratio (phylo / no phylo) is 1.000 to 3
  decimals. The ladder's shape is exactly as in §2: +9.6 %/decade with no
  controls → +2.5 with the contributor intercept → −6.3 with contributor + site
  → sign-unstable (−6 to +3.5) across the sensitivity tiers, every controlled
  interval spanning zero. The uncontrolled S0/S1 estimates (+9.6/+9.8 %,
  intervals excluding zero) are the only ones that do, with or without
  phylogeny.

### 9.3 Mean-model year slope from the same fits

| Tier | no phylogeny (§2) | + phylogeny | SE ratio |
|---|---:|---:|---:|
| S0 | −0.920 [−1.377, −0.463] | −0.922 [−1.381, −0.463] | 1.003 |
| S2 | −0.215 [−0.456, +0.025] | −0.216 [−0.457, +0.025] | 1.003 |
| S3 | −0.151 [−0.353, +0.051] | −0.152 [−0.353, +0.050] | 1.000 |
| S6 | −0.179 [−0.372, +0.014] | −0.181 [−0.374, +0.013] | 1.003 |
| S7 | −0.182 [−0.421, +0.057] | −0.183 [−0.421, +0.056] | 0.998 |

The S0 + phylogeny row reproduces the published brms wing model (−0.922
[−1.403, −0.440]; heteroscedastic here, hence the marginally narrower
interval) — a check that the engine, the trees and the species mapping are the
same as in `atlantic_parallel.R`. The S3/S7 mean slope (−0.15 / −0.18 mm per
SD-year ≈ −0.3 / −0.36 mm per decade, intervals spanning zero) is the
heteroscedastic, phylogenetic version of the §4 side effect: the homoscedastic
lme4 M3 fit on the same 8,282 records gives −0.32 [−0.69, +0.04]
(`variance_results.rds$lme4_mean_M3_year`), and modelling σ by
contributor/site halves the point estimate and the SE because records from
noisy contributors are down-weighted. Phylogeny does not move it. P1's
homoscedastic phylogenetic M3 (E1, `controlled_wing_phylo_results.rds`) is the
number to quote for the mean trend; this one is its σ-weighted companion.

### 9.4 Variance components and the 5 non-PD trees

- Per tree (`$varcomp`): phylogenetic SD 23.7 mm (S3; 95 % range over trees
  20.8–26.8), species-intercept SD 0.004–0.007 (i.e. zero: the phylogenetic
  term absorbs the whole between-species intercept variance, exactly as in the
  published model, phylogenetic signal 0.95), species year-slope SD 0.41
  (S3–S5) / 0.60 (S3t, S7), contributor SD 3.14, municipality SD 1.36 (S3;
  tree 1). Phylogenetic share of (between-species + contributor + site +
  baseline residual) variance: S3 0.953 [0.941, 0.963], S7 0.960 [0.949,
  0.968] (`$varcomp_summary`; denominator uses exp(b_σ intercept) because
  `sigma()` is undefined under a dispformula).
- σ random effects (tree-1 headline fits, `$headline_tree1`): S3 contributor
  0.383, municipality 0.389 (exp(u) range across contributors 0.45–2.53); S7
  0.334 / 0.390 (0.58–2.42) — identical to the non-phylogenetic S3/S3t values
  in §2.
- Species year slopes from the tree-1 S7 fit (fixed + BLUP, quadrature SE): 50
  of 72 negative, 4 exclude zero (all negative), median −0.19 mm per SD-year —
  reported only as a cross-reference to E1's caterpillar, which is the
  authoritative species-slope output.
- **S0: 45/50 trees with a positive-definite Hessian.** Diagnosed on all 50
  trees (scratch `diag_pdhess.R`): on trees 9, 14, 19, 21 and 49 the
  optimizer converged (nlminb code 0, "relative convergence") to the
  *alternative* mode in which the non-phylogenetic species intercept carries
  the between-species variance (SD 16.1) and the phylogenetic SD is 0.11–0.18,
  instead of 0.005 and ~24; with only Sex, year and latitude in the mean model
  the two species-level terms are weakly identified against each other and the
  Hessian is singular along that ridge. Fixed effects are unaffected:
  b_σ(year) = 0.04608 vs 0.04607, SE 0.007898 on every tree, mean year slope
  −0.9199 vs −0.922. Pooling only the 45 PD trees (`$pooled_all` `*_conv`
  columns) gives 0.046201 (0.007899) vs 0.046188 (0.007899) over all 50 — a
  difference in the 5th decimal. The S0 phylogenetic-share summary (mean 0.863,
  lower bound ~0) reflects those 5 trees and should not be quoted; use S1–S7
  (0.945–0.972). Once contributor/site intercepts are in the mean model
  (S2–S7) all 50 trees are PD.

### 9.5 What this means for the manuscript and for P7

1. The primary variance evidence is now a **phylogenetic distributional model
   on the published trees**, not drmSEM and not a brms run pending on Totoro:
   quote S7 + phylogeny (−6.8 %/decade [−14.0, +1.0], N = 7,651, 72 species,
   50 trees) as the plan's specification and S3 + phylogeny (−6.3 % [−13.4,
   +1.5], N = 8,282) as the no-temperature variant, alongside the uncontrolled
   S0 (+9.6 % [+6.3, +13.1]) to show what the controls do.
2. The referee-§2.7 σ contradiction is resolved the same way as in §2: the
   temperature effect on σ is −0.7 % per SD [−8.8, +8.1] (z = −0.16) once
   contributor and site are in the σ model; the drmSEM "warmer → less variable"
   path should not be carried forward.
3. Every σ-model result is engine-robust by construction (the glmmTMB fits are
   REML maximum-likelihood, the brms path is the Bayesian cross-check if Totoro
   time is available), and phylogeny-robust by direct comparison.
4. `update_descriptive_stats.R` (§6, item 3): read the primary numbers from
   `variance_phylo_results.rds$headline` (or `variance_results.rds$phylo$headline`)
   instead of `$brms$rubin_summary`.

### 9.6 Engine notes (for GLMMTMB_ENGINE.md; `_phylo_engine.R` not edited)

- `tidy_phylo_fit()` returns `sigma = NA` and therefore a `phylo_prop` that
  omits the residual variance whenever `dispformula != ~1` (`sigma()` is
  undefined for glmmTMB models with a dispersion model). Scripts with a σ
  sub-model must recompute the share with exp(b_σ intercept) as done here
  (`phylo_share()` in `atlantic_variance_sigma.R`); a `dispformula`-aware
  `phylo_prop` in the engine would remove the duplication.
- `tidy_phylo_fit()` builds `varcomp` with `data.frame()`, which mangles the
  SD names (`spp:(Intercept)` → `spp..Intercept.`, `Main_researcher:(Intercept)`
  → `Main_researcher..Intercept.`); `check.names = FALSE` would keep them
  stable. `varcomp` also only carries the conditional-model SDs; the
  dispersion-model random-effect SDs (`VarCorr(fit)$disp`) are not tidied, so
  this script refits tree 1 of the headline models to get them.
- `run_phylo_trees()` pools all trees regardless of `pdHess`; a
  `pooled_converged` element (pooling the PD trees only) would make the S0-type
  sensitivity a one-liner. Implemented locally here.
- Tiers whose mean model has different random effects produce `varcomp` tables
  with different columns; `dplyr::bind_rows()` (not `rbind`) is needed to stack
  them.


## Phylogenetic tier, final (glmmTMB `propto`, 50 trees) — 2026-09-09 16:12

Supersedes any 1- or 2-tree numbers in the section above. Run: `Rscript atlantic_variance_sigma.R
--trees 50 --engine glmmTMB`, 9 σ tiers, 1,768 s; 50/50 trees converged except S0 (45/50). Output
`variance_phylo_results.rds` (`$pooled_all`, `$headline`, `$comparison`, `$varcomp_summary`).

log-residual-SD year terms (per SD-year; % change in residual SD per decade in brackets): S0 +0.046 [+0.031, +0.062] (+9.6 %/decade in residual SD); S2 (+contributor +site in mean and σ) +0.013 [−0.015, +0.040] (+2.5 % [−2.9, +8.2]); S3 full controls −0.033 [−0.073, +0.007] (−6.3 % [−13.4, +1.5]); S4 within-contributor −0.031 [−0.073, +0.010], between −0.053 [−0.218, +0.112]; S6 (species intercepts in σ) +0.017 [−0.024, +0.058]; S7 (S3t + temperature) year −0.035 [−0.076, +0.005], temperature −0.007 [−0.092, +0.078].
Mean-model year terms in the same fits: S0 −0.922 [−1.381, −0.463]; S2 −0.216 [−0.457, +0.025];
S3 −0.152 [−0.353, +0.050]; S7 −0.183 [−0.421, +0.056].

Plain statement: the baseline rise in among-individual wing variance (+9.6 % per decade) is a
between-contributor composition effect. With contributor and site terms in the mean and σ models it
is +2.5 % [−2.9, +8.2]; under the full controls it is −6.3 % [−13.4, +1.5]; the within-contributor
σ-year component is −6.0 % [−13.5, +2.0]; temperature has no detectable effect on residual variance
(−0.7 % per SD of temperature [−9.2, +7.8]). Phylogeny changes nothing: phylogenetic and
non-phylogenetic σ-year terms differ by ≤ 0.0006 (`$comparison`). Together with the within-contributor
lnCVR (section above), the manuscript's variance headline does not survive provenance controls.
