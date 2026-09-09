# Revision notes — Phase 2 (multi-trait screen and wing–mass isometry), 2026-09-09

Companion to `REVISION_PLAN.md` §3 "Phase 2" and §1 "§2.5 Mass and allometry" /
"Bill width and other traits". Two new scripts, what was run, what was not, and
where every number lives. Style follows `LIVE_ONLY_MIGRATION.md` and
`REVISION_NOTES_P0.md`.

Executed locally with R 4.6.0: lme4 2.0.1, glmmTMB 1.1.14, brms 2.23.0 (rstan
2.32.7 backend; cmdstanr not installed here), prepR4pcm 0.5.0.9000, clootl
0.1.4, ape 5.8.1, MCMCglmm 2.36, phytools 2.5.2, posterior 1.7.0, dplyr 1.1.x.
**Everything in §2–§4 is a full run (lme4 / glmmTMB, seconds to a minute).
§5 is a brms SMOKE TEST only — its numbers are not results.** The 50-tree brms
run is for Totoro. **§8 (added 2026-09-09, afternoon) is the PHYLOGENETIC TIER: every
model of §2–§3 refitted with glmmTMB `propto` across the 50 published trees and
Rubin-pooled — a full run, quotable.**

Data: `data/derived/passer90.rda` as rebuilt by P0 (12,571 live known-sex
records, 73 species, 42 columns). `scaled_yr` SD = 5.0199 yr (attr
`scaled:scale`), so per-decade = β × 1.992. Early / late periods for
"spanning contributors" are ≤ 2006 / ≥ 2013 as in P0.

---

## 1. Scripts

| Script | What | Run time |
|---|---|---|
| `atlantic_multitrait.R` | lme4 (REML) screen of log body mass, bill width, bill length, tail, tarsus (+ wing as a cross-reference to Phase 1): M0 → M1a/M1b/M1 → M1cc → M2 → M3, Mundlak within/between contributor and municipality, mass capture-hour model. Writes `output/multitrait_results.rds` (`$table`, `$coverage`, `$hour`, `$isometry`, `$plan_comparison`, `$meta`) and `output/multitrait_table.md`. | 21 s |
| `atlantic_bivariate_wing_mass.R --fast-lme4` | Shared wing + mass subset (7,577 records): lme4 wing, log-mass and wing \| log-mass models (M0, M1, M1cc, M2), slope contrasts under independence, and a stacked glmmTMB joint model with a covariance-correct slope difference. Writes `output/bivariate_fast_lme4.rds`. | 22 s |
| `atlantic_bivariate_wing_mass.R --smoke` | brms bivariate `mvbind(wing, lnmass)` + `set_rescor(TRUE)` and the relative-wing model on 1 tree, chains = 2, iter = 400. Writes `output/bivariate_smoke_results.rds`, `output/models/bivariate_smoke.rda`. **Not a result.** | see §5 |
| `atlantic_multitrait.R --trees 50` | **Phylogenetic tier (glmmTMB `propto`, default engine)**: log mass, bill width, bill length, tail, tarsus at M0 / M1 / M3 + within-contributor Mundlak, and the mass capture-hour pair, across the 50 published trees via `_phylo_engine.R`. Writes `output/multitrait_phylo_results.rds` and appends the phylogenetic table to `multitrait_table.md`. §8. | 41.7 min |
| `atlantic_bivariate_wing_mass.R --trees 50` (= `--engine glmmTMB`, the default) | **Phylogenetic tier** on the shared subset: wing, log-mass and wing \| log-mass univariate models (M0, M1, M1cc, M2), delta-method slope contrasts (independence assumed), and a joint long-format model with two `propto` terms (covariance-correct contrast). Writes `output/bivariate_phylo_results.rds`. §8. | 37.6 min |
| `atlantic_bivariate_wing_mass.R --engine brms [--trees 50]` | Full brms run (bivariate `mvbind` + `rescor`, relative-wing model), Rubin-pooled across trees (default 50). Writes `output/bivariate_results.rds`, `output/models/bivariate_wing_mass_multiphylo.rda`. **Not run**; the Bayesian cross-check. | Totoro |

Model sequence (both scripts; `src` = `Main_researcher`, `site` = `Municipality`,
`ring` = `Ring × Binomial` with a unique level per unringed record):

```
M0   y ~ Sex + scaled_yr + scaled_lat + (1 + scaled_yr || spp)           # brms structure minus phylogeny
M1a  M0 + (1 | src)          M1b  M0 + (1 | site)          M1  M0 + (1 | src) + (1 | site)
M1cc M1 on the M2/M3 complete-case sample
M2   M1 + scaled_lon + scaled_alt + season
M3   M2 + (1 | ring)
MWsrc  M3 with scaled_yr -> yr_within_src + yr_mean_src        (Mundlak, contributor)
MWsite M3 with scaled_yr -> yr_within_site + yr_mean_site      (Mundlak, municipality)
```

---

## 2. Multi-trait table (generated: `output/multitrait_table.md`)

| Trait | n (% of 12,571) | contributors (spanning) | M0 baseline β_yr | M1 +src +site | M1 on M3 sample | M2 +lon +alt +season | M3 +ring | Within-contributor (Mundlak) | Between-contributor | Within-municipality |
|---|---:|---:|---|---|---|---|---|---|---|---|
| log body mass (ln g) | 11,256 (89.5 %) | 47 (6) | -0.0072 [-0.0165, 0.0021] t = -1.5 (n = 11,256) | -0.0056 [-0.0164, 0.0051] t = -1.0 (n = 11,256) | -0.0058 [-0.0157, 0.0041] t = -1.2 (n = 11,077) | -0.0051 [-0.0150, 0.0048] t = -1.0 (n = 11,077) | -0.0045 [-0.0138, 0.0047] t = -1.0 (n = 11,077) | -0.0031 [-0.0102, 0.0039] t = -0.9 (n = 11,077) | 0.0060 [-0.0303, 0.0424] t = 0.3 (n = 11,077) | -0.0034 [-0.0105, 0.0037] t = -0.9 (n = 11,077) (singular) |
| bill width (mm) | 3,209 (25.5 %) | 22 (2) | 0.24 [0.06, 0.42] t = 2.6 (n = 3,209) | -0.02 [-0.17, 0.13] t = -0.3 (n = 3,209) | -0.03 [-0.18, 0.13] t = -0.3 (n = 3,205) | -0.02 [-0.17, 0.14] t = -0.2 (n = 3,205) | -0.02 [-0.17, 0.14] t = -0.2 (n = 3,205) | 0.03 [-0.18, 0.23] t = 0.3 (n = 3,205) | 0.40 [-0.48, 1.29] t = 0.9 (n = 3,205) | -0.01 [-0.17, 0.15] t = -0.2 (n = 3,205) (singular) |
| bill length (mm) | 7,697 (61.2 %) | 42 (4) | -0.00 [-0.11, 0.10] t = -0.0 (n = 7,697) | 0.16 [-0.02, 0.34] t = 1.7 (n = 7,697) | 0.16 [-0.03, 0.35] t = 1.7 (n = 7,540) | 0.16 [-0.02, 0.35] t = 1.7 (n = 7,540) | 0.17 [-0.02, 0.36] t = 1.8 (n = 7,540) | 0.14 [-0.13, 0.41] t = 1.0 (n = 7,540) | 0.07 [-0.38, 0.52] t = 0.3 (n = 7,540) | 0.33 [0.11, 0.56] t = 2.9 (n = 7,540) (singular) |
| tail length (mm) | 8,872 (70.6 %) | 42 (6) | 0.28 [-0.27, 0.83] t = 1.0 (n = 8,872) | 0.30 [-0.19, 0.80] t = 1.2 (n = 8,872) | 0.22 [-0.32, 0.75] t = 0.8 (n = 8,679) | 0.17 [-0.36, 0.70] t = 0.6 (n = 8,679) | 0.25 [-0.28, 0.77] t = 0.9 (n = 8,679) | 0.24 [-0.21, 0.70] t = 1.0 (n = 8,679) | -0.19 [-1.33, 0.95] t = -0.3 (n = 8,679) | 0.24 [-0.19, 0.67] t = 1.1 (n = 8,679) |
| tarsus length (mm) | 4,361 (34.7 %) | 27 (3) | 0.42 [0.26, 0.59] t = 5.0 (n = 4,361) | 0.26 [-0.06, 0.58] t = 1.6 (n = 4,361) (singular) | 0.30 [-0.03, 0.63] t = 1.8 (n = 4,285) (singular) | 0.24 [-0.11, 0.58] t = 1.4 (n = 4,285) (singular) | 0.26 [-0.08, 0.59] t = 1.5 (n = 4,285) (singular) | 0.26 [-0.16, 0.68] t = 1.2 (n = 4,285) | 0.03 [-0.67, 0.72] t = 0.1 (n = 4,285) | 0.19 [-0.29, 0.68] t = 0.8 (n = 4,285) |
| wing length, coalesced (mm) [cross-reference; Phase 1 is authoritative] | 8,478 (67.4 %) | 42 (6) | -0.92 [-1.40, -0.44] t = -3.8 (n = 8,478) | -0.53 [-0.92, -0.13] t = -2.6 (n = 8,478) | -0.34 [-0.71, 0.02] t = -1.8 (n = 8,282) | -0.32 [-0.69, 0.04] t = -1.7 (n = 8,282) | -0.30 [-0.66, 0.06] t = -1.6 (n = 8,282) | 0.02 [-0.58, 0.62] t = 0.1 (n = 8,282) | -1.24 [-2.35, -0.12] t = -2.2 (n = 8,282) | -0.35 [-0.62, -0.08] t = -2.6 (n = 8,282) |

β_yr is per SD-year (SD = 5.0199 yr); multiply by 1.992 for per-decade units. Wald 95 % CI in brackets; lme4 REML, no phylogenetic term. "(singular)" = `isSingular()` TRUE (a variance component, usually the species year-slope or ring variance, estimated at zero); the fixed effects are still reported.

Intermediate single-control fits (in `$table` as `M1a_src`, `M1b_site`): bill width
+src only −0.24 [−0.39, −0.10] t = −3.4 (the plan's "flips sign"), +site only −0.07
(t = −0.9); wing +src −0.68 (t = −3.3), +site −0.78 (t = −3.8); tarsus +src +0.37
(t = 2.9), +site +0.13 (t = 0.8); bill length +src −0.07 (t = −0.9), +site +0.24
(t = 3.1).

### Body mass, capture hour (`multitrait_results.rds$hour`)

Records with body mass and a parseable clock time between 05:00 and 19:00
(complete for site and season): **n = 8,345** (8,391 mass records have any parsed
hour). `ln_body_mass ~ Sex + scaled_yr + scaled_lat + hour_num + season +
(1 + scaled_yr || spp) + (1 | src) + (1 | site)`:

- hour: **+0.404 % body mass per hour [0.291, 0.516], t = 7.02**
- year slope with hour in the model: −0.0063 [−0.0188, +0.0062], t = −0.99;
  without hour, same sample: −0.0075, t = −1.17.

### Contributor coverage per trait (`multitrait_results.rds$coverage`)

| Trait | n | contributors | spanning both periods (records) | who |
|---|---:|---:|---:|---|
| body mass | 11,256 | 47 | 6 (3,171) | A.Piratelli, A.Ross, C.Fontana, E.Carrano, L.Bugoni, M.Alves |
| **bill width** | 3,209 | 22 | **2 (620 of 3,209 = 19 %)** | **A.Piratelli, M.Alves** |
| bill length | 7,697 | 42 | 4 (1,837) | A.Piratelli, E.Carrano, L.Bugoni, M.Alves |
| tail | 8,872 | 42 | 6 (2,692) | as body mass |
| tarsus | 4,361 | 27 | 3 (530) | A.Piratelli, C.Fontana, L.Bugoni |
| wing | 8,478 | 42 | 6 (2,828) | as body mass |

Bill width: 945 early / 1,237 late records; the between-contributor Mundlak
component is +0.40 [−0.48, +1.29] and the within-contributor slope +0.03 [−0.18,
+0.23] — the baseline +0.24 is a contributor-composition effect with no
within-source support. **The Allen's-rule claim is not supported**; this
matches the plan.

### Isometry from the univariate lme4 fits (`multitrait_results.rds$isometry$lme4`; wing and mass from DIFFERENT samples)

Over the 23-year record (1995 → 2018), mass ∝ wing³ expectation vs observed:

| Model | wing change % [CI] | isometric mass expectation % [CI] | observed mass change % [CI] |
|---|---|---|---|
| M0 baseline | −5.9 [−9.0, −2.9] | −16.8 [−24.7, −8.3] | −3.3 [−7.3, +1.0] |
| M1 +src +site | −3.4 [−5.9, −0.9] | −9.9 [−16.8, −2.6] | −2.6 [−7.2, +2.4] |
| M3 (+lon +alt +season +ring) | −1.9 [−4.3, +0.4] | −5.7 [−12.3, +1.1] | −2.1 [−6.1, +2.2] |

The published (50-tree brms) baseline conversion in P0's `effect_scale.rds` is
−5.93 % wing → −16.7 % isometric expectation vs −3.3 % observed [−7.4, +1.0];
the lme4 baseline reproduces it. At M1 the intervals still separate (expected
−9.9 % vs observed upper +2.4 %, and the observed interval excludes −9.9); **at
M3 they overlap** (expected −5.7 [−12.3, +1.1] vs observed −2.1 [−6.1, +2.2]), so
on the fully controlled specification the departure from isometry is no longer
distinguishable from zero when wing and mass come from different records. The
shared-record test is §3.

---

## 3. Shared-record isometry test, fast mode (`output/bivariate_fast_lme4.rds`)

Shared subset: **7,577 records, 72 species** (one of the 73 has no record with
both traits), 42 contributors, 128 municipalities; mean wing 71.48 mm, mean mass
23.0 g. Complete-case (lon, alt, season) subset: 7,409.

### 3a. Univariate lme4 analogues (`$univariate`)

| Model | wing β_yr (mm/SD-yr) | log-mass β_yr | wing \| log-mass β_yr | allometric b(ln mass) |
|---|---|---|---|---|
| M0 baseline (n 7,577) | −0.83 [−1.33, −0.33] t = −3.3 | −0.0021 [−0.0102, +0.0061] t = −0.5 | −0.81 [−1.30, −0.33] t = −3.3 | 3.20 (se 0.33) |
| M1 +src +site (n 7,577) | **−0.55 [−0.97, −0.13] t = −2.6** | **−0.0027 [−0.0135, +0.0081] t = −0.5** | **−0.54 [−0.94, −0.15] t = −2.7** | 3.24 (se 0.32) |
| M1cc (M1 on n 7,409) | −0.35 [−0.75, +0.06] t = −1.7 | −0.0034 [−0.0118, +0.0049] t = −0.8 | −0.34 [−0.74, +0.05] t = −1.7 | 2.85 (se 0.33) |
| M2 +lon +alt +season (n 7,409) | −0.34 [−0.75, +0.06] t = −1.7 | −0.0030 [−0.0114, +0.0055] t = −0.7 | −0.34 [−0.73, +0.06] t = −1.7 | 2.85 (se 0.33) |

The M1 row reproduces REVISION_PLAN.md §0 exactly (wing −0.55 [−0.97, −0.14];
log-mass −0.0027 [−0.0135, +0.0081]; wing with log-mass −0.55 [−0.95, −0.15]).
With species-centred log mass as the covariate the year slopes are unchanged
(rows "wing | lnmass_c" in `$univariate`).

### 3b. Difference of year slopes, % per decade

`$contrasts_independent` (from the two univariate fits, treating the estimates
as independent) and `$stacked_glmmTMB` (joint fit, covariance-correct):

| Specification | wing %/decade | mass %/decade | **difference wing − mass** | isometry contrast mass − 3·wing (0 under isometry) | slope cov / cor |
|---|---|---|---|---|---|
| M0, independent | −2.31 [−3.69, −0.93] | −0.41 [−2.03, +1.21] | −1.90 [−4.03, +0.23] | (gap +6.4) | — |
| M1, independent | −1.53 [−2.69, −0.37] | −0.53 [−2.67, +1.60] | −1.00 [−3.43, +1.43] | (gap +4.0) | — |
| **M1, joint glmmTMB (n 7,577)** | **−1.49 [−2.64, −0.34]** | **−0.57 [−2.66, +1.53]** | **−0.92 [−3.27, +1.43], z = −0.77** | **+3.90 [−0.07, +7.88], z = 1.93** | 0.0061 / 0.039 |
| M2, independent (n 7,409) | −0.95 [−2.08, +0.18] | −0.59 [−2.25, +1.08] | −0.37 [−2.38, +1.65] | (gap +2.2) | — |
| **M2, joint glmmTMB (n 7,409)** | −0.88 [−1.93, +0.16] | −0.64 [−2.30, +1.02] | **−0.24 [−2.16, +1.68], z = −0.25** | **+2.01 [−1.47, +5.48], z = 1.13** | — |

Joint model: rows = record × trait, responses on the 100·ln scale (so slopes are
% per SD-year), trait-specific fixed effects, `diag(0 + trait + trait:scaled_yr |
spp)`, `diag(0 + trait | src)`, `diag(0 + trait | site)`, a shared `(1 | rec)`
intercept and `dispformula = ~ trait`. Converged with a positive-definite Hessian
in 6 s. Implied within-record correlation between wing and log mass ≈ 0.50 (M1)
— but note this parameterisation constrains the record-level covariance to a
single shared component; the brms `rescor` is the proper estimate. The
covariance between the two year slopes is tiny (correlation 0.04), so the
independence approximation in the first rows is adequate here.

**Reading, stated plainly.** On the same individuals with contributor and site
intercepts, the wing slope excludes zero and the mass slope does not — that
is what the manuscript currently says. But the *difference* between the two
year slopes, which is the quantity the referee asked for, does **not** exclude
zero in any specification (M1: −0.92 %/decade [−3.27, +1.43]). The isometry
contrast (mass declining less than wing³ predicts) is in the expected direction
in every specification and is marginal at M1 (z = 1.93, interval just
touching zero) but not at M2 (z = 1.13). The manuscript's "wings shorten while
mass does not" therefore cannot be presented as a demonstrated departure from
isometry; it is a consistent direction with an interval that includes zero,
driven by the wide interval on the mass slope (mass is ~2.5× more variable
than wing within species after controls: residual SD 3.7 % vs 2.2 %). The
brms bivariate posterior (Totoro) will give the definitive interval, but a
50-tree fit will not shrink the mass interval materially — the fast result is
the expected answer.

---

## 4. Comparison with the REVISION_PLAN.md §0 multi-trait table (`multitrait_results.rds$plan_comparison`)

Flag rule: coefficient differs by > 10 % or t by > 0.3.

| Trait | M0 (plan → here) | M1 +src +site (plan → here) | "fully controlled" (plan → here M3) | within-contributor (plan → here) |
|---|---|---|---|---|
| wing | −0.92 [t −3.8] → −0.92 [−3.8] ✔ | −0.53 [−2.6] → −0.53 [−2.6] ✔ | −0.33 [−1.9] → −0.30 [−1.6] ✔ | +0.02 [0.08] → +0.02 [0.06] ✔ |
| body mass | −0.0070 [−1.5] → −0.0072 [−1.5] ✔ | −0.0057 [−1.0] → −0.0056 [−1.0] ✔ | −0.0038 [−0.8] → −0.0045 [−1.0] ⚑ | −0.0035 [−0.7] → −0.0031 [−0.9] ⚑ |
| bill width | +0.23 [2.5] → +0.24 [2.6] ✔ | −0.04 [−0.5] → −0.02 [−0.3] ⚑ | −0.05 [−0.6] → −0.02 [−0.2] ⚑ | −0.05 [−0.7] → **+0.03 [+0.3]** ⚑ (sign) |
| bill length | +0.02 [0.4] → −0.00 [−0.0] ⚑ | +0.16 [1.7] → +0.16 [1.7] ✔ | +0.16 [1.7] → +0.17 [1.8] ✔ | +0.09 [0.8] → +0.14 [1.0] ⚑ |
| tail | +0.29 [1.0] → +0.28 [1.0] ✔ | +0.31 [1.2] → +0.30 [1.2] ✔ | +0.28 [1.1] → +0.25 [0.9] ⚑ | +0.12 [0.5] → +0.24 [1.0] ⚑ |
| tarsus | +0.42 [5.0] → +0.42 [5.0] ✔ | +0.24 [1.5] → +0.26 [1.6] ✔ | +0.22 [1.4] → +0.26 [1.5] ⚑ | +0.15 [0.9] → +0.26 [1.2] ⚑ |

- **M0 and M1 reproduce for every trait** (the bill-length M0 flag is −0.002 vs
  +0.02, both zero with t ≈ 0). Bill width "+contributor only" reproduces as
  −0.24 (t = −3.4) vs the plan's −0.28 (t = −3.8). Body-mass hour effect
  +0.404 %/h (t = 7.02, n = 8,345) vs the plan's +0.41 %/h (t = 7.25, n = 8,390;
  the plan evidently used all parsed hours — 8,391 here — rather than 05–19 h).
- The "fully controlled" and within-contributor columns differ in detail because
  the plan's scratch sequence added controls one at a time (wing-column proxy,
  then lon + alt, then season, then ring) and its Mundlak covariate set is not
  recorded; here M3 = M1 + lon + alt + season + ring on one complete-case sample
  and the Mundlak models use the M3 covariates. **No qualitative statement in
  the plan changes**: every within-contributor slope has an interval that
  includes zero, and every non-zero baseline trend (wing, bill width, tarsus)
  attenuates to a zero-including interval under M1. The bill-width
  within-contributor slope changes sign (+0.03 vs −0.05) — quote it as "≈ 0
  [−0.18, +0.23]", not as a negative value.
- **The wing attenuation from M1 (−0.53) to the plan's models 5–7 (−0.48 →
  −0.30 → −0.32) is mostly a change of SAMPLE, not of covariates.** `M1cc` —
  the M1 structure refitted on the 8,282 complete-case wing records — already
  gives −0.34 [−0.71, +0.02]; adding lon + alt + season on that same sample
  moves it only to −0.32, and ring to −0.30. The 196 records dropped for NA
  `Altitude` (190) or NA season (6) are mostly early: E. Carrano 1996–99
  (103), P. Serafini 2010–14 (40), C. Duca 2011–13 (23), M. Alves 1999 (8),
  L. Bugoni 1998/2004 (10); their species-sex-centred wing is **+0.45 mm** (mean
  69.2 vs 71.1 mm raw). Removing early, long-winged records flattens the trend
  by itself. Season IS a strong predictor of wing (M3: MAM +0.55, JJA +0.43,
  SON +0.57 mm vs DJF, t = 2.9–3.9), but on a fixed sample it barely moves the
  year slope. **P1 should fit the season/altitude steps on the M1 sample (e.g.
  an NA-altitude indicator or altitude imputed from the DEM the plan already
  mentions) before attributing the attenuation to season**, and P7 should not
  write "season attenuates the trend to −0.6 mm/decade" without that check.
- New observations, reported for completeness, not as claims: bill length is the
  only trait with a within-municipality slope excluding zero (+0.33 [0.11,
  0.56], t = 2.9; singular fit) and has a strong season effect (JJA −0.53 mm,
  t = −4.2); log mass increases with longitude (+0.063 per SD, t = 3.1);
  tail increases with altitude (+0.75 mm per SD, t = 3.1).

---

## 5. brms bivariate model — smoke test status (NOT results)

Two attempts, both local (laptop, rstan backend, 1 tree, chains = 2, iter = 400,
warmup = 200), while P1's `atlantic_parallel_controlled.R --smoke` Stan jobs were
running on the same machine:

1. **Full shared subset (7,577 records): did not finish within budget.** The
   bivariate model was still sampling after 25 min (tree retrieval and
   reconciliation took < 1 min; the rest was compile + sampling under
   `adapt_delta = 0.99`, `max_treedepth = 15`). Killed at 25 min per the ~20-min
   rule. No output written.
2. **`--smoke --subsample 800` (random 800 of the 7,577 records, 72 → fewer
   species): completed in 9.2 min, exit 0.** Bivariate fit 1 tree: max R-hat
   1.017, min bulk ESS 71, 0 divergences (fine for iter = 400); relative-wing
   fit 1.1 min, R-hat up to 1.08 (expected at 400 iterations). This verified
   the wiring end-to-end: the `mvbind` + `set_rescor(TRUE)` formula with
   per-response `gr(species_name, cov = A)`, `src` and `site` terms compiles;
   the parameter names used by the derived quantities exist
   (`b_wing_scaled_yr`, `b_lnmass_scaled_yr`, `rescor__wing__lnmass`,
   `sd_spp__wing_scaled_yr`, `sd_species_name__wing_Intercept`, `sigma_wing`,
   ...); `pool_rubin_draws()` and the posterior probabilities run; files are
   written to the namespaced smoke paths (`output/bivariate_smoke_results.rds`,
   `output/models/bivariate_smoke.rda`). **Its numbers mean nothing** (800
   records, 200 post-warmup draws per chain) and must not be quoted.
   The `--subsample` flag is refused without `--smoke`.

Defect found and fixed by the smoke test: the verbatim `pool_rubin()` from
`atlantic_parallel.R` returns NA standard errors when only one tree is fitted
(between-tree variance of a single value is NA). The script now routes m = 1
through `pool_rubin_draws()` (between-tree term = 0) via `pool_fixed()`, so the
`--trees 1` validation tier on Totoro returns intervals; for m > 1 the two
functions are numerically identical. The fix was unit-checked on the saved
smoke fit; the script was parse-checked but the 9-min smoke run was not
repeated after the edit.

Timing guidance for Totoro: the 7,577-record bivariate model did not complete
400 iterations × 2 chains sequentially in 25 min on a busy laptop. At the
production settings (4 chains × 4,000 iterations, cmdstanr, 4 cores per fit, 46
tree workers) budget several hours for the 50-tree run; run `--trees 1` first.

---

## 6. Decisions worth knowing

- **Bivariate structure.** `mvbind(wing, lnmass)` with `set_rescor(TRUE)`; each
  response has its own `(1 + scaled_yr || spp)`, `(1 | gr(species_name, cov = A))`,
  `(1 | src)`, `(1 | site)`. Group-level effects are **independent across
  responses** (no `|ID|` syntax) so the univariate Phase 1 / Phase 2 models are
  nested in it; only the residuals are correlated. Cross-response correlation of
  species slopes would be a separate question.
- **Priors.** Wing response as `atlantic_parallel.R` (Intercept N(71, 15)); log
  mass as `atlantic_parallel_mass.R` (Intercept N(3, 1.5)); b N(0, 10), sd and
  sigma half-Cauchy(0, 1) for both; `rescor ~ lkj(2)` (weakly regularising;
  brms default is lkj(1)).
- **Response names.** `conc.wing.length` and `ln_body_mass` are copied to `wing`
  and `lnmass` inside the script so brms parameter names are
  `b_wing_scaled_yr`, `b_lnmass_scaled_yr`, `rescor__wing__lnmass`.
- **Derived quantities** (posterior draws, then Rubin-pooled like any parameter):
  `wing_pct_decade = 100 · b_wing · (10/SD) / mean wing`;
  `mass_pct_decade = 100 · (exp(b_lnmass · 10/SD) − 1)`;
  `diff_pct_decade = wing − mass` (negative = wing shortens relative to mass);
  `iso_gap_pct_decade = mass_pct − 100·((1 + wing fraction)³ − 1)`;
  `iso_contrast_per_sdyr = b_lnmass − 3·b_wing/mean wing` (0 under isometry);
  plus over-record (23-yr) versions. `p_diff_negative` and `p_isogap_positive`
  are posterior probabilities averaged over trees. Intervals are estimate ±
  1.96 × Rubin SE as in `atlantic_parallel.R`; mean per-tree 2.5/97.5 %
  quantiles are stored alongside.
- **Tree cloud.** Identical retrieval / reconciliation code to
  `atlantic_parallel.R`; the only change is that when `AVESDATA_PATH` is unset
  the script uses `data/raw/AvesDataLite-main` if present (it is) instead of
  downloading. 100 trees pulled, `sample(…, N_TREES)` with `set.seed(20240101)`
  as the parent script. 0 unmatched species; 72 tips after alignment.
- **Smoke outputs are namespaced** (`bivariate_smoke_results.rds`,
  `models/bivariate_smoke.rda`) so a local smoke test can never overwrite a
  Totoro result.
- **glmmTMB joint model.** The textbook multivariate trick (`dispformula = ~0`
  with `us(0 + trait | rec)`) did not converge on these data (non-PD Hessian
  under REML and ML, raw or unit-SD responses); the shared record intercept
  with `dispformula = ~ trait` converges cleanly and reproduces the univariate
  slopes, at the cost of constraining the within-record covariance to one
  positive shared component. It is a screening tool; brms is the analysis.
- **Ring effect.** `interaction(Ring, Binomial)` per P0; unringed records get a
  unique level each (kept, uninformative). Several M3 / Mundlak fits are
  singular (ring or species-slope variance → 0); the fixed effects are still
  valid REML estimates and are flagged in `$table$singular`.
- **Complete cases.** M2/M3 drop records with NA Municipality, longitude,
  Altitude or season; `M1cc` isolates that (see §4, third bullet).

---

## 7. What downstream agents need

- **P1**: the complete-case finding in §4 (sample change, not season, drives most
  of the M1 → M2 attenuation); `multitrait_results.rds$table` has the wing rows
  (`trait == "wing"`) as an lme4 cross-check for M0–M3.
- **P6**: the wing-vs-mass slope-contrast figure should draw from the phylogenetic
  tier, `bivariate_phylo_results.rds$joint$M1_src_site$pooled_derived` (rows
  `wing_pct_per_sdyr`, `mass_pct_per_sdyr`, `diff_pct_per_sdyr`,
  `iso_contrast_pct_per_sdyr`; columns `pct_decade`, `pct_decade_lower`,
  `pct_decade_upper`), with `$contrasts_independent` for the M0–M2 ladder (§8.2).
  `bivariate_fast_lme4.rds$stacked_glmmTMB$M1_src_site` is the non-phylogenetic
  equivalent (identical to 2 decimals); `bivariate_results.rds$bivariate$derived`
  will be the brms cross-check if Totoro runs it.
  Bill-width before/after: `multitrait_results.rds$table` filtered to
  `trait == "Bill_width.mm."`, models `M0_baseline`, `M1a_src`, `M1_src_site`,
  `MWsrc_mundlak_contributor`.
- **P7**: quote the PHYLOGENETIC tier (§8; identical to the lme4 numbers below to
  two decimals, so the text can say "phylogenetic mixed models across 50 trees"):
  mass stable (M1 −0.0056 [−0.0164, +0.0051] per SD-yr = −1.1
  %/decade [−3.2, +1.0]; diurnal gain +0.40 %/h); bill width baseline +0.24
  [0.06, 0.42] → M1 −0.02 [−0.17, +0.13], 2 of 22 contributors span both
  periods (A. Piratelli, M. Alves; 620 records); tarsus +0.42 [0.26, 0.59] → M1
  +0.26 [−0.06, +0.58]. For isometry, do **not** claim a demonstrated departure:
  the slope difference on shared records is −0.9 %/decade [−3.3, +1.4] (M1) and
  the isometry contrast +3.9 [−0.1, +7.9]; report the direction and the
  interval. Replace "17 % expected vs 3 % observed" with the M1/M3 rows of §2
  and the shared-record contrast of §3.
- **Run order** (REVISION_PLAN.md §4, Step 5), all local now:
  `Rscript Analysis/scripts/atlantic_multitrait.R` (lme4 tier, 21 s) →
  `Rscript Analysis/scripts/atlantic_multitrait.R --trees 50` (phylogenetic tier,
  42 min) → `Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --fast-lme4`
  (22 s) → `Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --trees 50`
  (glmmTMB default, 38 min). The brms cross-check, if wanted on Totoro:
  `atlantic_bivariate_wing_mass.R --engine brms --trees 1` then `--trees 50`
  (hours).

---

## 8. Phylogenetic tier — glmmTMB `propto` across the 50 published trees (added 2026-09-09, afternoon)

**Why this exists.** The plan's phylogenetic models were brms 50-tree runs (hours
each, Totoro only). Williams, McGillycuddy, Drobniak, Bolker, Warton & Nakagawa
(2025, bioRxiv 10.64898/2025.12.20.695312) show that glmmTMB's `propto`
covariance structure fits the same phylogenetic GLMM in seconds with equivalent
estimates. `glmmtmb_validation_wing.R` confirmed this on our own data: the
published 50-tree brms wing model refitted in glmmTMB on the identical data and
per-tree covariance matrices gives year −0.922 [−1.400, −0.444] vs brms −0.922
[−1.403, −0.440], with Sex, latitude, species-slope SD, phylogenetic SD and sigma
agreeing to 2–3 decimals (phylo proportion 0.957 vs 0.95), in 24 s. Every model
of §2–§3 was therefore refitted with phylogeny through the shared engine
`scripts/_phylo_engine.R` (`phylo_species_name()`, `phylo_A_list()`,
`run_phylo_trees()`, `pool_rubin_df()`): same fixed and random structure as the
lme4 tier plus `propto(0 + species_name | g, A)`, REML, on the **identical 50
clootl trees of the published brms analysis** (correlation matrices cached from
`models/brm0_multiphylo.rda` to `data/derived/phylo_A_50trees.rds`), pooled with
Rubin's rules (interval = estimate ± 1.96 × Rubin SE; the `t`/`z` column is the
pooled z). glmmTMB 1.1.14, R 4.6.0, single core, laptop.

Files: `output/multitrait_phylo_results.rds` (`$generated` 2026-09-09 14:44:25;
`$table`, `$varcomp`, `$hour`, `$comparison_lme4`, `$timing`, `$per_tree`,
`$meta`) and `output/bivariate_phylo_results.rds` (`$generated` 14:40:17;
`$univariate{pooled, varcomp, per_tree}`, `$contrasts_independent`, `$joint`,
`$comparison_lme4`, `$timing`, `$definitions`). `multitrait_table.md` now carries
the phylogenetic table under the lme4 one. The lme4 tier (`multitrait_results.rds`)
was regenerated at 13:33 by the restructured script and is numerically identical
to the 11:54 file (`all.equal` on `$table`, `$hour`, `$isometry`, `$coverage`,
`$plan_comparison`; only `$meta$generated` differs).

### 8.1 Multi-trait screen with phylogeny (`multitrait_phylo_results.rds$table`)

Sample = the 12,571 live known-sex records, 73 species (72 with bill length /
tail, 71 tarsus, 68 bill width; `Herpsilochmus_sellowi` is not in the tree and has
0 records here). Wing is left to Phase 1 (`controlled_wing_phylo_results.rds`).

| Trait | M0 baseline β_yr | M1 +src +site | M3 +lon +alt +season +ring (complete cases) | Within-contributor (Mundlak, M3 covariates) | Between-contributor | Phylo proportion at M1 (phylo SD) |
|---|---|---|---|---|---|---|
| log body mass (ln g) | −0.0072 [−0.0165, +0.0021] z −1.5 (n 11,256) | −0.0056 [−0.0164, +0.0051] z −1.0 | −0.0046 [−0.0138, +0.0047] z −1.0 (n 11,077) | −0.0032 [−0.0104, +0.0041] z −0.9 | +0.0060 [−0.0305, +0.0424] z 0.3 | 0.954 (0.829) |
| bill width (mm) | **+0.24 [+0.06, +0.42] z 2.6** (n 3,209) | **−0.03 [−0.18, +0.12] z −0.4** | −0.03 [−0.18, +0.13] z −0.3 (n 3,205) | +0.02 [−0.19, +0.22] z 0.2 | +0.41 [−0.48, +1.30] z 0.9 | 0.539 (2.59) |
| bill length (mm) | −0.00 [−0.11, +0.11] z 0.0 (n 7,697) | +0.16 [−0.02, +0.34] z 1.8 (49/50 trees) | +0.17 [−0.02, +0.36] z 1.8 (n 7,540) | +0.15 [−0.12, +0.42] z 1.1 | +0.07 [−0.39, +0.52] z 0.3 | 0.600 (4.99) |
| tail (mm) | +0.27 [−0.27, +0.82] z 1.0 (n 8,872) | +0.30 [−0.19, +0.80] z 1.2 | +0.24 [−0.28, +0.76] z 0.9 (n 8,679) | +0.24 [−0.22, +0.70] z 1.0 | −0.20 [−1.34, +0.94] z −0.3 | 0.916 (28.0) |
| tarsus (mm) | **+0.42 [+0.24, +0.60] z 4.7** (n 4,361) | +0.26 [−0.06, +0.58] z 1.6 | +0.26 [−0.08, +0.60] z 1.5 (n 4,285) | +0.27 [−0.16, +0.70] z 1.2 | +0.02 [−0.67, +0.70] z 0.0 | 0.803 (7.82) |

Per SD-year (SD 5.0199); × 1.992 for per decade. In % per decade (`$table$pct_per_decade`),
M1: mass −1.1 [−3.2, +1.0]; bill width −0.8 [−4.6, +3.1]; bill length +2.5 [−0.3, +5.3];
tail +1.0 [−0.6, +2.7]; tarsus +2.2 [−0.6, +5.0].

Convergence (`$timing`; "converged" = positive-definite Hessian, `$sdr$pdHess`,
with finite SEs): **21 of 22 models 50/50 trees; bill length M1 49/50** (one
tree still without a positive-definite Hessian after two seeded restarts; it is
excluded from the Rubin pooling, its point estimate kept in `$per_tree`). The
first-pass counts before the seeded-restart repair were lower for several models
(tail M0 21/50, bill length M1 and M3 31/50, tarsus MWsrc 39/50, tail M1 36/50;
50/50 for every mass and bill-width model) — see "Two REML modes" below. Wall
time per model, all 50 trees: 15 s (bill width M0) to 303 s (bill length M3);
0.3–6.1 s per tree-fit; total 41.7 min.

**Body mass, capture hour, with phylogeny** (`$hour`; n = 8,345, 70 species,
35 contributors, 111 municipalities; 50/50 trees both models): **+0.404 % body
mass per hour [0.291, 0.517], z = 7.03**; year slope with hour in the model
−0.0063 [−0.0188, +0.0062] (z −0.99), without hour on the same sample −0.0076
(z −1.17). Identical to the lme4 tier (+0.404 [0.291, 0.516], t 7.02).

**Phylogenetic signal.** The phylogenetic proportion (phylo variance / total
species + residual variance) is 0.95 for log mass, 0.92 for tail, 0.80 for
tarsus, 0.60 for bill length and 0.54 for bill width — i.e. species means of
mass, tail and tarsus are almost entirely phylogenetically structured, bill
dimensions much less so. The non-phylogenetic species intercept SD is estimated
at ≈ 0 in every model once the phylogenetic term is present (as in the
validation), so the (1 | spp) intercept is redundant with `propto` here; the
species year-slope SD is retained and sits at the boundary (< 1e-3 × sigma; the
lme4 "singular" analogue) for bill length M1 and tarsus M1/M3, exactly the fits
lme4 flagged singular.

**Does phylogeny change anything? No.** `$comparison_lme4` matches 28 estimates
(5 traits × {M0, M1, M3, within, between} + the 3 hour-model terms) to the lme4
tier on the same sample: median |Δβ| / lme4 SE = **0.011**, maximum **0.106**
(bill width M1: −0.030 vs −0.022, both ≈ 0); median SE ratio phylo / lme4 =
**1.002** (range 0.993–1.047, the maximum being tarsus M0); **the zero-exclusion
verdict changes in 0 of 28**. Every statement of §2 and §4 stands verbatim under
phylogeny: mass is flat in every specification; bill width's baseline +0.24
collapses to −0.03 with contributor and site intercepts and the within-contributor
slope is +0.02 [−0.19, +0.22]; tarsus's +0.42 attenuates to +0.26 [−0.06,
+0.58]; bill length and tail never leave zero. This is expected: the
phylogenetic term re-partitions the between-species variance (almost all of the
iid species-intercept variance moves into the phylogenetic component) but the
year slope is identified within species, so its estimate and SE are essentially
untouched.

### 8.2 Shared-record isometry test with phylogeny (`bivariate_phylo_results.rds`)

Subset as §3: **7,577 records with both wing and mass, 72 species, 42
contributors, 128 municipalities** (complete-case 7,409); mean wing 71.48 mm,
mean mass 23.0 g. Every univariate model: 50/50 trees converged (M0 wing 29/50 on
the first pass, all 21 repaired by the seeded restart; every other model 50/50
first pass); 30–105 s per model (0.6–2.1 s per tree-fit). Total run incl. the
joint models 37.6 min.

**Univariate phylogenetic models** (`$univariate$pooled`):

| Model | wing β_yr (mm/SD-yr) | log-mass β_yr | wing \| log-mass β_yr | allometric b(ln mass) | phylo prop wing / mass |
|---|---|---|---|---|---|
| M0 baseline (n 7,577) | −0.83 [−1.32, −0.33] z −3.3 | −0.0021 [−0.0103, +0.0061] z −0.5 | −0.81 [−1.29, −0.33] z −3.3 | 3.21 [2.56, 3.86] | 0.957 / 0.966 |
| **M1 +src +site (n 7,577)** | **−0.55 [−0.97, −0.13] z −2.6** | **−0.0027 [−0.0135, +0.0081] z −0.5** | **−0.54 [−0.94, −0.14] z −2.6** | **3.25 [2.62, 3.89]** | 0.951 / 0.949 |
| M1cc (M1 on n 7,409) | −0.35 [−0.76, +0.07] z −1.6 | −0.0035 [−0.0118, +0.0049] z −0.8 | −0.34 [−0.74, +0.06] z −1.7 | 2.86 [2.21, 3.50] | 0.951 / 0.950 |
| M2 +lon +alt +season (n 7,409) | −0.34 [−0.75, +0.07] z −1.6 | −0.0030 [−0.0114, +0.0054] z −0.7 | −0.33 [−0.73, +0.06] z −1.6 | 2.86 [2.22, 3.50] | 0.951 / 0.949 |

Phylogenetic SD ≈ 22.4–22.9 mm for wing and 0.82–0.83 ln units for mass in every
specification. The relative-wing model says the same thing it did without
phylogeny: conditioning on the bird's own log mass (b ≈ 3.2 mm per ln-unit at
M1) leaves the wing year slope unchanged (−0.54 vs −0.55), so whatever wing
trend there is at M1 is not a by-product of a mass trend; at M2 both are
−0.33/−0.34 with intervals including zero.

**Difference of year slopes, % per decade.** Two estimators, both reported:

- `$contrasts_independent` — from the two univariate fits, **assuming
  Cov(b_wing, b_mass) = 0** (the two models are fitted separately, so the
  covariance is not estimated). Wing %/decade = 100·b·1.992/71.48 (linear);
  mass %/decade = 100·(exp(b·1.992) − 1) with delta-method SE
  100·exp(b·DEC)·DEC·se; SE of the difference = √(se_w² + se_m²). The joint model
  estimates the slope correlation at 0.04, so the independence assumption costs
  almost nothing here.
- `$joint` — a **joint long-format glmmTMB model** (rows = record × trait, both
  responses on the 100·ln scale so slopes are % per SD-year): trait-specific
  fixed effects, `diag(0 + trait:scaled_yr | spp)`, `diag(0 + trait | src)`,
  `diag(0 + trait | site)`, a shared record intercept `(1 | rec)` inducing the
  within-record correlation, `dispformula = ~ trait`, and **two phylogenetic
  terms** `propto(0 + species_name:is_w | g, A_w)` + `propto(0 + species_name:is_m | g, A_m)`
  (trait-specific phylogenetic species effects, each with its own variance, on
  the same tree). The slope difference and the isometry contrast use the
  estimated covariance of the two year slopes. **It converged: 50/50 trees
  usable and pooled for both M1 (first pass 31, 19 repaired by a seeded
  restart, 0 excluded; 679 s) and M2 (first pass 27, 23 repaired, 0 excluded;
  752 s).** So the covariance-correct estimator is available and is the one to
  quote; the independence estimator is the cross-check.

| Specification | wing %/decade | mass %/decade | **difference wing − mass** | isometry contrast mass − 3·wing (0 under isometry) | slope cor |
|---|---|---|---|---|---|
| M0, independent (n 7,577) | −2.30 [−3.69, −0.92] | −0.42 [−2.05, +1.20] | −1.88 [−4.01, +0.25] z −1.73 | +6.49 [+2.03, +10.95] z 2.85 | — |
| M1, independent (n 7,577) | −1.53 [−2.69, −0.36] | −0.54 [−2.68, +1.60] | −0.99 [−3.42, +1.45] z −0.79 | +4.04 [−0.06, +8.15] z 1.93 | (0) |
| **M1, joint, 50 trees (n 7,577)** | **−1.49 [−2.64, −0.33]** | **−0.58 [−2.66, +1.50]** | **−0.91 [−3.24, +1.43], z −0.76** | **+3.88 [−0.08, +7.84], z 1.92** | 0.039 |
| M1cc, independent (n 7,409) | −0.96 [−2.11, +0.18] | −0.69 [−2.35, +0.97] | −0.28 [−2.30, +1.74] z −0.27 | +2.21 [−1.62, +6.03] z 1.13 | — |
| M2, independent (n 7,409) | −0.95 [−2.09, +0.20] | −0.60 [−2.26, +1.07] | −0.35 [−2.38, +1.67] z −0.34 | +2.25 [−1.57, +6.07] z 1.15 | (0) |
| **M2, joint, 50 trees (n 7,409)** | −0.88 [−1.92, +0.15] | −0.66 [−2.31, +1.00] | **−0.23 [−2.14, +1.69], z −0.23** | **+1.99 [−1.47, +5.45], z 1.13** | 0.047 |

Joint-model variance components at M1 (100·ln scale): phylogenetic SD wing 31.6,
mass 82.0; contributor SD 3.3 / 11.1; site SD 2.4 / 4.4; record SD 2.9; residual
SD 2.16 (wing) / 3.74 (mass); implied within-record wing–mass correlation 0.50.
As in §3, the shared record intercept constrains the within-record covariance
to one positive component; the brms `rescor` model (`--engine brms`, not run)
would relax that, but the slope covariance it governs is tiny (r = 0.04), so the
interval on the difference would not move materially.

**Isometry statement, written honestly.** On the same individuals, with
contributor and municipality intercepts, species intercepts and slopes, and
phylogeny across 50 trees, wing length declines by 1.5 % per decade [0.3, 2.6]
while body mass changes by −0.6 % per decade [−2.7, +1.5]. **The difference
between the two year slopes is −0.9 % per decade with a 95 % interval of
[−3.2, +1.4] (z = −0.76); it does not exclude zero.** Under the isometric
expectation (mass ∝ wing³) the mass slope should be three times the wing slope
in % terms; the contrast mass − 3·wing is +3.9 % per decade [−0.1, +7.8]
(z = 1.92) at M1 — in the direction of "mass declines less than isometry
predicts", with an interval that just includes zero — and +2.0 [−1.5, +5.5]
(z = 1.13) once longitude, altitude and season are added (M2, complete cases).
Only the uncontrolled baseline (M0, no contributor or site terms) gives an
isometry contrast that excludes zero (+6.5 [+2.0, +10.9]), and that is the
specification the revision has already set aside. The manuscript's "wings
shorten while mass does not" is therefore a statement about two separate
tests (wing interval excludes zero at M1, mass interval does not), not a
demonstrated departure from isometry: the direct test of the difference is
compatible with no difference in every controlled specification, because the
mass slope is ~2.5× less precise than the wing slope (residual SD 3.7 % vs
2.2 %). Phylogeny does not change this; the lme4/glmmTMB tier of §3 gave
−0.92 [−3.27, +1.43] and +3.90 [−0.07, +7.88]. Report the direction and the
interval; do not write "departure from isometry demonstrated".

**lme4 tier vs phylogenetic tier on the shared subset** (`$comparison_lme4`):
16 matched univariate estimates (year slopes and allometric b, 4 specifications ×
3 responses), median |Δ| / lme4 SE = **0.010**, max **0.031**; SE ratio 1.000–1.014;
**verdict changes 0 of 16**. Slope differences: independent M1 lme4 −1.00 [−3.43,
+1.43] vs phylo −0.99 [−3.42, +1.45]; joint M1 lme4 −0.92 [−3.27, +1.43] vs phylo
−0.91 [−3.24, +1.43]; joint M2 −0.24 [−2.16, +1.68] vs −0.23 [−2.14, +1.69];
isometry contrast joint M1 +3.90 [−0.07, +7.88] vs +3.88 [−0.08, +7.84]. Identical
to two decimals throughout.

### 8.3 Engine notes (for GLMMTMB_ENGINE.md; `_phylo_engine.R` was not edited)

- **Two REML modes.** With both the iid species intercept `(1 | spp)` and the
  phylogenetic term (the published structure), the REML surface has two optima:
  the correct one (phylogenetic SD ≈ trait-scale, iid species SD ≈ 0) and a
  worse local one (phylogenetic SD ≈ 0, iid SD absorbing it, objective ≈ 24
  units higher, Hessian not positive-definite) that glmmTMB's default all-zero
  start reaches on some trees (e.g. 21 of 50 for M0 wing on the shared subset,
  29 of 50 for tail M0). The year slope differs by ~0.002 between modes, so the
  point estimates are safe, but the SEs of the bad mode are undefined and
  `pool_rubin_df()` would return NA. Both scripts wrap `run_phylo_trees()` in a
  `run_phylo_robust()` that refits every non-converged tree from the free
  variance parameters of a converged tree of the same model (then a generic
  seed), inserted into that tree's own theta vector (the `propto` entries encode
  A and are tree-specific), keeps the refit only if it converges and its
  objective is not worse (tolerance 0.01), excludes any tree still unusable
  from the pooling, and records `n_converged_first_pass`, `n_repaired`,
  `excluded_trees`. **Candidate addition to `_phylo_engine.R`** (identical code
  in both Phase-2 scripts; any other script combining `(1 | spp)` with the
  phylogenetic term will hit the same mode). Ideally the
  engine would also accept a `start` argument or drop `(1 | spp)` when the
  phylogenetic term is present, since that intercept is estimated at ~0 anyway.
- **Multiple `propto` terms.** `fit_phylo_glmmtmb()` appends exactly one
  `propto(0 + species_col | g, A)`. The joint model needs one per trait with
  matrix dimnames matching the interaction design columns
  (`species_name<tip>:is_w`); it is fitted directly in the script with
  `glmmTMB()` and pooled with the engine's `pool_rubin_df()`.
- **Non-phylogenetic species intercepts in the joint model** were dropped
  (`diag(0 + trait | spp)` gave a positive-definite Hessian but a non-invertible
  vcov; without them slopes identical to 4 decimals on a 1-tree probe).
- glmmTMB prints a "non-positive-definite Hessian" warning for every first-pass
  fit that later gets repaired; the bivariate run ended with "50 or more
  warnings" for that reason. Exit status 0 in both runs (scratch logs
  `bivariate_phylo_run3.log`, `multitrait_phylo_run3.log`).

### 8.4 What changes for downstream agents

- **Nothing in the numbers quoted in §7 changes; they can now be attributed to
  phylogenetic mixed models across 50 trees** (quote the `multitrait_phylo_results.rds`
  and `bivariate_phylo_results.rds` values; they equal the lme4 values to two
  decimals, differences are in the third decimal).
- **P7 isometry sentence**: use the joint M1 phylogenetic estimate — difference
  −0.9 % per decade [−3.2, +1.4]; isometry contrast +3.9 [−0.1, +7.8] — with
  the wording of §8.2. Delete "17 % expected vs 3 % observed".
- **P6**: `bivariate_phylo_results.rds$joint$M1_src_site$pooled_derived` for
  the slope-contrast panel; `$contrasts_independent` for the M0 → M2 ladder;
  `multitrait_phylo_results.rds$table` (`trait == "Bill_width.mm."`, models
  `M0_baseline`, `M1_src_site`, `MWsrc_mundlak_contributor`) for the bill-width
  before/after.
- **Synthesis / decision gate**: the phylogenetic tier of Phase 2 confirms the
  lme4 reading with no exceptions — mass flat, bill width and tarsus trends are
  contributor-composition effects, isometry not demonstrated on shared records.
- **brms path**: `atlantic_bivariate_wing_mass.R --engine brms --trees 50` remains
  the Bayesian cross-check for Totoro (rescor bivariate + relative-wing); it was
  not run and nothing in this section depends on it.
