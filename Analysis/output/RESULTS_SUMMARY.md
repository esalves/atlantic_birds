# Atlantic Forest passerines — curated results summary

**Updated:** 2026-06-18 (from the Totoro server run; source files in this folder).
This synthesises all current analyses. Per-run machine outputs are cited under
*Provenance*; this is the curated cross-analysis narrative.

## Dataset

Raw ATLANTIC_BIRD_TRAITS: 72,483 records → 60,066 passerine (83%) → 52,111 adult.
Analysis sample: **73 species, 12,571 records, 1995–2018** (wing: 8,478 records;
bill: 3,209). Species-level direction over the period: wing length decreased in
43 species / increased in 18; bill decreased in 13 / increased in 31.

## 1. Primary temporal trends (brms 50-tree multiphylo, Rubin-pooled)

| Trait | β per SD-year [95% CI] | Direction | Phylo signal |
|---|---|---|---|
| **Wing length** | **−0.92 [−1.40, −0.44]** | **decline** (CI excludes 0) | 0.95 [0.94, 0.96] |
| Bill | +0.24 [0.06, 0.43] | increase (excludes 0) | 0.41 [0.34, 0.48] |
| **log(body mass)** | **−0.007 [−0.017, 0.002]** | **no trend** (CI crosses 0) | 0.96 [0.95, 0.97] |

Sex and latitude (wing): male +2.64 [2.44, 2.85]; latitude −0.68 [−0.89, −0.46].
**Wing got shorter over time but body mass did not change** — a shape change
(relatively shorter wings at constant mass), with bill increasing.

## 2. Variance / morphological variability

- **Wing lnCVR meta-analysis = +0.183 [0.033, 0.334]** (k = 54, I² = 93.7%) →
  among-individual **variability in wing length increased** over the period.
- Bill lnCVR = +0.107 [−0.196, 0.411] (k = 36) → not significant.
- drmSEM σ-channel (full data) corroborates: `year → σ(wing)` = **+0.035,
  p < 1e-5** (variance rising), and `temperature → σ(wing)` = **−0.193,
  p < 1e-98** (warmer → less variable). The temperature→σ result is the single
  most robust effect across every model run.

## 3. Drivers — temperature and arthropods (drmSEM)

drmSEM was run in two scopes (separate caches/figures, `_arthro` / `_noarthro`):

| | Full-data model | Arthropod model |
|---|---|---|
| N / species / years | 7,810 / 72 / 1995–2018 | 1,697 / 68 / 1998–2008 |
| **Fisher's C** | 0.51, df 2, **p = 0.78 (clean)** | 13.5, df 4, **p = 0.009 (rejected)** |
| year → temperature | −0.005 (n.s.) | −0.048 (n.s.) |
| **temperature → wing (mean)** | −0.216 (p = 0.058) | −0.79 (p < 1e-3) |
| arthropod → wing (mean) | — | +0.07 (p = 0.56, n.s.) |
| year → σ(wing) | **+0.035 (p < 1e-5)** | −0.221 (p < 1e-9) |
| temperature → σ(wing) | **−0.193 (p < 1e-98)** | −0.040 (p = 0.16, n.s.) |
| year → wing (total indirect) | +0.001 [−0.001, 0.005] n.s. | −0.024 [−0.174, 0.123] n.s. |

Key points:
- **Temperature does not mediate the wing decline.** `year → temperature` ≈ 0 in
  both scopes, so the indirect `year → temperature → wing` effect is ≈ 0. The
  size decline seen in the brms model (§1) is a **direct temporal trend not
  explained by the measured climate or food mediators** (no contradiction — the
  brms model regresses wing directly on year, whereas drmSEM tests mediation).
- The mean Bergmann effect (`temperature → wing`) is significant only in the
  small arthropod subset (−0.79) and merely marginal on the full data (−0.22,
  p = 0.058) → **fragile / subset-dependent**, not a headline.
- **No food-limitation effect on body size** (`arthropod → wing` ≈ 0), despite a
  strong arthropod decline (below).
- The **arthropod-model DAG is now rejected** (p = 0.009; the `Sex ⊥ arthropod`
  claim fails on the smaller updated subset), so the full-data model is the
  trustworthy drmSEM. Phylogeny barely matters (λ = 1 on all trees; phylo ≈ i.i.d.).

**Arthropod abundance** itself declined steeply (PREDICTS, n = 63,150):
β = −1.40 [−1.46, −1.34] per SD-year → to **×0.25 [0.23, 0.26]** of baseline
(sensitivity variant: ×0.30). Arthropods crashed; birds' wings shrank; but the
two are not linked through the body-size path here.

## 4. Diet moderation of the wing decline

Diet × year interaction (65 species, 8,086 records; 30 invertebrate-leaning vs
35 "other"):
- Continuous: year main = −0.86 [−1.32, −0.39]; **year × diet_inv = +0.66
  [0.21, 1.11]** → more invertivorous species decline **less**.
- Categorical: year × Invertebrate = +1.02 [0.05, 2.00] (same direction).

So the wing decline is **strongest in non-invertivores** and buffered in
invertebrate eaters.

## 5. Sensitivity to inclusion thresholds

The wing decline is robust across min-records (30/40/50) and min-year-range
(5/8/10) thresholds **when the sex filter is applied** (β from −0.82 to −0.96,
all CIs exclude 0). **Without the sex filter** it weakens and is often
non-significant (β −0.26 to −0.43, several CIs cross 0; relaxed 158-species
model: −0.31 [−0.82, 0.19]). → The decline is **threshold-robust but sensitive
to sex filtering** — flag prominently; sex composition across years matters.

## Bottom line

Over 1995–2018, Atlantic Forest passerines show **shorter wings** (−0.92
SD-year, robust) at **unchanged body mass** and **slightly longer bills**, with
**rising among-individual variability** in wing length (lnCVR +0.18; drmSEM
`year → σ` +0.035). The decline is **not explained by temperature or arthropod
abundance** in the SEM (both mediation paths ≈ 0), is **buffered in
invertivores**, and is **sensitive to sex filtering**. The most robust climate
signal is distributional, not in the mean: **warmer → less variable wings**
(`temperature → σ` ≈ −0.19, p < 1e-98). Arthropods collapsed to ~25% of baseline
but that does not propagate to body size here.

## Provenance (this folder)

`descriptive_summary.rds` (primary trends, lnCVR), `mass_results.{rds,md}`,
`diet_interaction_results.rds`, `arthropod_estimates.rds`,
`sensitivity_thresholds_tier1/2.rds`, `drmsem_results_{arthro,noarthro}.{rds,md}`,
phylo/effects caches, and `../figures/drmsem_*_{arthro,noarthro}.png` (8 figures
regenerated 2026-06-18 from the saved results via
`../scripts/make_figures.R`, so they match the current result data). The earlier
drmSEM-only curated doc is retained in `../../archive/atlantic_drmsem_results.md`.

---

## Revision diagnostics, fast tier (2026-09)

**Added 2026-09-09** in response to the referee report (`REVISION_PLAN.md`; status
in `REVISION_STATUS.md`; per-phase notes `../scripts/REVISION_NOTES_P*.md`).
Everything in this section is a **fast-tier result: `lme4` / `glmmTMB` REML with
Wald ± 1.96 SE intervals and no phylogenetic term** (species intercepts absorb
it; the lme4 baseline reproduces the brms −0.92 to two decimals). §1–§5 above are
the June 2026 published-baseline narrative and are **qualified, not replaced**,
by this section.

**Phylogenetic tier, added same day (afternoon).** The 50-tree `brms` path
originally planned for this section's inferential numbers is **replaced by
default** with `glmmTMB`'s `propto` covariance structure (Williams,
McGillycuddy, Drobniak, Bolker, Warton & Nakagawa 2025), validated to reproduce
the published 50-tree `brms` wing model to 2–3 decimals in ~24 s instead of hours
(`glmmtmb_validation_wing.rds`; methods in `../scripts/GLMMTMB_ENGINE.md`; full
account in `REVISION_STATUS.md` §3.9). Phylogenetic-tier numbers are added
**alongside** the fast-tier lme4/glmmTMB numbers below, in-line, wherever a
50-tree phylogenetic run has completed — all five scripts (controlled wing, multi-trait,
bivariate wing-mass isometry, variance σ ladder, diet interaction) finished at 50 trees
on 2026-09-09 (the wing and variance ladders at 16:48 and 16:12). One `brms` cross-check (a Bayesian re-fit of
the primary controlled model, M3) remains planned and has not run.

Units as above: β per SD-year (SD = 5.020 yr; × 1.992 for per decade), mean wing
71.13 mm, live-only 73 species / 12,571 records / 8,478 wing records, 1995–2018.

### R1. Provenance audit (`audit_*.rds`, `effect_scale.rds`)

- **42 contributors** hold the wing records; **6 span both the early (≤ 2006) and
  late (≥ 2013) periods** (2,828 of 8,478 records). The wing column populated —
  a protocol proxy — flips **87.4 % right-wing early → 67.5 % unspecified late**.
- Of 139 named localities in the quartile comparison **6** occur in both periods
  (684 records); municipalities 15 of 99 (1,728 records); coordinate sites 4 of
  233. Later records are further north and higher (mean latitude −24.8° → −21.9°,
  median altitude 175 → 489 m).
- Individuals: 927 ringed birds have > 1 record (1,813 repeats); recapture share
  4.8 % → 12.3 %. Season shifts (DJF 12.6 % → 21.6 % of wing captures); sex ratio
  45.6 → 48.7 % female.
- Effect scale of the published wing trend: −0.92 per SD-year = **−1.83 mm/decade
  = −2.58 %/decade = −4.2 mm (−5.9 %, 95 % CI −9.0 to −2.8 %) over 1995–2018**.
  Mass −0.0073 → −1.44 %/decade, −3.3 % over record [−7.4, +1.0]. Isometric
  expectation for −5.9 % wing: −16.7 % mass.

### R2. Controlled wing trend (`controlled_wing_results.rds`; Phase 1, Tier 1)

Corrected 2026-09-09: an earlier version of this section quoted a 12:02 rds produced
with `Sex` silently dropped from every known-sex model (REVISION_REVIEW.md H1–H3); the
guard was fixed, `--fast-lme4` re-run, and every number below is from the regenerated
rds (`$generated` 12:35:53, `Sex` a fixed effect in all 35 known-sex specs), which
left the decision-gate reading unchanged. The 12:35 files are quotable at Tier-1 status.

| Model (cumulative) | β_year | 95 % CI | N | mm/decade |
|---|---:|---|---:|---:|
| M0 baseline (published structure) | −0.922 | [−1.400, −0.444] | 8,478 | −1.84 |
| M1 + (1\|contributor) + (1\|municipality) | −0.528 | [−0.923, −0.134] | 8,478 | −1.05 |
| M2 + wing column + longitude + altitude | −0.500 | [−0.894, −0.105] | 8,478 | −1.00 |
| **M3 + season + moult + (1\|individual)** | **−0.370** | **[−0.775, +0.034]** | 8,478 | **−0.74 [−1.54, +0.07]** |
| M3 on complete cases | −0.284 | [−0.643, +0.075] | 8,282 | −0.57 |
| M4 = M3 + contributor year slopes | −0.873 | [−1.750, +0.004] | 8,478 | −1.74 |
| M5 within-municipality year (Mundlak) | −0.735 | [−1.679, +0.209] | 8,478 | −1.46 |
| M5src within-contributor year (Mundlak) | −0.157 | [−0.797, +0.482] | 8,478 | −0.31 |
| M5src_cc within-contributor, complete cases | +0.028 | [−0.573, +0.628] | 8,282 | +0.06 |
| M6 unknown-sex records, M3 structure (no Sex term) | **−0.684** | **[−1.196, −0.172]** | 3,532 | −1.36 |
| M7 first wing captures only | −0.379 | [−0.781, +0.024] | 7,791 | −0.75 |
| M3, six long-running contributors only | +0.101 | [−0.517, +0.718] | 3,017 | +0.20 |
| E. Carrano alone (1,040 records, 1996–2017) | −0.029 | [−0.266, +0.208] | 1,040 | −0.06 |
| M3 without 6 anomalous contributor×site×year blocks | −0.428 | [−0.808, −0.048] | 8,285 | −0.85 |

Key points:
- **Contributor and municipality intercepts remove roughly half of the wing
  decline; the full control set leaves ~40 % of the baseline with an interval
  that includes zero** (`share_of_M0_remaining_in_M3` = 0.40; complete cases 0.31).
  Decision gate on this tier: **Scenario B** on both samples (REVISION_PLAN §2),
  provisional until the brms tiers run.
- The within-contributor slope overlaps zero in every specification (−0.41 to
  +0.03 across variants); the between-contributor component carries the signal
  (−1.29 [−2.41, −0.16]). The three largest spanning contributors show no trend
  (E. Carrano −0.008 [−0.25, +0.24]; A. Ross −0.29 [−1.28, +0.69]; M. Alves −0.07
  [−0.47, +0.32]); the two negative within-contributor series come from late blocks
  sitting 9–12 mm below species × sex means (C. Fontana São Francisco de Paula
  2016–17, A. Piratelli Itu 2016) — data-quality observations, not phenotypes.
- Much of the plan's "season attenuates the trend" step was a **sample change**:
  dropping 190 altitude-NA records (103 = E. Carrano 1996–99 Ilha Rasa) and 6
  undated records moves M1 from −0.53 to −0.34 by itself. Season is a strong wing
  predictor (M3_cc: MAM +0.59, JJA +0.40, SON +0.54 mm vs DJF) but barely moves the
  year slope on a fixed sample.
- **Two results cut the other way and must be reported**: the unknown-sex
  replication keeps a clear decline under full controls (M6 above; its baseline is
  flat, −0.016, so the controls *reveal* the trend there — the mirror image of the
  known-sex pattern), and removing the six sign-blind anomalous blocks makes M3
  exclude zero. Leave-one-contributor-out: M3 ranges −0.54 … −0.26, most negative
  without E. Carrano (4 of 42 refits exclude zero).
- Species slopes (M3, lme4 BLUP quadrature intervals — descriptive, not inference):
  49 of 72 negative, 9 intervals below zero and 2 above, median −0.56 mm/decade
  (M0: 56 of 72, 17 / 2, median −1.18).

**Phylogenetic tier (glmmTMB `propto`, 50 trees, `controlled_wing_phylo_results.rds`, 16:48; 39 specs, all 50/50 trees converged except one tree of M5_cc).** M0 −0.922 [−1.400, −0.444]; M1 −0.528 [−0.924, −0.133]; M2 −0.500 [−0.896, −0.104]; **M3 −0.371 [−0.777, +0.036]** (−0.74 mm/decade); M3_cc −0.285 [−0.648, +0.078] (−0.57 mm/decade); M4 −0.875 [−1.758, +0.007]; within-municipality (M5) −0.738 [−1.685, +0.209], rsTotal −0.402 [−0.816, +0.012]; within-contributor (M5src) −0.159 [−0.798, +0.481], between-contributor −1.290 [−2.421, −0.159]; within-contributor cc +0.026 [−0.578, +0.629]; M7 first captures −0.379 [−0.783, +0.025]; six long-running contributors +0.095 [−0.526, +0.715]; M3_noanom −0.432 [−0.815, −0.049]; M6 unknown-sex −0.690 [−1.206, −0.174]; M6_cc −0.587 [−1.124, −0.050]; phylo proportion 0.91–0.96. `$decision_gate$scenario` = **B** on both samples. Every year term is within 0.002 of the lme4 estimate above. Species slopes from the M3 fit (fixed + BLUP, quadrature intervals, not inference): 49 of 72 negative, 9 below zero, 2 above, median −0.54 mm/decade (M0: 54 / 14 / 2, −1.18).

### R3. Multi-trait table (`multitrait_results.rds`, `multitrait_table.md`; lme4 REML;
phylogenetic tier `multitrait_phylo_results.rds`, glmmTMB `propto`, **50 trees, complete**)

| Trait (n) | M0 baseline | + contributor + municipality | fully controlled (M3), lme4 | **fully controlled (M3), phylogenetic (50 trees)** | within-contributor |
|---|---|---|---|---|---|
| log body mass (11,256 / 11,077 at M3) | −0.0072 [−0.0165, +0.0021] | −0.0056 [−0.0164, +0.0051] | −0.0045 [−0.0138, +0.0047] | **−0.0046 [−0.0138, +0.0047], z −0.97, 50/50** | −0.0031 [−0.0102, +0.0039] |
| bill width (3,209 / 3,205) | **+0.24 [0.06, 0.42]** | **−0.02 [−0.17, +0.13]** | −0.02 [−0.17, +0.14] | **−0.027 [−0.182, +0.129], z −0.33, 50/50** | +0.03 [−0.18, +0.23] |
| bill length (7,697 / 7,540) | −0.00 [−0.11, +0.10] | +0.16 [−0.02, +0.34] | +0.17 [−0.02, +0.36] | **+0.171 [−0.018, +0.360], z 1.77, 50/50** | +0.14 [−0.13, +0.41] |
| tail length (8,872 / 8,679) | +0.28 [−0.27, +0.83] | +0.30 [−0.19, +0.80] | +0.25 [−0.28, +0.77] | **+0.243 [−0.278, +0.764], z 0.91, 50/50** | +0.24 [−0.21, +0.70] |
| tarsus length (4,361 / 4,285) | **+0.42 [0.26, 0.59]** | +0.26 [−0.06, +0.58] | +0.26 [−0.08, +0.59] | **+0.260 [−0.080, +0.599], z 1.50, 50/50** | +0.26 [−0.16, +0.68] |
| wing (8,478; cross-ref, R2 is authoritative) | −0.92 [−1.40, −0.44] | −0.53 [−0.92, −0.13] | −0.30 [−0.66, +0.06] | **−0.371 [−0.777, +0.036], z −1.79, 50/50** | +0.02 [−0.58, +0.62] |

Key points:
- **Every baseline trend that excluded zero (wing, bill width, tarsus) has a
  zero-including interval once contributor and municipality intercepts are added;
  every within-contributor slope includes zero.** Body mass is flat in every
  specification (diurnal gain **+0.404 % per hour [0.291, 0.516], t 7.02**, n 8,345;
  year slope with hour −0.0063, t −0.99).
- **The bill-width increase (§1 above, +0.24 [0.06, 0.43]) is a contributor-composition
  effect**: only **2 of 22** bill-width contributors span both periods (A. Piratelli,
  M. Alves; 620 of 3,209 records); with contributor alone the sign flips (−0.24,
  t −3.4). The Allen's-rule reading is not supported and is to be retracted.
- **Phylogeny changes nothing in this table**: the 50-tree phylogenetic tier
  (glmmTMB `propto`) reproduces every M3 lme4 estimate above to within one pooled
  SE (median |Δβ|/SE = 0.011, max 0.106, across 28 fixed-effect comparisons in
  `multitrait_phylo_results.rds$comparison_lme4`), with **0 of 28 significance
  verdicts changed**, and reproduces the diurnal mass gain (**+0.404 %/hr
  [0.292, 0.517], z 7.03**, phylogenetic tier, N 8,345, 50/50 converged; year
  with hour −0.0063, z −0.99). Phylogenetic signal in species-level means is
  substantial for wing (0.957, see R2/`GLMMTMB_ENGINE.md`) and body mass (0.949–0.950)
  but weaker for the other traits (bill width 0.52–0.54, bill length 0.59–0.63,
  tail 0.91–0.92, tarsus 0.81) — none of this changes the year-effect conclusions.
- **Isometry on shared records** (7,577 records with both traits, 72 species;
  lme4 tier `bivariate_fast_lme4.rds`; **phylogenetic tier `bivariate_phylo_results.rds`,
  glmmTMB, 50 trees, complete**): with contributor + site intercepts the wing
  slope excludes zero (lme4 −0.55 [−0.97, −0.13]; phylogenetic contrasts_independent
  estimator −1.53 %/decade [−2.64, −0.34], 50/50 converged) and the mass slope does
  not (lme4 −0.0027 [−0.0135, +0.0081]; phylogenetic −0.54 %/decade [−2.68, +1.60]) —
  the **difference of the two year slopes is −0.92 %/decade [−3.27, +1.43]** (lme4-tier
  joint glmmTMB, z −0.77) and **−0.99 %/decade [−3.42, +1.45], z −0.79** (phylogenetic
  tier, joint model with two trait-specific `propto` terms — the covariance-correct
  estimate; wing-mass slope correlation only 0.04–0.05, so it barely differs from the
  independence assumption). The isometry contrast (mass − 3·wing) is +3.90 [−0.07,
  +7.88] (lme4-tier, z 1.93) and **+4.04 %/decade [−0.06, +8.15], z 1.93** (phylogenetic
  tier) — **phylogeny does not change this conclusion**. The §1 statement "wings
  shorter at constant mass" is a consistent direction, **not a demonstrated departure
  from isometry**, with or without phylogeny; the mass interval is ~2.5× the wing
  interval (residual SD 3.7 % vs 2.2 %).

### R4. Variance (`variance_results.rds`; glmmTMB `dispformula`, metafor)

| Analysis | k / N | estimate | 95 % CI |
|---|---:|---:|---|
| lnCVR pooled per species (§2 definition, reproduced) | k 54 | +0.184 | [+0.034, +0.334], I² 93.6 % |
| lnCVR within sex (species × sex cells, n ≥ 5) | k 71 | +0.187 | [+0.048, +0.325] |
| **lnCVR within species × sex × contributor (n ≥ 5)** | k 15 | **−0.249** | [−0.515, +0.017] |
| lnCVR E. Carrano only (n ≥ 5) | k 7 | −0.429 | [−0.698, −0.159] |
| σ model S0: `sigma ~ year` (no controls) | N 8,282 | **+9.6 %/decade** | [+6.3, +13.1] |
| **σ model S3: + (1\|contributor) + (1\|site) in σ, M3 mean** | N 8,282 | **−6.3 %/decade** | [−13.5, +1.4] |
| σ model S6: S3 + (1\|species) in σ | N 8,282 | +3.4 %/decade | [−4.7, +12.2] |
| σ model S7: S3 + record temperature in σ (N 7,651) | | year −6.8 [−14.0, +1.0]; **tmean −0.7 % per SD [−8.7, +8.1]** | |

Key points:
- The CV rise is real as a description of the pooled data and is not a sex-ratio
  artefact (holds within sex), but **it vanishes, and turns negative, within the same
  measurer**, and the continuous σ model gives no consistent sign once contributor
  and site intercepts enter σ (−6 % to +3 %/decade, all intervals spanning zero).
- The mechanism the referee predicted is visible: contributors per species-period
  cell **3.7 → 5.5** (45 of 61 species gain contributors, 7 lose; Wilcoxon p 5e-8);
  the between-contributor share of within-cell variance rises from a median of
  0.25 to 0.45. Residual SD differs 0.45- to 2.5-fold between contributors.
- Temperature has **no effect on residual variance** once contributor and site are
  in σ (S7), so the §2 "warmer → less variable" drmSEM path (−0.193, no
  provenance terms) should not be carried forward as a headline; the Phase 3 brms
  distributional model on Totoro is the primary variance evidence.
- Side effect for the mean: the heteroscedastic S3 fit shrinks the controlled mean
  year slope on the same records from −0.32 [−0.69, +0.04] to −0.15 [−0.35, +0.05].

**Phylogenetic tier (glmmTMB `propto`, 50 trees, `variance_phylo_results.rds`, 16:12; 9 σ tiers, 50/50 trees converged except S0 45/50).** log-residual-SD year terms: S0 +0.046 [+0.031, +0.062] (+9.6 %/decade in residual SD); S2 (+contributor +site in mean and σ) +0.013 [−0.015, +0.040] (+2.5 % [−2.9, +8.2]); S3 full controls −0.033 [−0.073, +0.007] (−6.3 % [−13.4, +1.5]); S4 within-contributor −0.031 [−0.073, +0.010], between −0.053 [−0.218, +0.112]; S6 (species intercepts in σ) +0.017 [−0.024, +0.058]; S7 (S3t + temperature) year −0.035 [−0.076, +0.005], temperature −0.007 [−0.092, +0.078]. The phylogenetic and non-phylogenetic estimates differ by ≤ 0.0006 on every term (`$comparison`). Reading: the baseline rise in among-individual variance is a between-contributor composition effect; within contributor and site it is absent or reversed, and temperature has no detectable effect on residual variance.

### R5. ATLANTIC ANTS litter-ant richness index (`ants_results.rds`; Phase 4a, full run)

Region-matched replacement context for the trophic axis. 62,020 standardised
records → **855 campaign events** (393 Winkler / 462 pitfall; 61 contributor
files, 293 localities, every year 1994–2018). glmmTMB negative-binomial richness
per campaign with contributor + locality + site intercepts, Method × log(effort),
habitat; year on the bird scale.

| Model | N | β per SD-yr [95 % CI] | % per decade |
|---|---:|---|---|
| naive | 855 | −0.269 [−0.330, −0.208] | −41.5 |
| **pooled with provenance controls (headline)** | 855 | **−0.113 [−0.187, −0.039]** | **−20.1 [−31.1, −7.4]** |
| effort as offset | 855 | +0.041 [−0.051, +0.133] | +8.5 (sign flips) |
| within-contributor (Mundlak) | 855 | −0.126 [−0.203, −0.050] | −22.2 |
| leave PERD out (within) | 802 | −0.055 (SE 0.047) | −10.5, CI includes 0 |
| abundance (67 campaigns, 45 PERD; within) | 67 | −0.282 [−0.724, +0.161] | null |

Key points:
- With controls, litter-ant richness per standardised campaign declined
  **~20 %/decade (CI 7–31 %)**, half the naive slope — the contributor controls
  matter in the same direction as for the birds.
- **Not robust**: the offset effort model flips the sign (richness scales as
  effort^0.3–0.4, so the offset is the wrong form, but the sensitivity stands);
  about half the within-contributor signal comes from one programme (PERD) whose
  protocol changed within the series; the five ≥ 6-year programmes disagree in sign
  (+0.10, −0.39, −0.03, −0.47, +0.28); abundance is null. Effort and site fields
  mean different things in different contributor files (trap IDs, per-record codes,
  campaign totals).
- Wording rule: **"litter-ant richness per standardised sample", never prey biomass
  or availability**; context, not evidence of declining prey. The drmSEM arthropod
  node is dropped rather than rebuilt from it.

### R6. PREDICTS trend retired (Phase 4e)

The §3 arthropod result ("×0.25 of baseline", n = 63,150) cannot support a
temporal claim and is **retired**: `dat_no_grassland.rds` holds 56 studies of which
**53 are single-calendar-year** snapshots (three span two adjacent years); only
**17 studies / 6,270 rows** fall in Atlantic Forest ecoregions (1998–2009); there is
no within-study temporal information to model. The "shallower without ants"
sentence was a logic error (no-ants β −1.22 [−1.29, −1.15] vs −1.40 [−1.46, −1.34]
with ants: ants steepened the pooled slope). Fitted objects and data are archived
under `../../archive/predicts/` (README there); `arthropod_estimates.rds` and
`fig-arthropods.png` follow once `index.qmd` stops reading them. The drmSEM
(`atlantic_drmsem.R` v4) now refuses `INCLUDE_ARTHRO = TRUE` without an explicit
override, adds a contributor intercept on the wing node, and is supplementary /
exploratory with an authorship disclosure; its §3 numbers (2026-06-15 fit) are
superseded and await the v4 Totoro run.

### R7. Also from this session (bearing on §3–§4)

- **Climate at the sampled localities** (`climate_trends.rds`, 357 localities,
  year random effect): annual mean T **+0.246 °C/decade [+0.076, +0.415]**
  (+0.57 °C over 1995–2018); warm-quarter Tmax +0.235 [−0.014, +0.484]; SPEI-12
  −0.256 [−0.617, +0.105]; precipitation −44.7 mm/decade [−123.1, +33.7]. Warming
  supported, **no drying trend**. Record-level r(Year, temperature) = 0.30, mostly
  between-locality (0.32) rather than within (0.13): later sampling happened at
  warmer places.
- **Diet interaction** (§4; `diet_lme4_results.rds`): lme4 reproduces the published
  +0.66 (+0.665, t 2.93); family terms do not change it; **contributor + municipality
  intercepts halve it to +0.321 [−0.015, +0.658], t 1.87**. The substantive point
  (insectivores show no steeper decline; slope at Diet-Inv 100 ≈ 0 in every
  specification) survives; the size of the frugivore/omnivore contrast does not.
  **Phylogenetic tier (`diet_interaction_phylo.rds`, glmmTMB `propto`, 50 trees,
  complete):** Model A interaction (year × Diet-Inv) **+0.663 [+0.217, +1.108],
  z 2.92, 50/50 converged** — reproduces both the lme4 tier (+0.665) and the
  published `brms` fit (+0.660 [+0.215, +1.106]) to 3 decimals. With contributor +
  municipality intercepts (A_src_site): **+0.319 [−0.019, +0.657], z 1.85, 49/50
  converged** — matches the lme4 tier's +0.321 [−0.015, +0.658] and **confirms the
  halve-and-cross-zero result under phylogeny**. Categorical Model B interaction
  +1.024 [+0.056, +1.991] (published `brms` +1.025 [+0.046, +2.004]); with
  contributor + municipality +0.387 [−0.335, +1.109], z 1.05 (also crosses zero).
  Year slope at Diet-Inv p10/p50/p90 with provenance controls, phylogenetic tier:
  −1.82 / −1.07 / −0.14 mm/decade — the obligate-insectivore (p90) slope is
  indistinguishable from zero in every specification, phylogenetic or not.
  **Phylogeny changes nothing here.**

### Provenance (this section)

`audit_{sources,sites,individuals,season,traits}.rds`, `effect_scale.rds` (P0);
`controlled_wing_results.rds`, `controlled_wing_species_slopes.rds` (P1 — regenerated
12:35:53 after the H1 fix, `Sex` included; supersedes the 12:02 files); `multitrait_results.rds`, `multitrait_table.md`,
`bivariate_fast_lme4.rds` (P2); `variance_results.rds` (P3, mode `fast`);
`ants_results.{rds,md}`, `ants_annual_index.rds` (P4a); `../../archive/predicts/README.md`
(P4e); `climate_trends.rds`, `diet_lme4_results.rds`, `diet_distribution.rds` (P5);
`figure_data/*.csv` (P6). Smoke-test files (`bivariate_smoke_results.rds`,
`models/bivariate_smoke.rda`) are wiring checks and carry no results. Independent
re-run of every fast path: `../scripts/REVISION_REVIEW.md`.

Phylogenetic tier (added same day, afternoon; no separate reviewer re-run — see
`REVISION_STATUS.md` §3.9): `glmmtmb_validation_wing.rds` (engine validation);
`multitrait_phylo_results.rds`, `bivariate_phylo_results.rds`,
`diet_interaction_phylo.rds` (**complete, 50 trees**); `controlled_wing_phylo_results.rds`,
`controlled_wing_phylo_species_slopes.rds` (**complete, 50 trees**, 16:48), `variance_phylo_results.rds`
(**complete, 50 trees**, 16:12); methods note
`../scripts/GLMMTMB_ENGINE.md`.
