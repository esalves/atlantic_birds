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
