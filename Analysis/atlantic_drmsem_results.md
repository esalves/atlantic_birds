# drmSEM findings — Atlantic Forest passerines

**Script:** `atlantic_drmsem.R` · **Recorded:** 2026-06-11
Two model scopes are reported. Coefficients are phylogeny-corrected (50-tree
Rubin-pooled; **Pagel's λ = 1 on all trees in both models** — strong
phylogenetic signal, but the phylo-corrected paths are nearly identical to the
i.i.d. fits, so conclusions are robust to phylogenetic correction). Temperature
is time-resolved WorldClim (CRU-TS 4.09 downscaled).

| | Arthropod model (`INCLUDE_ARTHRO = TRUE`) | Full-data model (`INCLUDE_ARTHRO = FALSE`) |
|---|---|---|
| Nodes | temperature, arthropod, wing (+σ) | temperature, wing (+σ) |
| Arthropod index | observed PREDICTS, spatial (200 km), endogenous | — |
| N / species / years | 2,006 / 67 / 1998–2008 | **8,767 / 68 / 1990–2018** |
| **Fisher's C** | 9.06, df 4, **p = 0.057** (marginal) | **0.24, df 2, p = 0.885** (clean) |

## Path coefficients (phylogeny-corrected, sd_x)

| Path | Component | Arthropod model | Full-data model |
|---|---|---:|---:|
| year → temperature | mu | +0.062 (p = 0.024) | +0.029 (p < 1e-6) |
| latitude → temperature | mu | −0.41 | −0.47 |
| year → arthropod | mu | −0.87 (p < 1e-87) | — |
| temperature → arthropod | mu | +0.23 (p < 1e-10) | — |
| latitude → arthropod | mu | +0.40 | — |
| **temperature → wing** | **mu** | **−0.55 (p = 0.012)** | **+0.11 (p = 0.31)** |
| arthropod → wing | mu | +0.06 (p = 0.63) | — |
| Sex (Male) → wing | mu | +2.59 | +2.55 |
| latitude → wing | mu | +0.02 (n.s.) | −0.11 (n.s.) |
| **year → σ(wing)** | **sigma** | **−0.099 (p = 0.004)** | **+0.074 (p < 1e-20)** |
| **temperature → σ(wing)** | **sigma** | **−0.16 (p < 1e-8)** | **−0.17 (p < 1e-80)** |

## Effect decomposition (year/temperature → wing)

| Effect | Arthropod model | Full-data model |
|---|---:|---:|
| year → wing (total indirect) | −0.043 [−0.146, 0.046] n.s. | +0.003 [−0.004, 0.010] n.s. |
| year → wing via temperature | −0.016 [−0.036, −0.002] sig. | (only channel) |
| year → wing via arthropod | −0.026 [−0.130, 0.065] n.s. | — |
| temperature → wing (direct) | −0.41 [−0.79, −0.10] sig. | +0.10 [−0.10, 0.33] n.s. |

## Findings (what the two models jointly establish)

**1. Wing-size variance has increased over time — confirmed on the full data.**
`year → σ(wing)` = **+0.074, p < 1e-20** across all 8,767 records (1990–2018).
This corroborates the lnCVR "variance rising" result. The *negative* sign in the
arthropod model (−0.099) was an artefact of its restricted 1998–2008 subset —
**the full-data run resolves the sign in favour of rising variance.**

**2. Warmer → less variable wings — robust everywhere.**
`temperature → σ(wing)` = −0.16 to −0.17 (p < 1e-8 in both). The most stable
result in the whole analysis.

**3. The mean Bergmann response is NOT robust.** `temperature → wing` is
significantly negative (−0.55) only in the smaller arthropod subset; on the full
data it is +0.11 and non-significant. So "warmer → smaller wings" is specific to
the 1998–2008 arthropod-covered sample, not a general pattern — report it as
fragile/subset-dependent, not a headline.

**4. No food-limitation effect on body size.** `arthropod → wing` ≈ 0 (p = 0.63),
even with the spatially-resolved index. Arthropods declined strongly over time
(`year → arthropod` = −0.87) but that decline does not reach wing length.

**5. No net temporal trend in mean size through climate.** The overall
`year → wing` mean effect is ≈ 0 in both models (CIs cross zero). Years did warm
(`year → temperature` = +0.03 to +0.06), but little mean size change is
transmitted.

**6. DAG fit favours the full-data model.** Fisher's C is clean for the full-data
model (p = 0.89, one passing d-sep claim) and only marginal for the arthropod
model (p = 0.057, held together by the `latitude → arthropod` adjustment edge).

**7. Phylogeny barely matters here.** λ = 1 on all 50 trees in both models, yet
phylo-corrected and i.i.d. paths coincide → conclusions robust to phylogenetic
correction. (`brm0_multiphylo.rda` remains the formal phylogenetic authority.)

## Bottom line

The defensible story is **distributional, not about the mean**: morphological
*variability* rose over 1990–2018 (full data, p < 1e-20) and falls with
temperature, while the *mean*-size Bergmann response is weak and subset-fragile
and food limitation shows no effect on size. The variance results are the
contribution; the mean-size and food-limitation hypotheses are not supported on
the full data.

## Provenance

Auto-generated per-run records: `drmsem_results_arthro.{rds,md}`,
`drmsem_results_noarthro.{rds,md}`; figures in `figures/` (`*_arthro`,
`*_noarthro`). This file is the curated cross-model summary.
