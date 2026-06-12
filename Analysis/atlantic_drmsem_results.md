# drmSEM findings — Atlantic Forest passerines

**Script:** `atlantic_drmsem.R` · **Model tag:** `v3-arthro-spatArthro200km`
**Recorded:** 2026-06-11 (from the saved run caches + reported Fisher's C)

## Model specification

Distributional piecewise SEM, **arthropod model** (`INCLUDE_ARTHRO = TRUE`):

- **Nodes:** `scaled_tmean ~ scaled_yr + scaled_lat`; `arthro_obs_std ~ scaled_yr + scaled_tmean + scaled_lat`; `wing_length ~ scaled_tmean + arthro_obs_std + Sex + scaled_lat` with `sigma(wing_length) ~ scaled_yr + scaled_tmean`.
- **Arthropod index:** observed PREDICTS log(abundance/effort), **spatially resolved** as a locality × year mean within a 200 km radius (Atlantic-Forest-year fallback) — varies within year, mitigating pseudoreplication.
- **Temperature:** time-resolved WorldClim (CRU-TS 4.09 downscaled), per record × year.
- **Sample:** N ≈ 2,006 individual records, 67 species (arthropod-covered window).
- **Phylogeny:** `relmat()` on the wing-length node; Pagel's λ AIC-selected per tree; coefficients Rubin-pooled across 50 clootl trees. **λ = 1.0 for all 50 trees** (strong phylogenetic signal — the phylogeny-corrected coefficients below are authoritative).

## Global fit (main SEM)

**Fisher's C = 9.18 on 4 df, p = 0.0567.** The DAG is marginally **consistent** with the data (not rejected at α = 0.05) — a clean recovery from the earlier p ≈ 1e-59 rejection, achieved by spatially resolving the arthropod index and adding the `scaled_lat → arthro_obs_std` edge. Borderline; report as "consistent but marginal."

## Phylogeny-corrected path coefficients (50-tree Rubin-pooled, λ = 1)

| Path | Component | Estimate | 95% CI | p |
|---|---|---:|---|---:|
| scaled_yr → scaled_tmean | mu | +0.062 | [0.008, 0.115] | 0.024 |
| scaled_lat → scaled_tmean | mu | −0.406 | [−0.451, −0.361] | 8e-71 |
| scaled_yr → arthro_obs_std | mu | −0.896 | [−0.981, −0.811] | 2e-94 |
| scaled_tmean → arthro_obs_std | mu | +0.262 | [0.194, 0.330] | 5e-14 |
| scaled_lat → arthro_obs_std | mu | +0.403 | [0.343, 0.463] | 2e-39 |
| **scaled_tmean → wing_length** | **mu** | **−0.557** | **[−0.989, −0.126]** | **0.011** |
| arthro_obs_std → wing_length | mu | +0.062 | [−0.166, 0.290] | 0.59 |
| Sex (Male) → wing_length | mu | +2.587 | [2.132, 3.042] | 8e-29 |
| scaled_lat → wing_length | mu | +0.019 | [−0.419, 0.457] | 0.93 |
| **scaled_yr → sigma(wing_length)** | **sigma** | **−0.099** | **[−0.166, −0.031]** | **0.004** |
| **scaled_tmean → sigma(wing_length)** | **sigma** | **−0.160** | **[−0.212, −0.108]** | 2e-9 |

(mu estimates ≈ mm per 1 SD of predictor; sigma on the log scale.)

## Effect decomposition (main SEM, simulation-based)

| Effect | Estimate | 95% CI | Note |
|---|---:|---|---|
| scaled_yr → wing (both mediators, total indirect) | −0.043 | [−0.146, 0.046] | n.s. |
| scaled_yr → wing **via temperature** | **−0.016** | **[−0.036, −0.002]** | significant |
| scaled_yr → wing via arthropod | −0.029 | [−0.130, 0.065] | n.s. |
| scaled_tmean → wing (total) | −0.396 | [−0.748, −0.070] | significant |
| scaled_tmean → wing (**direct**) | −0.409 | [−0.795, −0.099] | significant |
| scaled_tmean → wing (indirect via arthropod) | +0.013 | [−0.032, 0.053] | n.s. |

Distribution-mediated components ≈ 0 throughout.

## Findings

1. **Temperature → smaller wings (Bergmann-consistent), acting directly.** Phylogeny-corrected `scaled_tmean → wing_length` = −0.56 mm/SD (p = 0.011); the temperature effect is essentially all direct (−0.41), not routed through arthropods. With warming in-sample (`year → temperature` = +0.06), the **year → temperature → wing** channel is small but significant (−0.016, CI excludes 0) — the temporal Bergmann mediation.

2. **No food-limitation effect on body size.** `arthro_obs_std → wing_length` = +0.06 (p = 0.59); the arthropod channel of the year→size effect is non-significant. The spatial index *strengthened* this null rather than reversing it, so it is not a level-mismatch artifact.

3. **Arthropods declined steeply but it is a dead-end path.** `year → arthropod` = −0.90 (p ≈ 1e-94); arthropods are higher in warmer (+0.26) and higher-latitude (+0.40) cells. None of this reaches wing length.

4. **Overall year → size effect is weak** (−0.043, CI crosses 0); only the temperature sub-channel is reliable.

5. **Variance — one robust, one fragile.** `temperature → sigma(wing)` = −0.16 (p ≈ 2e-9; warmer → less variable) holds across all runs. `year → sigma(wing)` = −0.099 (p = 0.004) is *negative* here (variability falling over time), **opposite to the lnCVR "variance rising" result** — this is sensitive to the arthropod-covered subset and must be checked on the full dataset (`INCLUDE_ARTHRO = FALSE`).

## Caveats

- **Borderline DAG fit** (p = 0.057). The `scaled_lat → arthro_obs_std` edge is partly a real spatial gradient (given the local index) and partly a sampling adjustment — see design decision (1).
- **Restricted sample** (N ≈ 2,006, arthropod-covered window). The variance-over-time result in particular is subset-dependent.
- **Authoritative phylogeny:** the primary brms analysis (`brm0_multiphylo.rda`) remains the formal authority; these are the drmSEM corroboration.
- **Pending:** the full-data temperature/variance model (`INCLUDE_ARTHRO = FALSE`, all ~8,700 records, 1990–2018) has not yet been run; it adjudicates the `year → sigma` sign.
