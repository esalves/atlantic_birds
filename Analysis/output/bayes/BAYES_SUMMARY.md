# Bayesian (Stan) tier: summary

Generated 2026-09-25 07:38:15. 35 fits. Posterior pooled over trees (mixture); REML = glmmTMB Rubin-pooled over the same trees.

## How to read the tables

- **P(<0)** is the posterior probability that the coefficient is negative: the share of posterior draws,
  pooled over all trees, that fall below zero. P(<0) = 0.96 means that, given the model, the priors and the
  data, there is a 96 % probability that the effect is negative (e.g. that wing length declined over time).
  The probability of a positive effect is 1 - P(<0). Values near 0.5 mean no evidence for either direction;
  values near 1 are strong evidence for a negative effect, values near 0 strong evidence for a positive one.
  It is not a p-value: it is a direct statement about the effect, not about data under a null hypothesis.
  For a roughly symmetric posterior, the 95 % CrI excludes zero exactly when P(<0) > 0.975 or < 0.025.
- **shift (post. SD)**: (posterior mean - REML estimate) / posterior SD. |shift| < 0.1 means the two engines agree.
- **SD ratio**: posterior SD / REML SE. Slightly above 1 is expected, because the posterior also carries the
  uncertainty in the variance components that REML plugs in as known.
- Year coefficients are per SD of calendar year (`scaled_yr`), on the response scale (mm for wing and bill,
  ln units for mass and the isometry contrast ln M - 3 ln L).

## Models

| model | code in source script | fit id | n | trees |
|---|---|---|---|---|
| Wing x diet: baseline | diet A_baseline | `atlantic_diet_interaction__002__conc.wing.length` | 8086 | 20 |
| Wing x diet: + contributor & municipality | diet A_src_site | `atlantic_diet_interaction__006__conc.wing.length` | 8086 | 20 |
| Wing x diet category: + contributor & municipality | diet B_src_site | `atlantic_diet_interaction__012__conc.wing.length` | 8086 | 20 |
| Body mass (ln): baseline | multitrait M0_baseline | `atlantic_multitrait__001__y` | 11256 | 20 |
| Body mass (ln): fully adjusted | multitrait M3_ring | `atlantic_multitrait__003__y` | 11077 | 20 |
| Body mass (ln): year within vs between contributors | multitrait MWsrc | `atlantic_multitrait__004__y` | 11077 | 20 |
| Bill width: baseline | multitrait M0_baseline | `atlantic_multitrait__005__y` | 3209 | 20 |
| Bill width: + contributor & municipality | multitrait M1_src_site | `atlantic_multitrait__006__y` | 3209 | 20 |
| Bill width: fully adjusted | multitrait M3_ring | `atlantic_multitrait__007__y` | 3205 | 20 |
| Bill width: year within vs between contributors | multitrait MWsrc | `atlantic_multitrait__008__y` | 3205 | 20 |
| Body mass (ln): capture-time sample | multitrait HOUR0 | `atlantic_multitrait__021__y` | 8345 | 20 |
| Body mass (ln): capture-time sample + capture hour | multitrait HOUR1 | `atlantic_multitrait__022__y` | 8345 | 20 |
| Wing: baseline | M0 | `atlantic_parallel_controlled__001__conc.wing.length__M0` | 8478 | 20 |
| Wing: fully adjusted | M3 | `atlantic_parallel_controlled__011__conc.wing.length__M3` | 8478 | 20 |
| Wing: fully adjusted [correlated phylo intercept-slope] | M3 | `atlantic_parallel_controlled__011__conc.wing.length__M3__pcor` | 8478 | 20 |
| Wing: fully adjusted + per-contributor trends | M4 | `atlantic_parallel_controlled__012__conc.wing.length__M4` | 8478 | 20 |
| Wing: year within vs between contributors | M5src | `atlantic_parallel_controlled__015__conc.wing.length__M5src` | 8478 | 20 |
| Wing: fully adjusted, first captures only | M7 | `atlantic_parallel_controlled__017__conc.wing.length__M7` | 7791 | 20 |
| Wing: fully adjusted, complete cases | M3_cc | `atlantic_parallel_controlled__023__conc.wing.length__M3_cc` | 8282 | 20 |
| Wing: year within vs between contributors, complete cases | M5src_cc | `atlantic_parallel_controlled__027__conc.wing.length__M5src_cc` | 8282 | 20 |
| Wing: fully adjusted, unknown-sex birds | M6 | `atlantic_parallel_controlled__038__conc.wing.length__M6` | 3532 | 20 |
| Wing (mean + SD model): baseline | S1 | `atlantic_variance_sigma__002__conc.wing.length__S1` | 8282 | 20 |
| Wing (mean + SD model): + contributor & municipality | S3 | `atlantic_variance_sigma__004__conc.wing.length__S3` | 8282 | 20 |
| Wing (mean + SD model): + local temperature | S7 | `atlantic_variance_sigma__009__conc.wing.length__S7` | 7651 | 20 |
| Wing vs temperature anomaly [+ municipality] | E3 | `referee_reruns__001__wing__site` | 7070 | 20 |
| Body mass (ln) vs temperature anomaly [+ municipality] | E3 | `referee_reruns__003__lnmass__site` | 7070 | 20 |
| Isometry vs temperature anomaly [+ municipality] | E3 | `referee_reruns__005__iso__site` | 7070 | 20 |
| Isometry vs temperature anomaly + year [+ municipality] | E3b | `referee_reruns__006__iso__site` | 7070 | 20 |
| Isometry vs year: parsimonious adjustment [+ municipality] | E4a | `referee_reruns__009__iso__site` | 7070 | 20 |
| Isometry vs year: fully adjusted [+ municipality] | E4b | `referee_reruns__010__iso__site` | 7070 | 20 |
| Isometry vs year: fully adjusted [+ municipality, correlated phylo intercept-slope] | E4b | `referee_reruns__010__iso__site__pcor` | 7070 | 20 |
| Wing (ln): fully adjusted, shared records [+ municipality] | E4c | `referee_reruns__011__lnwing__site` | 7070 | 20 |
| Body mass (ln): fully adjusted, shared records [+ municipality] | E4c | `referee_reruns__012__lnmass__site` | 7070 | 20 |
| Isometry vs detrended temperature anomaly [+ municipality] | E1b (phylo, M3) | `referee_reruns__013__iso__site` | 7070 | 20 |
| Isometry vs temperature anomaly + locality-year [+ municipality] | E2a (phylo, M3) | `referee_reruns__014__iso__site` | 7070 | 20 |

## Sampler diagnostics (worst tree per fit)

| model | trees | max R-hat | min bulk-ESS | divergences | median s/tree |
|---|---|---|---|---|---|
| Wing x diet: baseline | 20 | 1.021 | 233 | 0 | 997 |
| Wing x diet: + contributor & municipality | 20 | 1.017 | 312 | 0 | 1442 |
| Wing x diet category: + contributor & municipality | 20 | 1.013 | 247 | 0 | 2908 |
| Body mass (ln): baseline | 20 | 1.016 | 288 | 0 | 3549 |
| Body mass (ln): fully adjusted | 20 | 1.028 | 283 | 0 | 18411 |
| Body mass (ln): year within vs between contributors | 20 | 1.024 | 274 | 0 | 18639 |
| Bill width: baseline | 20 | 1.007 | 752 | 0 | 307 |
| Bill width: + contributor & municipality | 20 | 1.018 | 476 | 0 | 543 |
| Bill width: fully adjusted | 20 | 1.016 | 447 | 0 | 879 |
| Bill width: year within vs between contributors | 20 | 1.017 | 377 | 0 | 755 |
| Body mass (ln): capture-time sample | 20 | 1.022 | 271 | 0 | 3033 |
| Body mass (ln): capture-time sample + capture hour | 20 | 1.031 | 182 | 0 | 3413 |
| Wing: baseline | 20 | 1.035 | 159 | 0 | 988 |
| Wing: fully adjusted | 20 | 1.024 | 307 | 0 | 13135 |
| Wing: fully adjusted [correlated phylo intercept-slope] | 20 | 1.014 | 355 | 0 | 4050 |
| Wing: fully adjusted + per-contributor trends | 20 | 1.032 | 225 | 0 | 15879 |
| Wing: year within vs between contributors | 20 | 1.026 | 249 | 0 | 13866 |
| Wing: fully adjusted, first captures only | 20 | 1.016 | 293 | 0 | 1460 |
| Wing: fully adjusted, complete cases | 20 | 1.018 | 305 | 0 | 11973 |
| Wing: year within vs between contributors, complete cases | 20 | 1.024 | 254 | 0 | 11493 |
| Wing: fully adjusted, unknown-sex birds | 20 | 1.013 | 304 | 0 | 650 |
| Wing (mean + SD model): baseline | 20 | 1.024 | 270 | 0 | 1043 |
| Wing (mean + SD model): + contributor & municipality | 20 | 1.030 | 164 | 0 | 3896 |
| Wing (mean + SD model): + local temperature | 20 | 1.021 | 352 | 0 | 5836 |
| Wing vs temperature anomaly [+ municipality] | 20 | 1.015 | 284 | 0 | 8787 |
| Body mass (ln) vs temperature anomaly [+ municipality] | 20 | 1.019 | 273 | 0 | 8626 |
| Isometry vs temperature anomaly [+ municipality] | 20 | 1.028 | 155 | 0 | 1639 |
| Isometry vs temperature anomaly + year [+ municipality] | 20 | 1.031 | 174 | 0 | 3467 |
| Isometry vs year: parsimonious adjustment [+ municipality] | 20 | 1.053 | 81 | 0 | 959 |
| Isometry vs year: fully adjusted [+ municipality] | 20 | 1.032 | 154 | 0 | 1630 |
| Isometry vs year: fully adjusted [+ municipality, correlated phylo intercept-slope] | 20 | 1.046 | 107 | 0 | 1556 |
| Wing (ln): fully adjusted, shared records [+ municipality] | 20 | 1.016 | 342 | 0 | 1782 |
| Body mass (ln): fully adjusted, shared records [+ municipality] | 20 | 1.021 | 298 | 0 | 1741 |
| Isometry vs detrended temperature anomaly [+ municipality] | 20 | 1.057 | 57 | 0 | 3591 |
| Isometry vs temperature anomaly + locality-year [+ municipality] | 20 | 1.054 | 72 | 0 | 3866 |

## Time terms: posterior vs REML

| model | term | posterior mean [95% CrI] | P(<0) | REML [95% CI] | shift (post. SD) | SD ratio |
|---|---|---|---|---|---|---|
| Wing x diet: baseline | year | -0.860 [-1.337, -0.383] | 1.000 | -0.858 [-1.327, -0.389] | -0.01 | 1.02 |
| Wing x diet: baseline | year x diet (invertebrate share) | 0.661 [0.205, 1.120] | 0.003 | 0.663 [0.218, 1.108] | -0.01 | 1.03 |
| Wing x diet: + contributor & municipality | year | -0.465 [-0.875, -0.050] | 0.986 | -0.465 [-0.865, -0.066] | 0.00 | 1.03 |
| Wing x diet: + contributor & municipality | year x diet (invertebrate share) | 0.319 [-0.027, 0.670] | 0.036 | 0.321 [-0.017, 0.658] | -0.01 | 1.03 |
| Wing x diet category: + contributor & municipality | year | -0.663 [-1.199, -0.126] | 0.992 | -0.668 [-1.191, -0.146] | 0.02 | 1.03 |
| Wing x diet category: + contributor & municipality | b_scaled_yr:diet_cat2Invertebrate | 0.384 [-0.365, 1.135] | 0.154 | 0.396 [-0.327, 1.119] | -0.03 | 1.03 |
| Body mass (ln): baseline | year | -0.007 [-0.017, 0.002] | 0.934 | -0.007 [-0.017, 0.002] | -0.01 | 1.02 |
| Body mass (ln): fully adjusted | year | -0.005 [-0.014, 0.005] | 0.830 | -0.005 [-0.014, 0.005] | -0.01 | 1.02 |
| Body mass (ln): year within vs between contributors | year within contributor | -0.003 [-0.010, 0.004] | 0.798 | -0.003 [-0.010, 0.004] | 0.02 | 1.02 |
| Body mass (ln): year within vs between contributors | contributor mean year | 0.006 [-0.033, 0.044] | 0.374 | 0.006 [-0.030, 0.042] | -0.00 | 1.04 |
| Bill width: baseline | year | 0.245 [0.061, 0.433] | 0.005 | 0.244 [0.063, 0.425] | 0.01 | 1.03 |
| Bill width: + contributor & municipality | year | -0.028 [-0.185, 0.129] | 0.640 | -0.031 [-0.184, 0.123] | 0.04 | 1.02 |
| Bill width: fully adjusted | year | -0.024 [-0.181, 0.135] | 0.618 | -0.027 [-0.183, 0.129] | 0.04 | 1.01 |
| Bill width: year within vs between contributors | year within contributor | 0.020 [-0.190, 0.236] | 0.427 | 0.016 [-0.190, 0.223] | 0.04 | 1.03 |
| Bill width: year within vs between contributors | contributor mean year | 0.407 [-0.541, 1.362] | 0.192 | 0.408 [-0.482, 1.298] | -0.00 | 1.06 |
| Body mass (ln): capture-time sample | year | -0.008 [-0.020, 0.005] | 0.873 | -0.008 [-0.020, 0.005] | 0.00 | 1.02 |
| Body mass (ln): capture-time sample + capture hour | year | -0.006 [-0.019, 0.006] | 0.835 | -0.006 [-0.019, 0.006] | 0.01 | 1.01 |
| Body mass (ln): capture-time sample + capture hour | capture hour | 0.004 [0.003, 0.005] | 0.000 | 0.004 [0.003, 0.005] | -0.01 | 1.00 |
| Wing: baseline | year | -0.928 [-1.426, -0.432] | 1.000 | -0.922 [-1.400, -0.444] | -0.02 | 1.04 |
| Wing: fully adjusted | year | -0.369 [-0.780, 0.048] | 0.959 | -0.371 [-0.777, 0.035] | 0.01 | 1.02 |
| Wing: fully adjusted [correlated phylo intercept-slope] | year | -0.268 [-1.148, 0.679] | 0.763 | -0.371 [-0.777, 0.035] | 0.23 | 2.17 |
| Wing: fully adjusted + per-contributor trends | year | -0.874 [-1.809, 0.040] | 0.970 | -0.875 [-1.758, 0.007] | 0.00 | 1.04 |
| Wing: year within vs between contributors | year within contributor | -0.164 [-0.807, 0.491] | 0.692 | -0.159 [-0.798, 0.481] | -0.02 | 1.02 |
| Wing: year within vs between contributors | contributor mean year | -1.299 [-2.486, -0.097] | 0.983 | -1.289 [-2.421, -0.158] | -0.02 | 1.04 |
| Wing: fully adjusted, first captures only | year | -0.373 [-0.788, 0.040] | 0.962 | -0.379 [-0.783, 0.024] | 0.03 | 1.02 |
| Wing: fully adjusted, complete cases | year | -0.283 [-0.656, 0.091] | 0.932 | -0.285 [-0.648, 0.077] | 0.01 | 1.03 |
| Wing: year within vs between contributors, complete cases | year within contributor | 0.027 [-0.590, 0.665] | 0.471 | 0.025 [-0.578, 0.629] | 0.00 | 1.03 |
| Wing: year within vs between contributors, complete cases | contributor mean year | -1.206 [-2.448, 0.031] | 0.973 | -1.206 [-2.377, -0.034] | 0.00 | 1.05 |
| Wing: fully adjusted, unknown-sex birds | year | -0.698 [-1.229, -0.154] | 0.993 | -0.690 [-1.206, -0.173] | -0.03 | 1.04 |
| Wing (mean + SD model): baseline | year | -0.925 [-1.395, -0.452] | 1.000 | -0.923 [-1.382, -0.464] | -0.01 | 1.02 |
| Wing (mean + SD model): + contributor & municipality | year | -0.158 [-0.370, 0.056] | 0.930 | -0.152 [-0.353, 0.050] | -0.06 | 1.05 |
| Wing (mean + SD model): + local temperature | year | -0.184 [-0.439, 0.064] | 0.927 | -0.182 [-0.421, 0.056] | -0.01 | 1.05 |
| Wing vs temperature anomaly [+ municipality] | temperature anomaly | -0.007 [-0.417, 0.392] | 0.510 | 0.010 [-0.369, 0.390] | -0.09 | 1.07 |
| Body mass (ln) vs temperature anomaly [+ municipality] | temperature anomaly | 0.005 [-0.009, 0.020] | 0.241 | 0.005 [-0.008, 0.018] | 0.05 | 1.09 |
| Isometry vs temperature anomaly [+ municipality] | temperature anomaly | 0.012 [-0.008, 0.033] | 0.117 | 0.012 [-0.008, 0.031] | 0.06 | 1.05 |
| Isometry vs temperature anomaly + year [+ municipality] | temperature anomaly | 0.009 [-0.010, 0.029] | 0.175 | 0.009 [-0.010, 0.029] | -0.01 | 1.00 |
| Isometry vs temperature anomaly + year [+ municipality] | year | 0.011 [-0.003, 0.025] | 0.059 | 0.011 [-0.002, 0.025] | 0.00 | 1.03 |
| Isometry vs year: parsimonious adjustment [+ municipality] | year | 0.012 [-0.002, 0.026] | 0.046 | 0.012 [-0.002, 0.026] | -0.01 | 1.02 |
| Isometry vs year: fully adjusted [+ municipality] | year | 0.012 [-0.002, 0.026] | 0.050 | 0.012 [-0.002, 0.025] | -0.01 | 1.02 |
| Isometry vs year: fully adjusted [+ municipality, correlated phylo intercept-slope] | year | 0.009 [-0.017, 0.034] | 0.198 | 0.012 [-0.002, 0.025] | -0.19 | 1.81 |
| Wing (ln): fully adjusted, shared records [+ municipality] | year | -0.004 [-0.009, 0.002] | 0.909 | -0.004 [-0.009, 0.002] | 0.02 | 1.03 |
| Body mass (ln): fully adjusted, shared records [+ municipality] | year | -0.001 [-0.010, 0.007] | 0.631 | -0.001 [-0.010, 0.007] | -0.01 | 1.01 |

## Residual-SD terms (mean + SD models): posterior vs REML

Coefficients are on the log(residual SD) scale; `% SD change` = 100 (exp(b) - 1), the percentage change
in within-species residual SD per SD of year (or males vs females). P(<0) here is the probability that
the residual SD shrinks (for year: that wing length became less variable over time).

| model | term | posterior mean [95% CrI] | % SD change [95% CrI] | P(<0) | REML [95% CI] | shift (post. SD) | SD ratio |
|---|---|---|---|---|---|---|---|
| Wing (mean + SD model): baseline | year | 0.047 [0.032, 0.062] | +4.8 [+3.2, +6.4] | 0.000 | 0.047 [0.031, 0.062] | -0.01 | 1.00 |
| Wing (mean + SD model): baseline | male (vs female) | 0.022 [-0.010, 0.055] | +2.3 [-1.0, +5.7] | 0.092 | 0.022 [-0.010, 0.055] | -0.01 | 1.00 |
| Wing (mean + SD model): + contributor & municipality | year | -0.034 [-0.074, 0.006] | -3.3 [-7.1, +0.6] | 0.951 | -0.033 [-0.072, 0.007] | -0.07 | 1.00 |
| Wing (mean + SD model): + contributor & municipality | male (vs female) | -0.010 [-0.045, 0.025] | -1.0 [-4.4, +2.5] | 0.719 | -0.010 [-0.045, 0.025] | -0.00 | 1.00 |
| Wing (mean + SD model): + local temperature | year | -0.036 [-0.076, 0.005] | -3.5 [-7.3, +0.5] | 0.957 | -0.035 [-0.076, 0.005] | -0.02 | 1.01 |
| Wing (mean + SD model): + local temperature | local temperature | -0.005 [-0.091, 0.081] | -0.5 [-8.7, +8.4] | 0.549 | -0.007 [-0.092, 0.078] | 0.04 | 1.02 |
| Wing (mean + SD model): + local temperature | male (vs female) | -0.027 [-0.064, 0.009] | -2.7 [-6.2, +0.9] | 0.931 | -0.027 [-0.064, 0.009] | -0.01 | 0.99 |
