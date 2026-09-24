# Bayesian (Stan) tier: summary

Generated 2026-09-24 10:17:50. 16 fits. Posterior pooled over trees (mixture); REML = glmmTMB Rubin-pooled over the same trees.

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
| Body mass (ln): fully adjusted | multitrait M3_ring | `atlantic_multitrait__003__y` | 11077 | 20 |
| Bill width: fully adjusted | multitrait M3_ring | `atlantic_multitrait__007__y` | 3205 | 20 |
| Wing: baseline | M0 | `atlantic_parallel_controlled__001__conc.wing.length__M0` | 8478 | 20 |
| Wing: fully adjusted | M3 | `atlantic_parallel_controlled__011__conc.wing.length__M3` | 8478 | 20 |
| Wing: fully adjusted [correlated phylo intercept-slope] | M3 | `atlantic_parallel_controlled__011__conc.wing.length__M3__pcor` | 8478 | 20 |
| Wing: fully adjusted + per-contributor trends | M4 | `atlantic_parallel_controlled__012__conc.wing.length__M4` | 8478 | 20 |
| Wing: year within vs between contributors | M5src | `atlantic_parallel_controlled__015__conc.wing.length__M5src` | 8478 | 20 |
| Wing: fully adjusted, first captures only | M7 | `atlantic_parallel_controlled__017__conc.wing.length__M7` | 7791 | 20 |
| Wing: fully adjusted, unknown-sex birds | M6 | `atlantic_parallel_controlled__038__conc.wing.length__M6` | 3532 | 20 |
| Wing (mean + SD model): baseline | S1 | `atlantic_variance_sigma__002__conc.wing.length__S1` | 8282 | 20 |
| Wing (mean + SD model): + contributor & municipality | S3 | `atlantic_variance_sigma__004__conc.wing.length__S3` | 8282 | 20 |
| Isometry vs temperature anomaly [+ municipality] | E3 | `referee_reruns__005__iso__site` | 7070 | 20 |
| Isometry vs year: fully adjusted [+ municipality] | E4b | `referee_reruns__010__iso__site` | 7070 | 20 |
| Isometry vs year: fully adjusted [+ municipality, correlated phylo intercept-slope] | E4b | `referee_reruns__010__iso__site__pcor` | 7070 | 20 |

## Sampler diagnostics (worst tree per fit)

| model | trees | max R-hat | min bulk-ESS | divergences | median s/tree |
|---|---|---|---|---|---|
| Wing x diet: baseline | 20 | 1.021 | 233 | 0 | 997 |
| Wing x diet: + contributor & municipality | 20 | 1.017 | 312 | 0 | 1442 |
| Body mass (ln): fully adjusted | 20 | 1.028 | 283 | 0 | 18411 |
| Bill width: fully adjusted | 20 | 1.076 | 79 | 0 | 212 |
| Wing: baseline | 20 | 1.035 | 159 | 0 | 988 |
| Wing: fully adjusted | 20 | 1.024 | 307 | 0 | 13135 |
| Wing: fully adjusted [correlated phylo intercept-slope] | 20 | 1.014 | 355 | 0 | 4050 |
| Wing: fully adjusted + per-contributor trends | 20 | 1.032 | 225 | 0 | 15879 |
| Wing: year within vs between contributors | 20 | 1.026 | 249 | 0 | 13866 |
| Wing: fully adjusted, first captures only | 20 | 1.016 | 293 | 0 | 1460 |
| Wing: fully adjusted, unknown-sex birds | 20 | 1.013 | 304 | 0 | 650 |
| Wing (mean + SD model): baseline | 20 | 1.024 | 270 | 0 | 1043 |
| Wing (mean + SD model): + contributor & municipality | 20 | 1.030 | 164 | 0 | 3896 |
| Isometry vs temperature anomaly [+ municipality] | 20 | 1.028 | 155 | 0 | 1639 |
| Isometry vs year: fully adjusted [+ municipality] | 20 | 1.032 | 154 | 0 | 1630 |
| Isometry vs year: fully adjusted [+ municipality, correlated phylo intercept-slope] | 20 | 1.046 | 107 | 0 | 1556 |

## Time terms: posterior vs REML

| model | term | posterior mean [95% CrI] | P(<0) | REML [95% CI] | shift (post. SD) | SD ratio |
|---|---|---|---|---|---|---|
| Wing x diet: baseline | year | -0.860 [-1.337, -0.383] | 1.000 | -0.858 [-1.327, -0.389] | -0.01 | 1.02 |
| Wing x diet: baseline | year x diet (invertebrate share) | 0.661 [0.205, 1.120] | 0.003 | 0.663 [0.218, 1.108] | -0.01 | 1.03 |
| Wing x diet: + contributor & municipality | year | -0.465 [-0.875, -0.050] | 0.986 | -0.465 [-0.865, -0.066] | 0.00 | 1.03 |
| Wing x diet: + contributor & municipality | year x diet (invertebrate share) | 0.319 [-0.027, 0.670] | 0.036 | 0.321 [-0.017, 0.658] | -0.01 | 1.03 |
| Body mass (ln): fully adjusted | year | -0.005 [-0.014, 0.005] | 0.830 | -0.005 [-0.014, 0.005] | -0.01 | 1.02 |
| Bill width: fully adjusted | year | -0.024 [-0.183, 0.137] | 0.615 | -0.027 [-0.183, 0.129] | 0.04 | 1.02 |
| Wing: baseline | year | -0.928 [-1.426, -0.432] | 1.000 | -0.922 [-1.400, -0.444] | -0.02 | 1.04 |
| Wing: fully adjusted | year | -0.369 [-0.780, 0.048] | 0.959 | -0.371 [-0.777, 0.035] | 0.01 | 1.02 |
| Wing: fully adjusted [correlated phylo intercept-slope] | year | -0.268 [-1.148, 0.679] | 0.763 | -0.371 [-0.777, 0.035] | 0.23 | 2.17 |
| Wing: fully adjusted + per-contributor trends | year | -0.874 [-1.809, 0.040] | 0.970 | -0.875 [-1.758, 0.007] | 0.00 | 1.04 |
| Wing: year within vs between contributors | year within contributor | -0.164 [-0.807, 0.491] | 0.692 | -0.159 [-0.798, 0.481] | -0.02 | 1.02 |
| Wing: year within vs between contributors | contributor mean year | -1.299 [-2.486, -0.097] | 0.983 | -1.289 [-2.421, -0.158] | -0.02 | 1.04 |
| Wing: fully adjusted, first captures only | year | -0.373 [-0.788, 0.040] | 0.962 | -0.379 [-0.783, 0.024] | 0.03 | 1.02 |
| Wing: fully adjusted, unknown-sex birds | year | -0.698 [-1.229, -0.154] | 0.993 | -0.690 [-1.206, -0.173] | -0.03 | 1.04 |
| Wing (mean + SD model): baseline | year | -0.925 [-1.395, -0.452] | 1.000 | -0.923 [-1.382, -0.464] | -0.01 | 1.02 |
| Wing (mean + SD model): + contributor & municipality | year | -0.158 [-0.370, 0.056] | 0.930 | -0.152 [-0.353, 0.050] | -0.06 | 1.05 |
| Isometry vs temperature anomaly [+ municipality] | temperature anomaly | 0.012 [-0.008, 0.033] | 0.117 | 0.012 [-0.008, 0.031] | 0.06 | 1.05 |
| Isometry vs year: fully adjusted [+ municipality] | year | 0.012 [-0.002, 0.026] | 0.050 | 0.012 [-0.002, 0.025] | -0.01 | 1.02 |
| Isometry vs year: fully adjusted [+ municipality, correlated phylo intercept-slope] | year | 0.009 [-0.017, 0.034] | 0.198 | 0.012 [-0.002, 0.025] | -0.19 | 1.81 |

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
