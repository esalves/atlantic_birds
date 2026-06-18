# log(body mass) trend — 50-tree Rubin-pooled model (mirrors atlantic_parallel.R)
_generated 2026-06-15 15:34:26.340165 · 11256 body-mass records across 73 species (sample: 12571 records / 73 species), years 1995–2018_

## Pooled fixed effects (Wald–Rubin 95% CI)
```
                      par     estimate          se        lower       upper
b_Intercept   b_Intercept  2.838309326 0.397853652  2.058516168 3.618102485
b_SexMale       b_SexMale  0.013529575 0.002849060  0.007945417 0.019113733
b_scaled_yr   b_scaled_yr -0.007293741 0.004864934 -0.016829011 0.002241529
b_scaled_lat b_scaled_lat  0.004846753 0.002950338 -0.000935909 0.010629414
```

Phylogenetic signal (proportion of among-species variance): 0.96 [0.95, 0.97]

Model-free quartile comparison (late vs early period): 36 species decreased, 26 increased in mean log mass.

**Year effect on log(body mass): β = -0.007 [-0.017, 0.002] per SD-year.**
