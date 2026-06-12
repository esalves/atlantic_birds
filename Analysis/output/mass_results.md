# log(body mass) trend — 50-tree Rubin-pooled model (mirrors atlantic_parallel.R)
_generated 2026-06-12 16:44:10.745277 · 13518 body-mass records across 88 species (sample: 15332 records / 89 species), years 1990–2018_

## Pooled fixed effects (Wald–Rubin 95% CI)
```
                      par     estimate          se        lower       upper
b_Intercept   b_Intercept  2.891584068 0.398454517  2.110613215 3.672554922
b_SexMale       b_SexMale  0.012318612 0.002640120  0.007143978 0.017493247
b_scaled_yr   b_scaled_yr -0.001214742 0.003628972 -0.008327527 0.005898043
b_scaled_lat b_scaled_lat  0.007075987 0.002357862  0.002454578 0.011697396
```

Phylogenetic signal (proportion of among-species variance): 0.96 [0.95, 0.97]

Model-free quartile comparison (late vs early period): 42 species decreased, 38 increased in mean log mass.

**Year effect on log(body mass): β = -0.001 [-0.008, 0.006] per SD-year.**
