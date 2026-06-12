# drmSEM results — v3-noarthro
_generated 2026-06-11 18:59:44.345357 · N = 8767 records, 68 species, years 1990–2018_

## Fisher's C (global DAG fit)
```
   fisher_c df n_claims   p.value
1 0.2435746  2        1 0.8853367
```

## d-separation claims
```
                                            claim   x            y                 given df         LR
1 Sex _||_ scaled_tmean | {scaled_yr, scaled_lat} Sex scaled_tmean scaled_yr, scaled_lat  1 0.02079579
    p.value status
1 0.8853367     ok
```

## Path coefficients (raw)
```
          from           to component     link         term    estimate   std.error  statistic
1    scaled_yr scaled_tmean        mu identity    scaled_yr  0.02896548 0.005884933   4.921973
2   scaled_lat scaled_tmean        mu identity   scaled_lat -0.46521326 0.010290352 -45.208683
3 scaled_tmean  wing_length        mu identity scaled_tmean  0.11077474 0.108444372   1.021489
4          Sex  wing_length        mu identity      SexMale  2.55352258 0.108785712  23.472959
5   scaled_lat  wing_length        mu identity   scaled_lat -0.11390027 0.104238129  -1.092693
6    scaled_yr  wing_length     sigma      log    scaled_yr  0.07381578 0.007805659   9.456701
7 scaled_tmean  wing_length     sigma      log scaled_tmean -0.17414262 0.009107288 -19.121239
        p.value endogenous
1  8.567602e-07      FALSE
2  0.000000e+00      FALSE
3  3.070228e-01       TRUE
4 7.706204e-122      FALSE
5  2.745286e-01      FALSE
6  3.178131e-21      FALSE
7  1.680662e-81       TRUE
```

## Path coefficients (standardised, sd_x)
```
          from           to component     link         term    estimate   std.error  statistic
1    scaled_yr scaled_tmean        mu identity    scaled_yr  0.02896548 0.005884933   4.921973
2   scaled_lat scaled_tmean        mu identity   scaled_lat -0.46521326 0.010290352 -45.208683
3 scaled_tmean  wing_length        mu identity scaled_tmean  0.11077474 0.108444372   1.021489
4          Sex  wing_length        mu identity      SexMale  2.55352258 0.108785712  23.472959
5   scaled_lat  wing_length        mu identity   scaled_lat -0.11390027 0.104238129  -1.092693
6    scaled_yr  wing_length     sigma      log    scaled_yr  0.07381578 0.007805659   9.456701
7 scaled_tmean  wing_length     sigma      log scaled_tmean -0.17414262 0.009107288 -19.121239
        p.value endogenous std.estimate
1  8.567602e-07      FALSE   0.02978989
2  0.000000e+00      FALSE  -0.46375266
3  3.070228e-01       TRUE   0.10102579
4 7.706204e-122      FALSE   2.55352258
5  2.745286e-01      FALSE  -0.11354266
6  3.178131e-21      FALSE   0.07591670
7  1.680662e-81       TRUE  -0.15881686
```

## Phylogeny-corrected paths (Rubin-pooled across 50 trees)
Pagel's lambda (AIC-selected per tree), range: 1–1
```
                                               path    estimate          se       lower       upper
1         scaled_yr | scaled_tmean | mu | scaled_yr  0.02896548 0.005884933  0.01743101  0.04049995
2       scaled_lat | scaled_tmean | mu | scaled_lat -0.46521326 0.010290352 -0.48538235 -0.44504417
3    scaled_tmean | wing_length | mu | scaled_tmean  0.11456270 0.108421542 -0.09794352  0.32706892
4                  Sex | wing_length | mu | SexMale  2.55252752 0.108779627  2.33931945  2.76573559
5        scaled_lat | wing_length | mu | scaled_lat -0.11193447 0.104105038 -0.31598034  0.09211141
6       scaled_yr | wing_length | sigma | scaled_yr  0.07374582 0.007804585  0.05844883  0.08904280
7 scaled_tmean | wing_length | sigma | scaled_tmean -0.17428971 0.009106562 -0.19213857 -0.15644085
           z       p.value
1   4.921973  8.567602e-07
2 -45.208683  0.000000e+00
3   1.056641  2.906752e-01
4  23.465125 9.264843e-122
5  -1.075207  2.822821e-01
6   9.449038  3.419595e-21
7 -19.138915  1.197367e-81
```

## Effect: scaled_yr -> wing_length (all mediators)
```
      from          to      through target              quantity estimate conf.low conf.high
 scaled_yr wing_length scaled_tmean   mean            total_path   0.0034  -0.0038    0.0098
 scaled_yr wing_length scaled_tmean   mean                direct   0.0000   0.0000    0.0000
 scaled_yr wing_length scaled_tmean   mean              indirect   0.0034  -0.0038    0.0098
 scaled_yr wing_length scaled_tmean   mean         mean_mediated   0.0034  -0.0038    0.0097
 scaled_yr wing_length scaled_tmean   mean distribution_mediated   0.0000  -0.0002    0.0002
```

## Effect: scaled_tmean -> wing_length
```
         from          to through target              quantity estimate conf.low conf.high
 scaled_tmean wing_length           mean            total_path   0.0999  -0.0958    0.3258
 scaled_tmean wing_length           mean                direct   0.0999  -0.0958    0.3258
 scaled_tmean wing_length           mean              indirect   0.0000   0.0000    0.0000
 scaled_tmean wing_length           mean         mean_mediated   0.0000   0.0000    0.0000
 scaled_tmean wing_length           mean distribution_mediated   0.0000   0.0000    0.0000
```
