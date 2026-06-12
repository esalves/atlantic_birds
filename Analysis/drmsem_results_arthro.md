# drmSEM results — v3-arthro-spatArthro200km
_generated 2026-06-11 18:40:47.039257 · N = 2006 records, 67 species, years 1998–2008_

## Fisher's C (global DAG fit)
```
  fisher_c df n_claims    p.value
1 9.056557  4        2 0.05970126
```

## d-separation claims
```
                                                            claim   x              y
1                 Sex _||_ scaled_tmean | {scaled_yr, scaled_lat} Sex   scaled_tmean
2 Sex _||_ arthro_obs_std | {scaled_yr, scaled_tmean, scaled_lat} Sex arthro_obs_std
                                given df        LR    p.value status
1               scaled_yr, scaled_lat  1 0.5111002 0.47466297     ok
2 scaled_yr, scaled_tmean, scaled_lat  1 5.1873864 0.02275141     ok
```

## Path coefficients (raw)
```
             from             to component     link           term    estimate  std.error   statistic
1       scaled_yr   scaled_tmean        mu identity      scaled_yr  0.06152567 0.02719163   2.2626693
2      scaled_lat   scaled_tmean        mu identity     scaled_lat -0.40605738 0.02282456 -17.7903716
3       scaled_yr arthro_obs_std        mu identity      scaled_yr -0.86606415 0.04355845 -19.8828033
4    scaled_tmean arthro_obs_std        mu identity   scaled_tmean  0.22888083 0.03477478   6.5818050
5      scaled_lat arthro_obs_std        mu identity     scaled_lat  0.39506212 0.03071495  12.8622089
6    scaled_tmean    wing_length        mu identity   scaled_tmean -0.57860357 0.21991026  -2.6310895
7  arthro_obs_std    wing_length        mu identity arthro_obs_std  0.06167576 0.11694391   0.5273961
8             Sex    wing_length        mu identity        SexMale  2.59143777 0.23229215  11.1559423
9      scaled_lat    wing_length        mu identity     scaled_lat  0.03137900 0.22436790   0.1398551
10      scaled_yr    wing_length     sigma      log      scaled_yr -0.09596912 0.03434174  -2.7945329
11   scaled_tmean    wing_length     sigma      log   scaled_tmean -0.16061900 0.02662236  -6.0332367
        p.value endogenous
1  2.365608e-02      FALSE
2  8.392147e-71      FALSE
3  5.733891e-88      FALSE
4  4.647708e-11       TRUE
5  7.344739e-38      FALSE
6  8.511163e-03       TRUE
7  5.979186e-01       TRUE
8  6.697941e-29      FALSE
9  8.887745e-01      FALSE
10 5.197476e-03      FALSE
11 1.607077e-09       TRUE
```

## Path coefficients (standardised, sd_x)
```
             from             to component     link           term    estimate  std.error   statistic
1       scaled_yr   scaled_tmean        mu identity      scaled_yr  0.06152567 0.02719163   2.2626693
2      scaled_lat   scaled_tmean        mu identity     scaled_lat -0.40605738 0.02282456 -17.7903716
3       scaled_yr arthro_obs_std        mu identity      scaled_yr -0.86606415 0.04355845 -19.8828033
4    scaled_tmean arthro_obs_std        mu identity   scaled_tmean  0.22888083 0.03477478   6.5818050
5      scaled_lat arthro_obs_std        mu identity     scaled_lat  0.39506212 0.03071495  12.8622089
6    scaled_tmean    wing_length        mu identity   scaled_tmean -0.57860357 0.21991026  -2.6310895
7  arthro_obs_std    wing_length        mu identity arthro_obs_std  0.06167576 0.11694391   0.5273961
8             Sex    wing_length        mu identity        SexMale  2.59143777 0.23229215  11.1559423
9      scaled_lat    wing_length        mu identity     scaled_lat  0.03137900 0.22436790   0.1398551
10      scaled_yr    wing_length     sigma      log      scaled_yr -0.09596912 0.03434174  -2.7945329
11   scaled_tmean    wing_length     sigma      log   scaled_tmean -0.16061900 0.02662236  -6.0332367
        p.value endogenous std.estimate
1  2.365608e-02      FALSE   0.02949797
2  8.392147e-71      FALSE  -0.33761897
3  5.733891e-88      FALSE  -0.41522723
4  4.647708e-11       TRUE   0.16824093
5  7.344739e-38      FALSE   0.32847690
6  8.511163e-03       TRUE  -0.42530782
7  5.979186e-01       TRUE   0.06515708
8  6.697941e-29      FALSE   2.59143777
9  8.887745e-01      FALSE   0.02609027
10 5.197476e-03      FALSE  -0.04601159
11 1.607077e-09       TRUE  -0.11806445
```

## Phylogeny-corrected paths (Rubin-pooled across 50 trees)
Pagel's lambda (AIC-selected per tree), range: 1–1
```
                                                 path    estimate         se        lower       upper
1           scaled_yr | scaled_tmean | mu | scaled_yr  0.06152567 0.02719163  0.008230073  0.11482128
2         scaled_lat | scaled_tmean | mu | scaled_lat -0.40605738 0.02282456 -0.450793516 -0.36132125
3         scaled_yr | arthro_obs_std | mu | scaled_yr -0.86606415 0.04355845 -0.951438722 -0.78068959
4   scaled_tmean | arthro_obs_std | mu | scaled_tmean  0.22888083 0.03477478  0.160722257  0.29703940
5       scaled_lat | arthro_obs_std | mu | scaled_lat  0.39506212 0.03071495  0.334860818  0.45526343
6      scaled_tmean | wing_length | mu | scaled_tmean -0.55396924 0.21954037 -0.984268371 -0.12367011
7  arthro_obs_std | wing_length | mu | arthro_obs_std  0.05697613 0.11687841 -0.172105562  0.28605781
8                    Sex | wing_length | mu | SexMale  2.58799635 0.23217972  2.132924096  3.04306861
9          scaled_lat | wing_length | mu | scaled_lat  0.02183225 0.22321253 -0.415664313  0.45932880
10        scaled_yr | wing_length | sigma | scaled_yr -0.09863063 0.03434633 -0.165949438 -0.03131183
11  scaled_tmean | wing_length | sigma | scaled_tmean -0.15995423 0.02660941 -0.212108677 -0.10779978
              z      p.value
1    2.26266932 2.365608e-02
2  -17.79037164 8.392147e-71
3  -19.88280330 5.733891e-88
4    6.58180505 4.647708e-11
5   12.86220888 7.344739e-38
6   -2.52331372 1.162546e-02
7    0.48748203 6.259168e-01
8   11.14652186 7.446032e-29
9    0.09780923 9.220838e-01
10  -2.87165000 4.083349e-03
11  -6.01118941 1.841671e-09
```

## Effect: scaled_yr -> wing_length (all mediators)
```
      from          to                      through target              quantity estimate conf.low
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean            total_path  -0.0398  -0.1399
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean                direct   0.0000   0.0000
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean              indirect  -0.0398  -0.1399
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean         mean_mediated  -0.0397  -0.1395
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean distribution_mediated   0.0000  -0.0018
 conf.high
    0.0484
    0.0000
    0.0484
    0.0470
    0.0017
```

## Effect: scaled_yr -> wing_length via temperature
```
      from          to      through target              quantity estimate conf.low conf.high
 scaled_yr wing_length scaled_tmean   mean            total_path  -0.0164  -0.0357   -0.0022
 scaled_yr wing_length scaled_tmean   mean                direct   0.0000   0.0000    0.0000
 scaled_yr wing_length scaled_tmean   mean              indirect  -0.0164  -0.0357   -0.0022
 scaled_yr wing_length scaled_tmean   mean         mean_mediated  -0.0163  -0.0355   -0.0020
 scaled_yr wing_length scaled_tmean   mean distribution_mediated  -0.0001  -0.0016    0.0014
```

## Effect: scaled_yr -> wing_length via arthropod
```
      from          to        through target              quantity estimate conf.low conf.high
 scaled_yr wing_length arthro_obs_std   mean            total_path  -0.0258  -0.1245    0.0661
 scaled_yr wing_length arthro_obs_std   mean                direct   0.0000   0.0000    0.0000
 scaled_yr wing_length arthro_obs_std   mean              indirect  -0.0258  -0.1245    0.0661
 scaled_yr wing_length arthro_obs_std   mean         mean_mediated  -0.0258  -0.1254    0.0662
 scaled_yr wing_length arthro_obs_std   mean distribution_mediated   0.0000  -0.0005    0.0006
```

## Effect: scaled_tmean -> wing_length
```
         from          to        through target              quantity estimate conf.low conf.high
 scaled_tmean wing_length arthro_obs_std   mean            total_path  -0.3958  -0.7563   -0.0738
 scaled_tmean wing_length arthro_obs_std   mean                direct  -0.4059  -0.7929   -0.0973
 scaled_tmean wing_length arthro_obs_std   mean              indirect   0.0101  -0.0290    0.0446
 scaled_tmean wing_length arthro_obs_std   mean         mean_mediated   0.0101  -0.0290    0.0448
 scaled_tmean wing_length arthro_obs_std   mean distribution_mediated   0.0000  -0.0005    0.0006
```
