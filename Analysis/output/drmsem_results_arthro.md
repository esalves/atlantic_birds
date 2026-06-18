# drmSEM results — v3-arthro-spatArthro200km
_generated 2026-06-16 09:43:37.576859 · N = 1697 records, 68 species, years 1998–2008_

## Fisher's C (global DAG fit)
```
  fisher_c df n_claims     p.value
1 13.47291  4        2 0.009181997
```

## d-separation claims
```
                                                            claim   x
1                 Sex _||_ scaled_tmean | {scaled_yr, scaled_lat} Sex
2 Sex _||_ arthro_obs_std | {scaled_yr, scaled_tmean, scaled_lat} Sex
               y                               given df       LR     p.value
1   scaled_tmean               scaled_yr, scaled_lat  1 1.784354 0.181615546
2 arthro_obs_std scaled_yr, scaled_tmean, scaled_lat  1 7.396539 0.006534949
  status
1     ok
2     ok
```

## Path coefficients (raw)
```
             from             to component     link           term    estimate
1       scaled_yr   scaled_tmean        mu identity      scaled_yr -0.04770959
2      scaled_lat   scaled_tmean        mu identity     scaled_lat -0.47494959
3       scaled_yr arthro_obs_std        mu identity      scaled_yr -1.27359870
4    scaled_tmean arthro_obs_std        mu identity   scaled_tmean  0.07651912
5      scaled_lat arthro_obs_std        mu identity     scaled_lat  0.10024976
6    scaled_tmean    wing_length        mu identity   scaled_tmean -0.79336166
7  arthro_obs_std    wing_length        mu identity arthro_obs_std  0.07101351
8             Sex    wing_length        mu identity        SexMale  2.92396891
9      scaled_lat    wing_length        mu identity     scaled_lat -0.38888401
10      scaled_yr    wing_length     sigma      log      scaled_yr -0.22067212
11   scaled_tmean    wing_length     sigma      log   scaled_tmean -0.04016675
    std.error   statistic       p.value endogenous
1  0.03359681  -1.4200633  1.555893e-01      FALSE
2  0.03070100 -15.4701672  5.516359e-54      FALSE
3  0.04940115 -25.7807500 1.457860e-146      FALSE
4  0.03521361   2.1729982  2.978046e-02       TRUE
5  0.03973359   2.5230484  1.163424e-02      FALSE
6  0.21671620  -3.6608323  2.513972e-04       TRUE
7  0.12310400   0.5768579  5.640354e-01       TRUE
8  0.23076623  12.6706967  8.594762e-37      FALSE
9  0.26257844  -1.4810203  1.386012e-01      FALSE
10 0.03515970  -6.2762794  3.467709e-10      FALSE
11 0.02871764  -1.3986786  1.619094e-01       TRUE
```

## Path coefficients (standardised, sd_x)
```
             from             to component     link           term    estimate
1       scaled_yr   scaled_tmean        mu identity      scaled_yr -0.04770959
2      scaled_lat   scaled_tmean        mu identity     scaled_lat -0.47494959
3       scaled_yr arthro_obs_std        mu identity      scaled_yr -1.27359870
4    scaled_tmean arthro_obs_std        mu identity   scaled_tmean  0.07651912
5      scaled_lat arthro_obs_std        mu identity     scaled_lat  0.10024976
6    scaled_tmean    wing_length        mu identity   scaled_tmean -0.79336166
7  arthro_obs_std    wing_length        mu identity arthro_obs_std  0.07101351
8             Sex    wing_length        mu identity        SexMale  2.92396891
9      scaled_lat    wing_length        mu identity     scaled_lat -0.38888401
10      scaled_yr    wing_length     sigma      log      scaled_yr -0.22067212
11   scaled_tmean    wing_length     sigma      log   scaled_tmean -0.04016675
    std.error   statistic       p.value endogenous std.estimate
1  0.03359681  -1.4200633  1.555893e-01      FALSE  -0.02268844
2  0.03070100 -15.4701672  5.516359e-54      FALSE  -0.33252847
3  0.04940115 -25.7807500 1.457860e-146      FALSE  -0.60566358
4  0.03521361   2.1729982  2.978046e-02       TRUE   0.05485871
5  0.03973359   2.5230484  1.163424e-02      FALSE   0.07018829
6  0.21671620  -3.6608323  2.513972e-04       TRUE  -0.56878335
7  0.12310400   0.5768579  5.640354e-01       TRUE   0.07509210
8  0.23076623  12.6706967  8.594762e-37      FALSE   2.92396891
9  0.26257844  -1.4810203  1.386012e-01      FALSE  -0.27227101
10 0.03515970  -6.2762794  3.467709e-10      FALSE  -0.10494127
11 0.02871764  -1.3986786  1.619094e-01       TRUE  -0.02879667
```

## Phylogeny-corrected paths (Rubin-pooled across 50 trees)
Pagel's lambda (AIC-selected per tree), range: 1–1
```
                                                 path    estimate         se
1           scaled_yr | scaled_tmean | mu | scaled_yr -0.04770959 0.03359681
2         scaled_lat | scaled_tmean | mu | scaled_lat -0.47494959 0.03070100
3         scaled_yr | arthro_obs_std | mu | scaled_yr -1.27359870 0.04940115
4   scaled_tmean | arthro_obs_std | mu | scaled_tmean  0.07651912 0.03521361
5       scaled_lat | arthro_obs_std | mu | scaled_lat  0.10024976 0.03973359
6      scaled_tmean | wing_length | mu | scaled_tmean -0.78317519 0.21639867
7  arthro_obs_std | wing_length | mu | arthro_obs_std  0.06674253 0.12297297
8                    Sex | wing_length | mu | SexMale  2.92209188 0.23063417
9          scaled_lat | wing_length | mu | scaled_lat -0.42193520 0.26075080
10        scaled_yr | wing_length | sigma | scaled_yr -0.22163816 0.03514968
11  scaled_tmean | wing_length | sigma | scaled_tmean -0.04078359 0.02868758
          lower       upper           z       p.value
1  -0.113559337  0.01814015  -1.4200633  1.555893e-01
2  -0.535123547 -0.41477563 -15.4701672  5.516359e-54
3  -1.370424950 -1.17677244 -25.7807500 1.457860e-146
4   0.007500435  0.14553780   2.1729982  2.978046e-02
5   0.022371932  0.17812759   2.5230484  1.163424e-02
6  -1.207316579 -0.35903380  -3.6191313  2.955936e-04
7  -0.174284483  0.30776954   0.5427415  5.873078e-01
8   2.470048905  3.37413486  12.6698132  8.692113e-37
9  -0.933006781  0.08913637  -1.6181549  1.056292e-01
10 -0.290531536 -0.15274479  -6.3055526  2.871675e-10
11 -0.097011242  0.01544406  -1.4216463  1.551289e-01
```

## Effect: scaled_yr -> wing_length (all mediators)
```
      from          to                      through target
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean
 scaled_yr wing_length scaled_tmean, arthro_obs_std   mean
              quantity estimate conf.low conf.high
            total_path  -0.0236  -0.1737    0.1232
                direct   0.0000   0.0000    0.0000
              indirect  -0.0236  -0.1737    0.1232
         mean_mediated  -0.0235  -0.1741    0.1208
 distribution_mediated  -0.0001  -0.0020    0.0019
```

## Effect: scaled_yr -> wing_length via temperature
```
      from          to      through target              quantity estimate
 scaled_yr wing_length scaled_tmean   mean            total_path   0.0197
 scaled_yr wing_length scaled_tmean   mean                direct   0.0000
 scaled_yr wing_length scaled_tmean   mean              indirect   0.0197
 scaled_yr wing_length scaled_tmean   mean         mean_mediated   0.0196
 scaled_yr wing_length scaled_tmean   mean distribution_mediated   0.0000
 conf.low conf.high
  -0.0022    0.0479
   0.0000    0.0000
  -0.0022    0.0479
  -0.0026    0.0480
  -0.0019    0.0023
```

## Effect: scaled_yr -> wing_length via arthropod
```
      from          to        through target              quantity estimate
 scaled_yr wing_length arthro_obs_std   mean            total_path   -0.045
 scaled_yr wing_length arthro_obs_std   mean                direct    0.000
 scaled_yr wing_length arthro_obs_std   mean              indirect   -0.045
 scaled_yr wing_length arthro_obs_std   mean         mean_mediated   -0.045
 scaled_yr wing_length arthro_obs_std   mean distribution_mediated    0.000
 conf.low conf.high
  -0.1889    0.1039
   0.0000    0.0000
  -0.1889    0.1039
  -0.1895    0.1040
  -0.0007    0.0007
```

## Effect: scaled_tmean -> wing_length
```
         from          to        through target              quantity estimate
 scaled_tmean wing_length arthro_obs_std   mean            total_path  -0.5883
 scaled_tmean wing_length arthro_obs_std   mean                direct  -0.5922
 scaled_tmean wing_length arthro_obs_std   mean              indirect   0.0039
 scaled_tmean wing_length arthro_obs_std   mean         mean_mediated   0.0039
 scaled_tmean wing_length arthro_obs_std   mean distribution_mediated   0.0000
 conf.low conf.high
  -0.9393   -0.2695
  -0.9423   -0.2816
  -0.0103    0.0209
  -0.0102    0.0204
  -0.0007    0.0007
```
