# drmSEM results — v3-arthro-spatArthro200km
_generated 2026-06-12 09:13:18.990569 · N = 2149 records, 82 species, years 1998–2008_

## Fisher's C (global DAG fit)
```
  fisher_c df n_claims   p.value
1 7.668831  4        2 0.1044909
```

## d-separation claims
```
                                                            claim   x
1                 Sex _||_ scaled_tmean | {scaled_yr, scaled_lat} Sex
2 Sex _||_ arthro_obs_std | {scaled_yr, scaled_tmean, scaled_lat} Sex
               y                               given df        LR   p.value
1   scaled_tmean               scaled_yr, scaled_lat  1 0.5136052 0.4735823
2 arthro_obs_std scaled_yr, scaled_tmean, scaled_lat  1 3.9948579 0.0456393
  status
1     ok
2     ok
```

## Path coefficients (raw)
```
             from             to component     link           term    estimate
1       scaled_yr   scaled_tmean        mu identity      scaled_yr  0.06408012
2      scaled_lat   scaled_tmean        mu identity     scaled_lat -0.39220628
3       scaled_yr arthro_obs_std        mu identity      scaled_yr -0.86846043
4    scaled_tmean arthro_obs_std        mu identity   scaled_tmean  0.20289863
5      scaled_lat arthro_obs_std        mu identity     scaled_lat  0.39830573
6    scaled_tmean    wing_length        mu identity   scaled_tmean -0.63087904
7  arthro_obs_std    wing_length        mu identity arthro_obs_std  0.08529373
8             Sex    wing_length        mu identity        SexMale  2.58245537
9      scaled_lat    wing_length        mu identity     scaled_lat  0.03684786
10      scaled_yr    wing_length     sigma      log      scaled_yr -0.11262083
11   scaled_tmean    wing_length     sigma      log   scaled_tmean -0.12085624
    std.error   statistic      p.value endogenous
1  0.02560556   2.5025864 1.232895e-02      FALSE
2  0.02093813 -18.7316800 2.731587e-78      FALSE
3  0.04158411 -20.8844290 7.418615e-97      FALSE
4  0.03392106   5.9814951 2.210987e-09       TRUE
5  0.02930515  13.5916611 4.487783e-42      FALSE
6  0.22080218  -2.8572138 4.273779e-03       TRUE
7  0.11760303   0.7252681 4.682876e-01       TRUE
8  0.22897099  11.2785266 1.675284e-29      FALSE
9  0.22215814   0.1658632 8.682646e-01      FALSE
10 0.03414327  -3.2984783 9.721039e-04      FALSE
11 0.02523788  -4.7886846 1.678781e-06       TRUE
```

## Path coefficients (standardised, sd_x)
```
             from             to component     link           term    estimate
1       scaled_yr   scaled_tmean        mu identity      scaled_yr  0.06408012
2      scaled_lat   scaled_tmean        mu identity     scaled_lat -0.39220628
3       scaled_yr arthro_obs_std        mu identity      scaled_yr -0.86846043
4    scaled_tmean arthro_obs_std        mu identity   scaled_tmean  0.20289863
5      scaled_lat arthro_obs_std        mu identity     scaled_lat  0.39830573
6    scaled_tmean    wing_length        mu identity   scaled_tmean -0.63087904
7  arthro_obs_std    wing_length        mu identity arthro_obs_std  0.08529373
8             Sex    wing_length        mu identity        SexMale  2.58245537
9      scaled_lat    wing_length        mu identity     scaled_lat  0.03684786
10      scaled_yr    wing_length     sigma      log      scaled_yr -0.11262083
11   scaled_tmean    wing_length     sigma      log   scaled_tmean -0.12085624
    std.error   statistic      p.value endogenous std.estimate
1  0.02560556   2.5025864 1.232895e-02      FALSE   0.03094369
2  0.02093813 -18.7316800 2.731587e-78      FALSE  -0.32773784
3  0.04158411 -20.8844290 7.418615e-97      FALSE  -0.41937148
4  0.03392106   5.9814951 2.210987e-09       TRUE   0.14700079
5  0.02930515  13.5916611 4.487783e-42      FALSE   0.33283470
6  0.22080218  -2.8572138 4.273779e-03       TRUE  -0.45707416
7  0.11760303   0.7252681 4.682876e-01       TRUE   0.09000062
8  0.22897099  11.2785266 1.675284e-29      FALSE   2.58245537
9  0.22215814   0.1658632 8.682646e-01      FALSE   0.03079104
10 0.03414327  -3.2984783 9.721039e-04      FALSE  -0.05438355
11 0.02523788  -4.7886846 1.678781e-06       TRUE  -0.08756079
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
            total_path  -0.0562  -0.1553    0.0263
                direct   0.0000   0.0000    0.0000
              indirect  -0.0562  -0.1553    0.0263
         mean_mediated  -0.0563  -0.1544    0.0257
 distribution_mediated   0.0000  -0.0018    0.0018
```

## Effect: scaled_yr -> wing_length via temperature
```
      from          to      through target              quantity estimate
 scaled_yr wing_length scaled_tmean   mean            total_path  -0.0203
 scaled_yr wing_length scaled_tmean   mean                direct   0.0000
 scaled_yr wing_length scaled_tmean   mean              indirect  -0.0203
 scaled_yr wing_length scaled_tmean   mean         mean_mediated  -0.0202
 scaled_yr wing_length scaled_tmean   mean distribution_mediated  -0.0001
 conf.low conf.high
  -0.0417   -0.0033
   0.0000    0.0000
  -0.0417   -0.0033
  -0.0416   -0.0035
  -0.0015    0.0012
```

## Effect: scaled_yr -> wing_length via arthropod
```
      from          to        through target              quantity estimate
 scaled_yr wing_length arthro_obs_std   mean            total_path  -0.0331
 scaled_yr wing_length arthro_obs_std   mean                direct   0.0000
 scaled_yr wing_length arthro_obs_std   mean              indirect  -0.0331
 scaled_yr wing_length arthro_obs_std   mean         mean_mediated  -0.0331
 scaled_yr wing_length arthro_obs_std   mean distribution_mediated   0.0000
 conf.low conf.high
  -0.1238    0.0596
   0.0000    0.0000
  -0.1238    0.0596
  -0.1236    0.0599
  -0.0006    0.0006
```

## Effect: scaled_tmean -> wing_length
```
         from          to        through target              quantity estimate
 scaled_tmean wing_length arthro_obs_std   mean            total_path  -0.4485
 scaled_tmean wing_length arthro_obs_std   mean                direct  -0.4599
 scaled_tmean wing_length arthro_obs_std   mean              indirect   0.0114
 scaled_tmean wing_length arthro_obs_std   mean         mean_mediated   0.0113
 scaled_tmean wing_length arthro_obs_std   mean distribution_mediated   0.0000
 conf.low conf.high
  -0.7366   -0.1896
  -0.7438   -0.1973
  -0.0207    0.0477
  -0.0206    0.0474
  -0.0006    0.0006
```
