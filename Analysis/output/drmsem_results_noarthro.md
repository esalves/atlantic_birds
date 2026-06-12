# drmSEM results — v3-noarthro
_generated 2026-06-12 08:52:31.508181 · N = 9180 records, 86 species, years 1990–2018_

## Fisher's C (global DAG fit)
```
    fisher_c df n_claims   p.value
1 0.02112463  2        1 0.9894933
```

## d-separation claims
```
                                            claim   x            y                 given df
1 Sex _||_ scaled_tmean | {scaled_yr, scaled_lat} Sex scaled_tmean scaled_yr, scaled_lat  1
            LR   p.value status
1 0.0001734124 0.9894933     ok
```

## Path coefficients (raw)
```
          from           to component     link         term    estimate   std.error   statistic
1    scaled_yr scaled_tmean        mu identity    scaled_yr  0.02437741 0.005578465   4.3699141
2   scaled_lat scaled_tmean        mu identity   scaled_lat -0.45461494 0.009746521 -46.6438155
3 scaled_tmean  wing_length        mu identity scaled_tmean  0.12717097 0.109691547   1.1593507
4          Sex  wing_length        mu identity      SexMale  2.52283590 0.107263571  23.5199694
5   scaled_lat  wing_length        mu identity   scaled_lat -0.09784440 0.102811653  -0.9516859
6    scaled_yr  wing_length     sigma      log    scaled_yr  0.06842805 0.007620114   8.9799244
7 scaled_tmean  wing_length     sigma      log scaled_tmean -0.17441174 0.009029188 -19.3164357
        p.value endogenous
1  1.242954e-05      FALSE
2  0.000000e+00      FALSE
3  2.463133e-01       TRUE
4 2.548389e-122      FALSE
5  3.412563e-01      FALSE
6  2.709528e-19      FALSE
7  3.907154e-83       TRUE
```

## Path coefficients (standardised, sd_x)
```
          from           to component     link         term    estimate   std.error   statistic
1    scaled_yr scaled_tmean        mu identity    scaled_yr  0.02437741 0.005578465   4.3699141
2   scaled_lat scaled_tmean        mu identity   scaled_lat -0.45461494 0.009746521 -46.6438155
3 scaled_tmean  wing_length        mu identity scaled_tmean  0.12717097 0.109691547   1.1593507
4          Sex  wing_length        mu identity      SexMale  2.52283590 0.107263571  23.5199694
5   scaled_lat  wing_length        mu identity   scaled_lat -0.09784440 0.102811653  -0.9516859
6    scaled_yr  wing_length     sigma      log    scaled_yr  0.06842805 0.007620114   8.9799244
7 scaled_tmean  wing_length     sigma      log scaled_tmean -0.17441174 0.009029188 -19.3164357
        p.value endogenous std.estimate
1  1.242954e-05      FALSE   0.02533853
2  0.000000e+00      FALSE  -0.45381627
3  2.463133e-01       TRUE   0.11403803
4 2.548389e-122      FALSE   2.52283590
5  3.412563e-01      FALSE  -0.09767251
6  2.709528e-19      FALSE   0.07112592
7  3.907154e-83       TRUE  -0.15640024
```

## Phylogeny-corrected paths (Rubin-pooled across 50 trees)
Pagel's lambda (AIC-selected per tree), range: 1–1
```
                                               path    estimate          se       lower       upper
1         scaled_yr | scaled_tmean | mu | scaled_yr  0.02436719 0.005578865  0.01343261  0.03530176
2       scaled_lat | scaled_tmean | mu | scaled_lat -0.45418279 0.009747513 -0.47328792 -0.43507766
3    scaled_tmean | wing_length | mu | scaled_tmean  0.13073074 0.109669601 -0.08422167  0.34568316
4                  Sex | wing_length | mu | SexMale  2.52190021 0.107276598  2.31163807  2.73216234
5        scaled_lat | wing_length | mu | scaled_lat -0.09680337 0.102681925 -0.29805994  0.10445320
6       scaled_yr | wing_length | sigma | scaled_yr  0.06822036 0.007620490  0.05328420  0.08315652
7 scaled_tmean | wing_length | sigma | scaled_tmean -0.17430558 0.009028425 -0.19200129 -0.15660987
            z       p.value
1   4.3677680  1.255227e-05
2 -46.5947339  0.000000e+00
3   1.1920417  2.332449e-01
4  23.5083909 3.347485e-122
5  -0.9427499  3.458089e-01
6   8.9522276  3.483796e-19
7 -19.3063108  4.753380e-83
```

## Effect: scaled_yr -> wing_length (all mediators)
```
      from          to      through target              quantity estimate conf.low conf.high
 scaled_yr wing_length scaled_tmean   mean            total_path   0.0033  -0.0021    0.0115
 scaled_yr wing_length scaled_tmean   mean                direct   0.0000   0.0000    0.0000
 scaled_yr wing_length scaled_tmean   mean              indirect   0.0033  -0.0021    0.0115
 scaled_yr wing_length scaled_tmean   mean         mean_mediated   0.0033  -0.0020    0.0114
 scaled_yr wing_length scaled_tmean   mean distribution_mediated   0.0000  -0.0002    0.0002
```

## Effect: scaled_tmean -> wing_length
```
         from          to through target              quantity estimate conf.low conf.high
 scaled_tmean wing_length           mean            total_path   0.1128  -0.0832    0.3348
 scaled_tmean wing_length           mean                direct   0.1128  -0.0832    0.3348
 scaled_tmean wing_length           mean              indirect   0.0000   0.0000    0.0000
 scaled_tmean wing_length           mean         mean_mediated   0.0000   0.0000    0.0000
 scaled_tmean wing_length           mean distribution_mediated   0.0000   0.0000    0.0000
```
