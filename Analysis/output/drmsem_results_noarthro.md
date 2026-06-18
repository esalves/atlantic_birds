# drmSEM results — v3-noarthro
_generated 2026-06-15 17:50:34.128904 · N = 7810 records, 72 species, years 1995–2018_

## Fisher's C (global DAG fit)
```
   fisher_c df n_claims   p.value
1 0.5076503  2        1 0.7758275
```

## d-separation claims
```
                                            claim   x            y
1 Sex _||_ scaled_tmean | {scaled_yr, scaled_lat} Sex scaled_tmean
                  given df         LR   p.value status
1 scaled_yr, scaled_lat  1 0.08108868 0.7758275     ok
```

## Path coefficients (raw)
```
          from           to component     link         term     estimate
1    scaled_yr scaled_tmean        mu identity    scaled_yr -0.005244378
2   scaled_lat scaled_tmean        mu identity   scaled_lat -0.502278924
3 scaled_tmean  wing_length        mu identity scaled_tmean -0.215544939
4          Sex  wing_length        mu identity      SexMale  2.611727755
5   scaled_lat  wing_length        mu identity   scaled_lat -0.432180828
6    scaled_yr  wing_length     sigma      log    scaled_yr  0.034775356
7 scaled_tmean  wing_length     sigma      log scaled_tmean -0.192686715
    std.error   statistic       p.value endogenous
1 0.005896568  -0.8893949  3.737909e-01      FALSE
2 0.012006749 -41.8330490  0.000000e+00      FALSE
3 0.113904646  -1.8923279  5.844731e-02       TRUE
4 0.105435775  24.7707929 1.851169e-135      FALSE
5 0.118099416  -3.6594663  2.527411e-04      FALSE
6 0.007827440   4.4427498  8.881642e-06      FALSE
7 0.009129557 -21.1058112  7.033735e-99       TRUE
```

## Path coefficients (standardised, sd_x)
```
          from           to component     link         term     estimate
1    scaled_yr scaled_tmean        mu identity    scaled_yr -0.005244378
2   scaled_lat scaled_tmean        mu identity   scaled_lat -0.502278924
3 scaled_tmean  wing_length        mu identity scaled_tmean -0.215544939
4          Sex  wing_length        mu identity      SexMale  2.611727755
5   scaled_lat  wing_length        mu identity   scaled_lat -0.432180828
6    scaled_yr  wing_length     sigma      log    scaled_yr  0.034775356
7 scaled_tmean  wing_length     sigma      log scaled_tmean -0.192686715
    std.error   statistic       p.value endogenous std.estimate
1 0.005896568  -0.8893949  3.737909e-01      FALSE -0.005461886
2 0.012006749 -41.8330490  0.000000e+00      FALSE -0.477137498
3 0.113904646  -1.8923279  5.844731e-02       TRUE -0.193091757
4 0.105435775  24.7707929 1.851169e-135      FALSE  2.611727755
5 0.118099416  -3.6594663  2.527411e-04      FALSE -0.410548142
6 0.007827440   4.4427498  8.881642e-06      FALSE  0.036217650
7 0.009129557 -21.1058112  7.033735e-99       TRUE -0.172614660
```

## Phylogeny-corrected paths (Rubin-pooled across 50 trees)
Pagel's lambda (AIC-selected per tree), range: 1–1
```
                                               path     estimate          se
1         scaled_yr | scaled_tmean | mu | scaled_yr -0.005244378 0.005896568
2       scaled_lat | scaled_tmean | mu | scaled_lat -0.502278924 0.012006749
3    scaled_tmean | wing_length | mu | scaled_tmean -0.212251213 0.113863746
4                  Sex | wing_length | mu | SexMale  2.610330026 0.105430835
5        scaled_lat | wing_length | mu | scaled_lat -0.431312699 0.117853679
6       scaled_yr | wing_length | sigma | scaled_yr  0.034702393 0.007825895
7 scaled_tmean | wing_length | sigma | scaled_tmean -0.192724539 0.009127798
        lower        upper           z       p.value
1 -0.01680165  0.006312896  -0.8893949  3.737909e-01
2 -0.52581215 -0.478745696 -41.8330490  0.000000e+00
3 -0.43542415  0.010921728  -1.8640807  6.231037e-02
4  2.40368559  2.816974462  24.7586963 2.498954e-135
5 -0.66230591 -0.200319489  -3.6597305  2.524807e-04
6  0.01936364  0.050041147   4.4343035  9.237043e-06
7 -0.21061502 -0.174834055 -21.1140237  5.911898e-99
```

## Effect: scaled_yr -> wing_length (all mediators)
```
      from          to      through target              quantity estimate
 scaled_yr wing_length scaled_tmean   mean            total_path   0.0013
 scaled_yr wing_length scaled_tmean   mean                direct   0.0000
 scaled_yr wing_length scaled_tmean   mean              indirect   0.0013
 scaled_yr wing_length scaled_tmean   mean         mean_mediated   0.0013
 scaled_yr wing_length scaled_tmean   mean distribution_mediated   0.0000
 conf.low conf.high
  -0.0012    0.0049
   0.0000    0.0000
  -0.0012    0.0049
  -0.0010    0.0047
  -0.0003    0.0002
```

## Effect: scaled_tmean -> wing_length
```
         from          to through target              quantity estimate
 scaled_tmean wing_length           mean            total_path  -0.1939
 scaled_tmean wing_length           mean                direct  -0.1939
 scaled_tmean wing_length           mean              indirect   0.0000
 scaled_tmean wing_length           mean         mean_mediated   0.0000
 scaled_tmean wing_length           mean distribution_mediated   0.0000
 conf.low conf.high
  -0.3999    0.0463
  -0.3999    0.0463
   0.0000    0.0000
   0.0000    0.0000
   0.0000    0.0000
```
