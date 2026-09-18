# Referee re-runs (E1-E4, E6-E7)

Generated: 2026-09-18 17:00 | trees: 50 | glmmTMB 1.1.14

Commissioned from the co-author review of 2026-09-18. See `REVIEW_RESPONSE_PLAN.md`.

**Reproduction note.** The script that produced the reported thermal-coupling
estimates was not in the repository; block E0 refits that model from
`passer90_climate.rds` so the reported numbers have a reproducible source.

## E0. Published parsimonious thermal models (reproduction)

```
           model                   response term  estimate       se       z       p    n     lower     upper
 E0_parsimonious           wing length (mm)   dT -0.105300 0.201100 -0.5236 0.60060 7224 -0.499500 0.2889000
 E0_parsimonious              log body mass   dT  0.005740 0.006871  0.8355 0.40340 7224 -0.007726 0.0192100
 E0_parsimonious isometry contrast lnM-3lnL   dT  0.017460 0.010040  1.7380 0.08214 7224 -0.002225 0.0371400
 E0_parsimonious relative wing lnL-(1/3)lnM   dT -0.005818 0.003347 -1.7380 0.08214 7224 -0.012380 0.0007417

```

## E1. Anomaly with calendar year, and detrended anomaly (Mizuno c84)

```
                 model                   response      term  estimate       se       z         p    n     lower     upper
 E1a_anomaly_plus_year           wing length (mm)        dT  0.054130 0.198900  0.2722 0.7855000 7224 -0.335700  0.443900
 E1a_anomaly_plus_year           wing length (mm) scaled_yr -0.700300 0.225100 -3.1110 0.0018670 7224 -1.142000 -0.259000
 E1b_detrended_anomaly           wing length (mm)   dT_detr -0.060180 0.201500 -0.2987 0.7652000 7224 -0.455100  0.334700
 E1a_anomaly_plus_year              log body mass        dT  0.006797 0.006915  0.9831 0.3256000 7224 -0.006755  0.020350
 E1a_anomaly_plus_year              log body mass scaled_yr -0.003319 0.004461 -0.7441 0.4568000 7224 -0.012060  0.005423
 E1b_detrended_anomaly              log body mass   dT_detr  0.006036 0.006883  0.8769 0.3805000 7224 -0.007455  0.019530
 E1a_anomaly_plus_year isometry contrast lnM-3lnL        dT  0.011850 0.010070  1.1770 0.2394000 7224 -0.007890  0.031590
 E1a_anomaly_plus_year isometry contrast lnM-3lnL scaled_yr  0.024950 0.007096  3.5160 0.0004387 7224  0.011040  0.038850
 E1b_detrended_anomaly isometry contrast lnM-3lnL   dT_detr  0.015900 0.010060  1.5800 0.1140000 7224 -0.003818  0.035620
 E1a_anomaly_plus_year relative wing lnL-(1/3)lnM        dT -0.003949 0.003357 -1.1770 0.2394000 7224 -0.010530  0.002630
 E1a_anomaly_plus_year relative wing lnL-(1/3)lnM scaled_yr -0.008315 0.002365 -3.5160 0.0004388 7224 -0.012950 -0.003679
 E1b_detrended_anomaly relative wing lnL-(1/3)lnM   dT_detr -0.005300 0.003353 -1.5800 0.1140000 7224 -0.011870  0.001273

```

## E2. Locality-year clustering (Mizuno c85)

```
                        model                   response term  estimate       se       z      p    n     lower    upper
 E2a_locyear_random_intercept           wing length (mm)   dT -0.388000 0.411500 -0.9427 0.3458 7224 -1.195000 0.418600
 E2a_locyear_random_intercept              log body mass   dT  0.004861 0.010850  0.4480 0.6542 7224 -0.016410 0.026130
 E2a_locyear_random_intercept isometry contrast lnM-3lnL   dT  0.027260 0.018940  1.4390 0.1501 7224 -0.009867 0.064390
 E2a_locyear_random_intercept relative wing lnL-(1/3)lnM   dT -0.009087 0.006314 -1.4390 0.1501 7224 -0.021460 0.003289

```

## E3. Thermal models at M3 adjustment across trees (Mizuno c86)

```
                           model                   response      term  estimate       se       z         p    n     lower     upper
            E3_M3_adjusted_phylo           wing length (mm)        dT  0.086820 0.199600  0.4351 0.6635000 7070 -0.304300  0.477900
 E3b_M3_adjusted_phylo_plus_year           wing length (mm)        dT  0.231400 0.199400  1.1600 0.2459000 7070 -0.159500  0.622300
 E3b_M3_adjusted_phylo_plus_year           wing length (mm) scaled_yr -0.682500 0.209000 -3.2650 0.0010930 7070 -1.092000 -0.272800
            E3_M3_adjusted_phylo              log body mass        dT  0.004749 0.006857  0.6925 0.4886000 7070 -0.008692  0.018190
 E3b_M3_adjusted_phylo_plus_year              log body mass        dT  0.005666 0.006901  0.8212 0.4116000 7070 -0.007859  0.019190
 E3b_M3_adjusted_phylo_plus_year              log body mass scaled_yr -0.001205 0.003491 -0.3450 0.7301000 7070 -0.008047  0.005638
            E3_M3_adjusted_phylo isometry contrast lnM-3lnL        dT  0.009293 0.010130  0.9169 0.3592000 7070 -0.010570  0.029160
 E3b_M3_adjusted_phylo_plus_year isometry contrast lnM-3lnL        dT  0.004646 0.010180  0.4563 0.6482000 7070 -0.015310  0.024600
 E3b_M3_adjusted_phylo_plus_year isometry contrast lnM-3lnL scaled_yr  0.025200 0.006909  3.6470 0.0002653 7070  0.011660  0.038740
            E3_M3_adjusted_phylo relative wing lnL-(1/3)lnM        dT -0.003097 0.003378 -0.9169 0.3592000 7070 -0.009719  0.003524
 E3b_M3_adjusted_phylo_plus_year relative wing lnL-(1/3)lnM        dT -0.001548 0.003394 -0.4562 0.6482000 7070 -0.008201  0.005104
 E3b_M3_adjusted_phylo_plus_year relative wing lnL-(1/3)lnM scaled_yr -0.008400 0.002304 -3.6470 0.0002658 7070 -0.012920 -0.003885
 trees_used
         49
         50
         50
         50
         50
         50
         49
         50
         50
         48
         50
         50

```

## E4. Isometry contrast at M3 adjustment (Mizuno c81, c82)

```
                  model                   response      term  estimate       se       z         p    n     lower     upper trees_used
   E4a_iso_parsimonious isometry contrast lnM-3lnL scaled_yr  0.022380 0.007101  3.1520 0.0016210 7070  0.008464  0.036300         50
             E4b_iso_M3 isometry contrast lnM-3lnL scaled_yr  0.025320 0.006911  3.6640 0.0002484 7070  0.011780  0.038870         50
 E4c_separate_M3_slopes                     lnwing scaled_yr -0.009316 0.002791 -3.3380 0.0008428 7070 -0.014790 -0.003847         50
 E4c_separate_M3_slopes                     lnmass scaled_yr -0.001037 0.003484 -0.2977 0.7659000 7070 -0.007865  0.005791         50
 pct_decade pct_decade_lo pct_decade_hi
     4.4590         1.686         7.231
     5.0440         2.346         7.743
    -1.8390            NA            NA
    -0.2064            NA            NA

```

## E6/E7. Secular series: source calibration and age (Mizuno c41, c42)

```
                         model                                 sample                   term estimate      se       t     n mm_decade
        S0_asrun_status_offset      spanning cohort, Status intercept              scaled_yr -0.27800 0.14060 -1.9770 24580  -0.55380
      E6a_no_source_adjustment        spanning cohort, no Status term              scaled_yr -0.27610 0.13860 -1.9920 24580  -0.55010
 E6b_year_x_status_interaction                        spanning cohort              scaled_yr -0.51380 0.16840 -3.0510 24580  -1.02300
 E6b_year_x_status_interaction                        spanning cohort scaled_yr:Statusmuseum  0.36900 0.14410  2.5600 24580   0.73510
               E6c_museum_only  spanning cohort, museum only (n=3614)              scaled_yr -0.01615 0.08448 -0.1911  3315  -0.03216
                 E6d_live_only   spanning cohort, live only (n=21268)              scaled_yr -0.53620 0.26870 -1.9950 21270  -1.06800
             E7a_age_covariate         spanning cohort, Status + Age3              scaled_yr -0.25360 0.14050 -1.8050 24580  -0.50530
               E7b_adults_only spanning cohort, adults only (n=19245)              scaled_yr -0.21850 0.20750 -1.0530 18960  -0.43530
 pct_decade  pct_lo   pct_hi
   -0.70730 -1.4080 -0.00612
   -0.70250 -1.3940 -0.01138
   -1.30700      NA       NA
    0.93890      NA       NA
   -0.04108 -0.4624  0.38020
   -1.36400 -2.7040 -0.02420
   -0.64530 -1.3460  0.05551
   -0.55600 -1.5910  0.47880

```

