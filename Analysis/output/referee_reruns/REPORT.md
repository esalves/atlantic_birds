# Referee re-runs (E1-E4, E6-E7)

Generated: 2026-09-23 11:12 | trees: 50 | glmmTMB 1.1.14

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
                 model                   response      term  estimate       se       z        p    n      lower      upper
 E1a_anomaly_plus_year           wing length (mm)        dT -0.021950 0.195600 -0.1122 0.910600 7224 -0.4053000  0.3614000
 E1a_anomaly_plus_year           wing length (mm) scaled_yr -0.538200 0.219000 -2.4570 0.014020 7224 -0.9675000 -0.1088000
 E1b_detrended_anomaly           wing length (mm)   dT_detr -0.163000 0.195300 -0.8345 0.404000 7224 -0.5459000  0.2198000
 E1a_anomaly_plus_year              log body mass        dT  0.007267 0.006871  1.0580 0.290200 7224 -0.0062000  0.0207300
 E1a_anomaly_plus_year              log body mass scaled_yr -0.001889 0.005208 -0.3626 0.716900 7224 -0.0121000  0.0083190
 E1b_detrended_anomaly              log body mass   dT_detr  0.006403 0.006811  0.9402 0.347100 7224 -0.0069450  0.0197500
 E1a_anomaly_plus_year isometry contrast lnM-3lnL        dT  0.015710 0.009947  1.5790 0.114200 7224 -0.0037850  0.0352100
 E1a_anomaly_plus_year isometry contrast lnM-3lnL scaled_yr  0.019920 0.007211  2.7630 0.005734 7224  0.0057870  0.0340500
 E1b_detrended_anomaly isometry contrast lnM-3lnL   dT_detr  0.020140 0.009851  2.0450 0.040900 7224  0.0008333  0.0394500
 E1a_anomaly_plus_year relative wing lnL-(1/3)lnM        dT -0.005237 0.003316 -1.5790 0.114200 7224 -0.0117400  0.0012620
 E1a_anomaly_plus_year relative wing lnL-(1/3)lnM scaled_yr -0.006640 0.002404 -2.7630 0.005734 7224 -0.0113500 -0.0019290
 E1b_detrended_anomaly relative wing lnL-(1/3)lnM   dT_detr -0.006713 0.003284 -2.0450 0.040900 7224 -0.0131500 -0.0002778

```

## E2. Locality-year clustering (Mizuno c85)

```
                        model                   response term  estimate       se      z       p    n    lower     upper
 E2a_locyear_random_intercept           wing length (mm)   dT -0.329300 0.292000 -1.128 0.25940 7224 -0.90170 0.2430000
 E2a_locyear_random_intercept              log body mass   dT  0.005739 0.009051  0.634 0.52610 7224 -0.01200 0.0234800
 E2a_locyear_random_intercept isometry contrast lnM-3lnL   dT  0.025730 0.014540  1.770 0.07681 7224 -0.00277 0.0542400
 E2a_locyear_random_intercept relative wing lnL-(1/3)lnM   dT -0.008578 0.004848 -1.770 0.07681 7224 -0.01808 0.0009233

```

## E3. Thermal models at M3 adjustment across trees (Mizuno c86)

```
                           model                   response      term  estimate       se        z      p    n     lower     upper
            E3_M3_adjusted_phylo           wing length (mm)        dT  0.010430 0.193500  0.05389 0.9570 7070 -0.368900 0.3897000
 E3b_M3_adjusted_phylo_plus_year           wing length (mm)        dT  0.112000 0.195800  0.57170 0.5675 7070 -0.271900 0.4958000
 E3b_M3_adjusted_phylo_plus_year           wing length (mm) scaled_yr -0.293600 0.207700 -1.41400 0.1574 7070 -0.700600 0.1134000
            E3_M3_adjusted_phylo              log body mass        dT  0.004850 0.006797  0.71350 0.4755 7070 -0.008472 0.0181700
 E3b_M3_adjusted_phylo_plus_year              log body mass        dT  0.005960 0.006864  0.86830 0.3853 7070 -0.007494 0.0194100
 E3b_M3_adjusted_phylo_plus_year              log body mass scaled_yr -0.001788 0.004414 -0.40520 0.6854 7070 -0.010440 0.0068630
            E3_M3_adjusted_phylo isometry contrast lnM-3lnL        dT  0.011700 0.009931  1.17800 0.2388 7070 -0.007767 0.0311600
 E3b_M3_adjusted_phylo_plus_year isometry contrast lnM-3lnL        dT  0.009455 0.010050  0.94050 0.3470 7070 -0.010250 0.0291600
 E3b_M3_adjusted_phylo_plus_year isometry contrast lnM-3lnL scaled_yr  0.011170 0.006964  1.60400 0.1087 7070 -0.002479 0.0248200
            E3_M3_adjusted_phylo relative wing lnL-(1/3)lnM        dT -0.003899 0.003310 -1.17800 0.2388 7070 -0.010390 0.0025890
 E3b_M3_adjusted_phylo_plus_year relative wing lnL-(1/3)lnM        dT -0.003152 0.003351 -0.94050 0.3470 7070 -0.009720 0.0034160
 E3b_M3_adjusted_phylo_plus_year relative wing lnL-(1/3)lnM scaled_yr -0.003723 0.002321 -1.60400 0.1087 7070 -0.008273 0.0008265
 trees_used
         50
         50
         50
         38
         50
         50
         48
         50
         50
         47
         50
         50

```

## E4. Isometry contrast at M3 adjustment (Mizuno c81, c82)

```
                  model                   response      term  estimate       se       z       p    n     lower    upper trees_used
   E4a_iso_parsimonious isometry contrast lnM-3lnL scaled_yr  0.012330 0.007092  1.7390 0.08207 7070 -0.001569 0.026230         50
             E4b_iso_M3 isometry contrast lnM-3lnL scaled_yr  0.011700 0.006953  1.6820 0.09255 7070 -0.001932 0.025320         50
 E4c_separate_M3_slopes                     lnwing scaled_yr -0.003760 0.002702 -1.3910 0.16410 7070 -0.009056 0.001536         50
 E4c_separate_M3_slopes                     lnmass scaled_yr -0.001454 0.004396 -0.3306 0.74090 7070 -0.010070 0.007163         50
 pct_decade pct_decade_lo pct_decade_hi
     2.4570       -0.3125         5.226
     2.3300       -0.3849         5.045
    -0.7462            NA            NA
    -0.2891            NA            NA

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
   -1.36400 -2.7040 -0.02419
   -0.64530 -1.3460  0.05551
   -0.55600 -1.5910  0.47880

```

