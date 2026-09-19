# Secular threshold sensitivity

Generated: 2026-09-18 18:00

Follow-up to SCRIPT_CONSTANTS_AUDIT.md section B. Tests whether the secular
results depend on thresholds that were set without a stated rationale.

## S1. Species-inclusion rule

```
       analysis                             rule     n n_species estimate    se      t mm_decade pct_decade pct_lo   pct_hi
 T4_museum_only n>=10, span>=10 (as run, Tier 4)  3440       145  -0.1420 0.391 -0.363    -0.283     -0.306 -1.950  1.34000
 T4_museum_only n>=30, span>=15 (as run, Tier 5)  2180        42   0.0571 0.195  0.293     0.114      0.135 -0.769  1.04000
 T4_museum_only  n>=30, span>=5  (main analysis)  2180        42   0.0571 0.195  0.293     0.114      0.135 -0.769  1.04000
  T5_integrated n>=10, span>=10 (as run, Tier 4) 37700       322  -0.2610 0.132 -1.970    -0.520     -0.660 -1.320 -0.00473
  T5_integrated n>=30, span>=15 (as run, Tier 5) 34300       184  -0.1410 0.111 -1.270    -0.281     -0.362 -0.924  0.19900
  T5_integrated  n>=30, span>=5  (main analysis) 35800       205  -0.0536 0.154 -0.349    -0.107     -0.138 -0.916  0.63900
```

## S2. SD trim on the annual-mean series (breakpoint fixed at 1980)

```
 series trim n_records n_years pre_slope_decade   pre_t pre_p post_slope_decade  post_t post_p
   wing 3 SD     39000     121          0.04550  0.1610 0.872         -7.15e-02 -0.3270 0.7440
   wing 4 SD     39100     121          0.01960  0.0603 0.952         -4.42e-02 -0.1760 0.8610
   wing 5 SD     39100     121          0.01890  0.0458 0.964         -3.82e-02 -0.1200 0.9050
   wing none     39500     121          0.02230  0.0153 0.988         -3.46e-02 -0.0311 0.9750
   mass 3 SD     48600      50          0.00271  0.1020 0.919         -3.57e-03 -1.1300 0.2650
   mass 4 SD     49000      50         -0.01440 -0.5350 0.596         -2.27e-03 -0.7010 0.4860
   mass 5 SD     49200      50         -0.02640 -0.9650 0.340         -1.42e-03 -0.4240 0.6730
   mass none     49400      50         -0.02900 -1.0600 0.295         -6.99e-05 -0.0209 0.9830
    iso 3 SD     28800      41          8.50000  0.6100 0.545          1.79e+00  1.8400 0.0740
    iso 4 SD     29100      41          7.61000  0.4990 0.621          2.09e+00  1.9700 0.0566
    iso 5 SD     29200      41          6.65000  0.4270 0.671          2.38e+00  2.2000 0.0341
    iso none     29300      41          6.65000  0.4240 0.674          2.41e+00  2.2100 0.0332
```

## S3. Start year for the mass / allometry series

```
 series start_year n_records n_years pre_slope_decade  pre_t pre_p post_slope_decade post_t post_p
   mass       1900     49400      51          -0.0166 -0.646 0.521          -0.00221 -0.691 0.4930
   mass       1940     49400      51          -0.0166 -0.646 0.521          -0.00221 -0.691 0.4930
   mass       1960     49400      50          -0.0144 -0.535 0.596          -0.00227 -0.701 0.4860
   mass       1980     49200      37               NA     NA    NA          -0.00222 -0.621 0.5390
    iso       1900     29300      42           5.1800  0.398 0.692           2.13000  2.050 0.0470
    iso       1940     29300      42           5.1800  0.398 0.692           2.13000  2.050 0.0470
    iso       1960     29300      41           7.6100  0.499 0.621           2.09000  1.970 0.0566
    iso       1980     29100      31               NA     NA    NA           2.14000  1.740 0.0925
```

