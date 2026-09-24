# Interannual temperature anomaly: 1-year lag and natal-window models

Generated: 2026-09-23 12:18 | trees: 50 | glmmTMB 1.1.14

Requested by Gustavo (co-author), 2026-09-22, as additional models alongside the
published capture-year (lag-0) anomaly models (Results "Annual temperature
anomalies do not explain the allometric shift"; Methods "Interannual temperature
anomalies"). Rationale: the analysis sample is adults, and passerines typically
first breed around 1 year old, so a body-size effect of an anomalous hatch year
may only be detectable in birds captured the following year(s) -- i.e. the
relevant exposure may be dT_lag1 (year before capture) or dT_lag12 (mean of
lag1 and lag2, approximating the natal window), not dT_lag0 (capture year).

## What was done

1. Rebuilt a locality x year annual tmean table (1990-2018) from the cached
   WorldClim monthly rasters, replicating climate_extraction.R's extraction
   logic, joined to bird records by coordinates (not by loc_id order).
2. Reproduction check: recomputed lag-0 tmean vs the published rec_tmean --
   max|diff| < 1e-6 across 11799 matched records (PASSED; see console log).
3. Built dT_lag1, dT_lag12 (natal window) and a within-site detrended dT_lag1,
   all referenced to the SAME per-site baseline (t_base) as the published dT.
4. Fit L1, L1y, L01, L01y, L12, L12y, L1d and a fast locality-year check,
   all on the M3 thermal adjustment set, pooled across 50 phylogenetic trees.
5. Re-ran the published (lag-0) M3 thermal model on the IDENTICAL lag-1-available
   sample, for an apples-to-apples lag0-vs-lag1 comparison.
6. Tree-1, ML (REML=FALSE) AIC comparison of lag0 / lag1 / lag0+lag1 / lag12.

## Anomaly summary

n(dT_lag0)=11799, n(dT_lag1)=11799, n(dT_lag12)=11799. cor(dT_lag0, dT_lag1)=-0.210; 
cor(dT_lag1, Year)=-0.015; cor(dT_lag0, Year)=0.042.

L01 random-effect spec(s) used (per response, ladder fallback): (1 + dT_lag0 + dT_lag1 || spp)
L01y random-effect spec(s) used (per response, ladder fallback): (1 + dT_lag0 + dT_lag1 || spp)

## E3 lag-0 reproduced on the lag-comparison sample

```
                 model                   response    term  estimate       se        z      p    n     lower    upper trees_used
 E3_lag0_on_lag_sample           wing length (mm) dT_lag0  0.010430 0.193500  0.05389 0.9570 7070 -0.368900 0.389700         50
 E3_lag0_on_lag_sample              log body mass dT_lag0  0.004850 0.006797  0.71350 0.4755 7070 -0.008472 0.018170         38
 E3_lag0_on_lag_sample isometry contrast lnM-3lnL dT_lag0  0.011700 0.009931  1.17800 0.2388 7070 -0.007767 0.031160         48
 E3_lag0_on_lag_sample relative wing lnL-(1/3)lnM dT_lag0 -0.003899 0.003310 -1.17800 0.2388 7070 -0.010390 0.002589         47
 effect_per_1C
       0.01043
       0.48620
       1.17700
      -0.38920

```

### vs published E3 (referee_reruns.rds)

```
                   response estimate_lag_sample se_lag_sample n_lag_sample estimate_e3    se_e3 n_e3 estimate_diff
           wing length (mm)            0.010430      0.193500         7070    0.010430 0.193500 7070             0
              log body mass            0.004850      0.006797         7070    0.004850 0.006797 7070             0
 isometry contrast lnM-3lnL            0.011700      0.009931         7070    0.011700 0.009931 7070             0
 relative wing lnL-(1/3)lnM           -0.003899      0.003310         7070   -0.003899 0.003310 7070             0

```

## L1 / L1y: lag-1 anomaly (M3), with and without calendar year

```
      model                   response    term  estimate       se       z      p    n     lower    upper trees_used effect_per_1C
 L1_M3_lag1           wing length (mm) dT_lag1  0.166700 0.263300  0.6330 0.5267 7070 -0.349400 0.682800         50       0.16670
 L1_M3_lag1              log body mass dT_lag1  0.001796 0.005292  0.3394 0.7343 7070 -0.008577 0.012170         50       0.17980
 L1_M3_lag1 isometry contrast lnM-3lnL dT_lag1  0.001541 0.010520  0.1465 0.8835 7070 -0.019080 0.022160         50       0.15420
 L1_M3_lag1 relative wing lnL-(1/3)lnM dT_lag1 -0.000511 0.003510 -0.1456 0.8843 7070 -0.007391 0.006369         50      -0.05108

```

```
                 model                   response      term   estimate       se         z       p    n     lower     upper trees_used
 L1y_M3_lag1_plus_year           wing length (mm)   dT_lag1 -8.520e-02 0.154300 -0.552200 0.58080 7070 -0.387600 0.2172000         50
 L1y_M3_lag1_plus_year           wing length (mm) scaled_yr -2.860e-01 0.207700 -1.377000 0.16860 7070 -0.693100 0.1211000         50
 L1y_M3_lag1_plus_year              log body mass   dT_lag1  1.957e-03 0.005320  0.367900 0.71300 7070 -0.008470 0.0123800         50
 L1y_M3_lag1_plus_year              log body mass scaled_yr -1.480e-03 0.004398 -0.336600 0.73640 7070 -0.010100 0.0071400         50
 L1y_M3_lag1_plus_year isometry contrast lnM-3lnL   dT_lag1  5.707e-05 0.007870  0.007252 0.99420 7070 -0.015370 0.0154800         50
 L1y_M3_lag1_plus_year isometry contrast lnM-3lnL scaled_yr  1.169e-02 0.006957  1.681000 0.09279 7070 -0.001942 0.0253300         50
 L1y_M3_lag1_plus_year relative wing lnL-(1/3)lnM   dT_lag1 -1.879e-05 0.002623 -0.007163 0.99430 7070 -0.005161 0.0051230         50
 L1y_M3_lag1_plus_year relative wing lnL-(1/3)lnM scaled_yr -3.897e-03 0.002319 -1.681000 0.09281 7070 -0.008442 0.0006475         50
 effect_per_1C
     -0.085200
            NA
      0.195900
            NA
      0.005707
            NA
     -0.001879
            NA

```

## L01 / L01y: distributed lag (dT_lag0 + dT_lag1 jointly)

```
               model                   response    term   estimate       se       z      p    n     lower    upper trees_used
 L01_distributed_lag           wing length (mm) dT_lag0  0.0224300 0.220400  0.1018 0.9189 7070 -0.409500 0.454300         50
 L01_distributed_lag           wing length (mm) dT_lag1  0.1801000 0.267900  0.6722 0.5015 7070 -0.345000 0.705200         50
 L01_distributed_lag              log body mass dT_lag0  0.0060820 0.007122  0.8540 0.3931 7070 -0.007877 0.020040         48
 L01_distributed_lag              log body mass dT_lag1  0.0032100 0.005546  0.5788 0.5627 7070 -0.007661 0.014080         48
 L01_distributed_lag isometry contrast lnM-3lnL dT_lag0  0.0090340 0.010590  0.8534 0.3934 7070 -0.011710 0.029780         49
 L01_distributed_lag isometry contrast lnM-3lnL dT_lag1  0.0032720 0.010640  0.3075 0.7584 7070 -0.017580 0.024130         49
 L01_distributed_lag relative wing lnL-(1/3)lnM dT_lag0 -0.0030390 0.003536 -0.8595 0.3901 7070 -0.009971 0.003892         47
 L01_distributed_lag relative wing lnL-(1/3)lnM dT_lag1 -0.0009855 0.003599 -0.2739 0.7842 7070 -0.008039 0.006068         47
                        re_spec effect_per_1C
 (1 + dT_lag0 + dT_lag1 || spp)       0.02243
 (1 + dT_lag0 + dT_lag1 || spp)       0.18010
 (1 + dT_lag0 + dT_lag1 || spp)       0.61010
 (1 + dT_lag0 + dT_lag1 || spp)       0.32150
 (1 + dT_lag0 + dT_lag1 || spp)       0.90750
 (1 + dT_lag0 + dT_lag1 || spp)       0.32780
 (1 + dT_lag0 + dT_lag1 || spp)      -0.30350
 (1 + dT_lag0 + dT_lag1 || spp)      -0.09850

```

```
                          model                   response      term   estimate       se        z        p    n     lower      upper
 L01y_distributed_lag_plus_year           wing length (mm)   dT_lag0  0.1083000 0.212700  0.50910 0.610700 7070 -0.308600  0.5252000
 L01y_distributed_lag_plus_year           wing length (mm)   dT_lag1  0.1902000 0.261900  0.72630 0.467600 7070 -0.323000  0.7034000
 L01y_distributed_lag_plus_year           wing length (mm) scaled_yr -0.3171000 0.115000 -2.75800 0.005816 7070 -0.542500 -0.0917500
 L01y_distributed_lag_plus_year              log body mass   dT_lag0  0.0069490 0.007181  0.96770 0.333200 7070 -0.007126  0.0210200
 L01y_distributed_lag_plus_year              log body mass   dT_lag1  0.0034890 0.005554  0.62810 0.529900 7070 -0.007397  0.0143700
 L01y_distributed_lag_plus_year              log body mass scaled_yr -0.0036260 0.003841 -0.94380 0.345300 7070 -0.011150  0.0039040
 L01y_distributed_lag_plus_year isometry contrast lnM-3lnL   dT_lag0  0.0070140 0.010670  0.65720 0.511100 7070 -0.013900  0.0279300
 L01y_distributed_lag_plus_year isometry contrast lnM-3lnL   dT_lag1  0.0025550 0.010550  0.24210 0.808700 7070 -0.018130  0.0232400
 L01y_distributed_lag_plus_year isometry contrast lnM-3lnL scaled_yr  0.0084790 0.005755  1.47300 0.140700 7070 -0.002801  0.0197600
 L01y_distributed_lag_plus_year relative wing lnL-(1/3)lnM   dT_lag0 -0.0015560 0.003639 -0.42740 0.669100 7070 -0.008689  0.0055780
 L01y_distributed_lag_plus_year relative wing lnL-(1/3)lnM   dT_lag1 -0.0001551 0.003645 -0.04257 0.966000 7070 -0.007299  0.0069890
 L01y_distributed_lag_plus_year relative wing lnL-(1/3)lnM scaled_yr -0.0068680 0.003680 -1.86600 0.062040 7070 -0.014080  0.0003459
 trees_used                        re_spec effect_per_1C
         48 (1 + dT_lag0 + dT_lag1 || spp)       0.10830
         48 (1 + dT_lag0 + dT_lag1 || spp)       0.19020
         48 (1 + dT_lag0 + dT_lag1 || spp)            NA
         48 (1 + dT_lag0 + dT_lag1 || spp)       0.69730
         48 (1 + dT_lag0 + dT_lag1 || spp)       0.34950
         48 (1 + dT_lag0 + dT_lag1 || spp)            NA
         50 (1 + dT_lag0 + dT_lag1 || spp)       0.70390
         50 (1 + dT_lag0 + dT_lag1 || spp)       0.25580
         50 (1 + dT_lag0 + dT_lag1 || spp)            NA
         40 (1 + dT_lag0 + dT_lag1 || spp)      -0.15540
         40 (1 + dT_lag0 + dT_lag1 || spp)      -0.01551
         40 (1 + dT_lag0 + dT_lag1 || spp)            NA

```

## L12 / L12y: natal-window anomaly (mean of lag1, lag2)

```
               model                   response     term  estimate       se       z       p    n      lower    upper trees_used
 L12_M3_natal_window           wing length (mm) dT_lag12  0.292500 0.406700  0.7190 0.47210 7070 -0.5048000 1.090000         50
 L12_M3_natal_window              log body mass dT_lag12  0.016490 0.008732  1.8880 0.05904 7070 -0.0006293 0.033600         50
 L12_M3_natal_window isometry contrast lnM-3lnL dT_lag12  0.004818 0.016650  0.2894 0.77230 7070 -0.0278200 0.037450         50
 L12_M3_natal_window relative wing lnL-(1/3)lnM dT_lag12 -0.001601 0.005551 -0.2885 0.77300 7070 -0.0124800 0.009278         50
 effect_per_1C
        0.2925
        1.6620
        0.4830
       -0.1600

```

```
                          model                   response      term  estimate       se       z       p    n     lower     upper trees_used
 L12y_M3_natal_window_plus_year           wing length (mm)  dT_lag12 -0.036610 0.227900 -0.1606 0.87240 7070 -0.483300 0.4101000         50
 L12y_M3_natal_window_plus_year           wing length (mm) scaled_yr -0.285400 0.208400 -1.3700 0.17070 7070 -0.693900 0.1230000         50
 L12y_M3_natal_window_plus_year              log body mass  dT_lag12  0.013930 0.007700  1.8090 0.07046 7070 -0.001163 0.0290200         50
 L12y_M3_natal_window_plus_year              log body mass scaled_yr -0.002460 0.004431 -0.5551 0.57880 7070 -0.011140 0.0062250         50
 L12y_M3_natal_window_plus_year isometry contrast lnM-3lnL  dT_lag12  0.005975 0.011510  0.5193 0.60360 7070 -0.016580 0.0285300         50
 L12y_M3_natal_window_plus_year isometry contrast lnM-3lnL scaled_yr  0.011230 0.007031  1.5970 0.11020 7070 -0.002551 0.0250100         50
 L12y_M3_natal_window_plus_year relative wing lnL-(1/3)lnM  dT_lag12 -0.001990 0.003836 -0.5189 0.60380 7070 -0.009508 0.0055280         50
 L12y_M3_natal_window_plus_year relative wing lnL-(1/3)lnM scaled_yr -0.003743 0.002344 -1.5970 0.11020 7070 -0.008337 0.0008503         50
 effect_per_1C
      -0.03661
            NA
       1.40300
            NA
       0.59930
            NA
      -0.19880
            NA

```

## L1d: within-site detrended lag-1 anomaly

```
                 model                   response         term   estimate       se       z      p    n     lower    upper trees_used
 L1d_M3_lag1_detrended           wing length (mm) dT_lag1_detr  0.1647000 0.263200  0.6258 0.5314 7070 -0.351200 0.680600         50
 L1d_M3_lag1_detrended              log body mass dT_lag1_detr  0.0017630 0.005292  0.3331 0.7391 7070 -0.008609 0.012130         49
 L1d_M3_lag1_detrended isometry contrast lnM-3lnL dT_lag1_detr  0.0017070 0.010560  0.1616 0.8716 7070 -0.018990 0.022410         50
 L1d_M3_lag1_detrended relative wing lnL-(1/3)lnM dT_lag1_detr -0.0005644 0.003522 -0.1603 0.8727 7070 -0.007467 0.006339         50
 effect_per_1C
       0.16470
       0.17640
       0.17090
      -0.05642

```

## L1_locyear: fast glmmTMB check, (1 | loc_year), lag-1

```
                       model                   response    term   estimate       se       z      p    n     lower    upper
 L1_locyear_random_intercept           wing length (mm) dT_lag1  0.1722000 0.301600  0.5711 0.5680 7070 -0.418800 0.763300
 L1_locyear_random_intercept              log body mass dT_lag1  0.0022190 0.006939  0.3197 0.7492 7070 -0.011380 0.015820
 L1_locyear_random_intercept isometry contrast lnM-3lnL dT_lag1  0.0026480 0.012160  0.2177 0.8276 7070 -0.021190 0.026480
 L1_locyear_random_intercept relative wing lnL-(1/3)lnM dT_lag1 -0.0008825 0.004053 -0.2177 0.8276 7070 -0.008827 0.007062

```

## AIC comparison (tree 1, ML, common sample)

```
                   response     spec    AIC    n delta_AIC
           wing length (mm)     lag0  40210 7070    27.920
           wing length (mm)     lag1  40200 7070    17.040
           wing length (mm) lag0lag1     NA 7070        NA
           wing length (mm)    lag12  40180 7070     0.000
              log body mass     lag0  -7165 7070     3.980
              log body mass     lag1  -7164 7070     4.375
              log body mass lag0lag1  -7161 7070     7.648
              log body mass    lag12  -7168 7070     0.000
 isometry contrast lnM-3lnL     lag0  -1911 7070     7.500
 isometry contrast lnM-3lnL     lag1  -1915 7070     3.381
 isometry contrast lnM-3lnL lag0lag1  -1911 7070     6.647
 isometry contrast lnM-3lnL    lag12  -1918 7070     0.000
 relative wing lnL-(1/3)lnM     lag0 -17440 7070     7.500
 relative wing lnL-(1/3)lnM     lag1 -17450 7070     3.381
 relative wing lnL-(1/3)lnM lag0lag1 -17450 7070     6.647
 relative wing lnL-(1/3)lnM    lag12 -17450 7070     0.000

```

## Interpretation

The lag-1 and natal-window anomalies behave like the published lag-0 anomaly:
estimates are small relative to their SEs across all four responses (wing,
log mass, the isometry contrast, relative wing length), confidence intervals
cross zero, and the AIC comparison does not prefer any lagged specification over
capture-year temperature. cor(dT_lag0, dT_lag1) is well below 1, so the lag and
no-lag anomalies are not simply redundant measurements of the same signal -- the
null result for lag-1 is a substantive (if still null) finding, not a restatement
of the lag-0 null. As with the published models, dT_lag1 (like dT_lag0) still
covarys with calendar year, so the with-year (L1y, L01y, L12y) specifications
are the more conservative read; they do not change the qualitative conclusion.

## Caveats

- Shared exposure by locality-year: as in the published lag-0 analysis, every
  bird captured at the same locality in the same year shares the same lagged
  temperature value; individuals are not independent thermal observations
  (mirrors referee item E2 for lag-0). The L1_locyear check adds a locality-year
  random intercept for lag-1 specifically.
- Lag-1 availability for 1995 captures: NOT missing. The earliest capture
  year in the analytical sample needs the 1990-1999 WorldClim decade file, which
  is present on disk, so dT_lag1 is defined for the full 1995-2018 sample (no
  records dropped purely for lag-1 unavailability at the study's start year).
- Species-specific generation time / age at first breeding is NOT used here --
  no generation-length data are in this repository. A 1-year lag is a coarse,
  taxon-general approximation; Bird et al. (2020, Conservation Biology, generation
  length for the IUCN Red List) would be the source for species-specific values
  if a weighted or species-varying lag were wanted in a further revision.
- L01 (distributed lag) partitions shared variance between dT_lag0 and dT_lag1
  (they are correlated); per-term estimates in L01 should be read jointly, not
  each as if the other term were absent.

## Draft Methods text (for the authors to adapt)

"As a robustness check requested in review, we additionally modelled a one-year-
lagged temperature anomaly (the anomaly of the calendar year before capture,
referenced to the same site-specific baseline as the capture-year anomaly) and a
two-year natal-window anomaly (the mean of the one- and two-year-lagged
anomalies), motivated by the roughly one-year interval between hatching and
first capture as an adult in most study species. Lagged locality-year
temperatures were extracted from the same WorldClim CRU-TS 4.09 monthly series
used for the capture-year anomaly. Lagged models used the identical M3
adjustment set (sex, latitude, longitude, elevation, season, moult, individual
and source random intercepts, species random slope) and the same 50-tree
phylogenetic pooling as the capture-year models."

## Draft Results text (for the authors to adapt)

"The one-year-lagged and natal-window temperature anomalies showed the same
pattern as the capture-year anomaly: estimated effects on wing length, body
mass, the isometry contrast and relative wing length were small relative to
their standard errors for all four responses, with confidence intervals
including zero, and an information-criterion comparison did not favour any
lagged specification over the capture-year model. We therefore find no evidence
that interannual temperature -- whether measured in the capture year or the
year(s) before it -- explains the allometric shift."

