# ATLANTIC ANTS litter-ant richness index: year effects

Generated 2026-09-09 11:55:15. Campaign events = 855 (393 Winkler, 462 pitfall) from 61 contributor files, 293 ~11-km localities, 570 sites (55 resampled in >1 year); years 1994-2018. Family: nbinom2 (Poisson dispersion 1.66). Year is per SD of the bird model (5.020 yr); per-decade columns multiply by 1.992. Plot-date sensitivity: 3955 events (42 % with richness = 1).

| Level | Model | Term | N events | Estimate (per SD yr) | 95% CI | z | % per decade [95% CI] |
|---|---|---|---|---|---|---|---|
| campaign | naive: no random effects | scaled_year | 855 | -0.269 | [-0.330, -0.208] | -8.59 | -41.5 [-48.2, -33.9] |
| campaign | pooled: contributor + locality + site RE, effort covariate | scaled_year | 855 | -0.113 | [-0.187, -0.039] | -2.98 | -20.1 [-31.1, -7.4] |
| campaign | pooled: effort as offset | scaled_year | 855 | 0.041 | [-0.051, 0.133] | 0.87 | 8.5 [-9.7, 30.3] |
| campaign | Mundlak by contributor | year_within | 855 | -0.126 | [-0.203, -0.050] | -3.23 | -22.2 [-33.2, -9.4] |
| campaign | Mundlak by contributor | year_between | 855 | 0.060 | [-0.203, 0.322] | 0.45 | 12.7 [-33.2, 90.1] |
| campaign | Mundlak by site | year_within_site | 855 | -0.410 | [-0.631, -0.188] | -3.62 | -55.8 [-71.6, -31.2] |
| campaign | Mundlak by site | year_site_mean | 855 | -0.082 | [-0.158, -0.005] | -2.10 | -15.1 [-27.1, -1.1] |
| campaign | monitoring programmes (>= 6 yr): pooled | scaled_year | 251 | -0.187 | [-0.289, -0.084] | -3.57 | -31.1 [-43.8, -15.4] |
| campaign | monitoring programmes (>= 6 yr): Mundlak by contributor | year_within | 251 | -0.166 | [-0.268, -0.064] | -3.19 | -28.1 [-41.3, -12.0] |
| campaign | monitoring programmes (>= 6 yr): Mundlak by contributor | year_between | 251 | -0.453 | [-0.711, -0.195] | -3.45 | -59.5 [-75.7, -32.3] |
| campaign | without literature files: pooled | scaled_year | 704 | -0.150 | [-0.234, -0.066] | -3.51 | -25.8 [-37.2, -12.3] |
| campaign | without literature files: Mundlak by contributor | year_within | 704 | -0.166 | [-0.253, -0.079] | -3.74 | -28.2 [-39.6, -14.6] |
| campaign | without literature files: Mundlak by contributor | year_between | 704 | 0.039 | [-0.243, 0.322] | 0.27 | 8.2 [-38.4, 90.0] |
| campaign | constant protocol (modal effort per series): pooled | scaled_year | 655 | -0.229 | [-0.361, -0.096] | -3.39 | -36.6 [-51.3, -17.5] |
| campaign | constant protocol (modal effort per series): Mundlak by contributor | year_within | 655 | -0.290 | [-0.439, -0.142] | -3.84 | -43.9 [-58.3, -24.7] |
| campaign | constant protocol (modal effort per series): Mundlak by contributor | year_between | 655 | -0.029 | [-0.284, 0.227] | -0.22 | -5.6 [-43.3, 57.1] |
| campaign | abundance (secondary): pooled, offset | scaled_year | 67 | -0.142 | [-0.551, 0.266] | -0.68 | -24.7 [-66.6, 70.0] |
| campaign | abundance (secondary): Mundlak by contributor | year_within | 67 | -0.282 | [-0.724, 0.161] | -1.25 | -42.9 [-76.4, 37.9] |
| campaign | abundance (secondary): Mundlak by contributor | year_between | 67 | 0.338 | [-0.390, 1.066] | 0.91 | 95.9 [-54.1, 735.4] |
| plot | naive: no random effects | scaled_year | 3955 | -0.721 | [-0.772, -0.671] | -28.16 | -76.2 [-78.5, -73.7] |
| plot | pooled: contributor + locality + site RE, effort covariate | scaled_year | 3955 | -0.263 | [-0.322, -0.203] | -8.64 | -40.7 [-47.4, -33.3] |
| plot | pooled: effort as offset | scaled_year | 3955 | 0.465 | [0.400, 0.531] | 13.92 | 152.7 [121.8, 187.9] |
| plot | Mundlak by contributor | year_within | 3955 | -0.262 | [-0.323, -0.202] | -8.53 | -40.7 [-47.4, -33.1] |
| plot | Mundlak by contributor | year_between | 3955 | -0.273 | [-0.656, 0.111] | -1.40 | -41.9 [-73.0, 24.6] |
| plot | Mundlak by site | year_within_site | 3955 | -0.269 | [-0.512, -0.025] | -2.16 | -41.5 [-64.0, -4.9] |
| plot | Mundlak by site | year_site_mean | 3955 | -0.262 | [-0.324, -0.201] | -8.37 | -40.7 [-47.5, -33.0] |

Mundlak within - between (contributor, campaign level): -0.186 (SE 0.139, z = -1.34, p = 0.179).

## Monitoring programmes (>= 6 sampling years, literature files excluded): independent per-programme slopes

| Contributor file | N events | Years | Methods | Slope per SD yr (SE) | % per decade [95% CI] | Random-slope conditional |
|---|---|---|---|---|---|---|
| EQUIPE_LAMAT_UMC | 80 | 12 (2001-2017) | Pitfall/Winkler | 0.097 (0.172) | 21.2 [-38.1, 137.1] | -0.187 |
| SERVIO_RIBEIRO_PERD_2001_2017 | 53 | 9 (2001-2017) | Pitfall/Winkler | -0.385 (0.111) | -53.6 [-69.9, -28.4] | -0.187 |
| ELMO_and_DELABIE | 65 | 6 (1996-2002) | Winkler | -0.026 (0.062) | -5.1 [-25.4, 20.9] | -0.187 |
| IVAN_NASCIMENTO_UESC_CEPLAC_UESB | 27 | 6 (2003-2014) | Pitfall/Winkler | -0.475 (0.069) | -61.2 [-70.3, -49.2] | -0.187 |
| ROGERIO_ROSA_SILVA_BIOTA_FORMIGAS | 26 | 6 (1997-2003) | Winkler | 0.284 (0.134) | 76.0 [4.4, 196.7] | -0.187 |

## Leave-one-contributor-out (contributors with >= 20 events): Mundlak within slope

| Dropped | Events dropped | Within slope (SE) | % per decade |
|---|---|---|---|
| LITERATURE_2020 | 31 | -0.151 (0.038) | -26.0 |
| EQUIPE_LAMAT_UMC | 80 | -0.136 (0.039) | -23.7 |
| ELMO_and_DELABIE | 65 | -0.134 (0.042) | -23.4 |
| ROGERIO_ROSA_SILVA_BIOTA_FORMIGAS | 26 | -0.133 (0.040) | -23.3 |
| IVAN_NASCIMENTO_UESC_CEPLAC_UESB | 27 | -0.132 (0.041) | -23.1 |
| SERVIO_RIBEIRO_PERD_2001_2017_2018_PERD_Recente | 25 | -0.130 (0.039) | -22.8 |
| FREDERICO_NEVES_LEI_UFMG_PARNASC | 40 | -0.129 (0.040) | -22.6 |
| FLAVIO_RAMOS | 27 | -0.128 (0.039) | -22.4 |
| LASMAR_CJ_ETAL | 22 | -0.122 (0.038) | -21.6 |
| MAULYSSEA_LPPRADO | 43 | -0.121 (0.040) | -21.4 |
| LITERATURE_TEAM_PUBLISHED_DATA | 22 | -0.119 (0.040) | -21.2 |
| ROGERIO_ROSA_DA_SILVA_LITERATURE_DATA | 34 | -0.117 (0.041) | -20.8 |
| KOCH_SANTOS_and_DELABIE | 33 | -0.117 (0.039) | -20.7 |
| SERVIO_RIBEIRO_PERD_2001_2017 | 53 | -0.055 (0.047) | -10.5 |
