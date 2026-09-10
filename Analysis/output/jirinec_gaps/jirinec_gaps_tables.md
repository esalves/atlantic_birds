# Jirinec-gap analysis - numbers

Generated: 2026-09-10 11:42:01  | trees: 0  | runtime: 13.53 min

## Climate coverage

| n_records_total | n_records_climate | pct | n_mass_climate | n_wing_climate | n_spp_climate | n_loc_climate |
|---|---|---|---|---|---|---|
| 12571.0 | 11799.0 | 93.9 | 10874.0 | 7810.0 | 73.0 | 357.0 |

## G1 VIF (their max was 2.01)

| set | term | VIF |
|---|---|---|
| temperature | scaled_yr   | 2.293       |
| temperature | z_t_an_l0   | 1.497       |
| temperature | z_t_an_l1   | 1.356       |
| temperature | z_t_an_l2   | 1.347       |
| temperature | z_tmean_loc | 3.478       |
| temperature | scaled_lat  | 4.22        |
| precipitation | scaled_yr     | 1.109         |
| precipitation | z_p_an_l0     | 1.209         |
| precipitation | z_p_an_l1     | 1.218         |
| precipitation | z_p_an_l2     | 1.252         |
| precipitation | z_prec_loc    | 1.336         |
| precipitation | scaled_lat    | 1.443         |
| spei      | scaled_yr | 1.107     |
| spei      | z_s_an_l0 | 1.077     |
| spei      | z_s_an_l1 | 1.127     |
| spei      | z_s_an_l2 | 1.116     |
| spei       | z_spei_loc | 1.083      |
| spei       | scaled_lat | 1.223      |

## G1 year-slope correlation with climate anomalies

| term | r |
|---|---|
| r_yr_tanom_l0 | 0.402         |
| r_yr_panom_l0 | -0.024        |
| r_yr_sanom_l0 | -0.068        |

## G1 year slope, with and without lagged climate

| trait | model | estimate | lower | upper | pct_per_decade | pct_per_decade_lower | pct_per_decade_upper | n |
|---|---|---|---|---|---|---|---|---|
| ln_body_mass | M1_year_only | -0.004       | -0.015       | 0.006        | -0.862       | -2.869       | 1.186        | 10874        |
| ln_body_mass   | M1_temperature | -0.009         | -0.019         | 0.002          | -1.679         | -3.775         | 0.463          | 10874          |
| ln_body_mass     | M1_precipitation | -0.005           | -0.015           | 0.005            | -1.009           | -3.03            | 1.054            | 10874            |
| ln_body_mass | M1_spei      | -0.005       | -0.015       | 0.006        | -0.919       | -2.939       | 1.142        | 10874        |
| ln_body_mass | M3_year_only | -0.003       | -0.012       | 0.006        | -0.551       | -2.301       | 1.229        | 10711        |
| ln_body_mass   | M3_temperature | -0.007         | -0.017         | 0.002          | -1.418         | -3.269         | 0.469          | 10711          |
| ln_body_mass     | M3_precipitation | -0.004           | -0.013           | 0.005            | -0.758           | -2.531           | 1.047            | 10711            |
| ln_body_mass | M3_spei      | -0.003       | -0.012       | 0.006        | -0.596       | -2.363       | 1.203        | 10711        |
| wing         | M1_year_only | -0.601       | -1.019       | -0.182       | -1.667       | -2.828       | -0.505       | 7810         |
| wing           | M1_temperature | -0.509         | -0.94          | -0.078         | -1.413         | -2.609         | -0.218         | 7810           |
| wing             | M1_precipitation | -0.598           | -1.016           | -0.18            | -1.66            | -2.819           | -0.5             | 7810             |
| wing    | M1_spei | -0.626  | -1.046  | -0.205  | -1.736  | -2.904  | -0.569  | 7810    |
| wing         | M3_year_only | -0.354       | -0.748       | 0.041        | -0.981       | -2.076       | 0.114        | 7651         |
| wing           | M3_temperature | -0.283         | -0.691         | 0.125          | -0.786         | -1.918         | 0.346          | 7651           |
| wing             | M3_precipitation | -0.358           | -0.752           | 0.037            | -0.993           | -2.087           | 0.102            | 7651             |
| wing    | M3_spei | -0.381  | -0.779  | 0.018   | -1.056  | -2.16   | 0.049   | 7651    |

## G1 AIC across climate variants (ML refits)

| trait | tier | spec | AIC_ML | n | dAIC |
|---|---|---|---|---|---|
| ln_body_mass | M1           | temperature  | -11568.92    | 10874        | 0            |
| ln_body_mass | M1           | year_only    | -11568.81    | 10874        | 0.11         |
| ln_body_mass  | M1            | precipitation | -11566.68     | 10874         | 2.24          |
| ln_body_mass | M1           | spei         | -11564.01    | 10874        | 4.91         |
| ln_body_mass | M3           | temperature  | -11839.07    | 10711        | 0            |
| ln_body_mass | M3           | year_only    | -11838.72    | 10711        | 0.35         |
| ln_body_mass  | M3            | precipitation | -11836.24     | 10711         | 2.83          |
| ln_body_mass | M3           | spei         | -11835.68    | 10711        | 3.39         |
| wing     | M1       | spei     | 45379.01 | 7810     | 0        |
| wing      | M1        | year_only | 45383.26  | 7810      | 4.26      |
| wing          | M1            | precipitation | 45383.3       | 7810          | 4.3           |
| wing        | M1          | temperature | 45386.5     | 7810        | 7.49        |
| wing     | M3       | spei     | 44338.28 | 7651     | 0        |
| wing      | M3        | year_only | 44340.4   | 7651      | 2.12      |
| wing          | M3            | precipitation | 44342.48      | 7651          | 4.2           |
| wing        | M3          | temperature | 44345.06    | 7651        | 6.78        |

## G1 climate coefficients (all tiers)

| trait | model | term | estimate | se | lower | upper | t |
|---|---|---|---|---|---|---|---|
| ln_body_mass   | M1_temperature | z_t_an_l0      | 0.00178        | 0.00216        | -0.00245       | 0.006          | 0.82517        |
| ln_body_mass   | M1_temperature | z_t_an_l1      | 0.0042         | 0.00197        | 0.00035        | 0.00805        | 2.13633        |
| ln_body_mass   | M1_temperature | z_t_an_l2      | 0.0049         | 0.00209        | 8e-04          | 0.009          | 2.34295        |
| ln_body_mass   | M1_temperature | z_tmean_loc    | -0.00692       | 0.00764        | -0.02189       | 0.00805        | -0.90595       |
| ln_body_mass     | M1_precipitation | z_p_an_l0        | -8e-04           | 0.0018           | -0.00432         | 0.00272          | -0.44598         |
| ln_body_mass     | M1_precipitation | z_p_an_l1        | -0.0034          | 0.00205          | -0.00741         | 0.00061          | -1.66134         |
| ln_body_mass     | M1_precipitation | z_p_an_l2        | -0.00073         | 0.00181          | -0.00428         | 0.00281          | -0.40506         |
| ln_body_mass     | M1_precipitation | z_prec_loc       | 0.01105          | 0.00644          | -0.00158         | 0.02368          | 1.71503          |
| ln_body_mass | M1_spei      | z_s_an_l0    | -0.00176     | 0.00194      | -0.00557     | 0.00205      | -0.90563     |
| ln_body_mass | M1_spei      | z_s_an_l1    | -0.00176     | 0.00182      | -0.00532     | 0.0018       | -0.96856     |
| ln_body_mass | M1_spei      | z_s_an_l2    | 0.00069      | 0.0019       | -0.00304     | 0.00441      | 0.36191      |
| ln_body_mass | M1_spei      | z_spei_loc   | 0.008        | 0.00783      | -0.00734     | 0.02335      | 1.02222      |
| ln_body_mass   | M3_temperature | z_t_an_l0      | 0.00242        | 0.00213        | -0.00175       | 0.00659        | 1.13662        |
| ln_body_mass   | M3_temperature | z_t_an_l1      | 0.00441        | 0.00195        | 0.00059        | 0.00823        | 2.26339        |
| ln_body_mass   | M3_temperature | z_t_an_l2      | 0.00499        | 0.00208        | 9e-04          | 0.00907        | 2.393          |
| ln_body_mass   | M3_temperature | z_tmean_loc    | -0.01352       | 0.01663        | -0.04611       | 0.01907        | -0.81327       |
| ln_body_mass     | M3_precipitation | z_p_an_l0        | -0.00076         | 0.00177          | -0.00423         | 0.00271          | -0.43024         |
| ln_body_mass     | M3_precipitation | z_p_an_l1        | -0.00389         | 0.00203          | -0.00786         | 9e-05            | -1.91737         |
| ln_body_mass     | M3_precipitation | z_p_an_l2        | -0.0016          | 0.00175          | -0.00502         | 0.00182          | -0.91536         |
| ln_body_mass     | M3_precipitation | z_prec_loc       | 0.00805          | 0.0066           | -0.00489         | 0.021            | 1.21941          |
| ln_body_mass | M3_spei      | z_s_an_l0    | -0.00148     | 0.00194      | -0.00528     | 0.00231      | -0.76637     |
| ln_body_mass | M3_spei      | z_s_an_l1    | -0.00266     | 0.00179      | -0.00616     | 0.00085      | -1.48693     |
| ln_body_mass | M3_spei      | z_s_an_l2    | 0.00041      | 0.00184      | -0.00321     | 0.00402      | 0.22051      |
| ln_body_mass | M3_spei      | z_spei_loc   | 0.01044      | 0.00791      | -0.00506     | 0.02595      | 1.32003      |
| wing           | M1_temperature | z_t_an_l0      | -0.10026       | 0.07479        | -0.24684       | 0.04632        | -1.34062       |
| wing           | M1_temperature | z_t_an_l1      | -0.11492       | 0.06803        | -0.24826       | 0.01842        | -1.68923       |
| wing           | M1_temperature | z_t_an_l2      | -0.1025        | 0.07409        | -0.24771       | 0.04271        | -1.38348       |
| wing           | M1_temperature | z_tmean_loc    | -0.24454       | 0.27712        | -0.78769       | 0.29861        | -0.88245       |
| wing             | M1_precipitation | z_p_an_l0        | -0.01839         | 0.06303          | -0.14193         | 0.10514          | -0.29181         |
| wing             | M1_precipitation | z_p_an_l1        | 0.14751          | 0.07499          | 0.00053          | 0.29449          | 1.96701          |
| wing             | M1_precipitation | z_p_an_l2        | 0.16733          | 0.07406          | 0.02217          | 0.31249          | 2.2593           |
| wing             | M1_precipitation | z_prec_loc       | 0.22198          | 0.22007          | -0.20936         | 0.65332          | 1.00868          |
| wing      | M1_spei   | z_s_an_l0 | -0.06487  | 0.06923   | -0.20056  | 0.07081   | -0.93709  |
| wing      | M1_spei   | z_s_an_l1 | 0.21392   | 0.06593   | 0.08471   | 0.34314   | 3.24489   |
| wing      | M1_spei   | z_s_an_l2 | 0.04937   | 0.07006   | -0.08794  | 0.18668   | 0.70469   |
| wing       | M1_spei    | z_spei_loc | -0.04858   | 0.24937    | -0.53735   | 0.44018    | -0.19482   |
| wing           | M3_temperature | z_t_an_l0      | -0.05195       | 0.07579        | -0.2005        | 0.09661        | -0.68538       |
| wing           | M3_temperature | z_t_an_l1      | -0.11346       | 0.06833        | -0.24738       | 0.02047        | -1.66037       |
| wing           | M3_temperature | z_t_an_l2      | -0.08274       | 0.07515        | -0.23004       | 0.06456        | -1.101         |
| wing           | M3_temperature | z_tmean_loc    | -0.31093       | 0.56283        | -1.41408       | 0.79221        | -0.55245       |
| wing             | M3_precipitation | z_p_an_l0        | -0.01943         | 0.06384          | -0.14456         | 0.10569          | -0.30442         |
| wing             | M3_precipitation | z_p_an_l1        | 0.13234          | 0.07582          | -0.01626         | 0.28095          | 1.74554          |
| wing             | M3_precipitation | z_p_an_l2        | 0.12808          | 0.07501          | -0.01895         | 0.2751           | 1.70735          |
| wing             | M3_precipitation | z_prec_loc       | 0.2666           | 0.23917          | -0.20217         | 0.73536          | 1.1147           |
| wing      | M3_spei   | z_s_an_l0 | -0.08321  | 0.07085   | -0.22208  | 0.05566   | -1.17439  |
| wing      | M3_spei   | z_s_an_l1 | 0.18706   | 0.06671   | 0.05631   | 0.31782   | 2.80411   |
| wing      | M3_spei   | z_s_an_l2 | 0.05031   | 0.07068   | -0.08823  | 0.18884   | 0.71178   |
| wing       | M3_spei    | z_spei_loc | -0.10144   | 0.27691    | -0.64419   | 0.44131    | -0.36633   |

## G2 species-slope tallies

| trait | model | n_spp | n_negative_mean | n_positive_mean | n_CI_negative | n_CI_positive | n_CI_spans_zero | pct_CI_negative | median_pct_per_decade |
|---|---|---|---|---|---|---|---|---|---|
| ln_body_mass | M0           | 73           | 41           | 32           | 7            | 2            | 64           | 9.6          | -0.682       |
| ln_body_mass | M3           | 73           | 45           | 28           | 6            | 2            | 65           | 8.2          | -0.952       |
| wing   | M0     | 72     | 56     | 16     | 17     | 2      | 53     | 23.6   | -1.662 |
| wing   | M3     | 72     | 48     | 24     | 9      | 2      | 61     | 12.5   | -0.727 |

## G3 mass:wing ratio, pooled slope

| trait | model | estimate | lower | upper | pct_per_decade | pct_per_decade_lower | pct_per_decade_upper | n |
|---|---|---|---|---|---|---|---|---|
| ln_ratio | M0       | 0.0096   | 0.0019   | 0.0173   | 1.9343   | 0.3863   | 3.5062   | 7577     |
| ln_ratio | M1       | 0.0047   | -0.0035  | 0.0128   | 0.9337   | -0.6861  | 2.5799   | 7577     |
| ln_ratio | M3       | 0.0019   | -0.0064  | 0.0103   | 0.3832   | -1.2708  | 2.0649   | 7409     |
| raw_ratio | M0        | 0.0033    | 0.0012    | 0.0054    | 2.1215    | 0.7555    | 3.4875    | 7577      |
| raw_ratio | M1        | 9e-04     | -0.0016   | 0.0034    | 0.6031    | -1.0132   | 2.2195    | 7577      |
| raw_ratio | M3        | 2e-04     | -0.0024   | 0.0028    | 0.1298    | -1.5453   | 1.8049    | 7409      |

## G3 mass:wing species tally

| trait | model | n_spp | n_negative_mean | n_CI_negative | n_CI_positive | n_CI_spans_zero | pct_CI_negative |
|---|---|---|---|---|---|---|---|
| ln_ratio | M3       | 72       | 20       | 0        | 1        | 71       | 0        |
| raw_ratio | M3        | 72        | 15        | 0         | 0         | 72        | 0         |

## G4 phylogenetic signal in slopes

| metric | n_spp | moran_I | moran_expected | moran_p_median | moran_p_max | lambda | lambda_lo | lambda_hi | lambda_p_median |
|---|---|---|---|---|---|---|---|---|---|
| ln_body_mass (M0) | 73                | 0.0262            | -0.0139           | 0.0114            | 0.0536            | 0.0879            | 0.0879            | 0.0879            | 0.1433            |
| ln_body_mass (M3) | 73                | 2e-04             | -0.0139           | 0.3682            | 0.5746            | 0.0344            | 0.0344            | 0.0344            | 0.5198            |
| wing (M0) | 72        | -0.0107   | -0.0141   | 0.8394    | 0.9956    | 0.006     | 0.006     | 0.006     | 0.8889    |
| wing (M3) | 72        | -0.0185   | -0.0141   | 0.7701    | 0.9631    | 0         | 0         | 0         | 1         |
| ln(mass:wing) (M3) | 72                 | -0.0266            | -0.0141            | 0.446              | 0.6211             | 0                  | 0                  | 0                  | 1                  |

## G5 foraging-stratum interactions

| trait | model | term | estimate | se | lower | upper | t | n |
|---|---|---|---|---|---|---|---|---|
| ln_body_mass  | M1_continuous | scaled_yr     | -0.00616      | 0.00598       | -0.01787      | 0.00556       | -1.02997      | 10779         |
| ln_body_mass  | M1_continuous | au_std        | -0.08949      | 0.05882       | -0.20478      | 0.02581       | -1.5213       | 10779         |
| ln_body_mass     | M1_continuous    | scaled_yr:au_std | -0.00145         | 0.00471          | -0.01067         | 0.00778          | -0.307           | 10779            |
| ln_body_mass   | M1_categorical | scaled_yr      | -0.00393       | 0.00699        | -0.01764       | 0.00978        | -0.56219       | 10779          |
| ln_body_mass   | M1_categorical | stratumhigh    | -0.09453       | 0.14511        | -0.37895       | 0.18989        | -0.65142       | 10779          |
| ln_body_mass          | M1_categorical        | scaled_yr:stratumhigh | -0.00727              | 0.01119               | -0.0292               | 0.01466               | -0.65                 | 10779                 |
| ln_body_mass  | M3_continuous | scaled_yr     | -0.00473      | 0.00504       | -0.0146       | 0.00514       | -0.9386       | 10602         |
| ln_body_mass  | M3_continuous | au_std        | -0.08546      | 0.05865       | -0.20041      | 0.02949       | -1.45716      | 10602         |
| ln_body_mass     | M3_continuous    | scaled_yr:au_std | -0.00352         | 0.00375          | -0.01086         | 0.00382          | -0.94049         | 10602            |
| ln_body_mass   | M3_categorical | scaled_yr      | -0.0016        | 0.00575        | -0.01288       | 0.00968        | -0.27872       | 10602          |
| ln_body_mass   | M3_categorical | stratumhigh    | -0.08861       | 0.14454        | -0.37191       | 0.19468        | -0.61307       | 10602          |
| ln_body_mass          | M3_categorical        | scaled_yr:stratumhigh | -0.01088              | 0.00875               | -0.02804              | 0.00628               | -1.24308              | 10602                 |
| wing          | M1_continuous | scaled_yr     | -0.50723      | 0.2065        | -0.91196      | -0.1025       | -2.45637      | 8086          |
| wing          | M1_continuous | au_std        | -0.5105       | 1.79214       | -4.02309      | 3.00208       | -0.28486      | 8086          |
| wing             | M1_continuous    | scaled_yr:au_std | 0.17248          | 0.16388          | -0.14873         | 0.49369          | 1.05248          | 8086             |
| wing           | M1_categorical | scaled_yr      | -0.53015       | 0.24447        | -1.00931       | -0.05098       | -2.16855       | 8086           |
| wing           | M1_categorical | stratumhigh    | 0.54756        | 4.32563        | -7.93068       | 9.02579        | 0.12658        | 8086           |
| wing                  | M1_categorical        | scaled_yr:stratumhigh | 0.13784               | 0.39315               | -0.63273              | 0.90842               | 0.35061               | 8086                  |
| wing          | M3_continuous | scaled_yr     | -0.26807      | 0.18603       | -0.6327       | 0.09655       | -1.441        | 7892          |
| wing          | M3_continuous | au_std        | -0.52176      | 1.80608       | -4.06167      | 3.01815       | -0.28889      | 7892          |
| wing             | M3_continuous    | scaled_yr:au_std | 0.06911          | 0.14087          | -0.207           | 0.34521          | 0.49056          | 7892             |
| wing           | M3_categorical | scaled_yr      | -0.24716       | 0.21555        | -0.66963       | 0.17531        | -1.14666       | 7892           |
| wing           | M3_categorical | stratumhigh    | 0.55422        | 4.35941        | -7.99021       | 9.09866        | 0.12713        | 7892           |
| wing                  | M3_categorical        | scaled_yr:stratumhigh | -0.03198              | 0.33165               | -0.68201              | 0.61806               | -0.09641              | 7892                  |
| ln_ratio      | M1_continuous | scaled_yr     | 0.00458       | 0.00425       | -0.00375      | 0.01292       | 1.07751       | 7217          |
| ln_ratio      | M1_continuous | au_std        | -0.06705      | 0.03703       | -0.13964      | 0.00554       | -1.81053      | 7217          |
| ln_ratio         | M1_continuous    | scaled_yr:au_std | -0.0025          | 0.00264          | -0.00768         | 0.00267          | -0.94824         | 7217             |
| ln_ratio       | M1_categorical | scaled_yr      | 0.0056         | 0.00474        | -0.00368       | 0.01489        | 1.18231        | 7217           |
| ln_ratio       | M1_categorical | stratumhigh    | -0.07931       | 0.09104        | -0.25776       | 0.09913        | -0.87114       | 7217           |
| ln_ratio              | M1_categorical        | scaled_yr:stratumhigh | -0.00428              | 0.0063                | -0.01662              | 0.00807               | -0.67915              | 7217                  |
| ln_ratio      | M3_continuous | scaled_yr     | 0.00169       | 0.00422       | -0.00659      | 0.00996       | 0.39912       | 7051          |
| ln_ratio      | M3_continuous | au_std        | -0.0659       | 0.03702       | -0.13846      | 0.00665       | -1.78024      | 7051          |
| ln_ratio         | M3_continuous    | scaled_yr:au_std | -0.0051          | 0.00247          | -0.00994         | -0.00025         | -2.06181         | 7051             |
| ln_ratio       | M3_categorical | scaled_yr      | 0.00372        | 0.00473        | -0.00554       | 0.01299        | 0.78737        | 7051           |
| ln_ratio       | M3_categorical | stratumhigh    | -0.07736       | 0.09095        | -0.25563       | 0.10091        | -0.85057       | 7051           |
| ln_ratio              | M3_categorical        | scaled_yr:stratumhigh | -0.00836              | 0.0061                | -0.02031              | 0.0036                | -1.37018              | 7051                  |

## G5 species per dominant stratum

| dom | n |
|---|---|
| understory | 23         |
| ground | 21     |
| midhigh | 17      |
| canopy | 5      |

## G6 design-matched restricted tier

Criteria: span >= 8 yr, n >= 150; kept 8 contributors

| n_records | n_spp | n_mass | n_wing | n_src | yr_min | yr_max |
|---|---|---|---|---|---|---|
| 4114 | 73 | 3961 | 3497 | 8 | 1995 | 2017 |

| trait | model | estimate | lower | upper | pct_per_decade | pct_per_decade_lower | pct_per_decade_upper | n |
|---|---|---|---|---|---|---|---|---|
| ln_body_mass  | restricted_M0 | -0.0075       | -0.0184       | 0.0034        | -1.4818       | -3.6031       | 0.6862        | 3961          |
| ln_body_mass  | restricted_M1 | -0.0074       | -0.0196       | 0.0048        | -1.4582       | -3.8217       | 0.9634        | 3961          |
| ln_body_mass  | restricted_M3 | -0.0066       | -0.0154       | 0.0022        | -1.3037       | -3.0146       | 0.4373        | 3834          |
| wing          | restricted_M0 | -0.6191       | -1.3481       | 0.1099        | -1.7267       | -3.7598       | 0.3065        | 3497          |
| wing          | restricted_M1 | -0.1746       | -0.8155       | 0.4662        | -0.4871       | -2.2744       | 1.3002        | 3497          |
| wing          | restricted_M3 | -0.1056       | -0.5007       | 0.2896        | -0.2944       | -1.3965       | 0.8077        | 3377          |
| ln_ratio      | restricted_M0 | 0.0048        | -0.0031       | 0.0127        | 0.9656        | -0.6104       | 2.5665        | 3373          |
| ln_ratio      | restricted_M1 | 2e-04         | -0.0086       | 0.009         | 0.0405        | -1.7032       | 1.8151        | 3373          |
| ln_ratio      | restricted_M3 | -6e-04        | -0.0096       | 0.0085        | -0.1128       | -1.9026       | 1.7096        | 3256          |

## G6 contributor coverage (top 15 by span)

| src | n | n_mass | n_wing | yr_min | yr_max | span | n_years | n_site | n_spp |
|---|---|---|---|---|---|---|---|---|---|
| 31 | 539 | 533 | 433 | 1995 | 2017 | 22 | 23 | 28 | 41 |
| 15 | 1101 | 1075 | 1040 | 1996 | 2017 | 21 | 21 | 5 | 42 |
| 11 | 418 | 407 | 388 | 1999 | 2017 | 18 | 9 | 4 | 30 |
| 4 | 353 | 349 | 295 | 2001 | 2017 | 16 | 8 | 4 | 21 |
| 5 | 762 | 752 | 618 | 2001 | 2016 | 15 | 11 | 8 | 45 |
| 29 | 57 | 55 | 54 | 1998 | 2013 | 15 | 5 | 4 | 12 |
| 22 | 334 | 320 | 311 | 2008 | 2017 | 9 | 6 | 2 | 22 |
| 1 | 270 | 196 | 243 | 2000 | 2009 | 9 | 8 | 14 | 19 |
| 37 | 337 | 329 | 169 | 2009 | 2017 | 8 | 8 | 4 | 23 |
| 17 | 61 | 61 | 60 | 2007 | 2015 | 8 | 5 | 5 | 19 |
| 40 | 382 | 379 | 376 | 2011 | 2018 | 7 | 7 | 2 | 20 |
| 9 | 413 | 3 | 332 | 2009 | 2015 | 6 | 7 | 1 | 23 |
| 34 | 236 | 218 | 226 | 1999 | 2005 | 6 | 6 | 9 | 19 |
| 28 | 132 | 132 | 0 | 2011 | 2017 | 6 | 4 | 5 | 19 |
| 38 | 243 | 199 | 181 | 2010 | 2015 | 5 | 6 | 3 | 17 |

## G7 symmetric-coverage filter

| median_year | n_spp_total | n_spp_pass_mass | n_spp_fail_mass | n_spp_pass_wing | n_spp_fail_wing |
|---|---|---|---|---|---|
| 2010 | 73 | 67 | 6 | 62 | 11 |

| trait | model | estimate | lower | upper | pct_per_decade | pct_per_decade_lower | pct_per_decade_upper | n |
|---|---|---|---|---|---|---|---|---|
| ln_body_mass | symmetric_M0 | -0.0067      | -0.0153      | 0.002        | -1.3177      | -3.0041      | 0.3979       | 11086        |
| ln_body_mass | symmetric_M1 | -0.0066      | -0.0167      | 0.0034       | -1.3128      | -3.2638      | 0.6776       | 11086        |
| ln_body_mass | symmetric_M3 | -0.0039      | -0.0129      | 0.0052       | -0.7688      | -2.5378      | 1.0323       | 10911        |
| wing         | symmetric_M0 | -0.9366      | -1.4321      | -0.4412      | -2.6243      | -4.0124      | -1.2362      | 8268         |
| wing         | symmetric_M1 | -0.5066      | -0.9142      | -0.099       | -1.4193      | -2.5613      | -0.2774      | 8268         |
| wing         | symmetric_M3 | -0.3166      | -0.6541      | 0.0209       | -0.887       | -1.8327      | 0.0586       | 8077         |
| ln_ratio     | symmetric_M0 | 0.0099       | 0.0025       | 0.0174       | 2.0011       | 0.5022       | 3.5223       | 7433         |
| ln_ratio     | symmetric_M1 | 0.0025       | -0.0055      | 0.0106       | 0.5075       | -1.0923      | 2.1333       | 7433         |
| ln_ratio     | symmetric_M3 | 0.0017       | -0.0066      | 0.01         | 0.3401       | -1.3139      | 2.0218       | 7269         |

