# Revision notes — Phase 1 (artefact-controlled wing model, the decision gate), 2026-09-09

Companion to `REVISION_PLAN.md` §2–§3 "Phase 1" and to
`Analysis/scripts/atlantic_parallel_controlled.R`. Style follows `REVISION_NOTES_P0.md`.
**Corrected 2026-09-09 after the reviewer's H1–H3 (see §9, "Correction log").**

Executed locally: R 4.6.0, lme4 2.0.1, dplyr 1.2.1 (Tier 1, `--fast-lme4`, REML,
bobyqa, `calc.derivs = FALSE`). Every number below is read from
`Analysis/output/controlled_wing_results.rds` (`$generated` = 2026-09-09 12:35:53 MDT,
`$mode` = `fast-lme4`; slots `$before_after`, `$table`, `$jackknife`,
`$jackknife_summary`, `$spanning_contributor_slopes`, `$block_screen`,
`$anomalous_blocks`, `$decision_gate`) and `controlled_wing_species_slopes.rds`
(`$M3`, `$M3_cc`, `$M0`, `$summary`; same timestamp). **One specification throughout:
`Sex` is a fixed effect in every known-sex model** (`$specs$fixed` contains `Sex` for
35 of 39 specs; the 4 exceptions are the unknown-sex M6 models, which never had it).
**Tier-1 intervals are Wald (± 1.96 SE), there is no phylogenetic term, and the brms
tiers (`--trees`) remain the inferential result.** No brms fit from this script has
completed locally (§7).

Effect scale: SD(year) = 5.0199 (the SD `scaled_yr` was standardised with,
`attr(passer90, "scaling")$year`), mean wing 71.13 mm; mm/decade = β / 5.0199 × 10.

---

## 1. Two samples, and why

Altitude is NA for 190 wing records and the capture date for 6. The plan's ladder
(§0) dropped them (8,282 complete cases). **That deletion is not random:** 103 of
the 190 are E. Carrano's 1996–1999 records from Ilha Rasa (Guaraqueçaba) — the
early anchor of the one contributor who spans 1996–2017 with 1,040 wing records.
Dropping them alone moves the M1 estimate from −0.53 to −0.34 (M1 vs M1_cc below),
i.e. most of what the plan attributed to "season" (row 5 → 6) is the sample change.
The 4 undated C. Fontana 1999 records have the same kind of leverage: the year
effect is identified by a handful of records at the ends of a few contributors'
series.

| Sample | N wing | Definition |
|---|---:|---|
| **full** (primary, no suffix) | 8,478 | altitude imputed from the nearest coordinate site with a recorded altitude (190 records, 12 localities; `$sample_sizes$altitude_imputation`); season = "unrecorded" for the 6 undated records |
| **cc** (plan-comparable, `_cc`) | 8,282 | complete cases, as in REVISION_PLAN.md §0 |

Municipality-median imputation was rejected: E. Carrano's records are from an
island (Ilha Rasa) while the only Guaraqueçaba altitude on file (563 m) is R.
Bobato's inland farm. The imputation table (contributor, locality, imputed value,
distance) is in `$sample_sizes$altitude_imputation`. Altitude has a negligible
coefficient (M3: +0.21 [−0.16, +0.57] per SD) and municipality is a random intercept,
so the fill cannot drive the year effect; deleting the records can.

Other sample definitions: individual = Ring × Binomial (unringed birds are
singletons; 7,791 individuals, 373 with repeat wing records); moult = 3 levels
(no / Yes / unrecorded; 16.7 % unrecorded overall, so it is a level, not a filter);
first captures = first **wing** record per individual (7,791; P0's 7,734 requires the
first record of any kind to carry a wing); unknown-sex = `!known_sex` in
`passer90_allsex.rda` (3,532; the plan's 3,436 is the "Unknown" label only — same
estimate to two decimals).

## 2. Before/after table (scaled-year coefficient, mm per SD-year; lme4 REML, Sex included)

Full sample (primary); `$before_after`:

| Model | β_year | 95 % CI | t | N | mm/decade |
|---|---:|---|---:|---:|---:|
| M0 baseline (= brms structure) | −0.922 | [−1.400, −0.444] | −3.78 | 8,478 | −1.84 |
| M1a + (1\|contributor) | −0.678 | [−1.078, −0.279] | −3.33 | 8,478 | −1.35 |
| M1b + (1\|municipality) | −0.779 | [−1.180, −0.378] | −3.81 | 8,478 | −1.55 |
| **M1** + contributor + municipality | −0.528 | [−0.923, −0.134] | −2.63 | 8,478 | −1.05 |
| M2a + wing-column proxy | −0.500 | [−0.894, −0.106] | −2.49 | 8,478 | −1.00 |
| **M2** + longitude + altitude | −0.500 | [−0.894, −0.105] | −2.48 | 8,478 | −1.00 |
| M3a + season | −0.385 | [−0.794, +0.023] | −1.85 | 8,478 | −0.77 |
| M3b + (1\|individual) (no moult) | −0.365 | [−0.771, +0.040] | −1.77 | 8,478 | −0.73 |
| **M3** + moult (full controls) | **−0.370** | **[−0.775, +0.034]** | −1.79 | 8,478 | **−0.74 [−1.54, +0.07]** |
| M4 = M3 + contributor year slopes | −0.873 | [−1.750, +0.004] | −1.95 | 8,478 | −1.74 |
| M5 within-municipality year (REWB) | −0.735 | [−1.679, +0.209] | −1.53 | 8,478 | −1.46 |
| M5 between-municipality (site mean) | −0.274 | [−0.756, +0.208] | −1.11 | 8,478 | — |
| M5_rsTotal within-municipality (species slope on total year) | −0.402 | [−0.814, +0.011] | −1.91 | 8,478 | −0.80 |
| M5src within-contributor year (REWB) | −0.157 | [−0.797, +0.482] | −0.48 | 8,478 | −0.31 |
| M5src between-contributor | −1.287 | [−2.411, −0.163] | −2.24 | 8,478 | — |
| M5src_rsTotal within-contributor | −0.344 | [−0.751, +0.062] | −1.66 | 8,478 | −0.69 |
| **M6** unknown-sex, M3 structure (no Sex term) | **−0.684** | **[−1.196, −0.172]** | −2.62 | 3,532 | −1.36 |
| M6_m1 unknown-sex, M1 structure | −0.741 | [−1.253, −0.229] | −2.84 | 3,532 | −1.48 |
| M6_m0 unknown-sex, baseline | −0.016 | [−0.599, +0.567] | −0.05 | 3,532 | −0.03 |
| M7 first captures (M3, no ring term) | −0.379 | [−0.781, +0.024] | −1.84 | 7,791 | −0.75 |
| M3_longrun 6 long-running contributors | +0.101 | [−0.517, +0.718] | +0.32 | 3,017 | +0.20 |
| M_carrano E. Carrano alone (M3 minus src; ring and wing_col dropped, single level) | −0.029 | [−0.266, +0.208] | −0.24 | 1,040 | −0.06 |
| M1_noanom (6 anomalous blocks removed) | −0.495 | [−0.877, −0.112] | −2.53 | 8,285 | −0.99 |
| M3_noanom | −0.428 | [−0.808, −0.048] | −2.21 | 8,285 | −0.85 |
| M5src_noanom within-contributor (REWB) | −0.412 | [−0.992, +0.168] | −1.39 | 8,285 | −0.82 |
| M5src_rsTotal_noanom within-contributor | −0.374 | [−0.757, +0.009] | −1.91 | 8,285 | −0.74 |

Complete-case sample (plan-comparable):

| Model | β_year | 95 % CI | t | N | mm/decade |
|---|---:|---|---:|---:|---:|
| M0_cc | −0.936 | [−1.400, −0.473] | −3.96 | 8,282 | −1.87 |
| M1_cc | −0.342 | [−0.708, +0.023] | −1.83 | 8,282 | −0.68 |
| M2_cc | −0.318 | [−0.684, +0.049] | −1.70 | 8,282 | −0.63 |
| M3a_season_cc | −0.299 | [−0.665, +0.067] | −1.60 | 8,282 | −0.60 |
| M3b_ring_cc | −0.279 | [−0.640, +0.081] | −1.52 | 8,282 | −0.56 |
| **M3_cc** | **−0.284** | **[−0.643, +0.075]** | −1.55 | 8,282 | **−0.57 [−1.28, +0.15]** |
| M4_cc | −0.816 | [−1.652, +0.020] | −1.91 | 8,282 | −1.62 |
| M5_cc within-municipality (REWB) | −0.342 | [−0.609, −0.074] | −2.50 | 8,282 | −0.68 |
| M5_cc between-municipality (site mean) | −0.255 | [−0.751, +0.242] | −1.01 | 8,282 | — |
| M5_rsTotal_cc within-municipality | −0.292 | [−0.662, +0.078] | −1.55 | 8,282 | −0.58 |
| M5src_cc within-contributor (REWB) | +0.028 | [−0.573, +0.628] | +0.09 | 8,282 | +0.06 |
| M5src_cc between-contributor | −1.202 | [−2.369, −0.035] | −2.02 | 8,282 | — |
| M5src_rsTotal_cc within-contributor | −0.261 | [−0.621, +0.100] | −1.42 | 8,282 | −0.52 |
| M6_cc unknown-sex | −0.582 | [−1.115, −0.049] | −2.14 | 3,430 | −1.16 |
| M7_cc first captures | −0.295 | [−0.654, +0.064] | −1.61 | 7,595 | −0.59 |
| M3_longrun_cc | −0.024 | [−0.363, +0.315] | −0.14 | 2,897 | −0.05 |

No fit is singular and none carries a convergence message (`$before_after$singular`,
`$convergence_msg`). Fixed effects of the controls (`$table`, M3_cc): male +2.64 mm
(t 26.9); season MAM +0.59, JJA +0.40, SON +0.54 mm relative to DJF (t 2.6–4.0);
moulting −0.28 (t −2.3); moult unrecorded −0.34 (t −1.1). On the full sample (M3):
male +2.66, MAM +0.58, JJA +0.37, SON +0.50, moulting −0.34 (t −2.8); the 6 records
with season "unrecorded" carry a +9.9 mm coefficient (t 5.2), i.e. they are large
birds, not a season effect. `wing_col` adds nothing to the year effect once
contributor is in the model (M2a −0.500 vs M1 −0.528).

## 3. Comparison with REVISION_PLAN.md §0

The plan's scratch scripts (`quick_lmer.R`, `quick_lmer2.R`, in the session
scratchpad) were re-run and **reproduce the plan's table exactly** (−0.927, −0.682,
−0.783, −0.531, −0.503, −0.481, −0.315, −0.294, −0.786; within-site −0.342,
t = −2.5; within-contributor +0.023, t = 0.08; long-running −1.135; unknown-sex
−0.743 on N = 3,436). Differences between the plan and this script are therefore
definitional, and they matter ("Here" values from `$before_after`):

| Plan row | Plan | Here | Why |
|---|---:|---:|---|
| 0–4 (M0, M1a, M1b, M1, M2a) | −0.93 … −0.50 | −0.92, −0.68, −0.78, −0.53, −0.50 | same data; the plan's scratch scripts standardised year on the sample, hence the second-decimal differences |
| 5 "+ lon + altitude" −0.48 | N = 8,288 | M2 full −0.50 / M2_cc −0.32 | plan dropped 190 altitude-NA records here |
| 6 "+ season" −0.30 | N = 8,282 | M3a full −0.39 / cc −0.30 | **the drop 5→6 is the 6 undated records (4 of them C. Fontana 1999), not season**: M1_cc without season is already −0.34 |
| 7 "+ ring" −0.32 | | M3b_ring_cc −0.28; full −0.37 | plan used a 12-level month factor, here 4-level season |
| 8 contributor slopes −0.79 | | M4_cc −0.82 / full −0.87 | matches |
| Mundlak within-municipality −0.34 (t −2.5) | REWB spec | M5_cc −0.34 (t −2.5); full −0.74 [−1.68, +0.21] | matches on cc; on the full sample the REWB fit is poorly identified (wide; species slope SD on `yr_within_site` 3.6 vs 1.3 on total year) and the total-slope variant gives −0.40 [−0.81, +0.01] |
| Mundlak within-contributor +0.02 | REWB spec | M5src_cc +0.03; full −0.16 [−0.80, +0.48]; species slope on total year −0.26 (cc) / −0.34 (full) | **specification-sensitive**: which year variable carries the species random slope moves the point estimate by ~0.2–0.3; every variant overlaps zero |
| Long-running 6 contributors −1.14 [−2.51, 0.24] | M1 + contributor slopes, N 3,017 | M3_longrun +0.10 [−0.52, +0.72] (full), −0.02 [−0.36, +0.32] (cc) | plan fitted contributor slopes and no season/ring; with the M3 controls the long-running subset shows no trend |
| Unknown-sex −0.74 [−1.27, −0.21] | M1 structure, N 3,436 | M6_m1 −0.74 [−1.25, −0.23] N 3,532; M6 (M3 structure) −0.68 [−1.20, −0.17] | reproduces; survives the full controls |

## 4. Where the within-contributor signal comes from

Per-contributor year slopes (wing ~ Sex + year + lat + (1 + year || spp) + (1 | site),
fitted one contributor at a time for contributors with ≥ 8 sampling years;
`$spanning_contributor_slopes`; `Sex` is written into this formula directly, so these
fits were never affected by H1):

| Contributor | N | years | β_year | 95 % CI | species-sex-centred mean, early → late (mm) |
|---|---:|---|---:|---|---|
| E. Carrano | 1,040 | 1996–2017 | −0.008 | [−0.25, +0.24] | −0.4 (n 585) → −1.4 (n 205) |
| A. Ross | 618 | 2001–2016 | −0.29 | [−1.28, +0.69] | −0.1 (158) → +1.2 (148) |
| M. Alves | 433 | 1995–2017 | −0.07 | [−0.47, +0.32] | +1.6 (288) → +0.8 (55) |
| C. Fontana | 388 | 1999–2017 | −0.68 | [−2.24, +0.88] (singular fit) | −2.8 (247) → **−10.0** (141) |
| A. Piratelli | 295 | 2001–2017 | **−5.59** | [−8.43, −2.74] | +1.2 (118) → −2.4 (58) (Itu 2016 block: −12.0) |
| A. Bispo | 243 | 2000–2009 | +0.69 | [−0.82, +2.20] | +0.9 (98) → (no late data) |

The three largest spanning contributors show no within-contributor trend. The
negative within-contributor component comes from two contributors whose late
blocks sit 9–12 mm (13–17 % of a wing) below the species × sex means — C. Fontana,
São Francisco de Paula 2016 (n = 61, −9.3) and 2017 (n = 69, −10.9), and
A. Piratelli, Itu 2016 (n = 16, −12.0). Shifts of that size within one contributor
are protocol or data-entry changes, not phenotypes; they cannot be attributed to
biology and they are exactly what a "within-contributor" estimate is supposed to be
free of.

A sign-blind screen was therefore added (`$block_screen`, 355 contributor ×
municipality × year blocks; the screen applies to the 192 with n ≥ 10): blocks whose
species × sex-centred mean deviates by more than ± 6 mm (≈ 1.5 residual SD) are
flagged. It flags 6 of 192 blocks (193 records; `$anomalous_blocks`): the three above,
F. Santos Piraquara 2016 (n = 13, −6.7), and two **positive** blocks, C. Duca Guarapari
2014 (n = 21, +11.7) and M. Pizo Iracemápolis 2012 (n = 13, +7.1). Without them
(`_noanom` rows, N = 8,285) M3 is −0.43 [−0.81, −0.05] and the within-contributor
slope is −0.41 [−0.99, +0.17] (REWB) / −0.37 [−0.76, +0.01] (species slope on total
year). So removing the anomalous late blocks makes the pooled controlled estimate
*more* negative (two of the flagged blocks are positive and late), while the
within-contributor slope still overlaps zero. The screen is reported as a
sensitivity, not used to define the primary sample.

Leave-one-contributor-out jackknife (`$jackknife`, `$jackknife_summary`; 42 refits per
model): M1 ranges −0.753 … −0.432 around the full-sample −0.528 (all 42 refits exclude
zero); M3 −0.540 … −0.263 around −0.370 (4 of 42 exclude zero); M3_cc −0.511 … −0.190
around −0.284 (2 of 42). In every case the most negative estimate is obtained by
**dropping E. Carrano** — the contributor with the flattest series pulls the pooled
slope towards zero, as expected if the pooled slope is a between-contributor
contrast; the least negative by dropping L. Bugoni (M1, M3) or A. Ross (M3_cc).
(An earlier version of this note reported the M1 jackknife centred on −0.54 and
attributed the offset from the M1 point estimate to `finish()` re-deriving the Mundlak
means. That was wrong: the jackknife refits had `Sex` dropped (H1); with `Sex`
included the jackknife centre equals the point estimate exactly, as it should.)

## 5. Decision gate (Tier 1, provisional)

`$decision_gate` — `tier` = "1 (lme4 REML, Wald CI; provisional until the brms tiers
run)"; `rule` = "Scenario A requires the controlled year effect (M3) AND the
within-municipality year slope (M5, REWB) to exclude zero; otherwise B";
`scenario_full_sample` = **B**, `scenario_complete_cases` = **B**, `scenario` = **B**.
The rows the gate is read from:

| `$decision_gate` slot | estimate | 95 % CI | t | N |
|---|---:|---|---:|---:|
| `M0` | −0.922 | [−1.400, −0.444] | −3.78 | 8,478 |
| `M3` | −0.370 | [−0.775, +0.034] | −1.79 | 8,478 (−0.74 mm/decade) |
| `M3_cc` | −0.284 | [−0.643, +0.075] | −1.55 | 8,282 (−0.57 mm/decade) |
| `M5_within_site` | −0.735 | [−1.679, +0.209] | −1.53 | 8,478 |
| `M5_within_site_cc` | −0.342 | [−0.609, −0.074] | −2.50 | 8,282 |
| `M5_within_site_rsTotal` | −0.402 | [−0.814, +0.011] | −1.91 | 8,478 |
| `M5_between_site` | −0.274 | [−0.756, +0.208] | −1.11 | 8,478 |
| `M5src_within` (within-contributor) | −0.157 | [−0.797, +0.482] | −0.48 | 8,478 |
| `M5src_within_rsTotal` | −0.344 | [−0.751, +0.062] | −1.66 | 8,478 |
| `M5src_between` | −1.287 | [−2.411, −0.163] | −2.24 | 8,478 |
| `M6` | −0.684 | [−1.196, −0.172] | −2.62 | 3,532 |
| `M7` | −0.379 | [−0.781, +0.024] | −1.84 | 7,791 |
| `M_carrano` | −0.029 | [−0.266, +0.208] | −0.24 | 1,040 |
| `M3_noanom` | −0.428 | [−0.808, −0.048] | −2.21 | 8,285 |
| `M5src_within_noanom` | −0.412 | [−0.992, +0.168] | −1.39 | 8,285 |
| `M5src_within_rsTotal_noanom` | −0.374 | [−0.757, +0.009] | −1.91 | 8,285 |

`share_of_M0_remaining_in_M3` = 0.402; `share_of_M0_remaining_in_M3_cc` = 0.308.

- Full sample: M3 = −0.37 [−0.78, +0.03] (−0.74 mm/decade, −1.04 %/decade), 40 % of
  the baseline; M5 within-municipality −0.74 [−1.68, +0.21]; within-contributor
  (M5src) −0.16 [−0.80, +0.48]. → **Scenario B.**
- Complete cases: M3_cc = −0.28 [−0.64, +0.08] (−0.57 mm/decade), 31 % of the
  baseline; M5_cc within-municipality −0.34 [−0.61, −0.07] excludes zero but M3_cc
  does not; within-contributor +0.03 [−0.57, +0.63]. → **Scenario B.**

**The fast tier indicates Scenario B** in both samples: the controlled year
effect overlaps zero at roughly three-tenths to two-fifths of the published −0.92, and
the within-contributor slope is indistinguishable from zero in every
specification (−0.41 to +0.03), with the only large within-contributor declines
traceable to two contributors' anomalous late blocks. Two results cut the other way
and must be reported: (i) the unknown-sex replication (M6) keeps a clear negative
effect under the full controls, −0.68 [−1.20, −0.17] on 3,532 records — in the
unknown-sex data the baseline is flat (−0.02) and the controls *reveal* the decline,
the mirror image of the known-sex pattern, which itself argues that these estimates
are contrasts among contributors and sites rather than a stable phenotypic signal;
(ii) removing the six anomalous blocks moves M3 to −0.43 [−0.81, −0.05]. The
interval of the controlled estimate (−1.5 to +0.1 mm/decade) is compatible with the
published comparators as well as with zero; the brms tiers will narrow it but cannot
be expected to move the mean (REML and the published brms M0 agree to two decimals).

Species slopes (`controlled_wing_species_slopes.rds`; slope = fixed `scaled_yr` +
species BLUP, interval = slope ± 1.96 × √(SE_fixed² + conditional SD_BLUP²); these
are **lme4 BLUP quadrature intervals, descriptive, not inference** — the brms M3
posteriors replace them for the figure). Counts recomputed from the `$M3`, `$M3_cc`
and `$M0` tables (`lo`/`hi`; identical using `mm_per_decade_lo/hi`) and matching
`$summary`:

| Table | species | slope < 0 | interval entirely < 0 | interval entirely > 0 | median mm/decade |
|---|---:|---:|---:|---:|---:|
| `$M3` (full controls, full sample) | 72 | 49 | 9 | 2 | −0.56 |
| `$M3_cc` | 72 | 49 | 9 | 2 | −0.46 |
| `$M0` (baseline) | 72 | 56 | 17 | 2 | −1.18 |

M3 intervals below zero: *Tiaris fuliginosus*, *Molothrus bonariensis*, *Pyriglena
pernambucensis*, *Arremon taciturnus*, *Tangara seledon*, *Turdus flavipes*,
*Stephanophorus diadematus*, *Turdus amaurochalinus*, *Conopophaga melanops*; above
zero: *Chiroxiphia pareola*, *Drymophila malura*. SD of the species random year slope:
1.16 (M3), 0.86 (M3_cc), 1.73 (M0).

## 6. Numbers in REVISION_PLAN.md that need changing

- §0 table row 6 ("+ season") and the reading "season … attenuates the year
  effect": the attenuation between rows 5 and 6 is the deletion of 6 undated
  records (4 from C. Fontana 1999); season itself moves M2 → M3a by −0.02 (cc)
  / +0.12 (full). Season *is* a strong predictor of wing (+0.4–0.6 mm outside DJF)
  but it does not carry the year effect.
- "Within a contributor, there is no temporal trend at all" (+0.02, t = 0.08):
  true under the plan's specification on the complete-case sample; the honest
  statement is "the within-contributor slope overlaps zero in every specification
  (−0.41 to +0.03) and is specification-sensitive; the largest within-contributor
  declines are anomalous blocks of two contributors".
- "6 long-running contributors only … −1.14": with the M3 controls the subset gives
  +0.10 [−0.52, +0.72] (full) / −0.02 [−0.36, +0.32] (cc); E. Carrano alone −0.03
  [−0.27, +0.21].
- Unknown-sex N: 3,532 (`!known_sex`) rather than 3,436 (label "Unknown").
- Controlled M3 on the primary (full) sample: −0.37 [−0.78, +0.03], −0.74 mm/decade,
  not −0.33 [−0.69, 0.04] (which is the 12-level-month, complete-case variant).

## 7. What was and was not run

- **Run in full (2026-09-09 12:35, after the H1 fix):** Tier 1 (`--fast-lme4`) — the
  39 ladder/sensitivity fits in `$specs` (26 s of fitting in total, `$timings_sec`),
  the 6 per-contributor fits, the leave-one-contributor-out jackknife (3 × 42 = 126
  refits; the M3 refits with ~7,800 individual levels are the slow part), species
  slopes, block screen. The wall time of this run was not logged.
- **brms smoke test (`--smoke`: M3, 1 tree, chains = 2, iter = 400, warmup = 200,
  rstan backend): started, not completed.** The scratchpad log (`smoke_run.log`,
  10:57) shows tree retrieval, 73/73 species reconciled and the model announcement
  ("[brms] M3 (spec M3, N = 8478) on 1 tree(s): …"), then nothing — the process was
  killed before any Stan compilation or sampling output. `output/controlled_wing_brms_smoke.rds`
  and `output/models/controlled_smoke_M3.rda` **do not exist**. **The brms code path
  of this script is untested locally**; the first thing to run on Totoro is `--smoke`.
- **Not run:** Tier 2 (`--trees 1`, full sampler settings) and Tier 3
  (`--trees 50`, models M0/M3/M5/M6, `--sample full` by default). Run on Totoro.
  Expect one M3 fit (8,478 rows, ~7,800 individual levels, 73 phylo levels) to take
  substantially longer than the published M0; if it is prohibitive, drop `(1 | ind)`
  for the 50-tree run (M3b vs M3 shows the ring term changes the year estimate by
  < 0.01) and say so.

## 8. Notes for other agents

- P2/P3: use the `full` sample convention (nearest-site altitude, season
  "unrecorded") or state that you use complete cases; the two differ by 0.1–0.2
  in the wing year effect because of E. Carrano's 1996–99 Ilha Rasa records.
- P6: `controlled_wing_results.rds$before_after` (ladder figure),
  `$spanning_contributor_slopes` and `$jackknife` (supplementary figures),
  `controlled_wing_species_slopes.rds$M3` (caterpillar; replace with brms posteriors
  from `controlled_wing_brms_results.rds$species_slopes_M3` when available),
  `$anomalous_blocks` (supplementary table). Re-exported at 12:37 from the corrected
  12:35 files (Fig. 2: M0 −0.922, M3 −0.370; caterpillar 49 of 72 negative, 9/2 excluding zero).
- P7: the decision-gate statement in §5; the corrections in §6; the C. Fontana /
  A. Piratelli late-block anomaly belongs in the provenance Results, worded as a
  data-quality observation, not an accusation.
- P0 (if the pipeline is rebuilt): consider carrying `Locality`-level altitude or a
  DEM value so the nearest-site imputation becomes unnecessary; `Altitude` is NA
  for all of E. Carrano's Ilha Rasa records.

## 9. Correction log

**2026-09-09 (after `REVISION_REVIEW.md` §1, H1–H3).**

- *Defect (H1).* In `atlantic_parallel_controlled.R` the single-level guard that
  removes a factor from a formula when it has < 2 levels was written as
  `nlevels(dat[[fct]]) < 2`. `Sex` was a character column, `nlevels()` of a character
  vector is 0, so `Sex` was silently removed from every lme4 model on the known-sex
  frames (M0–M7, `_cc`, `_noanom`, `M_carrano`, all 126 jackknife refits) and the
  label "[Sex dropped: single level]" was appended to 35 of 39 specs. The
  per-contributor slopes (§4), which write `Sex` into their formula directly, and the
  unknown-sex M6 models (never had `Sex`) were not affected.
- *Consequence (H2, H3).* The first version of this note mixed two specifications:
  most rows came from a run made before the guard existed (with `Sex`), while the
  `_noanom`, `M_carrano`, `M6*` and `M3_longrun_cc` rows and the jackknife came from
  the final run (without `Sex`). The rds files then on disk (generated 12:02:47, the
  reviewer's re-run of the saved script) held the Sex-dropped fits throughout.
- *Fix.* `prep_frame()` now coerces `Sex = factor(Sex, levels = c("Female", "Male"))`
  (script L220), and the guard tests `length(unique(na.omit(dat[[fct]]))) < 2`, which
  is correct for character and factor columns (L413–419). `--fast-lme4` was re-run;
  the rds files were regenerated at 12:35:53 MDT. `$specs$fixed` contains `Sex` for
  35 of 39 specs; the only spec with a "dropped" label is `M_carrano` (ring term: no
  repeated individuals; `wing_col`: single level), which is legitimate.
- *Effect on the numbers.* Every table and sentence in §2–§5 was regenerated from the
  12:35 files under the one (Sex-included) specification. Compared with the
  Sex-dropped fits: M0 −0.922 (was −0.929), M1 −0.528 (−0.543), M3 −0.370 [−0.775,
  +0.034] (−0.387 [−0.793, +0.019]), M3_cc −0.284 (−0.308), M7 −0.379 [−0.781, +0.024]
  (−0.410 [−0.815, −0.006]), M3_noanom −0.428 [−0.808, −0.048] (−0.435 [−0.813,
  −0.056]), M_carrano −0.029 (−0.034), M3_longrun_cc −0.024 (−0.028); species slopes
  49 of 72 negative (was 53); M0 species slopes 56 negative, 17 intervals below zero
  (was 57, 18). M1_cc, M7, M5_rsTotal and M4 excluded zero without `Sex` and overlap
  zero with it. The gate reading (Scenario B on both samples) is unchanged.
- *Also corrected.* §4's attribution of the jackknife offset to `finish()` (the cause
  was the dropped `Sex`); §7's smoke-test sentence (the run never reached Stan
  output; no smoke files exist).


## 10. Phylogenetic tier (glmmTMB `propto`, 50 trees) — 2026-09-09 16:48

Engine: `_phylo_engine.R` (see `GLMMTMB_ENGINE.md`); the 50 species correlation matrices are the
published brms trees (`phylo_A_50trees.rds`). Run: `Rscript atlantic_parallel_controlled.R --trees 50
--engine glmmTMB`, 39 specs × 50 trees, 3,924 s (median 74 s, max 219 s per spec); 50/50 trees
converged for every spec except one tree of M5_cc. Output `controlled_wing_phylo_results.rds`
(`$before_after`, `$decision_gate`, `$varcomp`, `$comparison_tier1`, `$comparison_published_M0`)
and `controlled_wing_phylo_species_slopes.rds` (`$M3`, `$M3_cc`, `$M0`).

Year terms (per SD-year, Rubin-pooled Wald 95 % CI): M0 −0.922 [−1.400, −0.444]; M1 −0.528 [−0.924, −0.133]; M2 −0.500 [−0.896, −0.104]; **M3 −0.371 [−0.777, +0.036]** (−0.74 mm/decade); M3_cc −0.285 [−0.648, +0.078] (−0.57 mm/decade); M4 −0.875 [−1.758, +0.007]; within-municipality (M5) −0.738 [−1.685, +0.209], rsTotal −0.402 [−0.816, +0.012]; within-contributor (M5src) −0.159 [−0.798, +0.481], between-contributor −1.290 [−2.421, −0.159]; within-contributor cc +0.026 [−0.578, +0.629]; M7 first captures −0.379 [−0.783, +0.025]; six long-running contributors +0.095 [−0.526, +0.715]; M3_noanom −0.432 [−0.815, −0.049]; M6 unknown-sex −0.690 [−1.206, −0.174]; M6_cc −0.587 [−1.124, −0.050]; phylo proportion 0.91–0.96.

**Decision gate** (`$decision_gate`, same rule as Tier 1): **Scenario B** on the full sample and
on the complete-case sample. **Comparison with Tier 1**: no year term differs from the lme4
estimate by more than 0.002 (M0 −0.9221 vs −0.9217; M1 −0.5283 vs −0.5282; M3 −0.3709 vs −0.3704).
**Comparison with the published brms M0**: −0.922 [−1.400, −0.444] vs −0.922 [−1.403, −0.440].
**Species slopes** (fixed + BLUP, quadrature intervals; descriptive, not inference): M3 49 of 72
negative, 9 intervals below zero, 2 above, median −0.54 mm/decade; M3_cc 49 / 9 / 2, −0.46;
M0 54 / 14 / 2, −1.18. **Rows without a pooled CI**: M5_cc (one non-converged tree) and
M_carrano (single contributor; no between-contributor variance); `$errors` records no cause,
so the estimates are reported without intervals.

The brms path is retained as `--engine brms`. Its single-tree M0 run on Totoro (cmdstanr,
4 chains × 4,000 iterations, 55 min) agreed with glmmTMB; the ladder was stopped after M0 as
redundant. The `--smoke` statement in §7 is superseded: the brms path has now been exercised once
(M0) on the server, not locally.
