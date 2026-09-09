# Revision status log — referee revision, 2026-09-09

Status of the major revision planned in `REVISION_PLAN.md` (the authoritative plan)
after the first working day (agents P0–P7 plus an independent review). One entry per
phase: what was produced, the key numbers with their uncertainty and N, what has not
run yet, and what blocks the next step. Numbers are copied from the `.rds` outputs
named in each section, not from the agents' reports; where a report and the file on
disk disagree, the file is quoted and the disagreement is flagged.

Companion documents: `Analysis/scripts/REVISION_NOTES_P0.md` … `P6.md`,
`REVISION_NOTES_P4a.md`, `REVISION_NOTES_P4ef.md`, `Analysis/scripts/REVISION_REVIEW.md`
(reviewer), `Analysis/output/RESULTS_SUMMARY.md` § "Revision diagnostics, fast tier (2026-09)".

Conventions used throughout: β = scaled-year coefficient (mm per SD-year, SD = 5.0199 yr,
the `scale()` SD actually applied to `scaled_yr`; × 1.992 for per decade); mean wing
71.13 mm (n = 8,478); "early / late" = ≤ 2006 / ≥ 2013; sample = live-measured, adult,
known-sex Passeriformes inside the 20-km buffer, 73 species / 12,571 records / 8,478 wing
records; record spans 1995–2018 (24 calendar years; over-record conversions multiply by 23).

---

## 0. Where the revision stands (2026-09-09, end of day)

| | Status |
|---|---|
| Decision gate (REVISION_PLAN §2) | **Provisionally Scenario B** on the lme4 fast tier: the fully controlled wing year effect (M3) overlaps zero on both samples (full −0.370 [−0.775, +0.034]; complete-case −0.284 [−0.643, +0.075]) and the within-contributor slope overlaps zero in every specification (−0.412 to +0.028). **Not final**: the brms tiers (Tier 2 `--trees 1`, Tier 3 `--trees 50`) have not run. The Tier-1 bug (H1 below) was fixed and the tier re-run at 12:35; every Phase 1 number in this log is from the corrected rds, and the gate reading did not change (see §3.1). |
| Phases run in full locally | P0 (data + audit), P2 fast tier, P3 fast tier, P4a (complete — no Stan involved), P4ef (code + archive; no full SEM), P5 (climate in full; diet lme4 tier), P6 (figures from files on disk), P7 (mechanical manuscript edits; renders). |
| Phases partly run | P1: Tier 1 re-run in full at 12:35 after the H1 fix (`Sex` in 35 of 39 specs; the rds on disk is the corrected run); its brms `--smoke` was killed before any Stan output, so the brms path is untested. Every brms/Stan path (P1, P2, P3, P5 diet, drmSEM v4) untested or unfinished locally. |
| Reproducibility | Reviewer re-ran every fast path: all outputs identical to the agents' files (md5 / `all.equal` at 1e-6) except Phase 1 (H1–H3, since resolved: guard fixed, `--fast-lme4` re-run, rds `$generated` 12:35:53, notes rewritten from the rds under one specification). |
| Manuscript | `Manuscript/index.qmd` carries the mechanical corrections only (duration, MICE paragraph, ants logic, effect scale inline, "modest" removed); renders to HTML. The gate-dependent rewrite (title, abstract, Results order, PREDICTS deletion, drmSEM demotion, Allen's-rule retraction) has **not** started. |
| Git | Nothing committed by the agents. Working tree: 11 modified tracked files, 7 staged renames (`archive/predicts/`), ~50 untracked new outputs/notes. HEAD = `50aa4a3` (the plan). |

Two results that cut against the intended Scenario-B framing and must be carried into
the paper: (i) the unknown-sex replication keeps a clear negative wing trend under the
full controls (M6 −0.684 [−1.196, −0.172], N = 3,532), and (ii) removing six
sign-blind anomalous contributor × site × year blocks makes the controlled estimate
exclude zero (M3_noanom −0.428 [−0.808, −0.048], N = 8,285). Two results that cut
against the *published* manuscript beyond the wing trend: the wing-vs-mass slope
difference on shared records does not exclude zero (−0.92 %/decade [−3.27, +1.43]),
and the variance rise disappears within measurer (lnCVR within species × sex ×
contributor −0.249 [−0.515, +0.017], k = 15).

---

## 1. Reviewer's high-severity issues (verbatim from `Analysis/scripts/REVISION_REVIEW.md` §1)

All three are in Phase 1. Copied without alteration:

| # | File | Issue |
|---|---|---|
| H1 | `Analysis/scripts/atlantic_parallel_controlled.R` L412–417 | Single-level guard `nlevels(dat[[fct]]) < 2` is TRUE for the character column `Sex` (nlevels of a character vector is 0), so **`Sex` is dropped from every Tier-1 model on the known-sex frames** (M0–M7, `_cc`, `_noanom`, `M_carrano`, jackknife). M0 therefore no longer reproduces the published brms structure (−0.929 vs −0.922 with Sex). |
| H2 | `Analysis/scripts/REVISION_NOTES_P1.md` §2, §3, §5; P1 report | The before/after table quotes the **with-Sex** estimates (e.g. M0 −0.922, M1 −0.528, M3 −0.370 [−0.775, 0.034], M7 −0.379 [−0.781, 0.024]) — evidently from a run before the guard was added — while the saved script (and its own final log `scratchpad/fast_run.log`, 11:01) produces the **without-Sex** estimates (M0 −0.929, M1 −0.543, M3 −0.387 [−0.793, 0.019], M7 −0.410 [−0.815, −0.006]). The `_noanom`, `M_carrano`, `M6*` and `M3_longrun_cc` rows in the notes come from the final run and are without Sex. The notes are therefore a mixture of two specifications and do not match `controlled_wing_results.rds`. |
| H3 | `Analysis/output/controlled_wing_results.rds`, `controlled_wing_species_slopes.rds` | **Absent at review start** (the P1 run that produced the notes had not written them). Regenerated by the reviewer's re-run of the saved script, so the files now on disk carry the **Sex-dropped** estimates of H1. Do not quote or plot them; re-run after the fix. |

Independently re-checked for this log (2026-09-09, after the review): the script on disk
(`atlantic_parallel_controlled.R`, mtime 10:56) still has the guard at L413; the rds on disk
(generated 12:02:47, mode `fast-lme4`) has `Sex` in the fixed effects of **0 of 39** specs and
the label "[Sex dropped: single level]" on 35 (the other 4 are the unknown-sex M6 models,
which never had `Sex`). That was the state at ~12:10; superseded by the note below.

**Resolution (2026-09-09, after the review) — H1, H2, H3 resolved.** H1 fixed in
`atlantic_parallel_controlled.R`: `prep_frame()` now coerces
`Sex = factor(Sex, levels = c("Female", "Male"))` (L220) and the single-level guard tests
`length(unique(na.omit(dat[[fct]]))) < 2` (L413–419), which is correct for character and
factor columns. H2 and H3 resolved by re-running `--fast-lme4`: `controlled_wing_results.rds`
and `controlled_wing_species_slopes.rds` regenerated with `$generated` =
**2026-09-09 12:35:53 MDT**; `$specs$fixed` contains `Sex` for **35 of 39** specs (the 4
without are the unknown-sex M6 models), and the only spec carrying a "dropped" label is
`M_carrano` (ring term, `wing_col` — legitimate). `REVISION_NOTES_P1.md` rewritten from the
12:35 files under one specification (its §9 is the correction log); the `finish()`
attribution and the smoke-test sentence corrected. The same note is recorded under
REVISION_REVIEW.md §1. Still open: the brms path of the script is untested locally. P6
figures were re-exported from the corrected rds at 12:37 (Fig. 2 draws M0 −0.922, M3 −0.370;
caterpillar 49 of 72 negative, 9/2 intervals excluding zero).

Medium-severity items from the review (paraphrased, see REVISION_REVIEW.md §2 and §4):
P3 notes header/§8 imply a completed brms smoke test — the log shows it was killed at
iteration ~201/400 and `models/variance_sigma_SMOKE.rda` does not exist; P1 notes §7 said the
smoke test "compiled and sampling started" — `smoke_run.log` ends at the model announcement
with no Stan output (P1 §7 now says so); `index.qmd` bill-width sentence applies "the same standard" to an
interval that excludes zero and one that includes it (tagged `REVISE: gate`);
`REVISION_PLAN.md` §0 carries numbers that do not reproduce (§6 below);
`update_descriptive_stats.R` `sample.n.wing = n()` counts 2,630 NA-wing records (needs an
owner); smoke artefacts `bivariate_smoke*.{rds,rda}` are not git-ignored.

---

## 2. Phase-by-phase status

| Phase | Script(s) | Ran here? | Outputs on disk | Still to run |
|---|---|---|---|---|
| P0 audit + data | `rebuild_passer90_live.R` (modified), `audit_provenance.R` (new) | **Full** (< 1 min each) | `data/derived/passer90.rda` (rebuilt, 42 cols), `passer90_allsex.rda` (new), `output/audit_{sources,sites,individuals,season,traits}.rds`, `effect_scale.rds`, 3 figures | Nothing. `effect_scale.rds$comparators` are NA placeholders (`verified = FALSE`) to be filled from the papers. |
| P1 controlled wing (gate) | `atlantic_parallel_controlled.R` (new; H1 guard fixed) | Tier 1 (`--fast-lme4`) **re-run in full** after the H1 fix (12:35); `--smoke` killed before any Stan output (brms path untested) | `output/controlled_wing_results.rds`, `controlled_wing_species_slopes.rds` (`$generated` 12:35:53; `Sex` in 35 of 39 specs — quotable at Tier-1 status) | Totoro `--smoke` first (brms path never exercised), then `--trees 1`, then `--trees 50` (M0/M3/M5/M6). P6 re-export from the 12:35 files. |
| P2 multi-trait + isometry | `atlantic_multitrait.R`, `atlantic_bivariate_wing_mass.R` (new) | lme4/glmmTMB tiers **full**; brms smoke only on an 800-record subsample (9.2 min) | `output/multitrait_results.rds`, `multitrait_table.md`, `bivariate_fast_lme4.rds`; smoke files `bivariate_smoke_results.rds`, `models/bivariate_smoke.rda` (not results) | Totoro `atlantic_bivariate_wing_mass.R --trees 1` then `--trees 50` (hours). |
| P3 variance | `atlantic_variance_sigma.R` (modified) | `--fast` **full** (18 s); `--smoke` killed at iter ~201/400 | `output/variance_results.rds` (no `$brms`, no `$brms_smoke`), `figures/variance_summary.png`, `variance_lncvr_contributor.png` | Totoro `--smoke` then `--trees 50` (optionally `--trees 50 --no-tmean`). Notes §8 needs the smoke outcome written in. |
| P4a ATLANTIC ANTS | `atlantic_ants_index.R` (new) | **Full** (~25 s; glmmTMB, no Stan) — same script is the server version | `data/derived/ants_events.rds`, `ants_events_plotdate.rds`, `ants_records_std.rds`, `output/ants_results.{rds,md}`, `ants_annual_index.rds`, 2 figures, `data/raw/atlantic_ants/README.md` | Nothing computational. Phase 4d collaboration request (S. P. Ribeiro, PERD design) is what would make the index defensible. |
| P4e PREDICTS retire | `git mv` × 7 → `archive/predicts/` + README | Done (staged, not committed) | `archive/predicts/{4 .rda, 3 .rds, README.md}` | `arthropod_estimates.rds` and `Manuscript/images/fig-arthropods.png` still read by `index.qmd` L60/L196 — move after P7 deletes the references. |
| P4f drmSEM demote | `atlantic_drmsem.R` v3 → v4 (modified) | `--smoke` **end to end** (3.9 s, 1,500 records, outputs to tempdir); full SEM **not run** | none in repo (on-disk `drmsem_results_noarthro.rds` is the superseded June v3 fit) | Totoro `Rscript atlantic_drmsem.R` (add `--no-phylo` for a first pass); section-10 `relmat()` path unexercised in situ. |
| P5 climate + diet | `climate_extraction.R`, `atlantic_diet_interaction.R` (modified) | Climate **full** (~4 min, cached rasters); diet `--no-brms` **full** (10 s); diet brms **not run, no smoke** | `data/derived/passer90_climate.rds` (live-only, 12,571 rows), `output/climate_trends.rds`, `diet_distribution.rds`, `diet_lme4_results.rds`, 3 figures | Totoro `atlantic_diet_interaction.R --trees 10` (run `--smoke` first). TerraClimate SPEI route needs `ppt`/`pet` files (absent). |
| P6 figures | `Analysis/scripts/make_figures.R` (export mode), `Manuscript/make_figures.py` (modified) | **Full** from files on disk | `Manuscript/images/fig-map.png`, `fig-wingtrend.png` (redesigned), `fig-species-slopes.png`, `fig-s-records-by-source.png`, `fig-s-wingcol.png` (new); `Analysis/output/figure_data/*.csv` (11 files) | Re-export after the P1 fix (Fig. 2 and the caterpillar currently draw the Sex-dropped estimates). Auto-upgrade to brms when `controlled_wing_brms_results.rds` exists. |
| P7 manuscript (mechanical items) | `Manuscript/index.qmd`, `Manuscript/README.md` | **Full**; `quarto render --to html` OK | `_manuscript/index.html` (git-ignored) | Gate-dependent rewrite (Phase 7 proper) not started. |
| Review | `REVISION_REVIEW.md` | **Full** (all fast paths re-run) | see §1 | — |

---

## 3. Key numbers by phase (with intervals and N)

### 3.1 Phase 1 — controlled wing trend, Tier 1 (lme4 REML, Wald ± 1.96 SE, no phylogeny)

Source: `Analysis/output/controlled_wing_results.rds` (`$generated` 2026-09-09 12:35:53,
after the H1 fix; `Sex` a fixed effect in every known-sex model — one specification
throughout). The earlier version of this table (from the Sex-dropped rds of 12:02:47) is
superseded; the reviewer's verified with-Sex values equal the rds values below.

| Model | β | 95 % CI | t | N | mm/decade |
|---|---:|---|---:|---:|---:|
| M0 baseline (= published brms structure) | −0.922 | [−1.400, −0.444] | −3.78 | 8,478 | −1.84 |
| M1 + contributor + municipality | −0.528 | [−0.923, −0.134] | −2.63 | 8,478 | −1.05 |
| M2 + wing column + lon + alt | −0.500 | [−0.894, −0.105] | −2.48 | 8,478 | −1.00 |
| **M3 full controls (+ season + moult + individual)** | **−0.370** | **[−0.775, +0.034]** | −1.79 | 8,478 | **−0.74 [−1.54, +0.07]** |
| M3_cc (complete cases) | −0.284 | [−0.643, +0.075] | −1.55 | 8,282 | −0.57 |
| M1_cc | −0.342 | [−0.708, +0.023] | −1.83 | 8,282 | −0.68 |
| M4 = M3 + contributor year slopes | −0.873 | [−1.750, +0.004] | −1.95 | 8,478 | −1.74 |
| M5 within-municipality year (REWB) | −0.735 | [−1.679, +0.209] | −1.53 | 8,478 | −1.46 |
| M5_rsTotal within-municipality | −0.402 | [−0.814, +0.011] | −1.91 | 8,478 | −0.80 |
| M5_cc within-municipality | −0.342 | [−0.609, −0.074] | −2.50 | 8,282 | −0.68 |
| M5src within-contributor (REWB) | −0.157 | [−0.797, +0.482] | −0.48 | 8,478 | −0.31 |
| M5src_rsTotal within-contributor | −0.344 | [−0.751, +0.062] | −1.66 | 8,478 | −0.69 |
| M5src_cc within-contributor | +0.028 | [−0.573, +0.628] | +0.09 | 8,282 | +0.06 |
| M5src between-contributor | −1.287 | [−2.411, −0.163] | −2.24 | 8,478 | — |
| M6 unknown-sex, M3 structure (no Sex term) | −0.684 | [−1.196, −0.172] | −2.62 | 3,532 | −1.36 |
| M6_m0 unknown-sex baseline | −0.016 | [−0.599, +0.567] | −0.05 | 3,532 | −0.03 |
| M7 first wing captures | −0.379 | [−0.781, +0.024] | −1.84 | 7,791 | −0.75 |
| M3 on 6 long-running contributors | +0.101 | [−0.517, +0.718] | +0.32 | 3,017 | +0.20 |
| E. Carrano alone (1996–2017) | −0.029 | [−0.266, +0.208] | −0.24 | 1,040 | −0.06 |
| M3 without 6 anomalous blocks | −0.428 | [−0.808, −0.048] | −2.21 | 8,285 | −0.85 |

Decision-gate slots in the rds: `scenario_full_sample = "B"`, `scenario_complete_cases =
"B"`, `scenario = "B"`; `share_of_M0_remaining_in_M3 = 0.402` (cc 0.308). Jackknife (leave
one contributor out, 42 refits per model): M3 range −0.540 … −0.263, most negative when
dropping E. Carrano, 4 of 42 refits exclude zero; M3_cc −0.511 … −0.190, 2 of 42 exclude
zero; M1 −0.753 … −0.432, all 42 exclude zero. Species slopes (M3, lme4 BLUP ± quadrature
SE — descriptive, not inference; recomputed from `controlled_wing_species_slopes.rds$M3`):
72 species, 49 negative, 9 negative / 2 positive intervals exclude zero, median −0.56
mm/decade (M3_cc: 49, 9 / 2, −0.46; M0: 56 negative, 17 / 2, median −1.18).

Per-contributor year slopes (≥ 8 sampling years; `$spanning_contributor_slopes`):
E. Carrano −0.008 [−0.25, +0.24] (N 1,040); A. Ross −0.29 [−1.28, +0.69] (618);
M. Alves −0.07 [−0.47, +0.32] (433); C. Fontana −0.68 [−2.24, +0.88] (388);
A. Piratelli −5.59 [−8.43, −2.74] (295); A. Bispo +0.69 [−0.82, +2.20] (243).
Anomalous blocks (|species × sex-centred mean| > 6 mm, n ≥ 10; 6 of 192 blocks, 193 records):
C. Fontana São Francisco de Paula 2016 (n 61, −9.3 mm) and 2017 (69, −10.9), A. Piratelli
Itu 2016 (16, −12.0), F. Santos Piraquara 2016 (13, −6.7), C. Duca Guarapari 2014
(21, +11.7), M. Pizo Iracemápolis 2012 (13, +7.1). Data-quality observations, not biology.

What flipped between the Sex-dropped and the corrected specification (reviewer; historical,
the rds no longer carries the Sex-dropped fits): M1_cc, M7, M5_rsTotal and M4 excluded zero
without Sex and overlap zero with Sex. M3 overlapped zero either way, so the gate reading
did not change; secondary statements did. Definitional findings that stand regardless: the
plan's row-5→6 attenuation (−0.48 → −0.30, attributed to season) is caused by dropping 6
undated records (4 = C. Fontana 1999) plus the 190 altitude-NA records (103 = E. Carrano
1996–99 Ilha Rasa), not by season — M1 refitted on the 8,282 complete cases is already
−0.342 (M1_cc; P2's wing cross-reference agrees); season itself is a strong wing predictor
(M3_cc: MAM +0.59, JJA +0.40, SON +0.54 mm vs DJF, t 2.6–4.0) but barely moves the year
slope on a fixed sample.

### 3.2 Phase 0 — provenance audit and effect scale

Source: `audit_*.rds`, `effect_scale.rds`. 49 of 53 plan reference numbers reproduce.

- Contributors: 47 in the sample, **42 among wing records**, 0 NA; **6 of 42 span both
  periods** (E. Carrano, A. Ross, M. Alves, C. Fontana, A. Piratelli, L. Bugoni) =
  **2,828 of 8,478** wing records; long-running (≥ 8 yr, ≥ 200 rec) 6 contributors,
  3,017 records. Wing-column proxy **87.4 % right early → 67.5 % unspecified late**.
- Sites (5,240 quartile wing records): named localities **139, 6 shared, 684 records**
  (plan's 140 / 7 / 795 counted the NA level — corrected); municipalities 99, 15 shared,
  1,728 records; species × municipality 821, 61 shared; coordinate sites 93 early / 144
  late / 4 shared. Mean latitude −24.8° → −21.9°, median altitude 175 → 489 m;
  r(year, lat) 0.17, r(year, lon) 0.10, r(year, alt) 0.07.
- Individuals: Ring on 90.9 %; 36 ring strings on > 1 species → individual = Ring ×
  Binomial; 927 individuals with > 1 record (1,813 repeats); Recapture "Yes" 4.8 % →
  12.3 % (wing records, NA as no).
- Season: DJF share of wing captures 12.6 % → 21.6 %; Molt recorded 73 % (NA 43 % early
  vs 13 % late); sex ratio (wing) 45.6 % → 48.7 % female.
- Effect scale (SD actually used 5.020; plan's 5.05 = SD of the final sample, 0.5 %
  difference): published wing −0.92 [−1.40, −0.44] → **−0.183 mm/yr, −1.83 mm/decade,
  −2.58 %/decade, −4.2 mm (−5.9 %) over 23 yr [−9.0, −2.8 %]**; mass −0.0073
  [−0.0168, +0.0022] → −1.44 %/decade, −3.3 % over record [−7.4, +1.0] (plan's −7.6 % was
  a linear approximation — corrected); isometric expectation for −5.93 % wing = −16.7 % mass.
- Comparators (Jirinec 2021, Weeks 2020, Ryding 2024): NA placeholders, `verified = FALSE`.

### 3.3 Phase 2 — multi-trait screen and wing–mass isometry (lme4 / glmmTMB, no phylogeny)

Source: `multitrait_results.rds`, `multitrait_table.md`, `bivariate_fast_lme4.rds`.
(`Sex` is an explicit factor in these scripts; not affected by H1.)

| Trait (n) | M0 baseline β | M1 + contributor + municipality | M3 (+ lon, alt, season, ring; complete cases) | Within-contributor (Mundlak) |
|---|---|---|---|---|
| log body mass (11,256) | −0.0072 [−0.0165, +0.0021] t −1.5 | −0.0056 [−0.0164, +0.0051] t −1.0 | −0.0045 [−0.0138, +0.0047] t −1.0 (n 11,077) | −0.0031 [−0.0102, +0.0039] t −0.9 |
| bill width (3,209) | **+0.24 [0.06, 0.42] t 2.6** | **−0.02 [−0.17, +0.13] t −0.3** | −0.02 [−0.17, +0.14] t −0.2 (n 3,205) | +0.03 [−0.18, +0.23] t 0.3 |
| bill length (7,697) | −0.00 [−0.11, +0.10] t 0.0 | +0.16 [−0.02, +0.34] t 1.7 | +0.17 [−0.02, +0.36] t 1.8 (n 7,540) | +0.14 [−0.13, +0.41] t 1.0 |
| tail (8,872) | +0.28 [−0.27, +0.83] t 1.0 | +0.30 [−0.19, +0.80] t 1.2 | +0.25 [−0.28, +0.77] t 0.9 (n 8,679) | +0.24 [−0.21, +0.70] t 1.0 |
| tarsus (4,361) | **+0.42 [0.26, 0.59] t 5.0** | +0.26 [−0.06, +0.58] t 1.6 | +0.26 [−0.08, +0.59] t 1.5 (n 4,285) | +0.26 [−0.16, +0.68] t 1.2 |
| wing, cross-ref (8,478) | −0.92 [−1.40, −0.44] t −3.8 | −0.53 [−0.92, −0.13] t −2.6 | −0.30 [−0.66, +0.06] t −1.6 (n 8,282) | +0.02 [−0.58, +0.62] t 0.1 |

Bill width: 22 contributors, **2 span both periods** (A. Piratelli, M. Alves; 620 of
3,209 records); + contributor only −0.24 t −3.4; between-contributor Mundlak +0.40
[−0.48, +1.29]. Mass diurnal model (n 8,345, 05–19 h): **+0.404 % mass per hour
[0.291, 0.516], t 7.02**; year slope with hour −0.0063 t −0.99.

Shared-record isometry test (7,577 records, 72 species, 42 contributors, 128
municipalities): lme4 M1 wing −0.55 [−0.97, −0.13] t −2.6; log-mass −0.0027 [−0.0135,
+0.0081] t −0.5; wing | log-mass −0.54 [−0.94, −0.15], allometric b(ln mass) 3.24 (SE
0.32). Joint glmmTMB M1: wing −1.49 %/decade [−2.64, −0.34]; mass −0.57 [−2.66, +1.53];
**difference wing − mass −0.92 %/decade [−3.27, +1.43], z −0.77**; isometry contrast
(mass − 3·wing) +3.90 [−0.07, +7.88], z 1.93; slope correlation 0.04. M2 (+ lon, alt,
season; n 7,409): wing −0.88 [−1.93, +0.16]; difference −0.24 [−2.16, +1.68], z −0.25;
isometry contrast +2.01 [−1.47, +5.48], z 1.13. **The departure from isometry is
directionally consistent but does not exclude zero in any specification**; the mass
slope interval is ~2.5× wider than wing's (residual SD 3.7 % vs 2.2 %).

### 3.4 Phase 3 — variance (glmmTMB `dispformula`, metafor lnCVR)

Source: `variance_results.rds` (mode `fast`; no `$brms`).

| Analysis | k / N | Estimate | 95 % CI |
|---|---:|---:|---|
| lnCVR pooled per species (published definition, reproduced) | k 54 | +0.184 | [+0.034, +0.334], I² 93.6 % |
| lnCVR pooled, corrected cell n | k 54 | +0.177 | [+0.031, +0.323] |
| lnCVR within sex (species × sex, n ≥ 5) | k 71 | +0.187 | [+0.048, +0.325] |
| lnCVR within species × sex × contributor (n ≥ 5) | k 15 | **−0.249** | [−0.515, +0.017] |
| lnCVR E. Carrano only (n ≥ 5) | k 7 | −0.429 | [−0.698, −0.159] (plan's k 6 / −0.46 does not reproduce) |
| σ S0: `sigma ~ year`, baseline mean | N 8,282 | **+9.6 %/decade** | [+6.3, +13.1] |
| σ S3: M3 mean; `sigma ~ year + Sex + (1|contributor) + (1|site)` | N 8,282 | **−6.3 %/decade** | [−13.5, +1.4] |
| σ S6: S3 + species intercepts in σ | N 8,282 | +3.4 %/decade | [−4.7, +12.2] |
| σ S3t: S3 on tmean-complete records | N 7,651 | −6.9 %/decade | [−14.0, +0.9] |
| σ S7: S3t + `scaled_tmean` in σ | N 7,651 | year −6.8 [−14.0, +1.0]; **tmean −0.7 % per SD [−8.7, +8.1]** (z −0.16) | |
| two-stage log\|resid\| T1 (+ contributor + site) | N 8,282 | +6.3 %/decade | [−1.7, +14.9] |

Mechanics: contributors per species × period cell 3.68 → 5.50 (species in both periods
3.77 → 5.98; 45 gain / 7 lose, Wilcoxon p 5.4e-8); between-contributor share of
within-cell variance median 0.25 → 0.45. σ random-effect SDs in S3: contributor 0.383,
municipality 0.389 (log scale; `exp(u)` range 0.45–2.53). Side effect: the
heteroscedastic S3 fit shrinks the *mean* year slope on the same 8,282 records from −0.32
[−0.69, +0.04] (lme4 M3) to −0.15 [−0.35, +0.05]. `cor(scaled_yr, scaled_tmean)` 0.12.
**The variance-rise result does not survive the source structure.**

### 3.5 Phase 4a — ATLANTIC ANTS litter-ant richness index (glmmTMB nbinom2)

Source: `ants_results.rds`. 178,976 raw → 62,020 standardised records (9.5 % with
abundance) → **855 campaign events** (393 Winkler, 462 pitfall; 61 contributor files,
293 0.1° localities, 570 sites, 55 sites resampled in > 1 year; every year 1994–2018
covered, 1–17 contributors/yr). Year on the bird scale (SD 5.020 yr).

| Model | N | β per SD-yr [95 % CI] | % per decade [95 % CI] |
|---|---:|---|---|
| naive (no RE) | 855 | −0.269 [−0.330, −0.208] | −41.5 |
| **pooled: contributor + locality + site RE, Method × log(effort), habitat** | 855 | **−0.113 [−0.187, −0.039]**, z −2.99 | **−20.1 [−31.1, −7.4]** |
| effort as offset | 855 | +0.041 [−0.051, +0.133] | +8.5 [−9.7, +30.3] (sign flips) |
| Mundlak within-contributor / between | 855 | −0.126 [−0.203, −0.050] / +0.060 [−0.203, +0.322] | −22.2 [−33.2, −9.4] |
| Mundlak within-site (55 resampled sites) | 855 | −0.410 [−0.631, −0.188] | −55.8 [−71.6, −31.2] |
| monitoring programmes ≥ 6 yr, within | 251 | −0.166 [−0.268, −0.064] | −28.1 |
| constant protocol, within | 655 | −0.290 [−0.439, −0.142] | −43.9 |
| abundance (secondary; 67 campaigns, 45 PERD), within | 67 | −0.282 [−0.724, +0.161] | null |
| leave PERD out, within | 802 | −0.055 (SE 0.047) | −10.5, CI includes 0 |

Per-programme GLM slopes disagree in sign: LAMAT/UMC +0.097 (SE 0.172); PERD −0.385
(0.111); ELMO & DELABIE −0.026 (0.062); CEPLAC/UESC −0.475 (0.069); BIOTA FORMIGAS +0.284
(0.134). Effort/site fields are not standardised across contributor files (trap IDs,
per-record codes, campaign totals; PERD convention changes within the series). Random
slope SD across programmes collapsed to 0; log(effort) exponent 0.40 Winkler / 0.32 pitfall
(richness ∝ effort^0.3–0.4, so the offset is the wrong form). **Context only: "litter-ant
richness per standardised sample", never prey biomass/availability.**

### 3.6 Phase 4e/4f — PREDICTS retired, drmSEM demoted

From `archive/predicts/dat_no_grassland.rds`: 63,150 rows, 56 studies, 41 sources, 725
sites; **53 of 56** studies single calendar year (plan: 54; definition unrecorded);
Atlantic Forest ecoregions **17 studies / 6,270 rows**, 1998–2009 (plan: 12 / 1,652 — does
not reproduce). `arthropod_estimates.rds` (still in place): main β −1.40 [−1.46, −1.34]
(N 63,150), no-ants −1.22 [−1.29, −1.15] (N 53,139) — the "shallower without ants" sentence
was a logic error (ants steepened the pooled slope); corrected in `index.qmd` by P7 pending
deletion. drmSEM v4: `INCLUDE_ARTHRO = FALSE` with a `stop()` guard keyed on
`--arthro-override=PREDICTS_RETIRED`; wing node gains `(1 | Main_researcher)`; one d-sep
claim (Sex ⟂ scaled_tmean | {scaled_yr, scaled_lat}, Fisher's C on df 2); disclosure text
carried in results; MODEL_TAG `v4-noarthro-src` invalidates caches. The on-disk
`drmsem_results_noarthro.rds` (2026-06-15, N 7,810, 72 spp, C 0.51 df 2 P 0.78, σ year
+0.035 SE 0.008, σ tmean −0.193 SE 0.009) is the superseded v3 fit that `index.qmd` still reads.

### 3.7 Phase 5 — climate context and diet

Source: `climate_trends.rds` (year-RE pooled rows), `diet_lme4_results.rds`, `diet_distribution.rds`.

- Climate at 357 localities with raster data, 1995–2018, per decade: **annual mean T
  +0.246 °C [+0.076, +0.415], t 2.84** (+0.57 °C over the window [+0.18, +0.95]; 357/357
  per-locality slopes positive); warm-quarter Tmax +0.235 [−0.014, +0.484], t 1.85; SPEI-12
  −0.256 [−0.617, +0.105], t −1.39; annual precipitation −44.7 mm [−123.1, +33.7], t −1.12.
  **Warming supported; no supportable drying trend.** Naive locality-replicate model (t 47 /
  26 / −16 / −13) kept for comparison only. SPEI from WorldClim prec + Hargreaves–Samani
  PET, generalized-logistic L-moment fit (cross-check vs `SPEI::spei` r = 1.0000).
- `passer90_climate.rds` regenerated live-only (12,571 / 73 / 46 cols; `rec_tmean` NA for
  772 records at 98 coastal localities outside the WorldClim mask, median 10.8 km to the
  nearest valid cell — left NA). r(Year, rec_tmean) = **0.302** (n 11,799): within-locality
  0.131, between-locality 0.317 (record-weighted) / 0.250 — later sampling at warmer places;
  the plan's "uncorrelated" statement does not hold.
- Diet (65 species / 8,086 wing records; 7 species lack EltonTraits): Diet-Inv bimodal (20
  spp / 1,871 rec at 90–100; 11 / 1,520 at 0–10). lme4 analogue reproduces the published
  interaction (+0.665 SE 0.227 t 2.93 vs brms +0.660 [+0.215, +1.106]); family intercept /
  slope leave it at +0.664 (family slope variance 0); **contributor + municipality
  intercepts halve it: +0.321 (SE 0.172, t 1.87, CI −0.015 to +0.658)**, main year effect
  −0.856 → −0.465. Year slope at Diet-Inv p10 / p50 / p90 (10 / 50 / 100), mm/decade:
  baseline −3.56 / −2.01 / −0.08; with provenance controls −1.82 / −1.07 / −0.14 — the
  obligate-insectivore slope is ≈ 0 in every specification.

### 3.8 Phases 6 and 7

- Figures redrawn from files (P6): Fig. 1 early N 2,370 wing records / 93 sites / 63 spp /
  12 contributors, late 2,870 / 146 / 70 / 33, 4 shared sites (235 map sites vs 233 in
  `audit_sites.rds` — 2-site grouping-key gap inside P0's script, unresolved); Fig. 2 and the
  caterpillar were re-exported at 12:37 from the corrected 12:35 rds (M0 −0.922, M3 −0.370
  [−0.775, +0.034]; 49 of 72 species slopes negative, 9 below / 2 above zero); labels say
  "lme4 fast tier … brms pending".
- Manuscript (P7): duration 24 years everywhere; MICE paragraph replaced by an accurate
  complete-case statement; ants-logic sentences corrected; Fig. 3/7 captions aligned; Zenodo
  DOI marked placeholder; "modest" removed and replaced by inline conversions from
  `effect_scale.rds` (−1.8 mm/decade, −2.6 %; −4.2 mm, −5.9 % [−9.0, −2.8]); every
  gate-dependent sentence tagged `<!-- REVISE: gate -->`. Rendered HTML has 0 occurrences
  of "modest", "three decades", "28-year", "more than two decades", "mice".

---

## 4. What still runs on Totoro (in order)

Prerequisite on the laptop (done 2026-09-09 12:35+): H1 fixed, `--fast-lme4` re-run, `REVISION_NOTES_P1.md` rewritten from the rds.

```bash
cd Analysis/scripts
# Phase 1 (gate) — after the Sex-guard fix
Rscript atlantic_parallel_controlled.R --smoke          # brms code path never reached sampling locally
Rscript atlantic_parallel_controlled.R --trees 1        # Tier 2 validation
Rscript atlantic_parallel_controlled.R --trees 50       # Tier 3: M0 / M3 / M5 / M6 (drop (1|ind) if prohibitive; M3b vs M3 differ < 0.01)
# Phase 2
Rscript atlantic_bivariate_wing_mass.R --trees 1        # pool_fixed() m = 1 path; full 7,577-record smoke never finished locally
Rscript atlantic_bivariate_wing_mass.R --trees 50       # hours
# Phase 3
Rscript atlantic_variance_sigma.R --smoke               # killed at iter ~201/400 locally; unexercised end to end
Rscript atlantic_variance_sigma.R --trees 50            # optionally also --trees 50 --no-tmean (rename the first .rda before)
# Phase 5
Rscript atlantic_diet_interaction.R --smoke             # no local smoke at all
Rscript atlantic_diet_interaction.R --trees 10          # redraws Manuscript/images/diet_interaction_plot.png
# Phase 4f (climate_extraction.R already run locally; re-run only if rasters differ)
Rscript atlantic_drmsem.R --no-phylo                    # fast first pass of v4
Rscript atlantic_drmsem.R                               # full v4 incl. section-10 50-tree relmat refit (unexercised in situ)
# then
Rscript update_descriptive_stats.R                      # after the sample.n.wing fix (needs an owner)
cd ../../Manuscript && python3 make_figures.py && quarto render index.qmd
```

Expected outputs that do not yet exist: `output/controlled_wing_brms_results.rds`
(P6's figures auto-upgrade from it), `output/bivariate_results.rds`,
`variance_results.rds$brms`, `diet_interaction_results.rds$quantile_slopes`,
`drmsem_results_noarthro.rds` (v4).

---

## 5. Blockers and open items (with owner)

1. **H1 fix (P1) — done 2026-09-09 12:35.** `atlantic_parallel_controlled.R`: `Sex =
   factor(Sex, levels = c("Female", "Male"))` in `prep_frame()` (L220) and the guard now tests
   `length(unique(na.omit(dat[[fct]]))) < 2` (L413–419); `--fast-lme4` re-run (rds `$generated`
   12:35:53, `Sex` in 35 of 39 specs); every number in `REVISION_NOTES_P1.md` §2–§5
   regenerated from the rds under one specification; §4's `finish()` attribution (the −0.54
   jackknife centre was the dropped `Sex`) and §7's smoke sentence corrected (§9 = correction
   log). P6 re-export done (12:37). Still to do: P7 quoting Phase 1, Tier 2/3 comparison.
2. **P3 notes §8** — write in the smoke outcome (launched 11:23:20; killed at iteration
   ~201/400; no `.rda`, no `$brms_smoke`; brms path unexercised) and reword the header.
3. **`REVISION_PLAN.md` §0 corrections** (see §6) before P7 quotes from it.
4. **`update_descriptive_stats.R`** — `sample.n.wing = n()` counts 2,630 NA-wing records and
   reuses the wing n for the bill lnCVR (pooled lnCVR 0.184 → 0.177 with the fix; same
   conclusion). Needs an owner (Phase 8). Also add the lnCVR / σ / cell-heterogeneity slots
   listed in `REVISION_NOTES_P3.md` §6.
5. **P7 gate rewrite** — not started. Must: delete the PREDICTS Results section, Fig. 5 and
   the `arth`/`drm_a` reads (then the two remaining `git mv`), demote drmSEM to Supplement
   with disclosure and intervals, retract Allen's-rule bill width, present isometry as a
   consistent direction with a zero-including interval, use "litter-ant richness per
   standardised sample" wording with its sensitivities, switch diet to percentile slopes,
   update Fig. 1/2 captions and add the three new image files, fill `effect_scale.rds$comparators`.
6. **Housekeeping (Phase 8)** — `REPO_STRUCTURE.md`, `LIVE_ONLY_MIGRATION.md` L55,
   `RENDER_STEPS.md` L32/L59, `make_figures.py` docstring, `README.md` L146, `.gitignore`
   for `*smoke*.{rds,rda}`; `passer90_climate.rds` regeneration date; `Manuscript/README.md`
   export-step description.
7. **Data decisions pending** — nearest-valid-cell fill for the 772 coastal NA-tmean records
   (P5 declined; would change S3t/S7/brms N automatically); DEM or locality altitude to make
   P1's nearest-site imputation unnecessary; reconcile 235 vs 233 coordinate sites (P0).
8. **Collaboration request (4d)** — S. P. Ribeiro (raw PERD design) is what would make the
   ants index defensible; A. V. L. Freitas, J. H. C. Delabie, ICMBio Monitora as in the plan.
   GBIF occupancy (4b) and resource proxies (4c) not started.

---

## 6. REVISION_PLAN.md §0 numbers that did not reproduce

| Plan statement | Reproduced value | Source |
|---|---|---|
| 140 named localities, 7 shared, 795 records | **139 / 6 / 684** (NA `Locality` was counted as a level) | `audit_sites.rds` |
| Mass lower CI → −7.6 % over record | **−7.4 %** (exact `exp()`; −7.6 was linear) | `effect_scale.rds` |
| Row 5 → 6 attenuation "attributed to season" (−0.48 → −0.30) | Sample change (6 undated + 190 altitude-NA records); season moves the year slope by ≤ 0.02 on a fixed sample | `controlled_wing_results.rds`, `multitrait_results.rds` |
| "Within a contributor, there is no temporal trend at all" (+0.02, t 0.08) | Holds only under the REWB spec on complete cases; overlaps zero in every spec but ranges −0.41 … +0.03 | `controlled_wing_results.rds` |
| 6 long-running contributors −1.14 [−2.51, +0.24] | With M3 controls +0.08 [−0.55, +0.71] (full) / −0.03 (cc) | `controlled_wing_results.rds` |
| Unknown-sex N 3,436 | 3,532 with `!known_sex` (3,436 = "Unknown" label only); same estimate | `audit_sources.rds` |
| E. Carrano lnCVR k = 6, mean −0.46 | k = 7 (n ≥ 5), pooled −0.43 [−0.70, −0.16], mean −0.37, median −0.47 | `variance_results.rds` |
| Bill-width within-contributor −0.05 (t −0.7) | ≈ 0: +0.03 [−0.18, +0.23] (sign depends on covariate set) | `multitrait_results.rds` |
| PREDICTS 54 of 56 single-year; 12 AF studies / 1,652 rows | 53 of 56 (by `Sample_midpoint` year); 17 / 6,270 (five AF ecoregions) | `archive/predicts/README.md` |
| "Shallower without ants" | Logic error: −1.22 without ants vs −1.40 with — ants steepened the slope | `arthropod_estimates.rds` |
| SD(year) 5.05 | 5.020 is the SD `scaled_yr` was scaled with (5.047 = SD of the final sample) | `effect_scale.rds` |
| Contributors per cell 4.5 → 6.6 | 4.66 → 6.65 (all records, species in both periods); wing-only 3.68 → 5.50 | `audit_sources.rds` |
| Controlled M3 −0.33 [−0.69, +0.04] | that is the 12-level-month complete-case variant; primary full-sample M3 −0.370 [−0.775, +0.034] (rds of 12:35, Sex included) | `controlled_wing_results.rds` |

---

## 7. Working tree (2026-09-09 ~12:10; nothing committed by the revision agents)

- Modified, tracked: `Analysis/data/derived/passer90.rda` (188 kB → 608 kB, 42 columns),
  `Analysis/scripts/{rebuild_passer90_live.R, atlantic_drmsem.R, climate_extraction.R,
  atlantic_diet_interaction.R, make_figures.R}`, `Manuscript/{index.qmd, README.md,
  make_figures.py, images/fig-map.png, images/fig-wingtrend.png}`.
- Staged renames (P4ef): 7 files → `archive/predicts/`.
- Untracked, new: `REVISION_STATUS.md` (this file), `Analysis/scripts/REVISION_NOTES_P{0,1,2,3,4a,4ef,5,6}.md`,
  `REVISION_REVIEW.md`, `atlantic_parallel_controlled.R`, `audit_provenance.R`,
  `atlantic_multitrait.R`, `atlantic_bivariate_wing_mass.R`, `atlantic_variance_sigma.R`,
  `atlantic_ants_index.R`; `Analysis/output/{audit_*.rds, effect_scale.rds,
  controlled_wing_*.rds, multitrait_results.rds, multitrait_table.md, bivariate_fast_lme4.rds,
  bivariate_smoke_results.rds, variance_results.rds, ants_results.{rds,md}, ants_annual_index.rds,
  climate_trends.rds, diet_distribution.rds, diet_lme4_results.rds, figure_data/, models/bivariate_smoke.rda}`;
  `Analysis/data/derived/{passer90_allsex.rda, ants_*.rds}`; `Analysis/data/raw/atlantic_ants/`;
  `Analysis/figures/{audit_*, ants_*, climate_trends, diet_*, variance_*}.png`;
  `Manuscript/images/{fig-species-slopes, fig-s-records-by-source, fig-s-wingcol}.png`;
  `archive/predicts/README.md`.
- Raw data untouched; no existing fitted model object overwritten; git-ignored
  `passer90_climate.rds` regenerated (live-only).
