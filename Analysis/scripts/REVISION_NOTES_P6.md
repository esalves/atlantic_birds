# Revision notes — Phase 6 (figures), 2026-09-09

Companion to `REVISION_PLAN.md` §3 "Phase 6", `Manuscript/make_figures.py` and
`Analysis/scripts/make_figures.R`. Style follows `REVISION_NOTES_P0.md`.

Executed locally: python 3.11.1 (matplotlib 3.11.1, pandas 3.0.2, numpy 2.x;
system `python3`, the `.venv-massreport` interpreter hung on import and was not
used), R 4.6.0 (dplyr, data.table, ggplot2 4.0.3). Every figure listed in §2 was
regenerated from the files on disk; nothing here is a smoke test.

---

## 1. How the figures now get their numbers

`make_figures.py` used to read a 6-column `passer90_export.csv` and carry the
model estimate as a hard-coded constant (`WING_YR_BETA = -0.9216`). The
redesigned Fig. 1, the within-species-centred Fig. 2, the species-slope
caterpillar and the two provenance supplementary figures need columns and tables
that only exist in R objects (`passer90.rda`, `audit_*.rds`, `effect_scale.rds`,
`controlled_wing_results.rds`, `controlled_wing_species_slopes.rds`). `pyreadr`
is not installed, so:

- `Analysis/scripts/make_figures.R` gained an **`export` mode** that flattens
  those objects to plain CSV in **`Analysis/output/figure_data/`**
  (`fig_records.csv`, `fig_map_sites.csv`, `fig_contributors_per_year.csv`,
  `fig_records_per_contributor_year.csv`, `fig_contributor_table.csv`,
  `fig_wingcol_by_year.csv`, `fig_controlled_before_after.csv`,
  `fig_species_slopes.csv`, `fig_scalars.csv`). Missing inputs are skipped with
  a message; nothing is filled in by hand. The drmSEM figure code in that script
  is unchanged; with no argument the script now runs the export **and** both
  drmSEM modes, `Rscript make_figures.R export` runs the export only, and the
  old `noarthro` / `arthro` calls behave as before. The path block is now the
  standard `.find_analysis_dir()` + `raw_path()/derived_path()/out_path()/fig_path()`
  helpers instead of the earlier `.repo` locator.
- `make_figures.py` calls `Rscript Analysis/scripts/make_figures.R export` at
  start-up (skip with `--no-export`; falls back to the CSVs on disk when Rscript
  is absent), then reads only `figure_data/*.csv`. No model number is typed into
  the Python script any more. The lnCVR summary constants for Fig. 4 are the one
  exception and were left as they were (Fig. 4 is out of scope for this phase).
- `make_figures.py` writes back **`figure_data/fig_panel_numbers.csv`** (every
  N, count and coefficient printed on the new figures) and
  `figure_data/fig_wingtrend_annual_means.csv` (the annual means of Fig. 2b), so
  P7 can quote them without reading the figures. These are CSV rather than rds
  because they are produced from Python; the values themselves all trace back to
  the rds files listed above.

If the brms tiers are copied back from the server as
`output/controlled_wing_brms_results.rds`, the export writes
`fig_species_slopes_brms_M3.csv` / `fig_controlled_before_after_brms.csv`
(expects `$species_slopes_M3` with `spp, mm_per_decade, mm_per_decade_lo,
mm_per_decade_hi` and `$before_after` with the Tier-1 column names) and the
Python script switches to them and relabels the panels automatically; until
then every model line is labelled **"lme4 fast tier (Tier 1, REML, Wald 95 % CI;
brms phylogenetic fits pending)"**.

## 2. Figures regenerated (`Manuscript/images/`)

| File | Status | Content |
|---|---|---|
| `fig-map.png` (Fig. 1) | **redesigned** | early (1995–2006) vs late (2013–2018) panels; coordinate sites sampled in both periods drawn as red triangles with a black edge, period-unique sites as blue circles, marker area ∝ wing records; each panel labelled with N wing records, sites (shared), species, contributors. Data: `audit_sites.rds$map_sites` (P0) + `fig_records.csv`. |
| `fig-wingtrend.png` (Fig. 2) | **replaced** | (a) every wing record as its deviation from the species × sex mean, annual means ± 95 % CI, inset of wing records per year; (b) zoom on the annual means with the **M3** (full controls) year effect and its Wald band against the **M0** baseline (dashed). Lines are β·(Year − 2009.487)/5.0199, i.e. the marginal year effect through the sample mean year. Data: `controlled_wing_results.rds$before_after` (P1, lme4 fast tier). |
| `fig-species-slopes.png` (new) | **new** | caterpillar of the 72 species year slopes from M3 (fixed + random slope, ± 1.96 · √(SE²_fixed + SD²_BLUP)), in mm/decade, sorted; filled = interval excludes zero; each species' M0 slope as a hollow grey marker; pooled M3 effect as line + band; title carries the negative / excluding-zero counts. Data: `controlled_wing_species_slopes.rds` (P1). |
| `fig-s-records-by-source.png` (Supp.) | **new** | wing records per year stacked by the 12 largest contributors (+ "other (30)"); contributors spanning both quartile periods marked †; lower panel = contributors per year. Data: `audit_sources.rds` (P0). Python redraw of `Analysis/figures/audit_records_per_year_by_source.png`. |
| `fig-s-wingcol.png` (Supp.) | **new** | share of wing records by column populated (right / left / side-unspecified) per year, n above bars. Data: `audit_sources.rds$wingcol_by_year`. Python redraw of `audit_wingcol_by_year.png`. |
| `fig-trends.png` (Fig. 3) | unchanged | code untouched; regenerated from `passer90_export.csv` (identical input). |
| `fig-variability.png` (Fig. 4) | unchanged | code untouched; regenerated from `passer90_export.csv`. |

Not produced in this phase (inputs do not exist yet or belong to other agents):
the wing-vs-mass year-slope contrast (Phase 2 brms), the σ-model variance figure
(Phase 3), the comparator-magnitude table figure (comparators in
`effect_scale.rds` are unverified placeholders), climate, diet and bill-width
before/after figures (P5 / P2 draw their own into `Analysis/figures/`). Fig. 5
(PREDICTS) and Fig. 7 (drmSEM) are P4ef / P7 decisions and were not touched.

## 2b. Numbers the figures carry (from the rds files on disk, 2026-09-09 12:02 run)

The `controlled_wing_results.rds` present when the figures were drawn was
written by a **re-run** of `atlantic_parallel_controlled.R --fast-lme4` at
12:02 (the P1 report quoted an earlier run). The two runs differ in the third
decimal: M0 −0.929 [−1.406, −0.451] here vs −0.922 [−1.400, −0.444] in the P1
report; M3 −0.387 [−0.793, +0.019] = −0.77 mm/decade vs −0.370 [−0.775, +0.034]
= −0.74; M3_cc −0.308 vs −0.284; M5src within-contributor −0.172 vs −0.157; M6
−0.684 (identical). Species slopes (M3): 72 species, 53 negative, 9 negative and
2 positive intervals exclude zero, median −0.55 mm/decade (M0: 57 negative, 18/2,
median −1.18). The figures and `fig_panel_numbers.csv` reflect the files on
disk; if P1 reruns again, `python3 make_figures.py` picks the new values up.
Whoever reconciles the P1 numbers for the manuscript should read them from
`controlled_wing_results.rds`, not from the P1 report text. Fig. 1 panel
numbers: early 2,370 wing records / 93 sites / 63 species / 12 contributors;
late 2,870 / 146 / 70 / 33; 4 sites in both periods. *Tiaris fuliginosus*
(n = 29) has an M3 slope of −11.7 [−14.3, −9.1] mm/decade; its interval is
clipped at the axis and printed as text.

## 3. Decisions worth knowing

- **Centring for Fig. 2.** Records are centred on the species × sex mean, not
  the species mean, because Sex is a fixed effect in every model and the sex
  ratio shifts slightly between periods (45.6 → 48.7 % female, `audit_season.rds`).
  Annual means are simple means with ± 1.96 SE; they are descriptive and carry
  no provenance control — the point of the panel is to show how far the raw
  annual pattern is from the controlled slope, not to estimate it.
- **Which line is "the" fit.** The task asked for the M3 fast-tier line. M0 is
  drawn dashed alongside because the manuscript's argument is the before/after
  contrast; both come from the same `before_after` table. The M3 line is a
  marginal projection: M3 also contains season, moult, wing column, longitude,
  altitude and contributor/site/individual terms, so the line is "the year
  coefficient holding everything else fixed", which the caption must say.
- **Fig. 1 site counts.** `map_sites` has 235 coordinate sites (93 early, 146
  late, 4 shared). `audit_sites.rds$coordinate_site_wing` reports 233 / 93 / 144
  / 4 for what should be the same thing; the 2-site difference is inside P0's
  script (two grouping keys for the same lon/lat pairs), not in the figure. The
  figure prints what it draws (235 / 93 / 146 / 4); P0 may want to reconcile.
  Named-locality and municipality counts in the title (6 of 139, 15 of 99) are
  P0's corrected numbers, not the plan's 7 of 140.
- **Species slopes.** The fast-tier species slopes are lme4 BLUPs with a
  quadrature SE, not posteriors; they are drawn only because the brms posteriors
  do not exist yet, and the figure says so in its title. Do not quote the
  "intervals excluding zero" counts as inference.
- **Supplementary figures are redraws, not copies.** P0's ggplot versions in
  `Analysis/figures/` stay as the audit record; the Python versions read the same
  rds-derived tables so `make_figures.py` stays the single entry point for
  `Manuscript/images/`.

## 4. For other agents

- **P7 (`index.qmd`)**: new image files `fig-species-slopes.png`,
  `fig-s-records-by-source.png`, `fig-s-wingcol.png`; Fig. 1 and Fig. 2 captions
  must change (shared vs unique sites and per-panel N; within-species-centred
  deviations, M0 vs M3, lme4 fast tier pending brms; the M3 line is a marginal
  projection). Panel numbers are in `Analysis/output/figure_data/fig_panel_numbers.csv`.
  Fig. 2's current caption ("Shading shows the density of individual records
  (log scale)") is wrong for the new figure.
- **P1**: if `controlled_wing_brms_results.rds` gets a `$species_slopes_M3`
  table with `spp, mm_per_decade, mm_per_decade_lo, mm_per_decade_hi` and a
  `$before_after` table with the Tier-1 column names, the figures upgrade
  themselves on the next `python3 make_figures.py`.
- **`Manuscript/README.md`, `RENDER_STEPS.md`, top-level `README.md`** (P7): the
  figure list and the run line should mention the export step; `make_figures.py`
  now needs R on the path (or pre-existing CSVs) — not edited here.
- **`update_descriptive_stats.R`** still writes `passer90_export.csv` for Figs 3–4;
  once those figures are revised it can be retired in favour of `fig_records.csv`.
