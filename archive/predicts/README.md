# archive/predicts — retired PREDICTS arthropod trend (2026-09-09, revision Phase 4e)

Everything in this directory belonged to the "arthropod prey declined steeply"
result of the submitted manuscript (Results section, Fig. 5 `fig-arthropods`,
the `arthropod_estimates.rds` inline numbers). It was moved here with `git mv`
and is kept for provenance only. **Nothing in the active pipeline reads these
files any more**, and nothing from them is to be reported in the revision.

## Why it was retired (referee §2.6; REVISION_PLAN.md §1 "§2.6 PREDICTS", §3 Phase 4e)

The referee asked either for (a) a within-study temporal model of the PREDICTS
records or (b) removal of the arthropod trend. Option (a) is not estimable on
these data. Counts re-derived from `dat_no_grassland.rds` on 2026-09-09
(calendar year of `Sample_midpoint`; the plan quotes slightly different figures,
see the note below):

| Quantity | Value |
|---|---|
| Rows (abundance records, forest biomes, ants included) | 63,150 |
| Studies (`SS`) / sources / sites | 56 / 41 / 725 |
| Studies whose records fall in a single calendar year | **53 of 56** (the other three span 2005–2007 or 2006–2007) |
| Year range of all records | 1996–2011 |
| Records in Atlantic Forest ecoregions (mangroves excluded) | 6,270 rows, **17 studies**, sampled **1998–2009** |
| Largest single study | Cabra 2006 (Colombia), 11,310 rows |

So the fitted "temporal decline" (`brm_no_grass`: standardised year effect
−1.40, 95 % CrI −1.46 to −1.34; ×0.25 per SD-year) is a contrast **between**
independent single-year studies that differ in country, method (23 sampling
methods), taxon and effort, not a within-series trend. It cannot be separated
from study turnover — the same confound the referee identified in the bird data
— and no re-modelling fixes it. The "shallower without ants" sentence was also a
logic error (the ant-exclusion coefficient, −1.22, is *closer to zero*, i.e. a
weaker decline, so ants were steepening the pooled slope, not the reverse).

> Note on the plan's numbers. REVISION_PLAN.md §1 says "54 of 56 studies
> single-year; 12 Atlantic Forest studies (1,652 rows), 1998–2009". Re-deriving
> from the archived file gives 53 of 56 (by calendar year of the sample
> midpoint) and 17 studies / 6,270 rows in Atlantic Forest ecoregions
> (Serra do Mar, Bahia, Pernambuco, Alto Paraná; the plan's narrower study set
> could not be reproduced and its definition is not recorded). The 1998–2009
> span, the 56 studies, the 63,150 rows and the 11,310-row Cabra study do
> reproduce. Use the numbers in this table if any are quoted in the response
> letter, and say how they were defined.

The trophic axis of the paper is **replaced**, not deleted: Phase 4a builds a
region-matched litter-ant index from ATLANTIC ANTS (Silva et al. 2022, Ecology;
`Analysis/scripts/atlantic_ants_index.R`) and Phase 4b an annual GBIF
occupancy indicator for prey orders (`atlantic_gbif_occupancy.R`), both with
contributor, site, method and effort terms and framed as community/occupancy
indices, never as prey biomass. REVISION_PLAN.md §7 documents the search that
led to this choice.

## Files moved here (`git mv`, staged 2026-09-09)

| File | Was | Built by (legacy `atlantic_birds_ms.Rmd` chunk) | Content |
|---|---|---|---|
| `dat_test.rds` | `Analysis/data/derived/` | `subset_predicts` (eval = `REFRESH_PREDICTS`) | PREDICTS South America Insecta + Arachnida abundance records (object `dat_test`) |
| `dat_no_ants.rds` | `Analysis/data/derived/` | `abundance_dat_no_ants` | `dat_test` without Formicidae (56,895 rows; object `dat_no_ants`) |
| `dat_no_grassland.rds` | `Analysis/data/derived/` | `abundance_dat_no_grassland` | `dat_test` restricted to forest biomes (63,150 rows; object `dat_no_grassland`) — the "main" arthropod data |
| `brm_arthro_abund.rda` | `Analysis/output/models/` | `brm_arthro_abund` (eval = `REFIT_ARTHROPODS`) | negative-binomial brms fit on `dat_test` (object `brm1`) |
| `brm_arthro_abund_no_ants.rda` | `Analysis/output/models/` | `brm_arthro_abund_no_ants` | same, on `dat_no_ants` (object `brm_no_ants`) |
| `brm_no_grass.rda` | `Analysis/output/models/` | `brm_arthro_abund_no_grass` | same, on `dat_no_grassland` (object `brm_no_grass`) — the model behind Fig. 5 and `arthropod_estimates.rds$main` |
| `brm_no_ants_no_grass.rda` | `Analysis/output/models/` | (fitted on the server; loaded by `save_arthropod_estimates`) | forest biomes without ants — `arthropod_estimates.rds$sens` |

The raw PREDICTS download `Analysis/data/raw/predicts_extract.rds` is
git-ignored and was left where it is (raw data are immutable).

## Two files deliberately NOT moved yet

`Analysis/output/arthropod_estimates.rds` and `Manuscript/images/fig-arthropods.png`
are still read by `Manuscript/index.qmd` (line 60 `readRDS(...arthropod_estimates.rds)`
and the `![...](images/fig-arthropods.png){#fig-arthropods}` figure at line 176).
Moving them before the manuscript agent deletes that section would break
`quarto render`. Once `index.qmd` no longer references them, finish the archive
with:

```bash
git mv Analysis/output/arthropod_estimates.rds  archive/predicts/arthropod_estimates.rds
git mv Manuscript/images/fig-arthropods.png     archive/predicts/fig-arthropods.png
```

`arthropod_estimates.rds` holds: main β = −1.40 [−1.46, −1.34], ×0.25
[0.23, 0.26], N = 63,150; sens (no ants) β = −1.22 [−1.29, −1.15], ×0.30
[0.28, 0.32], N = 53,139.

## Legacy code that still references the old paths (left unedited on purpose)

These are all in the legacy R Markdown `Analysis/scripts/atlantic_birds_ms.Rmd`
(not part of the active pipeline; `REFRESH_PREDICTS` and `REFIT_ARTHROPODS` are
both `FALSE`, but the *loading* chunks are eval = TRUE and will now fail if the
Rmd is knitted):

| Chunk | Lines (2026-09-09) | Reads / writes |
|---|---|---|
| `load_predicts_main`, `subset_predicts` | 861–880 | `data/raw/predicts_extract.rds` → writes `dat_test.rds` |
| `load_predicts_brasil` | 882–893 | loads `data/derived/dat_test.rds` |
| `brm_arthro_abund`, `arthro_model_summary` | 905–936 | writes / loads `output/models/brm_arthro_abund.rda` |
| `abundance_dat_no_ants`, `brm_arthro_abund_no_ants`, `arthro_model_no_ants_summary` | 937–977 | writes `dat_no_ants.rds`; writes / loads `brm_arthro_abund_no_ants.rda` |
| `abundance_dat_no_grassland`, `brm_arthro_abund_no_grass`, `arthro_model_no_grass_summary` | 978–1021 | writes `dat_no_grassland.rds`; writes / loads `brm_no_grass.rda` |
| `save_arthropod_estimates` | 1024–1062 | loads `brm_no_grass.rda`, `brm_no_ants_no_grass.rda`, `dat_no_grassland.rds` → writes `output/arthropod_estimates.rds` |
| `compare_ant_models` | 1063–1102 | loads `brm_no_grass.rda`, `brm_arthro_abund_no_ants.rda` |
| `export_fig_arthropods` | 1104–1130 | loads `brm_no_grass.rda`, `dat_no_grassland.rds` → writes `../Manuscript/images/fig-arthropods.png` |

Other mentions that are documentation only and should be updated in Phase 8
housekeeping: `REPO_STRUCTURE.md` (old→new table), `Analysis/scripts/LIVE_ONLY_MIGRATION.md`
(line 55), `Analysis/output/RESULTS_SUMMARY.md` (§3 and the file list),
`Manuscript/RENDER_STEPS.md` (fig-arthropods paragraph), `Manuscript/make_figures.py`
(docstring), top-level `README.md` line 146 (`predicts_extract.rds`).

`Analysis/scripts/atlantic_drmsem.R` (Phase 4f) no longer loads
`dat_no_grassland.rds` unless the retired arthropod node is explicitly
re-enabled with `--arthro-override=PREDICTS_RETIRED`, in which case it looks
for the file here.
