# Repository structure

The repo was reorganised (2026-06) so that scripts, inputs, generated data,
results, and figures are cleanly separated. Large/regenerable data and superseded
work were moved out of the active tree.

```
atlantic_birds/
├── Analysis/
│   ├── scripts/          # all analysis code (.R, .Rmd, functions.R) + the Rmd's bib
│   ├── data/
│   │   ├── raw/          # immutable inputs (ABT csv, EltonTraits, worldclim/, terraclimate/, AvesDataLite-main/, predicts_extract.rds)
│   │   └── derived/      # generated datasets (passer90.rda, passer90_climate.rds, passer90_export.csv, dat_*.rds)
│   ├── output/           # results consumed downstream / by the manuscript
│   │   └── models/       # fitted brms model objects (*.rda)
│   └── figures/          # script-generated diagnostic figures (drmSEM, sensitivity)
├── Manuscript/           # Quarto manuscript: index.qmd, references.bib, make_figures.py,
│                         # _quarto.yml, images/, data/ (basemap). Build output
│                         # (_freeze/, _manuscript/, index_files/, .quarto/) is git-ignored.
├── archive/              # superseded / exploratory / scratch (kept for provenance, not in the pipeline)
│   ├── legacy_models/    # BirdTree-era models, loaded only by the Rmd's eval=FALSE chunks
│   ├── biotime/          # BioTIME exploration (not used by the manuscript)
│   ├── predicts/         # Retired PREDICTS arthropod models, data and estimates (Phase 4e)
│   ├── temperature_diagnostic/
│   ├── RegisteredReport/ # the registered report (separate, completed deliverable)
│   └── manuscript_agu_template/  # unused AGU Quarto extension + .cls/.sty (target is html+docx)
├── DataManagementPlan/   # project DMP (unchanged)
└── README.md, CODE_REVIEW.md, REPO_STRUCTURE.md, atlantic_birds.Rproj
```

## Path convention (scripts)

Every `.R` script resolves the `Analysis/` directory at runtime (works whether run
from the repo root, `Analysis/`, or `Analysis/scripts/`) and builds paths through
helpers — no hard-coded or `setwd()` paths:

```r
ANALYSIS_DIR <- .find_analysis_dir()                  # finds Analysis/ by locating data/derived + scripts
raw_path(...)      # Analysis/data/raw/...
derived_path(...)  # Analysis/data/derived/...
out_path(...)      # Analysis/output/...   (models: out_path("models", ...))
fig_path(...)      # Analysis/figures/...
script_path(...)   # Analysis/scripts/...
```

`atlantic_birds_ms.Rmd` sets `knitr::opts_knit$set(root.dir = <Analysis/>)` in its
`config` chunk, so its chunk paths are relative to `Analysis/` (`data/raw/…`,
`data/derived/…`, `output/…`, `output/models/…`, `../Manuscript/…`,
`../archive/legacy_models/…`).

`Manuscript/index.qmd` reads results from `../Analysis/output/…`; `make_figures.py`
reads `../Analysis/data/derived/passer90_export.csv`.

## Old → new locations (moved files)

| Old (flat `Analysis/`)                         | New                                            |
|------------------------------------------------|------------------------------------------------|
| `*.R`, `*.Rmd`, `functions.R`, `atlantic_bird_bibliography.bib` | `Analysis/scripts/`             |
| `ATLANTIC_BIRD_TRAITS_*.csv`, `BirdFuncDat.txt`, `worldclim/`, `terraclimate/`, `AvesDataLite-main/`, `AvesDataLite.zip`, `predicts_extract.rds` | `Analysis/data/raw/` |
| `passer90.rda`, `passer90_export.csv`, `passer90_climate.rds` | `Analysis/data/derived/` |
| `data_summary.rds`, `descriptive_summary.rds`, `diet_interaction_results.rds`, `sensitivity_thresholds_tier1.rds`, `drmsem_results_*.{rds,md}`, drmSEM caches | `Analysis/output/` |
| `brm0_multiphylo.rda`, `brm_bill_multiphylo.rda`, `brm_climate.rda` | `Analysis/output/models/` |
| `sensitivity_thresholds_tier1.png`, `drmsem_*_{arthro,noarthro}.png` | `Analysis/figures/` |

## Archived (moved out of the active pipeline)

- `dat_no_grassland.rds`, `dat_no_ants.rds`, `dat_test.rds`, `brm_no_grass.rda`, `brm_no_ants_no_grass.rda`, `brm_arthro_abund*.rda`, `arthropod_estimates.rds`, `fig-arthropods.png` → `archive/predicts/` (PREDICTS analyses retired in Phase 4e; studies lack multi-year temporal depth)
- `explore_biotime.R`, `biotime_brazil_arthropods.rds`, `biotime_plots/`, `biotime_cache/` → `archive/biotime/` (BioTIME exploration; not used by the manuscript)
- `brm.rda`, `brm0.rda`, `brm0.gape.rda`, `brm0.gape.model.rda`, `brm_imputed.rda`, `passer90_imp4a.rda`, `trees_passer.rda` → `archive/legacy_models/` (BirdTree-era; referenced only by the Rmd's `eval=FALSE` exploratory chunks)
- `fit_missing_model.R`, `atlantic_drmsem_results.md` (superseded by `drmsem_results_*.md`), `atlantic_birds_ms.html` (rendered output), the un-suffixed `drmsem_*` caches → `archive/`
- `temperature_diagnostic/` → `archive/temperature_diagnostic/`

## Deleted (regenerable junk)

- Stray `path/to/venv` directories (accidental `python -m venv path/to/venv`) at the repo root and `Manuscript/`
- Duplicate `AvesDataLite-main/` + `AvesDataLite.zip` at the repo root (canonical copy kept in `Analysis/data/raw/`)
- `Rplots.pdf` scratch files
