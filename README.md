# Morphological stability and subtle phenotypic shifts in Atlantic Forest birds across two decades of climate warming

Data, analysis code and Quarto source for a study of multi-decadal change in wing
length, body mass and bill width in Atlantic Forest passerines, using live-capture
records from ATLANTIC BIRD TRAITS (73 species, 12,571 records, 1995–2018).

**Authors:** Eduardo S. A. Santos, Ayumi Mizuno, Santiago Ortega, Gustavo Burin.
**Status:** manuscript in preparation (preprint forthcoming).

---

## Quick navigation

- **Manuscript source:** [`Manuscript/index.qmd`](Manuscript/index.qmd) and [`Manuscript/supplementary.qmd`](Manuscript/supplementary.qmd) (render to HTML and Word via Quarto)
- **Stan results used in the text:** [`Analysis/output/bayes/BAYES_SUMMARY.md`](Analysis/output/bayes/BAYES_SUMMARY.md)
- **Rendering instructions:** [`Manuscript/RENDER_STEPS.md`](Manuscript/RENDER_STEPS.md)
- **Repository layout and path conventions:** [`REPO_STRUCTURE.md`](REPO_STRUCTURE.md)
- **Computational environment:** [`renv.lock`](renv.lock) and [`Analysis/output/session_info.txt`](Analysis/output/session_info.txt)
- **Licences:** code MIT ([`LICENSE`](LICENSE)); data, results, figures and text CC BY 4.0 ([`LICENSE-DATA.md`](LICENSE-DATA.md)). How to cite: [`CITATION.cff`](CITATION.cff)

---

## Analytical approach in brief

- **Sample:** live-captured adults only (museum skins excluded because preparation shrinks them).
- **Main engine:** a custom Stan phylogenetic mixed model ([`Analysis/scripts/stan/phylo_lmm_marginal.stan`](Analysis/scripts/stan/phylo_lmm_marginal.stan), driven by [`_stan_engine.R`](Analysis/scripts/_stan_engine.R)). The individual random effect is integrated out analytically and the phylogenetic effect is fitted in the eigenbasis of the phylogenetic correlation matrix. Posteriors are pooled as a mixture over 20 trees from McTavish et al. (2025). Results are reported as posterior means, 95 % credible intervals and P(β < 0).
- **Frequentist check:** every model is also fitted with `glmmTMB` (REML, 50 trees, Rubin's rules; [`_phylo_engine.R`](Analysis/scripts/_phylo_engine.R), see [`GLMMTMB_ENGINE.md`](Analysis/scripts/GLMMTMB_ENGINE.md)). Supplementary Tables S1–S2 report this tier.
- **Sampling structure:** data contributor and municipality enter the models as random intercepts, because turnover in research groups and sites is confounded with calendar year in the compiled database.

---

## Repository structure

```text
atlantic_birds/
├── Analysis/
│   ├── data/
│   │   ├── raw/          # third-party inputs (ATLANTIC BIRD TRAITS csv, EltonTraits)
│   │   └── derived/      # analytical data sets (passer90*.rda, passer90_export.csv, phylo_A_50trees.rds)
│   ├── scripts/          # R analysis code; stan/ holds the Stan program
│   ├── output/           # result objects read by the manuscript (*.rds, *.md)
│   │   ├── bayes/        # pooled Stan results (per-tree fits are not tracked; see below)
│   │   ├── referee_reruns/, phylo_slopes/   # thermal, isometry and slope-structure analyses
│   │   └── figure_data/  # flat CSV tables read by the Python figure scripts
│   └── figures/          # diagnostic and supplementary figures
├── Manuscript/           # index.qmd, supplementary.qmd, references.bib, figure scripts, images/
├── archive/              # superseded and retired analyses (provenance only; not in the pipeline)
├── DataManagementPlan/
├── make_zenodo_archive.sh  # builds the curated Zenodo deposit from git HEAD
└── renv.lock, CITATION.cff, LICENSE, LICENSE-DATA.md
```

---

## Setup

- **R ≥ 4.5** with the packages in `renv.lock` (`renv::restore()`). `prepR4pcm` is installed from GitHub (`itchyshin/prepR4pcm`).
- **CmdStan ≥ 2.39** via `cmdstanr`. This is only needed to refit the Stan tier.
- **Python ≥ 3.11** for the figures: `pip install -r Manuscript/requirements.txt`
- **Quarto** ([quarto.org](https://quarto.org/)) to render the manuscript.

---

## Reproducing the work

### Option A: render the manuscript from the saved results (minutes)

The manuscript only reads saved result objects in `Analysis/output/`, so no model is refitted:

```bash
cd Manuscript
python3 make_figures.py            # main-text figures from Analysis/output/figure_data/
python3 make_isometry_figure.py    # wing–mass isometry figure
quarto render index.qmd
quarto render supplementary.qmd
```

### Option B: re-run the pipeline

Scripts are run from the repository root and find `Analysis/` on their own.

```bash
# 1. Analytical data set, audits and climate
Rscript Analysis/scripts/rebuild_passer90_live.R        # -> data/derived/passer90.rda
Rscript Analysis/scripts/audit_provenance.R             # contributor/site audits, effect_scale.rds
Rscript Analysis/scripts/climate_extraction.R           # needs a local WorldClim 2.1 download (see below)

# 2. glmmTMB phylogenetic tier (50 trees; minutes to ~1 h each)
Rscript Analysis/scripts/atlantic_parallel_controlled.R --trees 50   # wing length, adjustment ladder
Rscript Analysis/scripts/atlantic_parallel_mass.R
Rscript Analysis/scripts/atlantic_parallel_bill.R
Rscript Analysis/scripts/atlantic_multitrait.R --trees 50
Rscript Analysis/scripts/atlantic_bivariate_wing_mass.R --trees 50
Rscript Analysis/scripts/atlantic_diet_interaction.R
Rscript Analysis/scripts/atlantic_variance_sigma.R
Rscript Analysis/scripts/atlantic_phylo_slopes.R --trees 50
Rscript Analysis/scripts/referee_reruns.R               # thermal and isometry models
Rscript Analysis/scripts/referee_anomaly_lag.R          # lagged temperature anomalies
Rscript Analysis/scripts/update_descriptive_stats.R     # descriptive_summary.rds, passer90_export.csv

# 3. Stan tier (server scale: roughly 40 min to 6 h per model x tree)
#    Setting PHYLO_EXPORT_DIR before step 2 makes each glmmTMB model write a job file.
#    Each job is then refitted in Stan on each tree.
Analysis/scripts/run_bayes_totoro.sh Analysis/output/bayes/jobs_maintext.txt 1 20
Rscript Analysis/scripts/bayes_aggregate.R              # -> output/bayes/bayes_summary.rds
Rscript Analysis/scripts/bayes_derived.R                # -> output/bayes/bayes_derived.rds

# 4. Figures and manuscript
Rscript Analysis/scripts/make_figures.R export          # -> output/figure_data/*.csv
# then Option A
```

`run_bayes_totoro.sh` was written for the authors' Linux server. Adjust `CONC` and the
repository path at the top of the script for another machine. See the header of
[`bayes_fit_job.R`](Analysis/scripts/bayes_fit_job.R) for the per-fit options.

---

## What is not in the repository

- **Per-record WorldClim values** (`passer90_climate.rds`, `locality_year_tmean.rds`). WorldClim's terms do not allow redistribution. `climate_extraction.R` and `referee_anomaly_lag.R` rebuild these files from a WorldClim 2.1 download. The fitted climate trends (`climate_trends.rds`) are included.
- **Per-tree Stan fits** (`Analysis/output/bayes/fits/`, several GB). The pooled summaries the manuscript reads are included.
- **Large model objects** from the superseded brms engine (> 100 MB each).
- Internal planning and review notes. Some code comments still refer to these documents, but they are not needed to run anything.

---

## Data sources

Please cite each source in its own right. Full terms are in [`LICENSE-DATA.md`](LICENSE-DATA.md).

- **ATLANTIC BIRD TRAITS**: morphological records ([Rodrigues et al. 2019, *Ecology*](https://doi.org/10.1002/ecy.2647))
- **EltonTraits 1.0**: diet ([Wilman et al. 2014, *Ecology*](https://doi.org/10.1890/13-1917.1))
- **A complete and dynamic tree of birds**: phylogenies via `clootl` ([McTavish et al. 2025, *PNAS*](https://doi.org/10.1073/pnas.2409658122))
- **WorldClim 2.1**, downscaled from CRU-TS 4.09: temperature ([Fick & Hijmans 2017](https://doi.org/10.1002/joc.5086); [Harris et al. 2020](https://doi.org/10.1038/s41597-020-0453-3))

---

## Citation

Until the preprint is posted, please cite the repository using [`CITATION.cff`](CITATION.cff)
(GitHub's "Cite this repository" button). A Zenodo DOI will be added when the archive is minted.
