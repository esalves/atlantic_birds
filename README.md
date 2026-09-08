# Shrinking Body Size in Atlantic Forest Birds as a Response to Climate Change

This repository contains the data, analytical pipeline, and Quarto manuscript for the study of multi-decadal morphological changes (wing length, body mass, and bill width) in Atlantic Forest passerine birds in response to climate warming and food availability.

---

## Quick Navigation

- **Manuscript source:** [`Manuscript/index.qmd`](Manuscript/index.qmd) (renders to HTML and Word via Quarto)
- **Results synthesis:** [`Analysis/output/RESULTS_SUMMARY.md`](Analysis/output/RESULTS_SUMMARY.md) (comprehensive summary of all model findings)
- **Repository layout details:** [`REPO_STRUCTURE.md`](REPO_STRUCTURE.md)
- **Server execution order & live-only filter:** [`Analysis/scripts/LIVE_ONLY_MIGRATION.md`](Analysis/scripts/LIVE_ONLY_MIGRATION.md)
- **Manuscript rendering instructions:** [`Manuscript/RENDER_STEPS.md`](Manuscript/RENDER_STEPS.md)
- **R session information & versions:** [`Analysis/R_session_info.txt`](Analysis/R_session_info.txt)

---

## Repository Structure

The repository cleanly separates immutable inputs, derived data, analytical code, fitted models, and the publication manuscript:

```text
atlantic_birds/
├── Analysis/
│   ├── data/
│   │   ├── raw/             # Raw immutable datasets (ABT trait CSV, EltonTraits, PREDICTS, etc.)
│   │   └── derived/         # Generated analytical datasets (passer90.rda, passer90_climate.rds, etc.)
│   ├── scripts/             # Modular R analysis scripts and sampling configs
│   ├── output/              # Quoted summary statistics (*.rds, *.md)
│   │   └── models/          # Fitted brms Bayesian model objects across phylogenetic trees (*.rda)
│   ├── figures/             # Causal DAG diagrams, diagnostic and sensitivity plots
│   └── R_session_info.txt   # Exact R and package versions used
├── Manuscript/              # Publication source
│   ├── index.qmd            # Quarto manuscript (reads precomputed outputs; no refitting needed)
│   ├── make_figures.py      # Python script generating Figs 1–4 into images/
│   ├── requirements.txt     # Python dependencies for figure generation
│   └── references.bib       # BibTeX bibliography
├── archive/                 # Historical, exploratory, or superseded work (kept for provenance)
├── DataManagementPlan/      # Formal project data management plan
├── REPO_STRUCTURE.md        # Detailed directory mapping and path conventions
└── README.md                # This file
```

---

## Setup & Dependencies

The project uses **R (>= 4.4)** for data wrangling and Bayesian/SEM modeling, **Python (>= 3.9)** for geographical and trend figures, and **Quarto** for manuscript rendering.

### 1. R Environment & Packages

Install standard CRAN dependencies:
```R
install.packages(c(
  "tidyverse", "brms", "ape", "phytools", "MCMCglmm", "metafor",
  "data.table", "posterior", "future.apply", "mice", "readr", "ggplot2"
))
```

Install packages required for dynamic bird phylogenies (`clootl`, `prepR4pcm`) and distributional structural equation modeling (`drmTMB`, `drmSEM`):
```R
# Avian phylogeny retrieval
install.packages("clootl")

# Taxonomic reconciliation and tree auditing
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak("itchyshin/prepR4pcm")

# Distributional piecewise structural equation modeling (drmSEM)
pak::pak(c("drmTMB", "drmSEM"), repos = "https://itchyshin.r-universe.dev")
```

> **Note:** A snapshot of package versions and system configuration is recorded in [`Analysis/R_session_info.txt`](Analysis/R_session_info.txt).

### 2. Python Environment (for Manuscript Figures)

```bash
pip install -r Manuscript/requirements.txt
```

### 3. Quarto CLI

Install Quarto from [quarto.org](https://quarto.org/) to compile the manuscript.

---

## How to Reproduce the Work

### Option A: Fast Reproduction (Inspect Results & Render Manuscript)

Because fitting 50-tree Bayesian multi-level models in Stan takes substantial compute time (often hours or days on multi-core workstations), **all fitted model objects (`Analysis/output/models/`) and extracted statistics (`Analysis/output/`) are provided pre-computed in the repository**.

To inspect results or render the manuscript:

```bash
# 1. Regenerate manuscript Figures 1–4 (optional; pre-rendered images are included)
cd Manuscript
python3 make_figures.py

# 2. Render the manuscript to HTML and Word DOCX
quarto render index.qmd
```
The compiled outputs will be generated in `Manuscript/_manuscript/` (and `Manuscript/index.docx`).

---

### Option B: Full Pipeline Re-execution (From Raw Data to Models)

To re-run the entire pipeline from the raw datasets (`Status == 'live'`, 73 species, 12,571 records, 1995–2018):

```bash
cd Analysis/scripts

# 1. Rebuild the live-only analytical dataset
Rscript rebuild_passer90_live.R            # -> Analysis/data/derived/passer90.rda

# 2. Extract per-record time-resolved monthly climate (WorldClim)
Rscript climate_extraction.R               # -> Analysis/data/derived/passer90_climate.rds

# 3. Fit Bayesian phylogenetic models across 50 trees (Rubin's rules pooling)
Rscript atlantic_parallel.R                # Primary wing length model
Rscript atlantic_parallel_bill.R           # Bill width model
Rscript atlantic_parallel_mass.R           # Body mass model
Rscript atlantic_diet_interaction.R        # Diet x year interaction
Rscript atlantic_drmsem.R                  # Distributional piecewise SEM (drmSEM)
Rscript atlantic_sensitivity_thresholds.R  # Sensitivity across inclusion criteria

# 4. Update descriptive statistics and export table for figure generation
Rscript update_descriptive_stats.R         # -> Analysis/output/descriptive_summary.rds
                                           # -> Analysis/data/derived/passer90_export.csv

# 5. Refresh figures and render Quarto manuscript
cd ../../Manuscript
python3 make_figures.py
quarto render index.qmd
```

> **Hardware note:** Bayesian sampling configurations (`Analysis/scripts/_sampling_config.R`) automatically detect whether the script runs on a laptop or high-core server (Totoro), defaulting to safe single-worker memory usage on laptops to prevent out-of-memory errors.

---

## Data Provenance

- **`ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv`**: Comprehensive morphometric measurements for Atlantic Forest birds ([Rodrigues et al. 2019](https://doi.org/10.1002/ecy.2647)).
- **`BirdFuncDat.txt`**: Functional and dietary trait data from EltonTraits 1.0 ([Wilman et al. 2014](https://doi.org/10.1890/13-1917.1)).
- **`predicts_extract.rds`**: Invertebrate abundance data from the PREDICTS database ([Hudson et al. 2017](https://doi.org/10.1002/ece3.2579)).
- **Phylogenetic Trees**: Complete dynamic avian tree cloud ([McTavish et al. 2025](https://doi.org/10.1093/sysbio/syae054)) via `clootl` and `prepR4pcm`.

---

## Citation

If you use this codebase or data, please cite the associated manuscript (see [`Manuscript/index.qmd`](Manuscript/index.qmd) for author and citation details).