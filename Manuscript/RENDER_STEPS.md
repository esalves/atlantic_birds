# Rendering the manuscript

Steps to compile `index.qmd` so it reflects the current (live-only) analysis.
Run these on a machine with **R and Quarto installed** (e.g. your Mac or Totoro).

## Prerequisites

The manuscript reads pre-computed results from `../Analysis/output/*.rds` and
figures from `images/`. It does **not** refit any models, so rendering is fast.
Make sure the analysis outputs are current first (see
`../Analysis/scripts/LIVE_ONLY_MIGRATION.md` for the full server run order).

Current sample: **live-captured birds only, 73 species / 12,571 records, 1995–2018.**

## 1. Refresh the figures

Most figures are regenerated outside the render:

```bash
# Main-text figures (map, wing trend, per-species trends) and the
# supplementary lnCVR forests — Python
cd Manuscript
python3 make_figures.py        # reads ../Analysis/data/derived/passer90_export.csv
```

Two figures are produced from R, not by make_figures.py:

```bash
# Diet-interaction plot — REGENERATE (the committed copy may be stale)
Rscript ../Analysis/scripts/atlantic_diet_interaction.R   # writes images/diet_interaction_plot.png
```

- `fig-arthropods.png` and `arthropod_estimates.rds` were retired in Phase 4e
  and archived under `archive/predicts/` (PREDICTS studies lack multi-year temporal depth).
- The litter-ant (ATLANTIC ANTS) analysis, the century-scale museum/breakpoint
  series and the community-wide sampling expansion were dropped from both
  documents in the Sept 2026 revision. `fig-drmsem.png`, `fig-s-climate-*.png`,
  `fig-s-secular-trajectories.png` and `fig-s-model-comparison.png` are no
  longer referenced; the analysis scripts that produce them are still in
  `Analysis/scripts/`.

## 2. Render the document

```bash
cd Manuscript
quarto render        # _quarto.yml builds index.qmd → HTML + DOCX
```

Output goes to the Quarto manuscript output directory (`_manuscript/` /
`index.html` + `index.docx`).

## Notes

- All quantitative values in the text auto-populate via inline R from
  `descriptive_summary.rds`, `controlled_wing_phylo_results.rds`,
  `multitrait_phylo_results.rds`, `bivariate_phylo_results.rds`,
  `variance_phylo_results.rds`, `diet_interaction_phylo.rds`, and
  `climate_trends.rds`. If a number looks wrong,
  re-run the relevant analysis script, not the manuscript.
- `descriptive_summary.rds` now supplies the sample counts (`n_spp`, `n_rec`);
  `data_summary.rds` is no longer read by `index.qmd`.
- Quick consistency check after rendering: the abstract/intro should say
  73 species, 1995–2018, and the supplementary should contain exactly four
  notes (S1–S4), four figures (S1–S4) and three tables (S1–S3).
