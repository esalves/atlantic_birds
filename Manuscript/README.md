# Manuscript

Quarto source for *"Apparent morphological trends in Atlantic Forest birds
attenuate under observer and spatial controls."*

## Files

- `index.qmd` — the manuscript. Quoted statistics are read inline from the
  analysis outputs in `../Analysis/output/` (no model fitting here); figures are
  external PNGs in `images/`.
- `references.bib` — bibliography.
- `make_figures.py` — regenerates Figs 1–2, caterpillar, and supplementary figures
  into `images/` from `../Analysis/output/figure_data/` and `passer90_export.csv`.
- `data/south_america.geojson` — offline basemap for `make_figures.py`.
- `_quarto.yml` — renders `index.qmd` to HTML + Word.

## Build

```sh
python3 make_figures.py       # refresh figures (needs matplotlib, pandas, geopandas)
quarto render index.qmd       # → _manuscript/ (HTML + docx; git-ignored, regenerated)
```

See `../REPO_STRUCTURE.md` for how the manuscript connects to the analysis pipeline.
