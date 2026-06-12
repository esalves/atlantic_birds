# Manuscript

Quarto source for *"Shrinking body size and rising prey scarcity in Atlantic
Forest birds over three decades of climate change."*

## Files

- `index.qmd` — the manuscript. Quoted statistics are read inline from the
  analysis outputs in `../Analysis/output/` (no model fitting here); figures are
  external PNGs in `images/`.
- `references.bib` — bibliography.
- `make_figures.py` — regenerates Figs 1–4 (map, wing trend, per-species
  trajectories, lnCVR forest plots) into `images/` from
  `../Analysis/data/derived/passer90_export.csv`. The arthropod figure and the
  drmSEM figure are produced by the R analysis (see `../Analysis/scripts/`).
- `data/south_america.geojson` — offline basemap for `make_figures.py`.
- `_quarto.yml` — renders `index.qmd` to HTML + Word.

## Build

```sh
python make_figures.py        # refresh figures (needs the packages in requirements.txt)
quarto render index.qmd       # → _manuscript/ (HTML + docx; git-ignored, regenerated)
```

See `../REPO_STRUCTURE.md` for how the manuscript connects to the analysis pipeline.
