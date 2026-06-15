# Live-only analysis migration (2026-06)

## What changed and why

All main analyses now use **only live-measured birds** (`Status == "live"`),
dropping museum specimens. Museum skins shrink/shift on preservation (within
species: museum tail reads ~+5 %, body mass ~−1.8 %, tarsus ~+1 %; wing
essentially unaffected) **and** are concentrated in the early years (≈41 % of
records in the 1990s vs ≈6 % in 2015–2018), so they are a potential confound of
the temporal trends. The diagnosis is in `body_mass_descriptive.qmd`
(sections "museum vs live measurement confound" and the status-controlled shape
analysis).

The filter is applied **before** the `n >= 30` / `range >= 5` species filter, so
the sample stays internally consistent (every retained species has ≥ 30 live
records). Consequence:

| | species | records |
|---|---|---|
| old (all status) | 89 | 15,332 |
| **new (live only)** | **73** | **12,571** |

16 museum-dominated species drop below the 30-record threshold (e.g.
*Carpornis cucullata*, *Xiphorhynchus fuscus*, *Turdus albicollis*,
*Tangara palmarum*, *Lathrotriccus euleri*, *Sittasomus griseicapillus*).
`scaled_yr` / `scaled_lat` are **re-standardised on the live-only sample**, so
standardized coefficients are relative to the new mean/SD.

## Scripts edited (the `Status == "live"` filter lives at every raw-data entry point)

- `atlantic_birds_ms.Rmd` — `newdataframes` chunk (canonical builder of
  `passer90.rda`); `Status` added to the filter and kept in the `select()`.
- `atlantic_parallel_mass.R` — inline filter (re-reads the raw csv).
- `atlantic_sensitivity_thresholds.R` — `base_filter()` (re-reads the raw csv).
- `rebuild_passer90_live.R` — **new** helper that regenerates `passer90.rda`
  live-only without rendering the whole Rmd (validated to reproduce the old
  `passer90.rda` byte-for-byte when the live filter is off).

Scripts that simply `load(passer90.rda)` inherit the change automatically and were
**not** edited: `atlantic_parallel.R` (wing), `atlantic_parallel_bill.R`,
`atlantic_diet_interaction.R`, `climate_extraction.R`, `update_descriptive_stats.R`.
`atlantic_drmsem.R` inherits via `passer90_climate.rds`.

`body_mass_descriptive.qmd` is intentionally left on the full data — it is the
exploratory report that *diagnoses* the museum-vs-live effect and needs both.

## Stale artifacts — everything must be regenerated

`passer90.rda` has already been regenerated locally (live-only). **All other
generated outputs are still based on the old 89-species sample and must be
rebuilt on the server**: the fitted models in `output/models/*.rda`,
`passer90_climate.rds`, `passer90_mass.rda`, `descriptive_summary.rds`,
`passer90_export.csv`, the drmSEM caches, and the rendered manuscript.

The PREDICTS arthropod datasets (`dat_test.rds`, `dat_no_ants.rds`,
`dat_no_grassland.rds`) are **not** bird-measurement data and are unaffected.

## Run order on the server

```bash
cd atlantic_birds/Analysis/scripts

# 0. clear stale caches so cached chunks/results re-run on the new sample
rm -rf atlantic_birds_ms_cache atlantic_birds_ms_files
rm -f  ../output/drmsem_*cache*.rds        # drmSEM caches (also self-invalidate via signature)

# 1. data: regenerate the live-only analytical sample (fast, no models)
Rscript rebuild_passer90_live.R            # -> data/derived/passer90.rda (73 spp / 12,571)

# 2. climate: per-record WorldClim temps (needs WorldClim rasters; ~minutes)
Rscript climate_extraction.R               # -> data/derived/passer90_climate.rds

# 3. models (independent — can run in parallel; these are the heavy jobs)
Rscript atlantic_parallel.R                # wing   -> output/models/brm0_multiphylo.rda
Rscript atlantic_parallel_bill.R           # bill   -> output/models/brm_bill_multiphylo.rda
Rscript atlantic_parallel_mass.R           # mass   -> output/models/brm_mass_multiphylo.rda (re-reads csv)
Rscript atlantic_diet_interaction.R        # diet interaction
Rscript atlantic_drmsem.R                  # drmSEM (needs step 2)
Rscript atlantic_sensitivity_thresholds.R  # threshold grid (re-reads csv)

# 4. descriptive numbers + figure export (needs wing & bill models from step 3)
Rscript update_descriptive_stats.R         # -> descriptive_summary.rds, passer90_export.csv

# 5. manuscript (rebuilds passer90 itself via the edited chunk; uses model outputs)
quarto render atlantic_birds_ms.Rmd        # or rmarkdown::render(...)
```

Sanity check after step 1: the script prints
`passer90 (LIVE-only): 12571 records, 73 species`. The mass and sensitivity
scripts print their own live-only sample sizes (mass expects 12,571 / 73).
