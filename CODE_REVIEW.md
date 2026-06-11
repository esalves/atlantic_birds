# Code & manuscript review — Atlantic Forest birds

> ## Update — third pass: analysis upgrades from the lab's Bergmann/migration manuscript
>
> Source: *"Seasonal migration does not erase the interspecific Bergmann pattern
> in birds"* (Mizuno, Lundgren, Drobniak, Callaghan, **Santos**, Lagisz, Ortega,
> Lin & Nakagawa). Procedures from that paper that improve the present analysis,
> and what was changed:
>
> ### §A — Phylogeny (you asked to update it)
> - **Tree source updated.** That paper uses McTavish et al. (2025) *A complete
>   and dynamic tree of birds* via the **`clootl`** package (versioned, complete
>   avian tree). We replaced the Jetz et al. (2012) BirdTree/Hackett trees
>   (`read.nexus("Hackett_trees.nex")`, a file that wasn't even in the repo) with
>   clootl in both `atlantic_parallel.R` and the Rmd `phylo` chunk. **API note:**
>   `extractTree()` returns a *single* summary tree and takes `taxonomy_year` as a
>   *year* (2021–2024), not a version; the 100-tree dated cloud comes from
>   `sampleTrees()`, which first needs the AvesData repo downloaded once via
>   `get_avesdata_repo(path=".")` (unpacks to `AvesDataLite-main`). clootl uses
>   eBird scientific names **without underscores** (handled), and unmatched ABT
>   binomials are dropped — reconcile names first and check which species matched.
>   Record the tree version + taxonomy year used (`getCitations(tree)`).
> - **Phylogenetic uncertainty via Rubin's rules.** They sample 50 trees from a
>   pool of 100 and pool estimates across trees with **Rubin's rules**
>   (Nakagawa & de Villemereuil 2019) — pooled mean, within+between-tree variance,
>   Wald–Rubin 95% CI. `atlantic_parallel.R` now does this (`pool_rubin()`),
>   replacing the crude `combine_models()`/draw-concatenation.
> - **Optional:** add Pagel's-λ variance partitioning + parametric-bootstrap CIs
>   for the variance components (their `phylolm` step), as a variance-decomposition
>   check on a consensus tree rescaled to unit tip height. Note the trade-off:
>   `phylolm`+λ is cleaner for the *spatial* signal, but `brms` is kept here because
>   it supports the **random year slopes** that the *temporal* question needs.
>
> ### §B — Climate (you asked for better climate information)
> - **Their approach:** range/area-weighted species climate from **WorldClim v2.1**
>   rasters (temperature *and* precipitation), seasonal exposure windows.
> - **Key lesson for our temporal study:** WorldClim — and the
>   `Annual_mean_temperature`/`Annual_rainfall` fields in ATLANTIC BIRD TRAITS —
>   are long-term climatologies, **static in time** (confirmed: 97% of localities
>   have one temperature value across all sampling years). They cannot test change
>   over time. So we wired in **time-resolved monthly climate** (TerraClimate
>   1958–present, or CHELSA) extracted at each record's coordinates *and year* via
>   the new `Analysis/climate_extraction.R`, producing per-record temperature and
>   precipitation. The new Rmd `model_climate` chunk fits body size against actual
>   `scaled_tmean` + `scaled_ppt` (+ year retained), with random temperature
>   slopes — a genuine Bergmann/temperature-tracking test, and adds precipitation
>   as the second climatic axis.
>
> ### §D — Species/dataset/tree matching via prepR4pcm
> To make name matching robust (and auditable), both `atlantic_parallel.R` and
> the Rmd `phylo` chunk now use **prepR4pcm** (Nakagawa, Ortega, Mizuno, Santos
> et al. 2026; `pak::pak("itchyshin/prepR4pcm")`):
> - **Retrieval:** `pr_get_tree(spp, source = "clootl", n_tree = 100)` returns the
>   100-tree posterior (wraps clootl; `pr_cite_tree()` writes the provenance).
> - **Reconciliation:** `reconcile_tree(..., fuzzy = TRUE, resolve = "flag")` runs
>   the 4-stage cascade (exact → normalised → synonym → fuzzy); inspect with
>   `reconcile_summary()` / `reconcile_report()`; then `reconcile_apply(..., drop_unresolved = TRUE)`
>   returns aligned data + a pruned tree with **identical species sets** — the
>   precondition for any PCM. We prune all 100 trees to that set and subsample 50.
> - **Fixing misses:** `reconcile_override()` for one-off pairs;
>   `reconcile_crosswalk(crosswalk_birdlife_birdtree)` for the bundled
>   BirdLife↔BirdTree crosswalk; `reconcile_augment()` to graft a still-missing
>   species as sister to a congener — if used, fit models **with and without** the
>   grafts and report whether conclusions change.
> - **DONE — manual remaps removed.** The six `gsub()` name remaps in the Rmd
>   wrangling chunk (e.g. *Ceratopipra rubrocapilla* → *Pipra rubrocapilla*) were
>   made for the **old BirdTree/Jetz** taxonomy and have been deleted, so the
>   original ABT binomials now flow through to `reconcile_tree()`, which resolves
>   synonymy against the eBird/clootl tree in an auditable way. Their ABT names are
>   retained as a `historically_tricky` watch list in `atlantic_parallel.R` (with
>   the old BirdTree targets in comments, for reference only — they are **not**
>   assumed to be the eBird names). After running, confirm `reconcile_summary(rec)`
>   resolves all (or nearly all) ~68 species; for anything flagged, look up its
>   correct eBird tip and force it with `reconcile_override()` (or use
>   `reconcile_crosswalk(crosswalk_birdlife_birdtree)`).
> - **eBird genus changes fixed upfront (`ebird_synonyms`).** `clootl` errors
>   *during retrieval* on names absent from the eBird/Clements taxonomy (it does
>   not silently drop them on the `n_tree > 1`/`sampleTrees` path), so names must
>   be corrected before `pr_get_tree()` — reconciliation can't recover a species
>   that was never placed in the tree. Seven ABT binomials whose genera changed
>   since 2018 are remapped to current eBird names in both `atlantic_parallel.R`
>   and the Rmd: *Antilophia galeata*→*Chiroxiphia galeata*, *Tachyphonus
>   cristatus*→*Loriotus cristatus*, *Pyrrhocoma ruficeps*→*Thlypopsis pyrrhocoma*,
>   *Pyriglena pernambucensis*→*Pyriglena leuconota*, *Tangara sayaca*→*Thraupis
>   sayaca*, *Tangara cayana*→*Stilpnia cayana*, *Tiaris fuliginosus*→*Asemospiza
>   fuliginosa*. A `pr_get_tree(..., n_tree = 1)$unmatched` check (a `stopifnot`
>   in the script) guards the expensive 100-tree pull — extend `ebird_synonyms`
>   if it ever flags more.
>
> ### §C — Analytical framing borrowed
> - Separate **slope** change from **mean (intercept)** change (they cleanly
>   distinguished a migrant–resident slope difference from an average-size
>   difference). For us: distinguish "size declines over time/with temperature"
>   from "mean size differs among diet groups."
> - **z-standardise** body size and all continuous climate predictors on a common
>   scale; **AIC-based model selection** across candidate climate metrics
>   (lowest/mean/highest monthly temperature) and main-effects vs interaction models.
>
> ### To run (needs R + internet + downloads)
> 1. `install.packages("clootl")`; run `atlantic_parallel.R` → `brm0_multiphylo.rda`
>    + `rubin_summary`.
> 2. Download TerraClimate (or CHELSA) monthly rasters; run `climate_extraction.R`
>    → `passer90_climate.rds`; then the `model_climate` chunk → `brm_climate.rda`.
> 3. Refresh the manuscript Results and the phylogenetic-signal value (the current
>    numbers are from the old BirdTree/year-proxy pipeline; flagged `[AUTHOR ACTION]`).
>
> New references added to `Manuscript/references.bib`: `mctavish2025tree`,
> `miller2026clootl`, `nakagawa2019rubin`, `abatzoglou2018terraclimate`,
> `karger2017chelsa`.

---


> ## Update — second pass (executed)
>
> Working through the suggested order of work, the following were **done** using
> the raw `ATLANTIC_BIRD_TRAITS` CSV (present in the repo) and Python; only model
> re-fitting and the arthropod figure still need R + Stan + the missing data.
>
> 1. **Primary model decided & response harmonised** — the coalesced
>    wing-length model (`brm0`, N ≈ 9,500, decline excludes zero) is now the
>    primary inference in the manuscript; the 50-tree model is framed as a
>    robustness check. A new `Analysis/atlantic_parallel.R` re-fits the
>    phylogenetic-uncertainty model on the **coalesced** variable with modern
>    brms syntax (`data2` + `gr()`), ready to run on a Stan machine.
> 2. **Temperature/precipitation claim resolved (removed)** — inspection of the
>    data showed `Annual_mean_temperature` and `Annual_rainfall` are **static
>    per-locality climatologies** (97% of localities have one value across all
>    sampling years), so an in-sample "warming trend" would be a sampling
>    artifact. The manuscript now frames warming via the literature and adds a
>    "Climate context" Methods subsection explaining this.
> 3. **Models NOT re-fit** — needs R + Stan + `Hackett_trees.nex` (absent). Run
>    `atlantic_parallel.R`; save `brm0_multiphylo.rda`. See §2.6.
> 4. **Figures 1–3 regenerated from real data** by `Manuscript/make_figures.py`
>    (sampling map, wing/bill trends, corrected lnCVR forest plot). The Python
>    re-computation independently reproduced the pipeline: 69 species pass
>    filtering (68 after the Antilophia drop), wing 47 decreased, and a wing
>    lnCVR meta-analytic mean of 0.25 [0.11, 0.38] (vs 0.247 [0.10, 0.39] in R).
>    **Figure 4 (arthropods) is still a placeholder** — PREDICTS DB not in repo.
> 5. **Template cleanup** — `_quarto.yml` now renders only `index.qmd`; the
>    earthquake notebook, `la-palma*` assets, the stale `_manuscript/` build, and
>    a 229 MB virtualenv accidentally committed under `Manuscript/path/` could
>    **not be deleted from this sandbox** (mount is not unlink-able). Remove them
>    yourself:
>    ```
>    git rm -r Manuscript/notebooks Manuscript/path Manuscript/_manuscript \
>              Manuscript/images/la-palma-map.png Manuscript/la-palma.csv
>    echo "Manuscript/path/" >> .gitignore
>    ```
> 6. **Journal fit checked** — NClimate Letter limits are ~2,000 words, ≤5
>    display items, ≤30 references. The draft is ~1,310 main-text words, 4
>    figures, ~24 references — within limits. Output formats set to HTML + Word.
>
> **Remaining `[AUTHOR ACTION]` items:** re-run models (§2.6 / §3 step 3),
> regenerate Figure 4, add the repository archive DOI, and a final read-through.

---

Review of the analysis (`Analysis/atlantic_birds_ms.Rmd`) and the manuscript
(`Manuscript/index.qmd`), with the aim of getting the project ready for
submission to *Nature Climate Change* (Letter format). Items marked **[FIXED]**
were edited directly; items marked **[ACTION]** need a decision or a re-run that
requires R, the phylogeny, and the PREDICTS database (none available in this
environment).

The numbers used in the manuscript were taken from the rendered report
`Analysis/atlantic_birds_ms.html`, which contains the actual `brms`/`metafor`
output. No values were invented.

---

## 1. Scientific / statistical issues (highest priority)

### 1.1 The headline body-size result is model-dependent — state it honestly **[ACTION]**
The temporal decline in wing length is clear in the model fitted to the
coalesced wing-length variable (β = −0.66, 95% CrI −1.03 to −0.30; N = 9,387),
but its credible interval **overlaps zero** in the model that propagates
phylogenetic uncertainty over 50 trees and uses right-wing length only
(β = −0.43, 95% CrI −0.96 to 0.13; N = 4,629). These are the two natural
"primary" models and they disagree on whether the effect excludes zero.
Decide which is the primary inference and report the contrast openly. The
manuscript currently does this; confirm the framing is what you want before
submission. Avoid a flat "birds shrank" claim built on a single point estimate.

### 1.2 lnCVR confidence intervals were mis-scaled **[FIXED]**
The per-species forest-plot limits were computed as `estimate ± 1.96 * s2.lnCVR`,
i.e. 1.96 × the *sampling variance*. They should use the *standard error*,
`± 1.96 * sqrt(s2.lnCVR)`. Fixed for both wing and bill. The meta-analytic
estimates from `rma()` are unaffected, because `rma()` is correctly passed the
sampling variances via `vi`. Regenerate Figure 3 after re-running.

### 1.3 Temperature is asserted but never modelled **[ACTION]**
The introduction/abstract describe a warming, drying climate over the study
period, but `Annual_mean_temperature` is **not** a predictor in any body-size
model (the models use scaled *year*), and precipitation is not analysed at all.
Either (a) fit a trend of temperature (and a chosen precipitation variable)
versus year for the sampled localities and report it, or (b) soften the climate
framing to cite regional/global warming from the literature rather than an
in-sample trend. This also resolves the "XX years / which precipitation
variable" placeholders left in the old draft.

### 1.4 Two different wing-length responses across models **[ACTION]**
`brm0` and the imputed model use `conc.wing.length` (N = 9,387); the
phylogenetic-uncertainty model uses `Wing_length_right.mm.` (N = 4,629). This
inconsistency partly drives the difference in 1.1. Harmonise to one response
(the coalesced variable is recommended for power and consistency) and re-run the
phylogenetic-uncertainty model on it, or justify the difference explicitly.

### 1.5 Birds and prey are linked only by shared region/period **[ACTION]**
The trophic interpretation rests on parallel trends, not a joint model. The
strongest test — using predicted arthropod abundance per year as a time-varying
predictor in the wing-length model — was noted as an idea in the Rmd but never
run. Consider adding it; the manuscript flags it as the priority next step.

### 1.6 Latitude direction **[ACTION — check interpretation]**
`scaled_lat = scale(Latitude_decimal_degrees * -1)`, so larger values denote
more southerly (cooler) localities. The negative coefficient means wing length
is *smaller* toward the south — opposite to a naive spatial Bergmann prediction.
This is reported factually in the manuscript; decide how much to interpret it.

---

## 2. Reproducibility & deprecated API (mostly fixed)

### 2.1 `cov_ranef` removed from current brms **[FIXED]**
`cov_ranef = list(Binomial = A)` is deprecated and removed in recent `brms`.
Updated all model calls to pass the phylogenetic covariance via
`data2 = list(A = A)` and a `(1 | gr(Binomial, cov = A))` term.

### 2.2 Other deprecated calls **[FIXED]**
`marginal_effects()` → `conditional_effects()`; `stanplot()` → `mcmc_plot()`;
`posterior_samples()` → `as_draws_df()`; `stat_pointintervalh()` →
`stat_pointinterval()`; ggplot `size=` → `linewidth=` for lines/densities.

### 2.3 Hardcoded absolute paths **[FIXED]**
Removed `load("~/Google Drive/.../edu_joao/brm0.rda")` and the `/Volumes/...`
references in prose. The multi-tree model list should be saved in the repo as
`brm0_multiphylo.rda` (renamed to avoid clashing with the single-tree
`brm0.rda`). The PREDICTS path was already replaced with a relative placeholder.

### 2.4 Missing random seeds **[FIXED]**
Added `set.seed()` before `sample(trees, 50)` and a `seed=` argument to the
`mice()` call, so tree sampling and imputation are reproducible.

### 2.5 Fragile column indexing **[FIXED — comment added]**
`meth[c(12,13)] <- "norm"` breaks if column order changes. Added a comment
recommending name-based indexing; verify the intended variables before running.

### 2.6 Environment capture **[ACTION]**
Add a `sessionInfo()` / `renv.lock` and pin package versions (`brms`, `metafor`,
`mice`, `ape`, `MCMCglmm`, `tidyverse`). The saved `.rda` models were produced by
an older `brms` (output shows `Samples:` and `Eff.Sample`, i.e. pre-2.x naming);
re-running under current `brms` will rename these to `Draws:` / `Bulk_ESS` and
may shift estimates slightly.

---

## 3. Manuscript (`Manuscript/index.qmd`)

### 3.1 Rewritten from template **[FIXED]**
The file was a leftover "La Palma earthquake" Quarto template: placeholder
abstract, plain-language summary, key-points, and an embedded earthquake
notebook. It is now a complete *Nature Climate Change* Letter draft — lead
paragraph ("Here we show…"), Results-and-discussion woven with brief methods,
detailed Methods, data/code availability, and references — using the real
estimates from the rendered report.

### 3.2 Leftover template files to remove **[ACTION]**
`Manuscript/notebooks/data-screening.qmd` and `Manuscript/images/la-palma-map.png`
are earthquake-template artifacts and will break or pollute a manuscript render.
Remove them (and the stale `_manuscript/` build output) before submission.

### 3.3 Figures are placeholders **[ACTION]**
`Manuscript/images/fig-{map,trends,variability,arthropods}.png` are labelled
placeholders so the document renders. Regenerate the real figures from the named
Rmd chunks and overwrite them. Each caption names its source chunk.

### 3.4 Inline author-action notes **[ACTION]**
The qmd contains `[AUTHOR ACTION]` comments at each point that needs your input
(items 1.1, 1.3, 1.4, 2.6 above, plus the repository DOI). Search the file for
`AUTHOR ACTION` and resolve, then delete the leading note block.

### 3.5 Journal fit **[ACTION]**
Confirm format/length against current *Nature Climate Change* guidelines
(Letters are ~3,000 words, with Methods online). Adjust `_quarto.yml` output
formats as needed; `agu-pdf` is currently configured and may not match the
target.

---

## 4. Suggested order of work
1. Decide the primary body-size model and harmonise the response variable (1.1, 1.4).
2. Add or remove the temperature/precipitation trend claim (1.3).
3. Re-run models under current `brms`; save `brm0_multiphylo.rda` in-repo (2.1, 2.6).
4. Regenerate the four figures (3.3) and the corrected lnCVR forest plot (1.2).
5. Resolve `[AUTHOR ACTION]` notes, add repo DOI, remove template files (3.2, 3.4).
6. Final read-through against journal guidelines (3.5).
