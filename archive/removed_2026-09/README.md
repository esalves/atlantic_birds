# Analyses retired from the manuscript, September 2026

These analyses were part of earlier drafts and are **not** referenced by
`Manuscript/index.qmd` or `Manuscript/supplementary.qmd` any more. Nothing here
was found to be wrong — each was dropped for scope or because the evidence did
not support the weight the text put on it. The code and results are kept so the
decisions stay auditable and so the material can be restored if a reviewer asks
for it.

Retired on 2026-09-21, after comments from A. Mizuno.

## 1. Litter-ant / regional invertebrate context

| file | was |
| :--- | :--- |
| `scripts/atlantic_ants_index.R` | negative-binomial GLMMs of litter-ant richness (ATLANTIC ANTS, 855 campaigns, 1994–2018) |
| `output/ants_results.rds`, `ants_results.md`, `ants_annual_index.rds` | fitted year effects and the annual index |
| `figures/ants_events_by_year.png`, `ants_year_effect.png` | diagnostic and effect figures |

**Why dropped.** The analysis appeared as "exploratory regional context"
(Supplementary Note S7) behind the diet result. Litter-ant species richness is
not a measure of prey abundance or biomass, ground-dwelling ants are one
component of these birds' diets, and the estimate was highly sensitive to
sampling effort and site composition. It could not support the inference it was
being asked to contextualise. See also the earlier PREDICTS retirement in
`../predicts/`.

## 2. Distributional SEM (drmSEM)

| file | was |
| :--- | :--- |
| `scripts/atlantic_drmsem.R` | distributional causal DAG over warming, prey availability and morphology |
| `output/drmsem_results_{arthro,noarthro}.{rds,md}` | path coefficients and effect tables |
| `output/drmsem_phylo_prep_*.rds`, `drmsem_phylo_results_*.rds`, `drmsem_effects_cache_*.rds` | intermediate fits and caches |
| `figures/drmsem_*.png`, `images/fig-drmsem.png` | DAG, coefficient forest, sigma curves |

**Why dropped.** Supplementary Note S5 was exploratory and its prey node drew on
the litter-ant series retired above, so it went with it.

## 3. Century-scale museum series, climate breakpoint, sampling expansion

| file | was |
| :--- | :--- |
| `scripts/secular_museum_expansion.R` | museum study skins (n = 3,439, 1884–2017), the 1980 piecewise climate breakpoint, and the community-wide expansion under Jirinec et al. (2021) sampling criteria |
| `scripts/secular_threshold_sensitivity.R` | sensitivity of the secular result to filtering thresholds |
| `output/secular_expansion/` | `secular_expansion_results.rds` + `REPORT.md` |
| `figures/secular_expansion/`, `images/fig-s-climate-{morphometry,mass,allometry}.png`, `fig-s-secular-trajectories.png`, `fig-s-model-comparison.png` | the four supplementary figures |

**Why dropped.** The combined museum-plus-live series could not carry the claim
made of it: the apparent secular decline was carried entirely by the live
captures, vanished in museum specimens alone, and did not survive restriction to
adults. Museum preparation shrinks feathers and keratin and no sample was
measured by both routes, so there was no way to calibrate the offset. The
community-wide expansion (Supplementary Note S4 and old Table S1) went with it,
since both lived in the same Results and Discussion sections.

`images/fig-thermal-allometry.png` was already unreferenced before this pass and
is archived here for tidiness.

## Restoring any of this

`git mv` was used for tracked files, so `git log --follow` works on each path.
To bring one back, move it to its original location and re-add the section to
the relevant `.qmd`.

## Two loose ends left in the live tree, deliberately

- `Analysis/scripts/make_figures.R` still carries its drmSEM figure builders
  (modes `arthro` / `noarthro`). It guards every read with `file.exists()` and
  now prints `skip [mode]: drmsem_results_<mode>.rds not found`, so the export
  mode that `make_figures.py` depends on is unaffected.
- `Analysis/scripts/referee_reruns.R` still computes block **E6/E7** (the
  secular series) into `referee_reruns.rds`. It recomputes from source data and
  does not read anything archived here, so it runs unchanged. The manuscript no
  longer quotes those numbers, but the block is part of the written referee
  response and was left intact.
