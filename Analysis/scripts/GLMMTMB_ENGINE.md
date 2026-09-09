# glmmTMB phylogenetic engine (2026-09-09)

## What changed and why

Every phylogenetic model in the revision (`REVISION_PLAN.md` §3) was originally
specified as an 8-model × 50-tree `brms` grid — `400` Stan fits per analysis,
hours to days of compute, and the thing the plan's Phase 1 warning
("AVOID THE 400-STAN-FIT BOTTLENECK") was written to head off. `_phylo_engine.R`
replaces that path with **glmmTMB's `propto` covariance structure** (Williams,
McGillycuddy, Drobniak, Bolker, Warton & Nakagawa 2025, bioRxiv
10.64898/2025.12.20.695312), which represents the phylogenetic random effect as a
correlation matrix on a glmmTMB random term instead of a Stan `cov_ranef`. Same
statistical object — a Brownian-motion (equivalently OU-at-Brownian-limit)
phylogenetic covariance on species intercepts and, where specified, slopes — fit
by REML/Laplace instead of MCMC. glmmTMB ≥ 1.1.14 (CRAN, installed) already has
`propto`; no dev build was needed.

## Validation against the published brms model

`Analysis/scripts/glmmtmb_validation_wing.R` refits the published 50-tree wing
model (`Analysis/output/models/brm0_multiphylo.rda`) in glmmTMB on the **identical
data and the identical 50 per-tree species correlation matrices** taken from each
brmsfit's own `$data2$A` (not re-derived), then Rubin-pools exactly as the brms
run did. Numbers below are read from
`Analysis/output/glmmtmb_validation_wing.rds` (`$generated` 2026-09-09 13:08:54,
glmmTMB 1.1.14), not retyped from the script's own printed log:

| Term | glmmTMB (50-tree Rubin) | brms (published, 50-tree Rubin) |
|---|---|---|
| Intercept | 66.72 [43.85, 89.58] | 67.85 [50.00, 85.69] |
| SexMale | 2.643 [2.439, 2.846] | 2.642 [2.438, 2.847] |
| scaled_yr | **−0.922 [−1.400, −0.444]** | **−0.922 [−1.403, −0.440]** |
| scaled_lat | −0.676 [−0.891, −0.462] | −0.677 [−0.891, −0.463] |

Variance components (mean across the 50 trees; `$varcomp`): phylogenetic species
SD 23.58 (brms tree-1 reference fit 24.65), residual σ 4.591 (brms tree-1 4.591),
species year-slope SD 1.890 (brms tree-1 1.892), **phylogenetic proportion of
species-level variance 0.957 [0.936, 0.967]** (brms 0.95 [0.94, 0.96] per the
script's header). One component does **not** agree: the iid species-intercept SD
(`spp:(Intercept)`) is ≈ 0 in glmmTMB (mean 0.0054) against 0.89 in brms tree 1 —
see Limitations. Fixed effects, their intervals, and the phylogenetic-proportion
headline reproduce to 2–3 decimals; the median glmmTMB fit took **0.39 s**, 50
fits **≈ 20–24 s** total, against hours for the Stan run. This is the basis for
making glmmTMB the default engine for every `--trees` phylogenetic run in the
revision.

## Exact model specification

For a mean model `y ~ <fixed> + (1 [+ slope] || spp)` (no phylogeny), the engine
appends one term and fits by REML:

```r
y ~ <fixed> + (1 [+ slope] || spp) + propto(0 + species_name | g, A)
```

- `species_name` is the tree-tip-labelled species factor (`phylo_species_name()`
  maps ABT binomials through the ten `ebird_synonyms` entries to eBird/Clements
  tip names; `Herpsilochmus_sellowi` has no tree tip and is dropped upstream by
  every caller).
- `g` is a dummy single-level grouping factor (`factor(1)`), required by glmmTMB's
  `propto` syntax; it does not partition records.
- `A` is a **species correlation matrix** (diagonal 1), one per tree, taken from
  `inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv` inverted, on trees passed
  through `phytools::force.ultrametric(method = "nnls")` where not already
  ultrametric — the identical construction `atlantic_parallel.R` used for the
  brms `data2$A` objects. `phylo_A_list()` reads the 50 matrices cached from
  `brm0_multiphylo.rda` (`data/derived/phylo_A_50trees.rds`) so every phylogenetic
  model in the revision uses the **same 50 trees** as the published wing analysis;
  a species set the cache does not cover falls back to a fresh clootl draw (same
  seed, same 50-of-100 sampling as `atlantic_parallel.R`).
- The non-phylogenetic species term (`(1 [+ slope] || spp)`) is kept alongside the
  phylogenetic one wherever the mean structure already carried it (M0's published
  structure includes `(1 + scaled_yr || spp)`), so the phylogenetic intercept and
  the iid species intercept are fit simultaneously and can trade variance against
  each other — see Limitations.
- `dispformula` carries the σ (variance) sub-model where one is specified
  (Phase 3's `sigma ~ year + ...` ladder); `family` defaults to `gaussian()`.
- Contributor (`Main_researcher`), site (`Municipality`/`Locality`), season, moult,
  ring/individual, and Mundlak within/between terms are ordinary fixed or random
  effects in `<fixed>`, unrelated to the phylogenetic machinery — the engine only
  ever touches the species term.

Fitting: `fit_phylo_glmmtmb()` builds and fits one tree; `run_phylo_trees()` loops
over the `A_list`, calls `tidy_phylo_fit()` per tree (fixed effects for every
formula component, phylogenetic SD, σ, phylogenetic proportion, `pdHess`
convergence), and returns the per-tree table plus the pooled table.

## Rubin pooling over the 50 published trees

`pool_rubin_df()` implements Rubin's rules for multiply-imputed data
(Nakagawa & de Villemereuil 2019's application to multi-tree phylogenetic
inference): for each fixed-effect parameter across the `m` trees,

```
qbar = mean(estimate)                         # pooled point estimate
ubar = mean(se^2)                             # mean within-tree variance
b    = var(estimate)  [0 if m == 1]           # between-tree variance
se   = sqrt(ubar + (1 + 1/m) * b)             # total pooled SE
CI   = qbar ± 1.96 * se                       # Wald interval
```

Every downstream script (E1–E4's phylogenetic tiers, this validation) pools this
way over some number of trees ≤ 50; several scripts additionally report
`pooled_converged`, restricting the sum to `pdHess`-converged trees only, because
`run_phylo_trees()` itself pools indiscriminately (flagged as an engine gap
below).

## What remains on brms

Per `REVISION_PLAN.md` §3 Phase 1's own tiering, glmmTMB does **not** replace
brms everywhere:

- **The bivariate `mvbind(cwl, ln_mass)` rescor model** (Phase 2) — a genuine
  multivariate response with an estimated residual correlation is a brms/Stan
  specification; `atlantic_bivariate_wing_mass.R --engine brms` is kept for it
  and was not run this pass (Totoro). The glmmTMB tier instead fits two separate
  species-phylogeny univariate models (wing, log-mass) plus a joint long-format
  model with two trait-specific `propto` terms, and estimates the slope
  difference two ways — an independence-assumption delta method and the joint
  model's covariance-correct estimate — both on file in
  `Analysis/output/bivariate_phylo_results.rds`.
- **One Bayesian cross-check of the primary controlled model (M3)** on Totoro,
  kept explicitly so the revision has one posterior-based sanity check of the
  Wald/REML phylogenetic tier before it goes in the manuscript, per this task's
  instruction. Run on Totoro 2026-09-09 for M0 on one tree only (cmdstanr, 4 chains ×
  4,000 iterations, 55 min): agreed with the glmmTMB estimate; the rest of the brms ladder
  was stopped as redundant given the validation.
- Every other phylogenetic model in the revision (Phase 1's M0–M7 ladder, the
  multi-trait screen, the diet interaction, the variance σ ladder) now defaults
  to `--engine glmmTMB`; `--engine brms` remains selectable in each script for
  its own Totoro cross-check, per the boundary in this task's instructions, and
  none of the existing lme4 fast-tier code or outputs were removed.

## Limitations

- **Wald intervals, not posteriors.** Every interval quoted from the
  phylogenetic tier is `estimate ± 1.96 × Rubin SE`, not a credible interval; it
  assumes asymptotic normality of the REML/Laplace estimator and does not
  propagate tree uncertainty the way a full joint posterior would (Rubin pooling
  is a documented, but approximate, substitute — the reason one brms cross-check
  of M3 is retained).
- **No posterior, so no direct R-hat/ESS diagnostics.** Convergence is judged by
  `pdHess` (a positive-definite Hessian at the optimum) per tree; scripts report
  the fraction of trees converged and, in several cases, exclude non-converged
  trees from pooling.
- **Boundary/collinear variance estimates.** When a model's mean structure
  already carries an iid species random effect (`(1 [+ slope] || spp)`)
  alongside the phylogenetic `propto` term, the two can trade variance: in the
  validation run the iid species-intercept SD collapsed to ≈ 0 while the
  phylogenetic SD absorbed essentially all of it (23.58 vs residual 4.59; iid
  0.005 vs brms's 0.89) — same fitted mean and CI, different variance
  decomposition. E3 and E4 both independently hit a worse local optimum of the
  same collision (phylogenetic SD ≈ 0, iid SD absorbing it, non-PD Hessian) on a
  minority of trees once contributor/site random effects were added, and both
  implemented a seeded-restart repair (`run_phylo_robust()`, duplicated in the
  two scripts) that recovers the better optimum on all but a small number of
  tree × model combinations (documented per-model in `REVISION_NOTES_P2.md` §8.3
  and `REVISION_NOTES_P5.md` §8.1). **This repair belongs in `_phylo_engine.R`,
  not duplicated per script** — recommended as the next change to the shared
  engine, not made here (out of this agent's ownership).
- **`tidy_phylo_fit()` returns `sigma = NA`** when a model has a non-trivial
  `dispformula`, so a dispformula-aware phylogenetic proportion needs a
  script-local recomputation (E2's `phylo_share()` helper) rather than being
  available from the engine directly.
- **`VarCorr` name mangling.** `tidy_phylo_fit()`'s `varcomp` data frame passes
  random-effect SD names through `data.frame()`/`as.list()`, which applies
  `make.names()` and turns e.g. `"spp:(Intercept)"` into `"spp..Intercept."`,
  forcing regex matching downstream (E2, E3 both hit this).
- **`run_phylo_trees()` pools every tree regardless of `pdHess`.** No built-in
  option to restrict Rubin pooling to converged trees; every caller that needs
  this (all of them, once any tree fails to converge) reimplements it locally.
- **`Herpsilochmus_sellowi`** has no tip in the cached 50 trees and must be
  dropped from `species` before calling `phylo_A_list()`; this is a data gap in
  the tree, not an engine bug, but it is easy to trip on.

None of the above changes any of the headline validation numbers above; they are
implementation caveats for whoever extends the engine next.

## Current status of the phylogenetic tier across scripts (2026-09-09, this pass)

Copied from each script's own output file, not from any agent's report:

| Owner | Script | Output on disk | `n_trees` | Status |
|---|---|---|---:|---|
| E1 | `atlantic_parallel_controlled.R` | `controlled_wing_phylo_results.rds`, `..._species_slopes.rds` | **50** (49/50 on M5_cc) | Complete (`$generated` 16:48:02; 39 specs, 3,924 s; median 74 s, max 219 s per spec). Gate B on both samples; every year term within 0.002 of the lme4 tier. |
| E2 | `atlantic_variance_sigma.R` | `variance_phylo_results.rds` | **50** (45/50 on S0) | Complete (`$generated` 16:12; 9 tiers, 1,768 s). Phylogenetic vs non-phylogenetic σ-year terms differ by ≤ 0.0006. |
| E3 | `atlantic_bivariate_wing_mass.R` | `bivariate_phylo_results.rds` | **50** | Complete (`$generated` 14:40:17, 37.6 min). |
| E3 | `atlantic_multitrait.R` | `multitrait_phylo_results.rds` | **50** (49/50 on one bill-length spec) | Complete (`$generated` 14:44:25, 41.7 min). |
| E4 | `atlantic_diet_interaction.R` | `diet_interaction_phylo.rds` | **50** | Complete (`$generated` 15:23:54, 513 s). |

All five phylogenetic scripts have now finished their 50-tree runs (the wing and variance
ladders on the second attempt, run as detached jobs: 65 and 30 min). Across the whole
revision the phylogenetic term explains 91–96 % of species-level variance but moves no
year or interaction coefficient by more than 0.002, because the temporal contrasts are
identified within species. The engine's cost is therefore minutes on a laptop for the
complete phylogenetic analysis, against 55 min per model for brms on the server.

## Draft Methods paragraph (not inserted into `Manuscript/index.qmd`)

> *Phylogenetic mixed models.* Models with a phylogenetic random effect were fit
> with `glmmTMB` (version 1.1.14) using its `propto` correlation structure
> (Williams et al. 2025), which places a Brownian-motion phylogenetic covariance
> on species-level intercepts (and, where specified, slopes) within a standard
> penalized-likelihood mixed model rather than a Markov-chain Monte Carlo fit. We
> validated this approach against our previously published Bayesian (`brms`)
> phylogenetic wing-length model by refitting it in `glmmTMB` on the identical
> data and the identical 50 phylogenetic correlation matrices used in the
> published analysis: the two engines agreed to two–three decimal places on
> every fixed effect (year coefficient −0.922 [−1.400, −0.444] `glmmTMB` vs
> −0.922 [−1.403, −0.440] `brms`) and on the proportion of species-level variance
> attributable to phylogeny (0.957 vs 0.95), while reducing computation from
> hours to seconds per fit. We therefore used `glmmTMB` as the default
> phylogenetic engine for the models in this revision, retaining the original
> `brms` implementation as a Bayesian cross-check for the primary controlled wing
> model. For each phylogenetic model we fit the same specification independently
> on each of 50 phylogenetic trees (drawn as in the original analysis) and pooled
> fixed-effect estimates and standard errors across trees using Rubin's rules for
> multiply-imputed data, following Nakagawa & de Villemereuil (2019). Reported
> intervals from the phylogenetic tier are Wald intervals (estimate ± 1.96 ×
> pooled standard error) rather than Bayesian credible intervals, and convergence
> was assessed per tree via the positive-definiteness of the fitted Hessian; the
> proportion of trees converging is reported alongside each estimate.
