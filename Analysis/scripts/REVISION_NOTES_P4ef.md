# Revision notes — Phase 4e/4f (retire PREDICTS, demote drmSEM), 2026-09-09

Companion to `REVISION_PLAN.md` §1 "§2.6 PREDICTS", "§2.7 drmSEM" and §3
"Phase 4e/4f". Records exactly what was moved, what the drmSEM script now does,
what was and was not executed, and what the manuscript / figure / housekeeping
agents must change. Style follows `LIVE_ONLY_MIGRATION.md` and
`REVISION_NOTES_P0.md`.

Executed locally with R 4.6.0, drmTMB 0.1.4, drmSEM 0.5.0, prepR4pcm 0.5.0.9000,
brms 2.23.0, dplyr 1.1.x, ape. **No full SEM was run** (the task forbids it;
see §2.3 — the climate-file staleness guard that also blocked it in the first
pass no longer applies, because P5 has since regenerated the file). Every model
number quoted below is either (a) the archived June-2026 result, labelled as
such, or (b) a smoke test on a 1,500-record subset, labelled as such and not a
result. Every check was re-executed from scratch in a second pass on
2026-09-09 (§7); the earlier pass's claims were not taken on trust.

---

## 1. Phase 4e — PREDICTS trend retired

### 1.1 What moved (`git mv`, staged, NOT committed)

| From | To |
|---|---|
| `Analysis/output/models/brm_arthro_abund.rda` | `archive/predicts/brm_arthro_abund.rda` |
| `Analysis/output/models/brm_arthro_abund_no_ants.rda` | `archive/predicts/brm_arthro_abund_no_ants.rda` |
| `Analysis/output/models/brm_no_grass.rda` | `archive/predicts/brm_no_grass.rda` |
| `Analysis/output/models/brm_no_ants_no_grass.rda` | `archive/predicts/brm_no_ants_no_grass.rda` |
| `Analysis/data/derived/dat_test.rds` | `archive/predicts/dat_test.rds` |
| `Analysis/data/derived/dat_no_ants.rds` | `archive/predicts/dat_no_ants.rds` |
| `Analysis/data/derived/dat_no_grassland.rds` | `archive/predicts/dat_no_grassland.rds` |

Plus a new `archive/predicts/README.md` (why, counts re-derived from the
archived file, the legacy Rmd chunks that still point at the old paths, and the
two-file follow-up below). All seven files were tracked in git, so the moves are
renames in the index; nothing was deleted or rewritten.

### 1.2 What did NOT move, and why

`Analysis/output/arthropod_estimates.rds` and `Manuscript/images/fig-arthropods.png`
are **still read by `Manuscript/index.qmd`** (line 60
`arth <- readRDS("../Analysis/output/arthropod_estimates.rds")`; line 196 the
`{#fig-arthropods}` figure — line numbers from the 2026-09-09 ~11:30 working
copy; P7 is editing concurrently, so search for the anchors, not the numbers). The task condition was "move only if nothing
outside the arthropod pipeline reads them"; the manuscript does, and moving
them now would break `quarto render` for the P7 agent mid-edit. **P7: once the
arthropod Results section and Fig. 5 are deleted, run**

```bash
git mv Analysis/output/arthropod_estimates.rds  archive/predicts/arthropod_estimates.rds
git mv Manuscript/images/fig-arthropods.png     archive/predicts/fig-arthropods.png
```

(also listed in `archive/predicts/README.md`). The raw download
`Analysis/data/raw/predicts_extract.rds` is git-ignored raw data and stays put.

### 1.3 Grep of readers (2026-09-09, excluding `.git`, `_freeze`, `_manuscript`, `archive`)

| Reader | Files referenced | Action |
|---|---|---|
| `Manuscript/index.qmd` L60, L196 | `arthropod_estimates.rds`, `fig-arthropods.png` | P7 deletes (§3); then the two `git mv` above |
| `Analysis/scripts/atlantic_drmsem.R` | `dat_no_grassland.rds`, `brm_no_grass.rda` (commented) | guarded (§2) |
| `Analysis/scripts/atlantic_birds_ms.Rmd` chunks `load_predicts_brasil`, `arthro_model_summary`, `arthro_model_no_ants_summary`, `arthro_model_no_grass_summary`, `save_arthropod_estimates`, `compare_ant_models`, `export_fig_arthropods` (lines 882–1130) | all seven files + `predicts_extract.rds` | **legacy, left unedited** — listed with line numbers in `archive/predicts/README.md`. These chunks are `eval = TRUE` for loading, so knitting the Rmd now fails at chunk `load_predicts_brasil`; the Rmd is not in the pipeline |
| `Manuscript/make_figures.py` L10 (docstring), `Manuscript/RENDER_STEPS.md` L32/L59, `Analysis/output/RESULTS_SUMMARY.md` §3 + file list, `REPO_STRUCTURE.md` tables, `Analysis/scripts/LIVE_ONLY_MIGRATION.md` L55, `README.md` L146 | documentation mentions only | Phase 8 housekeeping (§5) |

### 1.4 Numbers re-derived from the archived data (for the README and response letter)

From `archive/predicts/dat_no_grassland.rds` (63,150 rows, 56 studies, 41
sources, 725 sites, forest biomes, ants included):

- **53 of 56 studies have all records in a single calendar year** (year of
  `Sample_midpoint`); the three multi-year studies span 2005–2007 (Raub, Brazil)
  and 2006–2007 (Poveda 1 & 2, Colombia). **The plan says 54 of 56** — the
  definition it used is not recorded; by `Sample_start_earliest`/`Sample_end_latest`
  year it is 31 of 56. Either way there is no within-study temporal information
  beyond two adjacent years.
- **Atlantic Forest ecoregions** (Serra do Mar, Bahia Coastal, Pernambuco
  Coastal/Interior, Alto Paraná; mangroves excluded): **17 studies, 6,270 rows,
  sampled 1998–2009.** The plan's "12 studies / 1,652 rows" did **not**
  reproduce; the 1998–2009 span did. Quote the reproduced numbers and say how
  they were defined, or leave counts out and keep the qualitative point.
- Year range of the file 1996–2011; largest study Cabra 2006 (Colombia)
  11,310 rows — reproduces.
- `arthropod_estimates.rds` (still in place): main β = −1.40 [−1.46, −1.34]
  (×0.25), N = 63,150; no-ants β = −1.22 [−1.29, −1.15] (×0.30), N = 53,139.
  The "shallower without ants" sentence is wrong on its own terms: −1.22 is a
  *weaker* decline than −1.40, so ants steepened the pooled slope.

---

## 2. Phase 4f — `atlantic_drmsem.R` (v3 → v4)

### 2.1 Changes

1. **Header rewritten** (WHAT / STATUS / DISCLOSURE / D-SEP CLAIMS / SIGMA
   WORDING / INPUTS / OUTPUTS / RUN / HISTORY): states that drmSEM is a
   supplementary exploratory analysis, not primary evidence; that drmSEM/drmTMB
   and prepR4pcm are co-authored by E.S.A.S. and are to be disclosed as
   unpublished development software (versions recorded); lists the
   d-separation claims to report (§2.4); records the σ-channel wording rule
   (year and temperature effects on residual SD are separable; not "linked to
   thermal extremes").
2. **`INCLUDE_ARTHRO` is FALSE and cannot be flipped by editing the file.**
   A guard `stop()`s if the constant is TRUE without the command-line flag
   `--arthro-override=PREDICTS_RETIRED` (verified: editing the line to TRUE
   halts with the message). The override exists only to reproduce the archived
   June-2026 numbers for the response letter and emits a warning.
3. **`dat_no_grassland.rds` load is guarded**: not loaded at all in the default
   mode (message: retired, archived); under the override it is looked for in
   `archive/predicts/` then `data/derived/`, and the script stops with a
   pointer to the README if neither exists. The commented section-9 code now
   points at `archive/predicts/brm_no_grass.rda`.
4. **Wing node gains `(1 | Main_researcher)`.** drmTMB accepts a second i.i.d.
   random intercept, and also `relmat(1 | sp_tip, K) + (1 | Main_researcher)`
   in one node (both tested 2026-09-09 on a 1,500-record subset; a random
   `ape::rtree` was used purely to exercise the relmat code path). drmSEM strips
   random-effect terms from the causal edge set, so the d-sep claim set is
   unchanged (one claim, df = 2, with or without the contributor term). The
   mean formula is built once by `wing_mu_formula()` / `wing_bf()` and used in
   section 4 (i.i.d.) and section 10 (relmat); `drmTMB::bf()` uses
   `substitute()`, so the built formula is spliced in with `bquote()`. If
   `Main_researcher` were missing the script fits species-only intercepts and
   records `src_intercept = FALSE` in the results. No contributor term is put in
   the σ submodel here (the Phase 3 brms σ model carries that).
5. **Climate-file staleness guard.** `passer90.rda` is always loaded;
   `passer90_climate.rds` is refused if it lacks `ID_ABT` (the pre-Phase-0
   file: 15,332 rows / 89 species incl. museum specimens — this is the file
   currently on disk locally). `--allow-stale-climate` (implied by `--smoke`)
   falls back to the static WorldClim climatology in `passer90.rda` with the
   existing warning. A current file is `semi_join`ed to `passer90` on `ID_ABT`
   (live-only guarantee) and any missing provenance column is attached.
   `temperature_source` is written to the results. **Second pass:** P5's
   regenerated `passer90_climate.rds` (written 2026-09-09 10:47; 12,571 rows /
   73 species / 46 columns incl. `ID_ABT`, `Main_researcher`, `rec_tmean`,
   `scaled_tmean`; `Status` all `live`) now passes the guard, and the script
   records `temperature_source = "worldclim_record_level"`. The static
   fallback is no longer exercised unless the file is removed.
6. **Flags**: `--smoke` (1,500-record subsample, `RUN_PHYLO = FALSE`,
   `EFFECT_B = EFFECT_NSIM = 20`, **all outputs redirected to a temp dir** so a
   smoke run can never overwrite `drmsem_results_noarthro.rds`, which the
   manuscript reads), `--no-phylo`, `--allow-stale-climate`,
   `--arthro-override=PREDICTS_RETIRED`.
7. `MODEL_TAG` moved to the end of section 3 and bumped to `v4-<mode>-src`
   (+`-SMOKE`), so all three caches (`drmsem_effects_cache_*`,
   `drmsem_phylo_prep_*`, `drmsem_phylo_results_*`) are invalidated on the next
   real run; the data signature now includes the contributor count.
   `scaled_yr` / `scaled_lat` are coerced from `scale()` matrices to numeric.
8. Results (`drmsem_results_<mode>.rds/.md`) now carry `status`, `smoke`,
   `src_intercept`, `n_contributors`, `temperature_source`, `wing_mu_formula`,
   `dsep_claims_to_report`, and `software` (versions + disclosure text) so P7
   can print the disclosure and the claim list inline. The `.md` header line
   flags smoke output as "NOT A RESULT".
9. `plot(sem_fit)` (side-effect drawing) only runs interactively, so
   `Rscript` no longer leaves a stray `Rplots.pdf`.

### 2.2 What was executed (all re-run in the second pass, 2026-09-09; logs in §7)

- `Rscript Analysis/scripts/atlantic_drmsem.R --smoke`: **ran end to end in
  3.9 s** (sections 1–8, 11, 12; section 10 skipped by design), twice — before
  and after the small `filter()` reorder in §7, with identical output. Uses the
  regenerated time-resolved climate file (12,571 live-only records loaded);
  PREDICTS not loaded (message: retired, archived); subsample **1,500 records /
  72 species / 39 contributors**, years 1995–2018; wing-node formula
  `wing_length ~ scaled_tmean + Sex + scaled_lat + (1 | spp) + (1 | Main_researcher)`;
  one d-sep claim (`Sex ⟂ scaled_tmean | {scaled_yr, scaled_lat}`), Fisher's C
  on df = 2; results `.rds`/`.md` and four figures written to
  `tempdir()/drmsem_smoke/`; **nothing in the repo written** (mtimes of all
  `Analysis/output/drmsem_*` and `Analysis/figures/drmsem_*` unchanged; no
  `Rplots.pdf`). **Smoke numbers are not results and are not quoted here.**
- `INCLUDE_ARTHRO <- TRUE` without the override (on a scratch copy of the
  script run from the repo root): stops with "INCLUDE_ARTHRO is retired
  (Phase 4e)…" before any data are loaded.
- Standalone drmTMB/drmSEM tests (scratch scripts, 1,500-record subsets of the
  Phase 0 `passer90.rda`, seed 20260909, not in the repo): (i) a wing node with
  `(1 | spp) + (1 | Main_researcher)` fits in drmTMB (convergence 0; both
  random-effect SDs estimated) and inside `drm_sem()`, where `dsep()` still
  returns exactly one claim on df = 2; (ii) `relmat(1 | sp_tip, K) +
  (1 | Main_researcher)` in one node also fits, alone and inside `drm_sem()`
  (a random `ape::rtree` was used purely to exercise the code path — its
  numbers mean nothing). Each test runs in ~1 s.
- `parse()` of the edited script: OK.

### 2.3 What was NOT executed and why

- **No full SEM** (task instruction). Note the change since the first pass:
  the climate-file guard no longer blocks it, so a plain
  `Rscript Analysis/scripts/atlantic_drmsem.R` **will now run the full
  supplementary model** (sections 1–12, including the 50-tree section 10) —
  this is the Totoro job. Not run locally. Section 10 was therefore exercised
  only through the standalone `relmat()` test above; the `bquote()`-spliced
  `wing_bf("relmat(1 | sp_tip, K = K)", env = environment())` call inside
  `fit_wing_node()` / `fit_phylo_sem()` has not been executed in situ.
- Server run order (REVISION_PLAN §4) is unchanged: `climate_extraction.R`
  first (done locally by P5), then `Rscript Analysis/scripts/atlantic_drmsem.R`
  (add `--no-phylo` for a fast first pass). Expected runtime is dominated by
  section 10. All three `drmsem_*_noarthro` caches self-invalidate (new
  `MODEL_TAG` `v4-noarthro-src`, new data signature).
- The June-2026 archived full-record result (N = 7,810, 72 spp, no contributor
  intercept): C = 0.51, df = 2, P = 0.78; σ paths year +0.035 (SE 0.008),
  temperature −0.193 (SE 0.009); year → tmean −0.005 (SE 0.006). These are
  **superseded** and must not be quoted once the v4 run exists; the
  `drmsem_results_noarthro.rds` currently on disk (2026-06-15) is that old fit.

### 2.4 d-separation claims to report (Supplement)

Full-record graph: exactly one missing edge → one claim on Fisher's C with
df = 2: **Sex ⟂ scaled_tmean | {scaled_yr, scaled_lat}**. Report C, df, P, the
claim's LR, N, species and contributors, and state that random-effect terms
(`(1|spp)`, `(1|Main_researcher)`, `relmat()`) are not causal edges. The
retired arthropod graph had two claims on df = 4 and was rejected (C = 13.5,
P = 0.009); it is not reported.

---

## 3. What the manuscript agent (P7, `Manuscript/index.qmd`) must change

Line numbers are from the 2026-09-09 ~11:30 working copy (P7 is editing the
file concurrently, so locate each item by the quoted anchor). Items already
fixed by P7 at that time are marked.

1. **Delete the arthropod Results section** `## Arthropod prey declined
   steeply over the same region and period` (L192; paragraph L194 with every
   `arth$...` inline chunk, including the "Excluding ants gave a shallower
   decline" logic error) and **Fig. 5**
   `![...](images/fig-arthropods.png){#fig-arthropods}` (L196). Remove
   `arth <- readRDS("../Analysis/output/arthropod_estimates.rds")` (L60).
   Replace with the Phase 4a/4b indices (P4a outputs) framed as
   community/occupancy indices with provenance controls, and state plainly
   that no multi-decadal arthropod monitoring exists for the biome
   (REVISION_PLAN §7). Then run the two `git mv` in §1.2.
2. **Drop the rejected arthropod SEM everywhere**: remove
   `drm_a <- readRDS("../Analysis/output/drmsem_results_arthro.rds")` (L64)
   and every `drm_a` inline value — the paragraph "To bring arthropod prey
   directly into the graph, we refitted it on the subset of years with
   arthropod coverage…" (L216); the sentence "Where arthropod prey can be
   brought into the model, its steep decline — estimated independently from
   the PREDICTS records above — shows no direct link to wing length" (L218);
   the Methods subsection `## Arthropod abundance` (L282–284); in
   `## Distributional causal model` the graph description "arthropod
   abundance as a function of year, temperature and latitude; … arthropod
   abundance, sex and latitude" (L288), the whole paragraph "The arthropod
   node used the observed log abundance per unit sampling effort…" (L290),
   and "To bring arthropod abundance directly into the graph we additionally
   fitted a variant… so we report its coefficients only as a qualitative
   check" (L292). "We fitted the graph at two scopes" (L292) becomes one scope.
3. **Discussion sentence relying on the rejected graph** (L234): "Here, by
   integrating PREDICTS arthropod surveys into an explicit distributional
   piecewise structural equation model, we demonstrated that despite a steep
   regional collapse in invertebrate prey, arthropod availability had no
   direct causal connection to wing length (path = …, P …)" — delete, together
   with the rest of that paragraph's "arthropod collapse" reasoning; the diet
   interaction null and the Phase 4a/4b indices carry the trophic test. Also:
   the caution paragraph "We caution that our analytical framework links avian
   morphology and invertebrate prey…" (L238) — rewrite around the new indices;
   "and the PREDICTS database [@hudson2017database] are publicly…" in
   `## Data and code availability` (L296) — replace with ATLANTIC ANTS / GBIF;
   the PREDICTS mention in the Introduction (L154). (The first-pass anchor
   "arthropod prey abundance fell steeply … robust to how ants were treated"
   is no longer in the file — already removed.)
4. **Demote drmSEM to the Supplement** (Results `## A causal model separates
   the prey and temperature pathways`, L208–220 → Supplementary text; Methods
   `## Distributional causal model {#sec-drmsem}`, L286–292 → Supplementary
   methods): retitle without "prey"; say it is supplementary and exploratory;
   drop the P-values (`.dp(...)` calls at L214) and report intervals
   (`drmsem_results_noarthro.rds$paths_raw` has `estimate`/`std.error`); state
   the single d-sep claim (§2.4; also in `results$dsep_claims_to_report`); add
   the **disclosure** that drmSEM/drmTMB and prepR4pcm are co-authored by
   E.S.A.S. and unpublished (text in `results$software$disclosure`, versions
   in `results$software`), plus a one-sentence validation note (drmSEM i.i.d.
   paths agree in sign and magnitude with the brms fits on the same data).
   Replace "Species entered as i.i.d. random intercepts" (L288) with species
   **and contributor (`Main_researcher`)** random intercepts on the wing node.
5. **Resolve the σ contradiction.** Already done at L232 ("the two
   associations are separable: variability rose with year, not with
   temperature") — keep. Still to do at L214: delete "This temperature →
   variance association is the single most robust signal across every drmSEM
   fit" (there is now one fit) and "corroborating the lnCVR result under an
   explicit causal model" → "consistent in sign with the brms σ model (Phase
   3), which is the primary variance evidence". "linked to thermal extremes"
   is no longer in the file.
6. **Fig. 7 (`fig-drmsem`) caption** (L220): move to Supplement; drop
   "of the hypothesised causal graph" (→ "exploratory graph"); the sentence
   "the arthropod-coverage variant (whose graph was rejected by the
   directed-separation test) is not drawn" goes (nothing to allude to);
   replace "solid = P < 0.05, dashed = not significant" with an
   interval-based convention (or say lines are coloured by standardised
   coefficient with 95 % intervals in Table S-x); add that species and
   contributor random intercepts are part of the wing node but are not causal
   edges and are not drawn; rephrase "rises with year and falls in warmer
   years" as two independent associations. Numbers in the caption must come
   from the v4 run (`drm_f` re-read after Totoro).
7. Title (L2 "rising prey scarcity"), abstract (L49 prey-collapse sentence),
   keywords (L46 "trophic interactions"): remove or reframe as
   tested-and-rejected, per REVISION_PLAN §3.1.

## 4. What the figures agent (P6) must change

- `Analysis/scripts/make_figures.R` regenerates `fig-drmsem.png` from
  `drmsem_results_noarthro.rds` (mode `noarthro`) and also iterates over
  `arthro`; the `arthro` mode should be dropped or made conditional on the file
  existing (it will not be regenerated). The DAG builder's node list still
  contains `arthro_obs_std`; harmless (nodes absent from the paths are not
  drawn) but the `p.value < 0.05` solid/dashed encoding must go if the caption
  switches to intervals. Move the figure to the Supplement.
- `Manuscript/RENDER_STEPS.md` paragraph on `fig-arthropods.png` (L32) is
  obsolete.

## 5. Housekeeping (Phase 8) — documentation mentions to update

`REPO_STRUCTURE.md` (old→new table rows for `dat_*.rds`, `brm_no_grass*.rda`,
`brm_arthro_abund*.rda`; add `archive/predicts/`), `Analysis/scripts/LIVE_ONLY_MIGRATION.md`
L55, `Analysis/output/RESULTS_SUMMARY.md` §3 and file list,
`Manuscript/make_figures.py` docstring L10, `README.md` L146
(`predicts_extract.rds` description), `REVISION_PLAN.md` §4 run order comment
("INCLUDE_ARTHRO = FALSE" is now the enforced default).

## 6. Inputs the drmSEM run needs from other agents

- **P5**: regenerated `Analysis/data/derived/passer90_climate.rds` from the
  Phase 0 `passer90.rda` — **delivered** (2026-09-09 10:47; 12,571 rows / 73
  spp, carries `ID_ABT`, `Main_researcher`, `rec_tmean`, `scaled_tmean`); the
  guard in `atlantic_drmsem.R` now passes and the smoke run used it.
- **P4a** (optional, future hook only): a locality × year ants index would
  enter as a new node joined on `Municipality × Year`; not wired in because the
  P4a outputs were not final. The PREDICTS switch is not the place for it.

## 7. Verification log — second pass (2026-09-09, R 4.6.0, drmTMB 0.1.4, drmSEM 0.5.0, prepR4pcm 0.5.0.9000, brms 2.23.0)

All commands run from the repo root; scratch files live in the session
scratchpad, not the repo.

| Check | Command | Result |
|---|---|---|
| Two i.i.d. intercepts on the wing node | `Rscript <scratch>/test_drmtmb_2re.R` | exit 0, 1.1 s; drmTMB convergence 0, `sd(1|spp)` and `sd(1|Main_researcher)` both estimated; `dsep()` = 1 claim, df 2 |
| `relmat()` + contributor intercept | `Rscript <scratch>/test_relmat_2re.R` | exit 0, 1.1 s; single node and `drm_sem()` both fit (random tree — numbers meaningless) |
| Script parses | `Rscript -e 'parse("Analysis/scripts/atlantic_drmsem.R")'` | OK |
| Smoke, end to end | `Rscript Analysis/scripts/atlantic_drmsem.R --smoke` | exit 0, 3.9 s; time-resolved temperature used; 1,500 rec / 72 spp / 39 contributors; outputs in `tempdir()`; repo untouched |
| `INCLUDE_ARTHRO` guard | scratch copy with `INCLUDE_ARTHRO <- TRUE`, run with `--smoke` and no override | halts: "INCLUDE_ARTHRO is retired (Phase 4e)…" |
| Default (no flags) | not run | would now fit the full SEM (guard passes) — forbidden here; Totoro job |
| `git status` | — | 7 renames staged (`R  … -> archive/predicts/…`), `archive/predicts/README.md` untracked, `atlantic_drmsem.R` modified; nothing committed |

**Small code fix in this pass** (`atlantic_drmsem.R`, section 3 "SEM dataset"
block): the `mutate(scaled_yr = as.numeric(...), scaled_lat = as.numeric(...))`
now precedes the `filter(!is.na(scaled_lat), ...)`, because dplyr ≥ 1.1
deprecates filtering on one-column `scale()` matrices (it was a warning, not an
error; the sample is identical — the smoke Fisher's C is 3.412 on df 2 before
and after).

**PREDICTS counts in `archive/predicts/README.md` re-derived** from
`archive/predicts/dat_no_grassland.rds` in this pass: 63,150 rows; 56 studies /
41 sources / 725 sites / 23 methods; 53 of 56 studies single calendar year by
`Sample_midpoint` (multi-year: Raub 2005–2007, Poveda 1 & 2 2006–2007);
31 of 56 by `Sample_start_earliest`/`Sample_end_latest` year; midpoint years
1996–2011; Atlantic Forest ecoregions (Serra do Mar, Bahia Coastal, Pernambuco
Coastal + Interior, Alto Paraná; mangroves excluded) 17 studies / 6,270 rows /
1998–2009; largest study Cabra 2006, 11,310 rows. All match the README. The
plan's "54 of 56" and "12 studies / 1,652 rows" still do not reproduce under
any definition tried.
