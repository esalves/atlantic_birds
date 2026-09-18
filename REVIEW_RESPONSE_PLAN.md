# Co-author review — consolidated action plan

Sources (read 2026-09-18):

| Source | Form | Volume |
|---|---|---|
| `Authors Reviews Comments.docx` | Unsigned PT-BR letter (content points to the co-author involved at the initial design stage) | 1 general letter |
| `atlantic_birds_JSOR_260916.docx` | J. Santiago Ortega Ramírez | 4 general paragraphs (tracked insertions) + 24 margin comments |
| `Ed_MS_AM.docx` | Ayumi Mizuno | 53 margin comments, no tracked edits |

No tracked deletions anywhere; Ortega's "edits" are inserted commentary paragraphs, not replacement text.

---

## Decisions taken (2026-09-18)

1. **Biology first.** Model ladder demoted; Results lead with conclusions.
2. **Dispersion model becomes the primary variance analysis** (Mizuno E10); lnCVR becomes the complementary species-level summary.
3. Ants, century-scale, and detailed model tables demoted to Supplementary, retained only to substantiate main-text claims.
4. **Analyses to run: E1–E4 and E6–E7.** (E5, E8–E12 not commissioned.)
5. **João C.T. Menezes moved from the author list to Acknowledgements**, per his request.

## Progress

- [x] **F. Authorship** — Menezes removed from YAML author block; new `# Acknowledgements` section added before References.
- [x] **C1** — all 4 "credible interval" → "confidence interval"; "credible declines/increases" reworded.
- [x] **C2** — "observer" → "contributor" throughout (~35 sites). "inter-observer measurement error" retained in 2 places where it refers to the literature concept about actual measurers. Methods §"Contributor turnover and research-group random effects" rewritten to define the term once, state that `Main_researcher` is a proxy, and say explicitly that it absorbs protocol/locality/species-composition differences as well as measurement conventions (Mizuno c65, c66). Discussion model description corrected: contributor + individual identity are random intercepts, lat/long/elev/season/moult are fixed effects (Mizuno c57).
- [x] **B. Overclaiming** — abstract, plain-language summary, aims 3–4, 4 section headings, thermal Results + figure caption, bill-width Results, body-mass Results, variance Results, diet Results, ant paragraph, and 8 Discussion passages (Mizuno c18, c20, c30, c35, c43, c46, c47, c48, c51, c53, c54, c56, c92; Ortega c40, c47, c48, c51, c52, c59, c60, c63).
- [x] Manuscript re-renders cleanly (`quarto render index.qmd`); no inline R expressions or citations lost.

- [x] **D12** — "resident" dropped everywhere it described our sample (abstract, Introduction ×2, Discussion); Methods §4.1 now states that migratory status is not recorded in the compilation, was not used as a filter, and that the sample may include partial or altitudinal migrants.
- [x] **D5** — the glmmTMB-vs-brms validation claim removed entirely; the `propto` rationale now stands on computational grounds alone, with no equivalence claim. No Supplementary section needed.
- [x] **A1/A2 Introduction rebuilt** — now: (1) warming alters morphology via several pathways, predictions exist, but responses are heterogeneous [@radchuk2019adaptive added to `references.bib`]; (2) the evidence base is temperate/Northern-Hemisphere skewed, Jirinec as the tropical benchmark, framing question stated; (3) Atlantic Forest as the contrasting test system, mechanisms deferred to Discussion; (4) variability and diet motivated explicitly (Mizuno c5); (5) compilation hazard in two sentences (Ortega c27); (6) one broad aim with sub-analyses in place of five numbered questions.
- [x] **A3 Results reordered** to Ortega's sequence: warming → wing → body mass → wing–mass decoupling → temperature associations → bill width → variability → diet → robustness.
- [x] **A4 Model-ladder prose removed** — wing, mass and bill-width sections now lead with the conclusion; Ortega's drafted replacement prose (c35, c40, c55, c60, c63) used.
- [x] **A5 Demotions** — `tbl-controlled-wing` and `tbl-mass-allometry` moved to `supplementary.qmd` as Tables S2/S3 (the `load-estimates` chunk was copied across so they still build from the same `.rds` files); lnCVR forest plot → Supplementary Fig. S10 (new Note S6); litter-ant analysis → new Note S7, one sentence left in the main text.
- [x] **Dispersion model is now primary** for variance (your decision 2); lnCVR reframed as the complementary species-level summary in both Results and Methods.
- [x] Discussion summary reordered to match the Results and reworded to the softened claims; opens with the framing question.
- [x] `supplementary.qmd`: Menezes removed from its author list; "observer" terminology aligned.
- [x] Both documents re-render cleanly; no unresolved cross-references, no inline-R leaks.

### E-cluster (commissioned: E1–E4, E6–E7)

- [x] `Analysis/scripts/referee_reruns.R` written — all six items in one script, writing `output/referee_reruns/{referee_reruns.rds,REPORT.md}`. Runs in ~30 min on 50 trees; `REFEREE_N_TREES=2` for a smoke test.
- [!] **Reproducibility gap found.** The script that produced the reported thermal-coupling estimates (β = +0.0212, t = 2.17, p = 0.030 on `passer90`) is **not in the repository** — only its figure (`figures/secular_expansion/manuscript_thermal_coupling_allometry.png`) and its numbers in `output/secular_expansion/REPORT.md` §5.5. `secular_museum_expansion.R` produces the *secular* coupling figures but never fits the `passer90` model. Block E0 of the new script refits it from `passer90_climate.rds` so the number has a reproducible source.
- [!] **E0 does not reproduce the reported estimate.** Closest reconstruction gives β = +0.0175, t = 1.74, **p = 0.082** on n = 7,224 (the figure's report says n = 7,577; 772 records lack `rec_tmean`, so that model's climate linkage differs from `passer90_climate.rds`). Needs resolution before the thermal claim can stand as published.
- [!] **clubSandwich CR2 unavailable** for these models (non-nested random effects). E2 therefore uses a locality-year random intercept, which addresses the same shared-exposure concern within the mixed-model framework.

### E-cluster results (50 trees, `output/referee_reruns/referee_reruns.rds`)

**E4 — isometry: HOLDS, slightly stronger under stricter adjustment.**
Contrast fitted directly as a response (exact SE; answers Mizuno c82 — no delta-method covariance needed).
| spec | β/SD-yr | 95% CI | p | %/decade |
|---|---|---|---|---|
| parsimonious | +0.0224 | 0.0085–0.0363 | 0.0016 | +4.46 |
| **full M3** | **+0.0253** | 0.0118–0.0389 | **0.00025** | **+5.04** (2.35–7.74) |
| ln wing alone, M3 | −0.0093 | −0.0148–−0.0038 | 0.0008 | −1.84 |
| ln mass alone, M3 | −0.0010 | −0.0079–+0.0058 | 0.77 | −0.21 |

**E1–E3 — thermal association: DOES NOT HOLD.** Anomaly (`dT`) on the isometry contrast:
parsimonious z = 1.74 (p = 0.082) → +year z = 1.18 (p = 0.239) → detrended z = 1.58 (p = 0.114)
→ locality-year clustered z = 1.44 (p = 0.150) → **M3-adjusted z = 0.92 (p = 0.359)** → M3 + year z = 0.46 (p = 0.648).
In every combined model `scaled_yr` carries the signal (iso: z = 3.65, p = 0.00027) while `dT` does not.
Same pattern for wing, mass and relative wing. 48–50 of 50 trees converged per model.

**E6 — secular source bridging: NOT ROBUST.** Spanning cohort (n = 24,583):
| spec | mm/SD-yr | t | %/decade |
|---|---|---|---|
| as run (Status intercept) | −0.278 | −1.98 | −0.71 |
| no source adjustment | −0.276 | −1.99 | −0.70 |
| **year × Status interaction** | — | **+2.56** | — |
| **museum only** (n = 3,315) | **−0.016** | **−0.19** | **−0.04** |
| live only (n = 21,268) | −0.536 | −2.00 | −1.36 |
The decline is carried entirely by the live records; the museum series alone is flat. A Status *intercept* does not fix this because the *slopes* differ.

**E7 — age composition: NOT ROBUST.** P(Adult) ~ Year: β = −0.0122, z = −9.12, **p = 7.3×10⁻²⁰** — the proportion of adults falls steadily over the series. Adding Age as a covariate: −0.254 (t = −1.81). **Adults only (n = 18,957): −0.219, t = −1.05, CI spans zero.**

### Reproducibility finding: the thermal figure has no generating code

- `manuscript_thermal_coupling_allometry.png` was committed in `bad2941` alongside `secular_museum_expansion.R`, which **does not produce it** (neither the current nor the committed version). `git log --all -S 'manuscript_thermal_coupling'` returns only that commit, i.e. the binary. The committed `secular_expansion_results.rds$climate_coupling$direct_temp_model` is a different model (`wing ~ af_temp + Status`). Figure mtime Sep 10 15:52 vs all sibling figures Sep 11 15:38.
- The figure's own data ARE reproducible: it plots 24 annual means weighted by birds/year, and the reconstructed series matches its axes (anomaly −0.34→+0.26 vs plotted −0.35→+0.28; iso −8.6→+9.7% vs plotted −10→+10).
- The reported numbers (β = +0.0212, t = 2.17, p = 0.030, n = 7,577) do not follow from it. Local anomaly on annual means: t = 0.19–1.03. GISS regional anomaly: t = 2.09–2.48 but β = +7.9 to +9.6 %/°C, ~4× the reported value. Individual-level mixed model, local anomaly: β = +0.0200, z = 1.98, p = 0.048 — closest, still not it. n = 7,577 is impossible (only 7,224 shared records have temperature); p = 0.030 is the *normal* p for 2.17 (a t with 22 df gives 0.041).
- `climate_extraction.R` carries its own caveat: *"capture-year temperature cannot cause a fixed adult wing length; the locality trends are CONTEXT for the study region and period, not a per-individual exposure."*

### `rec_tmean` provenance (answered)

WorldClim **historical monthly weather** (CRU-TS 4.09 downscaled with WorldClim 2.1, monthly tmin/tmax 1950–2024, 2.5 arc-min) — the time-resolved product, explicitly not the v2.1 normals. `rec_tmean` = annual mean of (tmax+tmin)/2 in the **capture year** at the record's coordinates. The 772 NAs are a coastal land-mask artefact (98 of 455 localities, 91 around Guarapari/ES) with no nearest-cell fill applied by design.

### Structural changes made after the E-cluster (2026-09-18)

All numbers below are pulled into the text by inline R from `Analysis/output/referee_reruns/referee_reruns.rds` via accessors in the `load-estimates` chunk (`rr_therm`, `rr_iso`, `rr_sec`). Nothing is typed by hand.

- **Thermal claim removed and replaced by the null result.** The Results section is now "Annual temperature anomalies do not explain the allometric shift": one short section reporting that the anomaly effect is indistinguishable from zero under M3 adjustment, and under three further specifications, while calendar year carries the signal. `fig-thermal-allometry` is dropped from the main text (main-text figures 7 → 6); the figure file itself has no generating code and should not be republished. Abstract, Discussion point 4 and the Methods section updated to match. Methods now describe what was actually run and state the exposure-window objection.
- **New Supplementary Note S8** tabulates all four thermal specifications (E1/E2/E3) as Supplementary Table S4, built from the saved object, with the convergence count per model and a note that CR2 is unavailable for non-nested random effects.
- **Secular section rewritten.** Keeps only the claim that survives — no morphological trend before 1980, so the contemporary pattern is not the tail of a long trajectory — and reports the failure explicitly: decline carried by live captures, absent in museum specimens alone, year x source interaction non-zero, and gone when restricted to adults. Discussion paragraph rewritten to match; the "emerged with post-1980 warming" claim is withdrawn.
- **Isometry strengthened.** The Results now add the M3 confirmation (+5.0%/decade, 95% CI 2.3-7.7, p < 0.001) and note that fitting the contrast directly yields its SE exactly, which answers Mizuno c82.
- Both documents render clean: no unresolved cross-references, no inline-R leaks, no NA leaks.

### Still open

- **D1–D4, D6–D11, D13–D15** — the remaining Methods clarifications (11 items), all text-only.
- **E1–E4, E6–E7** re-runs.
- The Discussion body sections still follow the old order internally; worth a pass once the E-cluster numbers are final.

---

## A. Structural / strategic — all three reviewers converge

- **A1. Rebuild the Introduction around one question.** Current intro develops Bergmann, Allen, nutritional limitation, the Amazon comparison, fragmentation, topographic refugia, and trait-compilation bias. Ortega (c15, c16, c24, c25, c27) wants: morphological responses to warming are *heterogeneous*, and the evidence base is geographically biased (Radchuk et al. 2019) → the Atlantic Forest is a test of whether the Amazonian pattern (Jirinec et al. 2021) generalises. Mechanisms move to Discussion. The methodological-bias paragraph shrinks to 1–2 sentences.
- **A2. Replace the five numbered aims with one broad aim** plus sub-analyses underneath (Ortega c28). Mizuno c5 separately notes aims 5 (variability, diet) are the least motivated in the Introduction.
- **A3. Reorder Results** per Ortega's explicit sequence: (1) regional warming; (2) weak/heterogeneous wing change; (3) body-mass stability; (4) wing–mass coordination vs Amazonia; (5) temperature associations; (6) bill width + variability; (7) diet; (8) historical + broad-sampling as robustness. Discussion follows the same hierarchy.
- **A4. Stop narrating the model ladder.** Ortega c34, c35, c52, c54, c60 and Mizuno c7, c8 all say the same thing: report the biological conclusion, keep M0→M3 in a table (supplementary). Ortega supplied ready-to-use replacement prose for wing (c35), variance (c55), diet (c60), ants (c63), and a summary sentence for wing heterogeneity (c40).
- **A5. Demote to Supplementary:** detailed model tables (Mizuno c13), Figure 7 (Ortega c57), the litter-ant prey analysis (Ortega c62, Mizuno c37 — and note Mizuno c38: ant *richness* is not prey availability). Century-scale/museum section must not compete with the 24-year story — keep it narrow and framed as robustness (Ortega c65, Mizuno c40).
- **A6. PT-BR letter:** asks the framing question directly — is the paper about the biology or the modelling? Recommends biology-first, most models to Supplementary, a forest plot showing M0→M3→Mundlak progression, and a **separate methods paper** as a worked example of which factors produce spurious trends. **Also asks to be removed from the author list.**

## B. Overclaiming / causal language — Ortega + Mizuno agree

- "directly elicit", "actively modulates", "directly drive", "directly modulate", "thermal stress" → association wording (Ortega c47, Mizuno c4, c25, c26). Affects abstract (l.71), aim 4 (l.262), Results §3.x (l.388, l.392), Discussion point 5 (l.439).
- **Headings:** "Direct thermal sensitivity and temperature-driven allometric modulation" → "Associations between annual temperature anomalies and morphology" (Ortega c48). "Bill width invariance and testing Allen's rule" → "Bill width shows little temporal change" (Ortega c51). "Trophic guild buffering…" → "Dietary specialization and morphological change" (Ortega c59). Discussion "preservation of allometric scaling" → "departure from geometric isometry" / "mass–wing decoupling" (Mizuno c51).
- Frame bill width as *little evidence for Allen-type change*, not a test/rejection of Allen's rule — we measure temporal change within populations, not the ecogeographical pattern (Ortega c52).
- Soften: "trajectories are negative" when the M3 CI overlaps zero (Mizuno c18); body mass "stable" → "little or no clear evidence of directional change" (Mizuno c20); "does not drive morphological sensitivity" (Mizuno c35); "artifact of later contributors" (Mizuno c30); "confirms that phenotypic stability reflects a genuine biological characteristic" (Mizuno c48, l.455); "provides compelling evidence … functionally linked" (Mizuno c53, l.473); "confirms our primary conclusions do not depend on filtering" — body mass actually changes in the expanded sample (Mizuno c54, l.475); "intermediate between temperate and equatorial extremes" not established by a modest decline alone (Mizuno c43); the ~60% attenuation cannot be attributed specifically to observer + spatial turnover (Mizuno c56).
- **Amazon comparison:** separate the observed cross-study contrast from candidate explanations; the two studies differ in design, coverage, composition and measurement structure, and our data cannot discriminate among fragmentation / refugia / plasticity (Mizuno c46). Diurnal mass gain (+0.40%/h) was interpreted in Results as a feeding rhythm and cannot then be used as evidence of thermal plasticity (Mizuno c47).

## C. Terminology — mechanical, do first

- **C1.** "credible interval" → "confidence interval" (4 occurrences: l.311, l.317, l.414, l.465 incl. "credible declines/increases"). Main inference is glmmTMB REML + Rubin's rules; brms is validation only (Mizuno c15).
- **C2.** `Main_researcher` is the PI / research-group coordinator, **not** the measurer (`Measurer` is 100% NA). Pick one term — "contributor" — and drop "observer", "within-observer", "observer control", "observer turnover" (~35 occurrences in `index.qmd`) (Mizuno c65, c66, c57). Also fix the M3 description: lat/long/elev/season/moult are fixed effects; `Main_researcher` and `ID_individual` are random intercepts (Mizuno c57).

## D. Methods clarifications — text only, no re-runs

| # | Comment | What to add |
|---|---|---|
| D1 | Mizuno c71 | How missing covariates enter main M3 vs `M3_cc` (n = 8,282); exact difference |
| D2 | Mizuno c72 | `M6`: state the exact formula modification (Sex dropped?) |
| D3 | Mizuno c73 | `M3_noanom`: how contributor × site × year block deviation was computed; why 6 mm |
| D4 | Mizuno c76 | Confirm phylogenetic covariance is on species intercepts only, slopes independent |
| D5 | Mizuno c77 | Document the glmmTMB-vs-brms validation (models, parameters, priors, agreement metric) in Supplementary — the "0.001 units across 50 trees" claim currently has no backing |
| D6 | Mizuno c79 | Capture-hour model: full M3 + hour, or year + hour only? |
| D7 | Mizuno c82 | How the mass–wing slope covariance for the delta-method SE was obtained (joint vs separate fits) |
| D8 | Mizuno c88 | Breakpoint: estimated per response, or 1980 imposed from the climate series? |
| D9 | Mizuno c90 | lnCVR within-contributor subset: report species, contributors, k, overlap with pooled |
| D10 | Mizuno c91 | Why 1995–2006 / 2013–2018, and why 2007–2012 is unused |
| D11 | Mizuno c61 | Justify the 20-km biome buffer |
| D12 | Mizuno c62 | **No residency/migratory filter is described in §4.1** — either add it or stop calling them "resident Atlantic Forest passerines" |
| D13 | Mizuno c28 | Fig. 6: does the line come from the main mixed model? |
| D14 | Mizuno c3 | Why bill width is the focal Allen's-rule appendage |
| D15 | Mizuno c92 | Within-contributor comparison removes between-contributor pooling — it does not isolate "genuine biological variance" |

## E. New analyses / re-runs — the costly set (needs your call)

| # | Comment | Request |
|---|---|---|
| E1 | Mizuno c84 | Temperature models: add calendar year simultaneously, or use detrended anomalies — anomaly still covaries with year |
| E2 | Mizuno c85 | Cluster SEs by locality-year — anomaly is shared by all birds at a site-year |
| E3 | Mizuno c86 | Re-fit temperature models at the M3 adjustment level (or justify the parsimonious set) |
| E4 | Mizuno c81 | Same for the bivariate isometry model — why omit longitude, elevation, season, moult? |
| E5 | Mizuno c31 | Test bill width directly against local annual temperature anomaly |
| E6 | Mizuno c41 | Secular series: source (museum vs live) calibration or museum-only sensitivity — §4.1 excludes skins precisely because of shrinkage |
| E7 | Mizuno c42 | Secular series uses all age classes; control for age |
| E8 | Mizuno c68 | Within-locality sensitivity: locality-specific year slopes, or synthesise within-locality slopes |
| E9 | Mizuno c94 | Add phylogenetic non-independence to the lnCVR meta-analysis |
| E10 | Mizuno c93 | Promote the continuous location-scale dispersion model to the primary variance analysis; lnCVR becomes complementary |
| E11 | Mizuno c96 | Ant model: replace the `log(effort)` offset with an estimated (possibly nonlinear) effort covariate — only if the ant analysis survives A5 |
| E12 | Mizuno c22 | Isometry: justify wing as *the* linear dimension, or use a structural size measure (tarsus is in the database) — rejecting mass = 3 × wing shows mass–wing decoupling, not necessarily whole-body stoutening |

## F. Administrative

- Authorship: the PT-BR letter asks to be taken off the author list. Your call.
- Decide whether the methods material becomes a companion paper (both the PT-BR letter and Ortega's framing point that way).

---

## Decision gate — needed before revision starts

1. **Framing:** biology-first with the model ladder in Supplementary (all three reviewers' recommendation) — confirm?
2. **Scope of demotion:** ants → Supplementary; century-scale → narrow robustness section; variance → which model is primary (E10)?
3. **Which of E1–E12 to run.** E1–E4 are the ones a journal referee is most likely to ask for anyway. E6/E7 are the ones that could change a reported result.
4. **Authorship.**
