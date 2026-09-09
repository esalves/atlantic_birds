# Revision notes — Phase 4a (ATLANTIC ANTS litter-ant richness index), 2026-09-09

Companion to `REVISION_PLAN.md` §3 Phase 4a and §7. Script:
`atlantic_ants_index.R`. Everything below is from the FULL run (glmmTMB, no
Stan; ~25 s on a laptop, R 4.6.0, data.table 1.18.4, glmmTMB 1.1.14, ggplot2
4.0.3, patchwork 1.3.2). Nothing here is a smoke test. All numbers are stored
in `Analysis/output/ants_results.rds` (`$year_effects`, `$programme_slopes`,
`$loo_within`, `$annual_index`, `$contributors`, `$events_campaign`, ...) and
tabulated in `Analysis/output/ants_results.md`.

---

## 1. What this index is, and what it is not

**It is** the number of ant species (distinct `Genus + Species`; contributor
morphospecies codes count as species within a contributor file) recorded in one
standardised sampling campaign of the ATLANTIC ANTS compilation (Silva et al.
2022, *Ecology* 103:e3580, CC-BY), adjusted for the number of Winkler samples or
pitfall traps in that campaign, with random intercepts for contributor file,
~11-km locality and site, and a coarse habitat class. It is a **litter/ground-ant
community richness index per standardised sample** for the Brazilian Atlantic
Forest, 1994–2018.

**It is not** prey biomass, prey availability for insectivorous birds, or
arthropod abundance. Ants are one taxon; richness per sample is not abundance;
only 9.5 % of the standardised records (67 of 855 campaigns, 45 of them from one
programme) carry counts. Any sentence in the manuscript that uses this index
must say "litter-ant richness per standardised sample", never "prey".

**It is also not a monitoring series.** Like PREDICTS, ATLANTIC ANTS is a
compilation of independent studies. The difference — and the only reason it is
usable here — is that it is region-matched and carries contributor, site,
method, effort and habitat fields that allow the controls the referee asked for.
The results below show that those controls change the estimate by a factor of
two and that the sign depends on how effort is treated. Read §4 before quoting
anything.

## 2. Data and filters (counts from the run)

| Step | Records |
|---|---|
| Raw file | 178,976 (57 fields, 113 contributor files) |
| Country == BRAZIL | 165,657 |
| `ma.lmt == 1` (authors' Atlantic Forest limit) | 148,837 |
| inside bbox lat −31…−5, lon −57…−34 (guard; 1,498 dropped) | 147,339 |
| `Exclude != 1` | 145,983 |
| `Start.year` 1990–2019 | 93,461 |
| Method == Winkler with `Winkler.Number` > 0 | 44,391 of 46,264 Winkler records (96 %) |
| Method == Pitfall with `Pitfall.Number` > 0 | 17,934 of 23,488 pitfall records (76 %) |
| campaigns spanning ≥ 2 calendar years dropped | 305 |
| **standardised records** | **62,020** (9.5 % with abundance) |

Mixed-method records ("Pitfall and Winkler", baited/vegetation pitfalls, etc.)
are excluded. Reading: `data.table::fread(sep = "\t", quote = "", encoding =
"Latin-1")`; `Pitfall.Number` / `Winkler.Number` / `Total.Ant.Abundance` are
character in the file (`"5to15"` token, thousands separators) and are coerced.
`ma.lmt` was used instead of the bbox alone because the data authors delimit
the biome themselves; the bbox removed 1,498 further records (coastal islands,
far-west Paraná).

Habitat: the free-text `Habitat.Type` (200+ levels) was collapsed by keyword
into forest / open / mixed / other / unknown (campaign events: 736 / 78 / 10 /
21 / 10). `Disturbance` was 49 % NA and was not used.

## 3. The sampling-event problem (the non-obvious decision of this phase)

The plan assumed events could be defined as contributor × locality × year ×
method with `Winkler.Number` / `Pitfall.Number` as effort. Profiling showed
that the effort and site fields do **not** mean the same thing across
contributor files:

1. **Effort field is a trap identifier in the Servio Ribeiro PERD/PEIT pitfall
   files.** Within one plot and year, `Pitfall.Number` takes 19–26 distinct
   values (1…30) and the same species is repeated under different values
   (e.g. *Pheidole* sp04 with 30, 25, 24, 21). Those records are species ×
   trap. 2,010 records / 42 plot-dates are affected (PERD 6 % of plot-dates,
   PEIT 2 %, LAMAT Winkler 17 %, LITERATURE_TEAM 7 %). Taking the maximum per
   plot-date recovers the number of traps (30 in 2001–02, ~5 in 2014–16).
2. **`ID.codLoc` is a per-record code, not a sampling area, in several large
   files.** At the plot-date level GRACIELI_ARAUJO has 724 "plot-dates" of
   which 99.7 % contain exactly one species with 20 traps each; PEMS 154
   (98.7 %, 600 traps); ERICK_VILLARREAL 50 (100 %, 320 traps); UFV_LABECOL
   256 (76 %). Pooled to the campaign these become 46 species per 20-trap
   campaign etc., which is plausible.
3. **Effort is a campaign total in others.** PERD Winkler 2005–06: 80 plot
   codes at two coordinate pairs sharing one `Winkler.Number = 360`
   ("all leaf-litter inside the 1-m² plot", i.e. 360 plots ≈ 80 codes).

The metadata (aants_metadata.pdf, p. 77–78) define both effort fields as the
number of samples for the sampled site and `ID.codLoc` as the "name or code of
the sampling area provided by the reference paper". The **campaign** —
contributor file × Method × Start date × coordinates rounded to 0.01° ×
`Habitat.Type`, all `ID.codLoc` pooled — is therefore the primary unit. Effort
is built in two steps: per plot-date `max(effort field)`; per campaign the
shared value if all plot-dates report the same number (a campaign total),
otherwise the sum of plot-date efforts. This yields **855 campaign events**
(393 Winkler, 462 pitfall; 61 contributor files, 293 ~11-km localities, 570
sites, 55 sites sampled in more than one year; every year 1994–2018 has at
least one event, 1–17 contributors per year; the plan's "~700 events"
corresponds to the 651 contributor × method × coordinate × year cells).
The plot-date level (3,955 events, **42 % with richness = 1**) is kept only as
a sensitivity and is flagged as an unreliable unit in the figure and table.

Consequence for inference: a reporting convention that is constant within a
contributor is absorbed by the contributor intercept and cannot bias the
within-contributor slope. A convention that **changes within a contributor**
can, and this is exactly the situation in the largest programme (PERD: trap-ID
pitfalls with 30 traps per plot in 2001–02 → 360-sample Winkler campaigns in
2005–06 → 5-trap pitfall plots in 2014–17, at partly different coordinates).
Hence the constant-protocol check in §4.

Not verifiable without the data owners: which convention each file used. This
is the reason Phase 4d (collaboration request to S. P. Ribeiro for the raw PERD
design) matters more than any further modelling here.

## 4. Results (campaign level, N = 855; year per SD of the wing model = 5.020 yr)

Family: negative binomial 2 (Poisson dispersion 1.66; ΔAIC 608 in favour of
NB; θ = 5.8). Random-intercept SDs (pooled model): contributor 0.98, locality
0.30, site 0.00 (site variance collapsed to zero — the contributor and
locality terms absorb it). log(effort) coefficient 0.40 (SE 0.05) for Winkler,
0.32 for pitfall: richness scales as effort^0.3–0.4, **not** proportionally.

| Model | Term | N | Estimate per SD yr [95 % CI] | % per decade [95 % CI] |
|---|---|---|---|---|
| naive (no random effects) | year | 855 | −0.269 [−0.330, −0.208] | −41.5 [−48.2, −33.9] |
| **pooled: contributor + locality + site RE, effort covariate** | year | 855 | **−0.113 [−0.187, −0.039]**, z = −2.99 | **−20.1 [−31.1, −7.4]** |
| pooled: effort as offset | year | 855 | +0.041 [−0.051, 0.133] | +8.5 [−9.7, +30.3] |
| Mundlak by contributor | within | 855 | −0.126 [−0.203, −0.050] | −22.2 [−33.2, −9.4] |
|  | between | 855 | +0.060 [−0.203, 0.322] | +12.7 [−33.2, +90.1] |
| Mundlak by site (55 resampled sites, 206 events) | within-site | 855 | −0.410 [−0.631, −0.188] | −55.8 [−71.6, −31.2] |
| monitoring programmes ≥ 6 yr (5 files) | pooled year | 251 | −0.187 [−0.289, −0.084] | −31.1 [−43.8, −15.4] |
|  | within | 251 | −0.166 [−0.268, −0.064] | −28.1 [−41.3, −12.0] |
| without literature-compilation files | within | 704 | −0.166 [−0.253, −0.079] | −28.2 [−39.6, −14.6] |
| constant protocol (modal effort per contributor × method series) | within | 655 | −0.290 [−0.439, −0.142] | −43.9 [−58.3, −24.7] |
| abundance (secondary; 67 events, 45 PERD) | pooled, offset | 67 | −0.142 [−0.551, 0.266] | −24.7 [−66.6, +70.0] |
|  | within | 67 | −0.282 [−0.724, 0.161] | −42.9 [−76.4, +37.9] |
| plot-date level (unreliable unit) | pooled | 3,955 | −0.263 [−0.322, −0.203] | −40.7 [−47.4, −33.3] |

Mundlak within − between = −0.186 (SE 0.139, p = 0.18): no evidence that
contributor composition drives the pooled slope, i.e. the decline is
within-contributor.

**Per-programme slopes (independent GLMs, ≥ 6 sampling years, literature files
excluded):** LAMAT/UMC (80 events, 2001–2017) +0.097 (SE 0.172); PERD (53,
2001–2017) −0.385 (0.111); ELMO & DELABIE (65, 1996–2002) −0.026 (0.062);
CEPLAC/UESC (27, 2003–2014) −0.475 (0.069); BIOTA FORMIGAS (26, 1997–2003)
+0.284 (0.134). Two positive, one flat, two strongly negative. The random-slope
model (diag covariance) converged with slope SD = 0, so its conditional slopes
are all equal to the fixed slope (−0.187) and carry no information; the GLM
slopes are the ones to quote. LITERATURE_TEAM_PUBLISHED_DATA (10 "years",
+1.10) and the two EQUIPE_COLETA_PUBLISHED_DATA files are compilations of
different source papers and were excluded from the programme set (they remain
in the pooled models as contributors).

**Leave-one-contributor-out (Mundlak within slope, contributors ≥ 20 events):**
range −0.151 (drop LITERATURE_2020) to **−0.055 (SE 0.047; −10.5 %/decade,
CI ≈ [−0.15, +0.04]) when PERD is dropped**. PERD supplies 53 of 855 events but
about half of the within-contributor signal.

**Annual adjusted index** (`ants_annual_index.rds`, year-factor model,
reference 2015 = year with most contributors, 17): no monotone pattern —
1997–2012 sit at +0.2 to +0.5 log units above 2015 with wide intervals, 2013 is
a dip (−0.53 [−0.79, −0.28], 6 contributors), 2014–2017 ≈ 0, 2018 (4 events)
+1.28. Years 1994–1996 have 1–4 events from a single contributor.

## 5. Honest reading

1. With contributor, locality, method, effort and habitat controls, litter-ant
   richness per standardised campaign declined by about 20 % per decade
   (−20 %, 95 % CI −31 to −7); the within-contributor estimate is the same
   (−22 %). The naive compilation trend (−42 %) is twice as steep — the
   contributor controls matter, in the same direction as for the birds.
2. The estimate is **not robust to the effort assumption**: treating effort as
   an offset (richness ∝ samples) flips the sign (+8.5 %, CI includes 0). The
   offset is the wrong model (fitted exponent 0.3–0.4, and effort per campaign
   fell over time, r(year, log effort) = −0.23), but a referee will ask.
3. It is **half-driven by one programme (PERD)** whose protocol changed within
   the series in ways the data cannot fully encode (§3). Dropping PERD leaves
   −0.055 (SE 0.047). The constant-protocol subset is steeper (−0.29), but it
   retains only 23 multi-year series, most of 2–3 years.
4. The five long programmes disagree in sign. The within-site estimate
   (−56 %/decade) rests on 55 resampled sites dominated by PARNASC 2014–15,
   MAULYSSEA 2013–16 and PERD.
5. Abundance (the quantity closest to what the manuscript once called "prey")
   is available for 67 campaigns, 45 from PERD, and is null with an interval
   from −67 % to +70 % per decade.

Recommended manuscript wording (for P7): "The only region-matched arthropod
compilation with temporal coverage (ATLANTIC ANTS) shows a decline in
litter-ant richness per standardised sample of roughly 20 % per decade after
contributor, site, method and effort controls (95 % CI 7–31 %), but the
estimate depends on how sampling effort is modelled, is driven to a large
extent by a single long-term programme whose protocol changed over time, and
the five multi-year programmes disagree in sign. We therefore treat it as
context, not as evidence of declining prey." The arthropod → wing causal path
stays out of the title and abstract, as the plan requires.

## 6. Outputs and downstream use

- `Analysis/data/derived/ants_events.rds` — 855 campaign events (`event`,
  `fl.nm`, `contributor`, `literature_file`, `Method`, `year`, `month`,
  `lat/lon`, `locality` (0.1° cell), `site`, `n_plot_codes`, `habitat3`,
  `effort`, `effort_is_id`, `richness`, `abundance`, `scaled_year` on the bird
  scale, `year_within`/`year_between`, `long_programme`, `modal_effort`).
- `Analysis/data/derived/ants_events_plotdate.rds` — 3,955 plot-date events
  (sensitivity only; 42 % richness = 1).
- `Analysis/data/derived/ants_records_std.rds` — the 62,020 standardised
  records with both event keys, `effort_raw`, `effort`, `effort_is_id`.
- `Analysis/output/ants_results.rds` — everything quoted here, plus the fitted
  glmmTMB objects (`$models`) and per-model convergence flags (all TRUE).
- `Analysis/output/ants_results.md`, `Analysis/output/ants_annual_index.rds`.
- `Analysis/figures/ants_events_by_year.png` (events per year by contributor
  and method), `Analysis/figures/ants_year_effect.png` (forest plot of all
  year effects in % per decade; within-programme fits for the 5 programmes).

For Phase 4f (drmSEM `INCLUDE_ARTHRO` node): a locality × year arthropod
covariate can be built from `ants_events.rds` (residual richness after the
pooled model, averaged by `locality` × `year`), but only 293 localities × 25
years are covered, mostly one campaign each; the bird localities and the ant
localities overlap little. My recommendation to the P4ef agent is to drop the
node rather than build it on this (documented in `REVISION_PLAN.md` as
acceptable). `ants_annual_index.rds` gives a biome-level annual index if one is
wanted for a supplementary figure — with 1–17 contributors per year it should
not be given a causal role.

## 7. Things I did not do / open

- Did not attempt to match contributor files to the source papers to recover
  the true per-plot design; that is the Phase 4d collaboration request.
- Did not model species accumulation explicitly (e.g. rarefaction); the
  log-linear effort covariate with method interaction is the pragmatic
  substitute, and the offset variant shows the sensitivity.
- The `Manuscript/index.qmd` text and `update_descriptive_stats.R` are P7's
  and P0's files; the numbers to read are in `ants_results.rds$year_effects`
  (row `model == "pooled: contributor + locality + site RE, effort covariate"`
  for the headline, `"Mundlak by contributor"` for within/between,
  `$loo_within` for the PERD sensitivity).
- The raw zip is git-ignored (`*.zip`); the README in
  `Analysis/data/raw/atlantic_ants/` gives the re-download command.
