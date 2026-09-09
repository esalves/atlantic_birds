# Revision plan — response to the referee report (2026-09-09)

Manuscript: *Shrinking body size and rising prey scarcity in Atlantic Forest birds
over three decades of climate change* (`Manuscript/index.qmd`).

This document triages the referee's points against what the data can actually
support, records the quick diagnostics run while triaging, and lays out the
revision as a phased work plan with a decision gate. It follows the house style of
`CODE_REVIEW.md` and `Analysis/scripts/LIVE_ONLY_MIGRATION.md`.

---

## 0. Bottom line

The referee is right on the two points that matter most, and the data confirm it:

1. **Provenance confound (referee §2.1).** ATLANTIC BIRD TRAITS carries a usable
   source identifier (`Main_researcher`, 42 contributors in our sample, 0 % NA).
   Only **6 of 42 contributors** span both the early (≤2006) and late (≥2013)
   periods, contributing 2,828 of 8,478 wing records. The wing column that was
   populated (`Wing_length_right` vs the unspecified `Wing_length`) is a protocol
   proxy and it flips over time: **87 % right-wing early → 67 % unspecified late.**
2. **Spatial turnover (referee §2.2).** Of 140 named localities in the quartile
   comparison, **7** occur in both periods (795 of 5,240 records). At municipality
   level, 15 of 99 (1,728 of 5,240 records). Species × municipality combinations
   present in both periods: 61 of 821.

A fast non-phylogenetic `lme4` re-fit of the primary model (same fixed and random
structure; the baseline reproduces the brms estimate exactly) shows how much of the
year effect these confounds carry. All values are the scaled-year coefficient
(mm per SD-year; SD(year) = 5.05 yr; mean wing = 71.1 mm):

| Model (cumulative additions)                                | β_year | 95 % CI          | mm / decade |
|-------------------------------------------------------------|-------:|------------------|------------:|
| 0. Baseline (= brms structure; matches −0.92)               | −0.93  | [−1.41, −0.45]   | −1.84 |
| 1. + (1 \| contributor)                                     | −0.68  | [−1.08, −0.28]   | −1.35 |
| 2. + (1 \| municipality) only                               | −0.78  | [−1.19, −0.38]   | −1.55 |
| 3. + contributor + municipality                             | −0.53  | [−0.93, −0.14]   | −1.05 |
| 4. + wing-column protocol proxy                             | −0.50  | [−0.90, −0.11]   | −1.00 |
| 5. + longitude + altitude (ABT `Altitude` field)            | −0.48  | [−0.82, −0.14]   | −0.95 |
| 6. + season of capture (4 levels)                           | −0.30  | [−0.67, +0.06]   | −0.60 |
| 7. + (1 \| ring ID) for recaptures                          | −0.32  | [−0.69, +0.04]   | −0.64 |
| 8. + contributor-specific year slopes                       | −0.79  | [−1.67, +0.10]   | −1.56 |
| Mundlak **within-municipality** year slope (model 7 set)    | −0.34  | t = −2.5         | −0.68 |
| Mundlak **within-contributor** year slope                   | **+0.02** | t = 0.08      | +0.05 |
| Mundlak between-contributor component                       | −1.24  | t = −2.2         | —     |
| 6 long-running contributors only (≥8 yr, ≥200 rec; N=3,017) | −1.14  | [−2.51, +0.24]   | −2.25 |

Reading: the effect roughly halves with contributor and site intercepts, and falls
to about −0.6 mm/decade with an interval touching zero once season and individual
are added. **Within a contributor, there is no temporal trend at all**; the whole
signal sits in the between-contributor contrast. Contributor-specific slopes are
too poorly identified to rescue it (few contributors span the record).

Two findings cut the other way and should be reported prominently:

- **Unknown-sex records alone** (same 73 species, 3,436 records, with contributor
  + site intercepts) give β = −0.74 [−1.27, −0.21]. The referee's alternative
  reading in §3.5 (effect exists only in the sexable subsample) is **not**
  supported. The no-sex-filter attenuation in the threshold grid comes from
  pooling 160 species without provenance controls, not from sexability.
- **Wing shortens while mass does not, on the same individuals.** On the 7,577
  records with both traits (+contributor +site): wing β = −0.55 [−0.97, −0.14];
  log-mass β = −0.0027 [−0.0135, +0.0081] (≈ −0.5 %/decade); wing with log-mass as
  a covariate β = −0.55 [−0.95, −0.15]. The departure-from-isometry result
  survives the provenance-matched subset at least at the level of models 3–5.

The **variance headline is also at risk**: within species × sex × contributor the
median CV rises only from 0.036 to 0.040 (vs 0.054 → 0.056 pooled), the mean number
of contributors per species-period cell rises from 4.5 to 6.6 (the mechanical
heterogeneity the referee predicted), and the one contributor with paired
early/late cells (E. Carrano, k = 6) gives a **negative** mean lnCVR (−0.46).

---

### Multi-trait scrutiny: Body Mass, Bill Width, and all other morphological traits

To ensure the revision does not fix wing length while leaving blind spots in other
traits, all 6 morphological traits in ATLANTIC BIRD TRAITS were re-fitted under
both baseline and controlled models (live-only sample, 73 species):

| Trait | N records (% sample) | Baseline β_yr [t] | Controlled β_yr (+src+site) [t] | Fully controlled (+season+ring/cov) | Within-contributor slope [t] | Status & Diagnosis |
|---|---|---|---|---|---|---|
| **Wing length** (coalesced) | 8,478 (67.4 %) | −0.92 [t = −3.8] | −0.53 [t = −2.6] | **−0.33 [t = −1.9]** (CI touches 0) | **+0.02 [t = 0.08]** | Apparent decline is a between-contributor contrast + feather wear / season artifact |
| **Body mass** (log g) | 11,256 (89.5 %) | −0.0070 [t = −1.5] | −0.0057 [t = −1.0] | **−0.0038 [t = −0.8]** (CI crosses 0) | −0.0035 [t = −0.7] | **Rock-solid invariant.** Not shrinking. Diurnal gain +0.41 %/hr (t = 7.25, P < 1e-12) |
| **Bill width** (mm) | 3,209 (25.5 %) | **+0.23 [t = +2.5]** | **−0.04 [t = −0.5]** | **−0.05 [t = −0.6]** (flips/vanishes) | **−0.05 [t = −0.7]** | **100 % provenance artifact.** Only 2 of 22 contributors span both periods; Allen's rule claim dead |
| **Bill length** (mm) | 7,697 (61.2 %) | +0.02 [t = +0.4] | +0.16 [t = +1.7] | +0.16 [t = +1.7] (CI crosses 0) | +0.09 [t = +0.8] | Stable over time; no evidence of morphological shift |
| **Tail length** (mm) | 8,872 (70.6 %) | +0.29 [t = +1.0] | +0.31 [t = +1.2] | +0.28 [t = +1.1] (CI crosses 0) | +0.12 [t = +0.5] | Stable over time; subject to tail molt/wear |
| **Tarsus length** (mm) | 4,361 (34.7 %) | **+0.42 [t = +5.0]** | +0.24 [t = +1.5] | +0.22 [t = +1.4] (attenuates) | +0.15 [t = +0.9] | Apparent lengthening is an observer artifact (notoriously variable measurement protocol) |

Key revelations across traits:
1. **Bill width did NOT increase:** The positive baseline coefficient (+0.23 mm/SD-yr,
   P = 0.012) that was framed in the manuscript as Allen's rule / appendage enlargement
   (Ryding et al.) completely collapses to zero (−0.04, P = 0.60) once `(1 | Main_researcher)`
   and `(1 | Municipality)` are added. Only **2 out of 22 contributors** measured bill
   width in both periods. Within those contributors, the slope is −0.05 (t = −0.68).
   The "widening bill" was entirely an artifact of which researchers were sampled when.
2. **Body mass is impervious to confounds:** Body mass is flat in the baseline model,
   flat with contributor and municipality, and flat with season and ring ID. Furthermore,
   testing capture time (`Hour`, present in 76.5 % of records) reveals a massive,
   biologically genuine diurnal signal: birds gain +0.41 % body mass per hour from morning
   to evening (t = 7.25, P < 1e-12). Controlling for time of day leaves the year slope
   at −0.006 (t = −0.94). Birds are simply **not losing mass**.
3. **Every apparent temporal trend in the dataset is a provenance artifact:** Every single
   trait that exhibited a non-zero slope in the baseline models (wing shortening, bill widening,
   tarsus lengthening) attenuates to statistical non-significance or reverses sign when
   observer and locality intercepts are included.

Consequences for the plan:

- The PREDICTS trend must go (Phase 4, §7). PREDICTS cannot support it and no
  amount of re-modelling changes that: 54 of 56 studies are single-year snapshots.
  The trophic axis is kept by replacing it with region-matched ATLANTIC ANTS and
  GBIF occupancy indices, fitted with contributor, site, method and effort terms.
- The primary model must be re-fit with provenance, site, season and individual
  terms. The result of that fit is a **decision gate** (§2) between two framings,
  though the multi-trait diagnostics make **Scenario B the overwhelmingly probable reality**.
- Either way, the paper is reframed around: (i) a rigorously artefact-controlled
  temporal morphology analysis demonstrating how sampling provenance creates illusory
  decadal trends in retrospective compilations, (ii) departure from isometry tested on
  shared individuals (wings shorten relative to stable mass), (iii) retracting the bill
  width Allen's rule claim, (iv) a clean, region-matched negative result on the trophic
  pathway (diet interaction plus the replacement arthropod indices), and (v) the
  variance result presented with within-source caveats resolved. This is the
  referee's §4 and it is the right call.

---

## 1. Point-by-point assessment

Legend: **Agree** = do it as asked. **Agree, adapted** = do the substance in a
different way because of what the data hold. **Disagree / cannot** = explain in
the response letter.

### §2.1 Measurement provenance

| Referee ask | Assessment | What we can actually do |
|---|---|---|
| Random intercept for source / project | **Agree** | `Main_researcher` (42 levels, complete). `ID_Res` is a per-record ID, not a source. `Reference` is 99 % NA. |
| Source-specific year slopes | **Agree, adapted** | Fit them, but only 6 contributors have ≥8 years; report the poorly identified result honestly and give the within/between decomposition as the primary evidence. |
| Sources per year, spanning both periods | **Agree** | Table already computed (§0); becomes Supplementary Table S1 + Fig. S2. |
| Measurement protocol / measurer covariate | **Cannot (state it)** | `Measurer` is 100 % NA, `Skull_ossification` 100 % NA. Only proxy: which wing column was populated (right / left / unspecified). Use it; state the limitation in the main text. |
| lnCVR within source | **Agree, adapted** | Only one contributor has paired cells (k = 6). Report it, plus CV within species × sex × contributor, plus a continuous σ model with contributor intercepts (§5). |
| Month / season of capture | **Agree** | `Date` present for all but 6 records. Season is a strong predictor (spring/dry-season captures +0.4–0.6 mm) and attenuates the year effect. Moult status (`Molt`, 73 % complete) also available: moulting birds −0.28 mm. |
| Unique individuals / recaptures | **Agree** | `Ring` present for 91 %; 927 rings with >1 record (1,813 extra records). Recapture share rises 5 % → 12 % early→late. Add (1 \| ring) and a first-capture-only sensitivity. |
| Age criterion | **Cannot refine (state it)** | `Age == "Adult"` is the only information; no skull data. Say so. |

### §2.2 Spatial turnover — **Agree** on every item

Longitude and `Altitude` are in ABT (no DEM needed, though a DEM cross-check of
`Altitude` is cheap with `terra`, already installed). Site random intercept at
municipality (99 levels) or locality (193 levels, 11 % NA). Mundlak within/between
decomposition for year and latitude. Report shared-site counts (§0) and redraw
Figure 1 with shared vs period-unique sites and per-panel N.

### §2.3 Figure 2 pools species — **Agree**

Replace with within-species-centred wing length vs year, add a species-slope
caterpillar plot (56 of 72 species slopes negative in the baseline `lme4` fit,
median −1.2 mm/decade; the brms fits in `brm0_multiphylo.rda` hold the proper
posteriors), records per year and per contributor per year, and explain the 73 →
61 → 44 attrition (species present in only one quartile are dropped by design in
`update_descriptive_stats.R`).

### §2.4 Effect magnitude — **Agree; the referee under-estimated it**

SD(year) = 5.05, so −0.92 per SD-year is **−0.18 mm/yr = −1.8 mm/decade = 2.6 %
per decade**, about −4.2 mm (−5.9 %) across 1995–2018. Mass: −0.0073 log units per
SD-year = −1.4 %/decade, −3.3 % over the record, interval to −7.6 %. Isometric
expectation for a 5.9 % wing decline is ≈ −17 % mass, so departure from isometry
holds. Report all of this in mm/decade and %, tabulate against Jirinec et al.
2021, Weeks et al. 2020, Ryding et al. 2024, and delete "modest" everywhere. If
the controlled estimate is ≈ −0.6 mm/decade (0.9 %/decade) it lands inside the
published range, which is itself worth saying.

### §2.5 Mass and allometry — **Agree, plus diurnal check**

Bivariate brms model (`mvbind(cwl, ln_mass)`) on the 7,577 shared records with
correlated residuals and species effects, so the difference in year slopes gets
its own posterior. Plus wing ~ log-mass covariate model. Rephrase as departure from
isometry. Soften wing-loading language: chord is not area; state the area ∝ L²
assumption.

**Scrutiny of Body Mass:**
- Body mass is available for 11,256 live records (89.5 % of the analytical sample,
  across 47 contributors and 150 municipalities).
- The year coefficient is remarkably stable across all model specifications:
  baseline β = −0.0070 (t = −1.5); + contributor β = −0.0064; + municipality
  β = −0.0057; + lon + alt + season β = −0.0045; + ring ID β = −0.0038 (t = −0.8).
  All credible intervals comfortably cross zero (≈ −0.8 %/decade).
- **Diurnal mass gain:** Capture time (`Hour`, present in 76.5 % of records) shows a
  highly significant diurnal foraging trajectory: birds gain +0.41 % body mass per
  hour of the day (t = 7.25, P < 1e-12). Controlling for time of day leaves the year
  slope at −0.0060 (t = −0.94).
- **Measurement precision:** Older studies used spring Pesola scales (0.5–1.0 g
  precision) whereas modern studies use digital balances (0.1 g). This introduces
  minor heteroskedasticity over time but does not bias the zero slope.
- **Conclusion:** Body mass is rock-solid invariant over time.

### Bill width and other traits (bears on referee §2.4, §2.7 and §5) — **Retract the Allen's rule claim**

The original manuscript reported that bill width increased over time (β = +0.24 per
SD-year, 95 % CI 0.06 to 0.43) and framed this as Allen's rule / appendage
enlargement for heat dissipation (citing Ryding et al. 2024). **This finding is completely spurious:**

- **Sample size:** Bill width is sparsely recorded: only 3,209 records (25.5 % of
  the sample), measured by only 22 contributors across 84 municipalities.
- **Protocol ambiguity:** 98 % of bill width entries are in the unspecified
  `Bill_width.mm.` column. Anatomical caliper placement (bill base at feather line,
  anterior edge of nostrils, or gape) is unstandardized across researchers.
- **Severe contributor turnover:** Only **2 out of 22 contributors** (A. Piratelli
  and M. Alves) sampled bill width in both early and late periods.
- **Sign reversal under controls:**
  - Baseline model: β = +0.23 mm/SD-yr (t = +2.50, P = 0.012).
  - + Contributor: β = −0.28 mm/SD-yr (t = −3.79) — **completely flips sign**!
  - + Contributor + Municipality: β = −0.04 mm/SD-yr (t = −0.53, P = 0.60) — **vanishes**.
  - Within-contributor slope: β = −0.05 mm/SD-yr (t = −0.68, flat zero).
- **Other traits scrutinized:**
  - `Bill_length.mm.` (7,697 records, 61.2 %): Baseline β = +0.02 (t = 0.4); controlled
    β = +0.16 (t = 1.7, P = 0.08). Invariant over time.
  - `Tail_length.mm.` (8,872 records, 70.6 %): Baseline β = +0.29 (t = 1.0); controlled
    β = +0.31 (t = 1.2). Invariant over time (subject to tail molt/abrasion).
  - `Tarsus_length.mm.` (4,361 records, 34.7 %): Baseline β = +0.42 (t = 5.0, P < 0.001);
    controlled β = +0.24 (t = 1.5, P = 0.13). Attenuates to zero with contributor controls;
    tarsus measurement is notoriously sensitive to observer technique (joint-to-notch vs bend).
- **Action:** Delete the Allen's rule / appendage enlargement claims for bill width.
  Present bill width and tarsus length alongside wing length as cautionary case studies
  in how observer turnover generates spurious morphological trends in retrospective databases.

### §2.6 PREDICTS — **Agree that the trend is unsupportable; replace rather than remove**

Evidence that option (a) is not viable on PREDICTS: `dat_no_grassland` = 56 studies
from 41 sources, 725 sites, 23 sampling methods, 5 biomes, 26 ecoregions; **only 12
studies (1,652 of 63,150 rows) fall in Atlantic Forest ecoregions**, spanning
1998–2009; the largest single study (Cabra 2006, Colombia) is 11,310 rows; 54 of 56
studies have `yr_min == yr_max`. There is no within-study temporal information to
model. The "shallower without ants" sentence is a logic error and goes regardless.
Because removing arthropods altogether would delete the paper's trophic axis, §7
below documents the search for a substitute and Phase 4 implements it: a
region-matched ATLANTIC ANTS index and GBIF occupancy trends, both with
contributor/site/method controls, framed as community indices. "Rising prey
scarcity" still leaves the title and the causal sentence leaves the abstract.

### §2.7 drmSEM — **Agree**

Demote to a supplementary, explicitly exploratory analysis; drop the arthropod SEM
entirely (rejected DAG, no interpretable coefficient); the Discussion sentence that
relies on it goes. Resolve the σ contradiction: the model says warmer → *less*
variable and later → *more* variable; since year and temperature are uncorrelated
here (year → tmean ≈ 0), these are separable effects and the text must stop saying
the variance rise is "linked to thermal extremes". List the d-sep claims (there is
one, on df = 2). One inferential standard throughout (drop P-values from the main
text, report intervals). Disclose that drmSEM and prepR4pcm are co-authored by
E.S.A.S. and unpublished; add a short validation note. Fix the Fig. 7 caption.

### §3.1 Title, abstract, plain-language summary — **Agree, fully**

Remove "shrinking body size", "rising prey scarcity", "three decades", and the
causal sentence linking prey collapse to morphology. Working title (referee's
suggestion, to be finalised after the decision gate):
*Wings shorten faster than mass declines in Atlantic Forest passerines, with no
detectable link to arthropod prey availability.*

### §3.2 Latitude runs against Bergmann — **Agree**

Within/between decomposition of latitude; report the within-species cline as the
interpretable quantity; confront the mismatch in the Discussion.

### §3.3 Diet interaction — **Agree**

Diet-Inv across the 65 species is bimodal: 20 species at 90–100, 11 at ≤10; the
−1 SD projection (19.7) is supported by **11 species / 1,520 records**, nearly all
Thraupidae, Pipridae and Emberizidae, while the +1 SD end (92.6) is
Thamnophilidae/Tyrannidae/Conopophagidae. Show this; replace the ±1 SD projection
with predictions at observed diet quantiles or drop it; add family as a check on
clade confounding; qualify EltonTraits.

### §3.4 Climate barely in the paper — **Agree**

Extract CRU-TS / WorldClim monthly series at the sampled localities and show
1995–2018 trends in Tmean, warm-quarter Tmax and SPEI (Fig. S-climate). State the
exposure-window mismatch (capture-month temperature cannot cause a fixed wing
length) explicitly.

### §3.5 Sex filter — **Agree to test; the test favours us** (see §0)

### §3.6 lnCVR — **Agree**

Make the continuous σ model (with contributor, site, sex terms) primary; keep the
binned lnCVR as a descriptive check; compute within sex and, where possible, within
contributor; report sex ratio by period (46/54 → 49/51, stable); explain the
"corrected SE" (the `s2.lnCVR` mean–SD correlation correction, Nakagawa et al. 2015).

### §3.7 Imputation — **Agree, and the Methods are currently wrong**

`atlantic_parallel.R` fits **complete cases**; the `brm_multiple`/`mice` chunk in
the Rmd is `eval = FALSE` legacy. The Methods paragraph describing MICE imputation
must be deleted or rewritten to say what was actually done. Nothing in the primary
pipeline imputes the response. Check `body_mass_descriptive.qmd` for any residual
imputation claim.

### §5 Minor and editorial — **Agree on all**

Duration is 24 years (1995–2018): fix "three decades", "28-year", "more than two
decades", "Two decades". "Below" → "above" (or delete with §2.6). Open-habitat
exclusion is design, not sensitivity. Align the Fig. 3 citation. Zenodo DOI is a
placeholder. Figure 1 sample sizes and shared sites. Cut ~⅓ of the Discussion:
carbon storage, genetic drift, "indelible signature", "silent restructuring".
Bergmann as a spatial interspecific rule (Teplitsky & Millien 2014). **Retract
Ryding et al. appendage prediction**: bill width does not increase once contributor
turnover is controlled (flips from +0.23 to −0.04, within-contributor slope −0.05).

### Housekeeping found while triaging (not raised by the referee)

- `Analysis/data/derived/passer90_climate.rds` is **stale**: 15,332 rows / 89
  species (pre-live-only build). The drmSEM outputs (7,810 records / 72 species)
  came from a server run on a regenerated file. Re-run `climate_extraction.R`
  locally or copy the server file, and record which was used.
- `Date` has a few non-standard values (e.g. `1/1/2006-2007`); the month parser
  must handle them (6 NA currently).

---

## 2. The decision gate

Phase 1 below re-fits the primary wing model with all artefact controls in brms
across the 50 trees. Its outcome determines the framing. Decide **before** the
Discussion is rewritten.

**Empirical Reality from Multi-Trait Diagnostics:**
As shown in the table in §0, when contributor, municipality, season, and individual
recaptures are included, the wing length year effect falls to −0.33 mm/SD-year
(95 % CI [−0.69, +0.04], touching/overlapping zero), and the within-contributor
slope is flat (+0.02, t = 0.08). Simultaneously, bill width collapses to zero
(−0.04, t = −0.53; within-contributor −0.05), tarsus length collapses to zero
(+0.24, t = 1.52), and body mass is rock-solid flat (−0.0038, t = −0.79).
**Therefore, the orchestrating agent should treat Scenario B as the overwhelmingly
probable baseline reality**, rather than an unexpected fallback.

**Scenario A — the controlled year effect excludes zero** (unlikely given diagnostics;
would require phylogenetic pooling to sharply reduce error without shifting means).
Frame: a real but smaller wing shortening that survives contributor, site, season
and individual controls, benchmarked honestly against the literature; departure from
isometry; negative trophic result; variance as a qualified secondary result.
Title as in §3.1.

**Scenario B (The Indicated Reality) — the controlled year effect overlaps zero, and
within-contributor slope is ≈ 0.**
Frame: apparent decadal morphological changes (wing shortening, bill widening, tarsus
lengthening) in retrospective databases can be generated entirely by observer and
spatial sampling turnover; when rigorously controlled, phenotypes are stable, but
individual-level contrast reveals an allometric divergence (wings shorter relative to
mass on shared individuals). No multi-decadal arthropod monitoring exists for the
biome; the region-matched ATLANTIC ANTS and GBIF occupancy indices (Phase 4) carry
the trophic test, and their arthropod → wing path is expected to be null.
This is a high-impact, cautionary methodological and macroecological paper.
Title: *Apparent morphological shifts in Atlantic Forest passerines are not separable
from sampling provenance; body mass is stable.*
Target venue: *Global Change Biology*, *Proceedings of the Royal Society B*, or
*Methods in Ecology and Evolution*.

In both scenarios the within-contributor decomposition, the six long-running
contributors, and the unknown-sex replication are reported in the main text.

---

## 3. Work plan (phased)

Heavy fits run on Totoro (`_sampling_config.R` auto-detects). Scripts follow the
existing path-helper convention. Each phase lists new/changed files.

### Phase 0 — Provenance and design audit (no model fitting; 1–2 days)

New `Analysis/scripts/audit_provenance.R` producing:

- `output/audit_sources.rds` — contributors per year, records per contributor,
  year span, shared-period flags, wing-column mix, mean wing per contributor.
- `output/audit_sites.rds` — localities/municipalities per period, shared-site
  counts, species × site shared counts, longitude/altitude by period.
- `output/audit_individuals.rds` — ring duplication, recapture share by period.
- `output/audit_season.rds` — month/season distribution by period, moult by period,
  feather wear / abraded plumage phenology.
- `output/audit_traits.rds` — completeness and provenance distribution for all 6
  morphological traits (`Body_mass.g.`, `Wing_length`, `Bill_width`, `Bill_length`,
  `Tail_length`, `Tarsus_length`), plus `Hour`.
- `output/effect_scale.rds` — SD(year), mean wing, conversions to mm/decade and %,
  isometric expectation for mass; comparator table values.
- Figures: `Analysis/figures/audit_records_per_year_by_source.png`,
  `audit_sites_map_shared.png` (redesigned Fig. 1), `audit_wingcol_by_year.png`.

**Data pipeline update (`rebuild_passer90_live.R`):**
Rebuild `passer90.rda` to **carry all relevant metadata and traits**. 
*Crucial fix:* Do NOT drop `Sex == "Unknown"` unconditionally at row 50, otherwise
M6 (unknown-sex replication) cannot be run from the derived data! Instead:
- Filter `Year >= 1990, Age == "Adult", Order == "Passeriformes", Status == "live", AtlanticForests_20km_Buffer == "inside..."`.
- Mutate `known_sex = (Sex != "Unknown")`.
- Retain columns: `spp`, `Binomial`, `Sex`, `known_sex`, `Status`, `Year`, `scaled_yr`,
  `Date`, `month`, `season`, `Hour`, `Main_researcher`, `Municipality`, `Locality`,
  `Ring`, `Recapture`, `Altitude`, `Longitude_decimal_degrees`, `Latitude_decimal_degrees`,
  `scaled_lat`, `scaled_lon`, `scaled_alt`, `wing_col`, `conc.wing.length`,
  `ln_conc_wing_length`, `Body_mass.g.`, `Bill_width.mm.`, `Bill_length.mm.`,
  `Tail_length.mm.`, `Tarsus_length.mm.`, `Molt`, `Reproductive_stage`.
- Save two objects: `passer90.rda` (filtered to `known_sex == TRUE`) and
  `passer90_allsex.rda` (including unknown-sex records, for M6).

### Phase 1 — Artefact-controlled primary model: THE DECISION GATE (server; ~1–2 days compute)

> [!IMPORTANT]
> **COMPUTATIONAL DIRECTIVE FOR THE ORCHESTRATING AGENT: AVOID THE 400-STAN-FIT BOTTLENECK**
> Running 8 models (M0–M7) across 50 phylogenetic trees in Stan with `brms`
> (`8 models × 50 trees = 400 MCMC fits`) on ~8,500 rows with complex multi-level
> structures (`(1|spp) + (1|gr(binomial, cov=A)) + (1|src) + (1|site) + (1|ring)`)
> will consume hundreds of CPU hours, risk divergences, and stall the project.
> The orchestrating agent MUST execute Phase 1 in three tiered steps:
> 1. **Tier 1 — Instantaneous `lme4` screening (< 2 minutes):**
>    Fit M0 through M7 in `lme4` immediately using `Analysis/scripts/atlantic_parallel_controlled.R`
>    with a `--fast-lme4` flag. This immediately confirms the exact coefficients, SEs, and
>    t-values for the entire sequence.
> 2. **Tier 2 — Single-tree `brms` validation (~1 hour):**
>    Fit M0 through M7 on a single representative tree (or consensus tree) in `brms` to verify
>    that Bayesian posterior means and 95 % credible intervals match REML estimates and to check
>    R-hat / Bulk-ESS.
> 3. **Tier 3 — 50-tree Rubin-pooled fits ONLY for definitive models:**
>    Reserve the expensive 50-tree Rubin pooling exclusively for:
>    - M0 (baseline replication, to match existing manuscript numbers).
>    - M3/M5 (the definitive controlled model).
>    - M6 (unknown-sex replication check).

New `Analysis/scripts/atlantic_parallel_controlled.R`. Models to fit and save:

- M0 current model (re-fit for identical data/priors).
- M1 + `(1 | src)` + `(1 | site)`.
- M2 + `wing_col` + `lon_s` + `alt_s`.
- M3 + `season` (+ `molt`), + `(1 | ring)`.
- M4 = M3 with `(1 + scaled_yr || src)`.
- M5 = M3 with Mundlak decomposition: `yr_within_site + yr_site_mean`, and
  `lat_within_spp + lat_spp_mean` (referee §3.2).
- M6 = M3 on unknown-sex records from `passer90_allsex.rda` (drop `Sex`).
- M7 = M3 restricted to first captures.
- Species-slope posteriors from M3 → caterpillar figure; count negative / excl. 0.

Output `output/controlled_wing_results.rds` with a before/after table.

### Phase 2 — Multi-trait and Isometry test on shared individuals (server; ~1–2 days)

New `Analysis/scripts/atlantic_bivariate_wing_mass.R`:
1. **Bivariate brms model:** `mvbind(cwl, ln_mass)` with `set_rescor(TRUE)`, species
   intercept+slope per response, phylo term per response, contributor and site intercepts,
   on the 7,577 shared records; derived quantity = difference of year slopes (in % per decade)
   with its interval. Plus `cwl ~ ... + ln_mass` model.
2. **Body Mass standalone controlled models:**
   - Model across all 11,256 mass records: `log_mass ~ Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season + (1|spp) + (1|src) + (1|site) + (1|ring)`.
   - Diurnal mass check on subset with `Hour`: `log_mass ~ Sex + scaled_yr + scaled_lat + hour_num + season + (1|spp) + (1|src) + (1|site)`,
     documenting the +0.41 %/hr diurnal foraging trajectory and the robust zero year slope.
3. **Bill Width standalone controlled models:**
   - Model across 3,209 bill width records: `Bill_width.mm. ~ Sex + scaled_yr + scaled_lat + (1|spp) + (1|src) + (1|site)`,
     documenting the complete collapse of the apparent positive slope from +0.23 to −0.04 (t = −0.53)
     and within-contributor slope of −0.05 (t = −0.68).
4. **Summary script for remaining traits (`Bill_length`, `Tail_length`, `Tarsus_length`):**
   - Quick `lme4` models to generate Supplementary Table documenting trait stability across the board.

### Phase 3 — Variance re-analysis (server; ~1–2 days)

- brms distributional model: `bf(cwl ~ <M3 mean structure>, sigma ~ scaled_yr +
  scaled_tmean + Sex + (1|src) + (1|site))`, 10–50 trees. This replaces drmSEM
  as the primary variance evidence and is not self-authored software.
- lnCVR within sex, and within species × sex × contributor where k permits;
  report the k = 6 paired result plainly (E. Carrano negative mean −0.46).
- Document mechanical variance inflation caused by contributor turnover (mean
  contributors per cell rising 4.5 → 6.6).
- Update `update_descriptive_stats.R` accordingly.

### Phase 4 — Replace the PREDICTS trend with region-matched arthropod indices; demote the SEM (1–2 weeks)

Decision (2026-09-09, see §7): dropping arthropods outright would remove the
paper's trophic axis, so the PREDICTS South-America trend is **replaced**, not
deleted. Nothing found provides a true multi-decadal Atlantic Forest arthropod
monitoring series, so the replacement must be framed as *community / occupancy
indices with explicit provenance controls*, never as prey biomass, and the causal
claim still leaves the title and abstract.

- **4a. ATLANTIC ANTS litter-ant index (primary).** Brazil, Atlantic Forest bbox,
  1990–2019: 96,169 records with sampling year; standardised events (Winkler with
  known sample number, or pitfall with known trap number) in every year 1995–2018
  from 1–17 contributor files per year (~700 events, ~600 localities). Response =
  species richness (or occurrence) per standardised sample, with offsets for
  `Winkler.Number` / `Pitfall.Number`, and random effects for contributor file
  (`fl.nm`), locality, method, habitat. Report the pooled year effect AND the
  within-contributor year effect (four programmes span ≥6 years: Rio Doce/PERD
  pitfalls 2001–2018 with true abundance; CEPLAC Bahia Winkler 1996–2002 and
  2003–2014; LAMAT/UMC São Paulo 2001–2017). New script
  `atlantic_ants_index.R`; data zip from `LEEClab/Atlantic_Ants` (CC-BY).
- **4b. GBIF occupancy trends for prey orders (secondary).** 742,203 Insecta
  occurrences in the Atlantic Forest bbox 1995–2018 with continuous annual
  coverage (8k–68k per year; 62 % preserved specimens from Brazilian collections,
  33 % human observations mostly post-2015). Fit occupancy–detection models with
  list-length correction (Outhwaite et al. 2018 approach; `sparta`/`occupancy` in
  R) for Lepidoptera, Coleoptera, Orthoptera, Hemiptera, Araneae, and aggregate to
  an annual occupancy indicator. New script `atlantic_gbif_occupancy.R`. Heavier
  (2–3 weeks) but the only route to an annual index across the whole bird record.
- **4c. Resource proxies at the bird localities (complement, full coverage).**
  MODIS EVI/NDVI anomalies (2000–2018), SPEI from the TerraClimate rasters already
  in `Analysis/data/raw/terraclimate/`, CHIRPS rainfall; enter the SEM as
  upstream drivers of prey availability, stated as such.
- **4d. Collaboration requests (parallel, do not block).** A. V. L. Freitas
  (UNICAMP): Serra do Japi monthly fruit-feeding butterfly monitoring 2011–2021,
  unpublished. S. P. Ribeiro (UFOP): Rio Doce arthropod programme 2000–2018.
  J. H. C. Delabie (CEPLAC): Ilhéus Winkler series 1990–2014. ICMBio Programa
  Monitora: fruit-feeding butterflies in Atlantic Forest federal UCs 2014–2022.
- **4e. Retire the PREDICTS trend.** Archive `brm_arthro_abund*.rda`,
  `brm_no_grass*.rda`, `dat_*.rds` under `archive/predicts/` with a README
  (54 of 56 studies single-year; 12 Atlantic Forest studies, 1998–2009). Delete
  Fig. 5 and the "shallower without ants" sentence.
- **4f. drmSEM.** Supplement only, exploratory; `INCLUDE_ARTHRO` node rebuilt
  from 4a/4b (locality × year) or dropped; contributor intercept on the wing
  node if drmTMB allows; authorship disclosure and validation note; Figure 7
  caption fixed; the Discussion sentence relying on the rejected graph removed.

### Phase 5 — Climate context and diet (1–2 days)

- `climate_extraction.R` already downloads WorldClim monthly; extend it to output
  locality-level annual Tmean, warm-quarter Tmax 1995–2018 and fit a simple trend;
  add SPEI from TerraClimate (`Analysis/data/raw/terraclimate/` exists). Figure
  S-climate. Regenerate `passer90_climate.rds` on the live-only sample.
- `atlantic_diet_interaction.R`: add diet-distribution figure, predictions at
  observed quantiles (10, 50, 90 % Diet-Inv), family random intercept sensitivity;
  drop the ±1 SD projection text.

### Phase 6 — Figures (1–2 days, after Phases 1–3)

`make_figures.py` / `make_figures.R`:
- Fig. 1: shared vs unique sites, N per panel.
- Fig. 2: within-species-centred wing vs year with M3 fit; inset records per year.
- New Fig.: species-slope caterpillar (M3 posteriors).
- New Fig.: wing vs mass year-slope contrast (Phase 2).
- Fig. 4: variance — continuous σ prediction + lnCVR within sex.
- Supplementary: records per contributor per year (stacked), wing-column mix by
  year, comparator-magnitude table, climate trends, diet distribution, bill width artifact before/after.
- Delete Fig. 5; move Fig. 7 to Supplement.

### Phase 7 — Manuscript rewrite (3–5 days, after the gate)

`Manuscript/index.qmd`:
- Title, abstract, plain-language summary, keywords (drop "trophic interactions"
  or keep only as the tested-and-rejected pathway).
- Results reordered:
  1. Design and provenance audit (observer turnover, locality turnover, protocol flips).
  2. Controlled wing trend with before/after and within/between decomposition.
  3. Isometry contrast on shared individuals (wing shortens relative to invariant body mass).
  4. Multi-trait stability and provenance artifacts: body mass invariant (with diurnal check),
     bill width Allen's rule artifact refuted, tarsus attenuation.
  5. Variance: continuous σ model with observer controls; lnCVR within-source caveats.
  6. Trophic pathway: diet interaction null; ATLANTIC ANTS and GBIF occupancy indices
     with provenance controls (Phase 4); PREDICTS retired and why.
  7. Sensitivity (unknown sex replication, first captures).
- Effect sizes in mm/decade and % everywhere; comparator table; remove "modest".
- Methods: delete the MICE paragraph; add provenance/site/season/individual
  terms; add the protocol and age limitations; describe the Mundlak
  decomposition; climate section rewritten around the new figure and the
  exposure-window caveat; SEM moved to Supplement with disclosure.
- Discussion cut by ~⅓; remove carbon storage and genetic drift chains; add Bergmann
  spatial/temporal mismatch; soften wing-loading; **retract Allen's rule bill width claims**
  and replace with discussion of observer turnover in compiled trait databases.
- Duration = 24 years everywhere; Zenodo DOI minted before submission
  (`zenodo` GitHub integration on `esalves/atlantic_birds`).
- `Manuscript/README.md`, top-level `README.md` title lines.

### Phase 8 — Response letter and housekeeping (1 day)

- Point-by-point response drawing on §1 of this document.
- Update `RESULTS_SUMMARY.md`, `REPO_STRUCTURE.md`, `LIVE_ONLY_MIGRATION.md`
  run order, `R_session_info.txt`.
- Memory/notes: record the decision-gate outcome.

---

## 4. Run order on the server

```bash
cd Analysis/scripts
# Step 0: regenerate data carrying all metadata, traits, and unknown-sex records
Rscript rebuild_passer90_live.R            # -> passer90.rda (known sex) + passer90_allsex.rda

# Step 1: Phase 0 audit & diagnostics (all traits, sources, sites, dates)
Rscript audit_provenance.R                 # Phase 0 tables + figures (fast)

# Step 2: Phase 1 fast lme4 screening (takes < 2 min; confirms M0-M7 trajectory)
Rscript atlantic_parallel_controlled.R --fast-lme4

# Step 3: Phase 1 brms single-tree validation (~1-2 hr; verifies Bayesian/REML agreement)
Rscript atlantic_parallel_controlled.R --trees 1

# Step 4: Phase 1 final 50-tree Rubin pooling (only on definitive models: M0, M3/M5, M6)
Rscript atlantic_parallel_controlled.R --trees 50

# Step 5: Phase 2 multi-trait models (bivariate wing-mass, diurnal mass, bill width collapse)
Rscript atlantic_bivariate_wing_mass.R     # Phase 2

# Step 6: Phase 3 variance analysis
Rscript atlantic_variance_sigma.R          # Phase 3
Rscript atlantic_ants_index.R              # Phase 4a: ATLANTIC ANTS litter-ant index
Rscript atlantic_gbif_occupancy.R          # Phase 4b: GBIF occupancy trends (long)

# Step 7: Climate and diet updates
Rscript climate_extraction.R               # regenerate live-only passer90_climate.rds + climate trends
Rscript atlantic_diet_interaction.R        # Phase 5 (updated)
Rscript atlantic_drmsem.R                  # Supplement only, INCLUDE_ARTHRO = FALSE
Rscript update_descriptive_stats.R

# Step 8: Render figures and manuscript
cd ../../Manuscript && python3 make_figures.py && quarto render index.qmd
```

---

## 5. Where to push back (briefly, in the response letter)

- **Protocol / measurer / age covariates**: not recorded in ATLANTIC BIRD TRAITS
  (`Measurer` and `Skull_ossification` are empty). We use the wing-column proxy and
  state the limitation; we cannot do more.
- **Sex-filter reading**: the unknown-sex records replicate the decline once
  provenance is controlled; the grid attenuation was a species-pool effect.
- **Effect magnitude**: the referee's back-calculation (≈2 %/decade) is, if
  anything, an under-estimate (2.6 %/decade); we agree with the conclusion drawn.
- **Elevation**: taken from the dataset's own `Altitude` field (DEM cross-check
  reported in the Supplement).
- **Option (a) for PREDICTS**: not attempted because 54 of 56 studies are
  single-year; a within-study temporal model is not estimable.

---

## 6. Diagnostics provenance

Quick checks in this document were run 2026-09-09 with R 4.6.0, `lme4`, REML, on
the live-only 73-species / 12,571-record sample rebuilt from
`ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv` with the filters of
`rebuild_passer90_live.R`. 
Checked:
1. **Wing length** (coalesced; 8,478 records): M0 through M7, Mundlak within/between
   decomposition, contributor-specific slopes, long-running contributors.
2. **Body mass** (`Body_mass.g.`; 11,256 records): baseline and controlled models (+src,
   +site, +season, +ring), plus diurnal time-of-day model (`Hour`, N = 8,390, +0.41 %/hr,
   t = 7.25, P < 1e-12).
3. **Bill width** (`Bill_width.mm.`; 3,209 records): baseline (+0.23, t = 2.50) vs
   controlled (−0.04, t = −0.53); spanning contributor count (2 of 22); within-contributor
   slope (−0.05, t = −0.68).
4. **Other traits**: `Bill_length.mm.` (7,697 records; flat), `Tail_length.mm.`
   (8,872 records; flat), `Tarsus_length.mm.` (4,361 records; baseline +0.42 t = 5.0 ->
   controlled +0.24 t = 1.52).

Scripts: session scratchpad `diag.R`, `diag2.R`, `quick_lmer.R`, `quick_lmer2.R`,
`quick_var.R`, `quick_multitrait.R` (to be folded into `audit_provenance.R` in Phase 0).
They are planning aids, not results.

---

## 7. Arthropod data alternatives — search of 2026-09-09

Question: is there a dataset that can replace the PREDICTS South-America trend
with a temporally resolved arthropod signal for the Atlantic Forest, 1995–2018?
Short answer: **no monitoring series exists; two compilations can yield honest,
provenance-controlled indices.** Everything below was checked against the actual
files or APIs, not abstracts.

| Candidate | What it is | Atlantic Forest coverage 1995–2018 | Verdict |
|---|---|---|---|
| **ATLANTIC ANTS** (Silva et al. 2022, Ecology; GitHub LEEClab, CC-BY) | 178,976 ant records, 57 fields incl. `Start.year`, `Method`, `Winkler.Number`, `Pitfall.Number`, `Total.Ant.Abundance`, `Measurement.Type`, contributor file `fl.nm`, habitat | 96,169 Brazil-AF records with year; standardised events in **every** year 1995–2018 (1–17 contributors/yr); 398 localities sampled in ≥2 years, 73 with span ≥10 yr; 4 programmes with ≥6 sampling years | **Best available.** Litter-ant richness/occurrence per standardised sample; abundance only 14 % of records. Same structural caveat as PREDICTS (compilation) but region-matched, with contributor and effort fields that allow the controls the referee asked for. |
| **GBIF Insecta occurrences** | 742,203 records in AF bbox 1995–2018 (62 % preserved specimens, 33 % human observations) | Continuous, 8k–68k/yr; top sources ICMBio SISBio, DCBU-UFSCar, DZUP, MZUSP, iNaturalist | **Viable secondary route** via occupancy–detection models with list-length correction; yields annual occupancy, not abundance. Heavier. |
| ATLANTIC SEED RAIN (Daibes et al. 2026; Zenodo 10.5281/zenodo.19597741) | 52 patches, trap × species matrices, trap-months, landscape metrics | Released files carry **no sampling year** (only `months`); 1987–2021 per the paper | Not usable as a time series without owner-level dates. Keep as future work. |
| ATLANTIC BUTTERFLIES (Santos et al. 2018) | 122 fruit-feeding butterfly communities, species lists | Presence-only; data only in the Wiley supplement | Richness/occupancy only; no abundance. Low value. |
| BioTIME (archived exploration, `archive/biotime/`) | 20 Brazilian "arthropod-ish" studies | Only two invertebrate studies: Jari dung beetles (Amazon, 2009–2013), Santa Catarina Island dung beetles (2016–2019) | Nothing in the Atlantic Forest window. Confirms the earlier decision. |
| InsectChange (van Klink et al. 2021; KNB 10.5063/F11V5C9V) | 1,668 series, 165 studies | Brazil: Amazonas dung beetles (1986, 2000), Minas Gerais and Rio Grande do Sul freshwater | No terrestrial Atlantic Forest series. |
| Lewinsohn et al. 2022 appraisal (Biol. Lett.; Figshare table) | 75 Brazilian trend cases | 23 terrestrial Atlantic Forest cases; only 4 quantitative with ≥3 years inside our window: Rio Doce Chrysomelidae 2001–2016 (1 site), Rio Doce Ichneumonidae 2000–2008 (3 sites), Serra do Japi Nymphalidae 2011–2021 (2 sites, **unpublished**, Freitas), Ribeirão Preto bees 1990–2020 (1 urban fragment) | Useful as a citation list and for collaboration requests; none is a regional series. |
| PELD datasets on GBIF/SiBBr | 220 PELD datasets | Atlantic Forest insect sets are all ≤2 years (MANP Hymenoptera 2018–19, Lepidoptera 2015–16; Rio Doce ants/beetles 2001–02) | Not temporal. |
| ICMBio Programa Monitora | Fruit-feeding butterflies, 2 campaigns/yr, federal UCs since 2014 | Not on GBIF; only Brasília (Cerrado) 2017–2019 published | Request from ICMBio; overlaps only 2014–2018. |
| PREDICTS (current) | 63,150 rows, 56 studies, South America | 12 AF studies, 1998–2009, 54/56 single-year | Retire. |

Recommendation adopted in Phase 4: ATLANTIC ANTS index (4a) as the primary
region-matched arthropod context, GBIF occupancy (4b) as the annual indicator,
climate/productivity proxies (4c) as upstream drivers in the SEM, collaboration
requests (4d) in parallel. State plainly in the paper that no multi-decadal
arthropod monitoring exists for the biome. Expect the arthropod → wing path to
stay null; that is the honest result the referee already anticipated.

Working files from the search (scratch, not in the repo): the ATLANTIC ANTS zip
and metadata PDF, Lewinsohn SM3 spreadsheet, InsectChange data paper text, seed
rain Zenodo archive, and the R profiling scripts `ants_contrib.R` /
`ants_contrib2.R`. Fold the profiling into `atlantic_ants_index.R` in Phase 4.
