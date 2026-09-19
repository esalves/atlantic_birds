# Audit: hardcoded constants in the analysis scripts

**Why this exists.** While answering referee comment c88 we found that
`secular_museum_expansion.R` sets the 1980 climate breakpoint by hand
(`bp_val <- 1980`, line 429) while the manuscript reported it as the result of a
maximum-likelihood grid search over 1900–2005 yielding "1980 ± 3 years". Neither
the search range nor the ± 3 has any basis in the code or outputs, and the grid
the script does run (1940–1995) favours 1956 by 18.4 AIC units. This audit checks
every analysis script for constants of the same kind: values that gate an
inference but are asserted rather than derived or cited.

**Method.** Two scans over `Analysis/scripts/*.R`: (1) every assignment of a bare
numeric literal to a name (26 found); (2) every inline numeric threshold inside
`filter()`/`subset()`/`if()`/`which()` and every `k * sd()` trim (131 found).
Each hit was then read in context and classified. Seeds, worker counts, tree
counts, plotting limits and optimizer start values are excluded as
non-inferential.

---

## A. Confirmed problem (corrected in the manuscript)

| Location | Constant | Finding |
|---|---|---|
| `secular_museum_expansion.R:429` | `bp_val <- 1980` | Set by hand on the strength of a code comment ("Practical climate regime break ~ 1976-1980"). The script's own grid search (`search_temp_bp`, range **1940–1995**) returns **1956** (AIC −111.5); 1980 scores −93.1, i.e. ΔAIC = 18.4. The manuscript claimed a search over **1900–2005** identifying "1980 ± 3 years". The range was wrong, the result was not what the search returned, and no source exists anywhere for "± 3 years". |

Corrected in commit `8fca277`: the text now states 1980 is an a priori boundary,
reports the grid result against it, and frames the pre/post split as a comparison
of two fixed periods. Values are pulled from the stored grid.

## B. Undocumented thresholds that change which data enter a model

These are not necessarily wrong, but none carries a stated rationale, and several
are inconsistent with the main analysis without explanation.

| Location | Constant | Concern |
|---|---|---|
| `secular_museum_expansion.R:234` | Tier 4 museum series: `n >= 10, span >= 10` | Main analysis uses `n >= 30, span >= 5`. No rationale for the change; not disclosed in the manuscript. |
| `secular_museum_expansion.R:251` | Tier 5 integrated series: `n >= 30, span >= 15` | A third, different rule. No rationale; not disclosed. |
| `secular_museum_expansion.R:407, 455, 526, 602` | `abs(x - mean) <= 4 * sd` trims | Four separate 4-SD trims on wing, mass and the isometry contrast before the annual means that the secular figures and their regressions use. No rationale for 4 SD; the effect on the reported slopes is untested. |
| `secular_museum_expansion.R:520, 590` | `Year >= 1960` | The mass and allometry secular series start in 1960 while the wing series starts in 1880. No rationale given. |
| `secular_museum_expansion.R:169` | `med_yr <- 2010` | Splits Tier 2 into halves for a coverage filter. Arbitrary; low impact (it only gates `n1 >= 5, n2 >= 5`), but undocumented. |
| `atlantic_drmsem.R:226` | `ARTHRO_RADIUS_KM <- 200L` | 200 km search radius for nearby PREDICTS sites. Commented as to purpose but not justified as to value. (drmSEM is exploratory and not in the manuscript.) |

## C. Documented and defensible

| Location | Constant | Basis |
|---|---|---|
| `atlantic_parallel_controlled.R:319` | `ANOM_MM <- 6` | Comment gives the reasoning: ≈1.5 residual SD of the M3 fit, ≈8% of a wing, and identifies the two contributor blocks that motivated it. Now also stated in the manuscript Methods. Note the "1.5 residual SD" equivalence is asserted in prose, not computed in code. |
| `rebuild_passer90_live.R:191` | `n >= 30, range >= 5` | The primary species-inclusion rule; stated in the manuscript Methods. |
| `atlantic_parallel_controlled.R:307` | `n_years >= 8, n_wing >= 200` | Defines "long-running contributors"; used only for a sensitivity tier and described in the script header. |
| `atlantic_multitrait.R:299` | `hour_num >= 5, hour_num <= 19` | Daylight capture window; now stated in the manuscript Methods. |
| `audit_provenance.R:79` | `EARLY_MAX <- 2006`, `LATE_MIN <- 2013` | Outer quartiles of the sampling record; now stated in the manuscript Methods. |
| `climate_extraction.R:134` | `SPEI_SCALE <- 12L` | SPEI-12 is the standard 12-month accumulation window. |
| `secular_museum_expansion.R:74`, `jirinec_gap_analysis.R:93`, `referee_reruns.R:67` | `SD_YR <- 5.019925`, `CTR_YR <- 2009.487` | **Verified against the data**: `attr(passer90$scaled_yr, "scaled:scale")` = 5.019925 and `"scaled:center"` = 2009.4865. Correct today, but hardcoded in three scripts, so a rebuild of `passer90` would silently desynchronise them. Recommend reading the attribute instead. |

## C2. Sensitivity results for section B (run 2026-09-18)

`Analysis/scripts/secular_threshold_sensitivity.R` refits the secular models
under the alternatives. Full tables in
`Analysis/output/referee_reruns/SECULAR_THRESHOLD_SENSITIVITY.md`.

**S1 — species-inclusion rule.** The museum-only series is null under all three
rules (−0.31%/decade, t = −0.36 as run; +0.14%/decade, t = +0.29 under either
n >= 30 rule), so the conclusion we kept ("no trend before 1980") is robust. The
**integrated live + museum series is not**: −0.66%/decade (t = −1.97, CI
excluding zero) under the Tier 4 rule, −0.36%/decade (t = −1.27) under the Tier 5
rule actually used, and −0.14%/decade (t = −0.35) under the main analysis rule.
Three undocumented rules give three materially different answers.

**S2 — 4-SD trims.** No material effect on wing or mass (null under every trim).
The isometry contrast's post-1980 slope moves across the 0.05 boundary with the
trim: p = 0.074 (3 SD), 0.057 (4 SD, as run), 0.034 (5 SD), 0.033 (untrimmed).
More trimming makes it less significant, so the as-run choice was conservative
rather than favourable, but the result is not stable.

**S3 — start year.** Body mass is insensitive (null from 1900, 1940, 1960 and
1980 starts). The isometry post-1980 slope again straddles the boundary:
p = 0.047 (1900 or 1940 start), 0.057 (1960, as run), 0.093 (1980 start).

Taken together, S2 and S3 put the post-1980 isometry result anywhere in
p = 0.033-0.093 across defensible analytic choices. That claim has already been
withdrawn from the manuscript on other grounds (referee item E6); these results
confirm it could not have been defended.

## D. Recommendations

1. `secular_museum_expansion.R` needs the Tier 4 / Tier 5 thresholds, the four
   4-SD trims and the 1960 start either justified in the script or aligned with
   the main analysis, and disclosed wherever the corresponding results are
   reported. This script is the one that produced both the fabricated breakpoint
   and the unreproducible thermal figure; it warrants the closest reading.
2. Hardcoded `SD_YR` / `CTR_YR`: DONE. All three scripts
   (`jirinec_gap_analysis.R`, `secular_museum_expansion.R`, `referee_reruns.R`)
   now read `attr(passer90$scaled_yr, ...)` with a `stopifnot` guard.
3. Sensitivity checks: DONE, see section C2. The integrated secular series and
   the post-1980 isometry slope are both threshold-dependent; the museum-only
   null and the body-mass null are robust.

## E. Scope note

This audit covers numeric constants. It does not verify that every reported
number has generating code — that is a separate problem, already found twice
(the thermal-coupling estimate and the 1980 breakpoint), both in
`secular_museum_expansion.R` and both introduced in commit `bad2941`.
