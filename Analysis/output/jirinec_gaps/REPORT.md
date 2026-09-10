# What happens if we run the analyses Jirinec et al. (2021) ran and we didn't

Exploratory side-report, 2026-09-10. **No manuscript file was modified.**
Numbers: `jirinec_gaps_tables.md`. Code: `Analysis/scripts/jirinec_gap_analysis.R`.
Tier: `lme4` (no phylogenetic term) — see caveats in `README.md`.

## Headline

None of the seven gaps overturns anything in the manuscript. Every headline null
survives. Three of them (G1, G2, G3) would materially strengthen the paper and
are worth adding; one (G4) is a clean supporting result; three (G5, G6, G7) are
supplementary robustness.

The one genuinely new substantive finding is in G1: the Atlantic Forest mass–climate
relationship is not merely absent, it is **sign-reversed** relative to Amazonia.
Birds here are slightly *heavier* after warm years, where Jirinec's birds were lighter.

---

## G1 — Lagged climate covariates in the mean models

**Verdict: add to the paper. This is the gap that mattered.**

Coverage: 11,799 of 12,571 records (93.9%) carry complete lag-0/1/2 climate;
mass n = 10,874, wing n = 7,810, 357 localities.

**Design note.** Jirinec had one site, so their climate coefficients are pure
within-site temporal anomalies by construction. Ours span 357 localities across
~15° of latitude, where raw temperature is dominated by the spatial gradient
(Bergmann-in-space). Every climate variable was therefore split into a
between-locality mean (1990–2018) and a within-locality anomaly, and the temporal
response is read off the **within** terms. Their seasonal lags are not
transferable — captures here run year-round across a biome with no common dry
season — so lags are annual (0 = capture year, 1 = previous, 2 = two years back).
Collinearity is acceptable (max VIF 4.2 on latitude; year VIF 2.3 in the
temperature set, r = 0.40 between year and the lag-0 temperature anomaly).

### Does conditioning on climate remove the year trend? No — it deepens it slightly.

| Model (mass, M3 tier) | %/decade | 95% CI |
|---|---|---|
| year only | −0.55 | −2.30 to +1.23 |
| + lagged temperature | −1.42 | −3.27 to +0.47 |
| + lagged precipitation | −0.76 | −2.53 to +1.05 |
| + lagged SPEI | −0.60 | −2.36 to +1.20 |

Jirinec's central mechanistic result was that the mass trend *vanished* under
lagged temperature. Ours moves the other way and stays firmly across zero.

### Do birds respond to interannual climate at all? Weakly, and in the opposite direction.

Mass, M3 tier, per 1 SD anomaly (SD ≈ 0.37–0.39 °C, ≈ 250–300 mm, ≈ 0.86 SPEI):

| Term | β (ln g) | 95% CI | t |
|---|---|---|---|
| temperature, lag 0 | +0.0024 | −0.0018 to +0.0066 | 1.14 |
| temperature, lag 1 | **+0.0044** | +0.0006 to +0.0082 | 2.26 |
| temperature, lag 2 | **+0.0050** | +0.0009 to +0.0091 | 2.39 |
| precipitation, lag 1 | −0.0039 | −0.0079 to +0.0001 | −1.92 |

Warmer (and, consistently, drier) previous years predict **heavier** birds —
roughly +1.2% body mass per °C of prior-year anomaly. Jirinec found the reverse:
lighter birds after hot, dry dry-seasons, strongest at lag 2.

### Model fit: climate buys nothing for mass here

ΔAIC (ML refits, mass, M3): temperature 0.00, year-only 0.35, precipitation 2.83,
SPEI 3.39. In Jirinec's data precipitation beat temperature by **27.6 AIC** for
mass; here the climate variants are indistinguishable from a model with no
climate at all.

For **wing**, SPEI is the best variant (ΔAIC 2.1 over year-only at M3), driven by
a positive lag-1 SPEI coefficient (+0.19 mm/SD, 95% CI +0.06 to +0.32): wetter
previous years, marginally longer wings. The wing year slope is essentially
unchanged by any climate variant (−0.98 → −0.79 %/decade).

**Why this matters.** The manuscript's null is currently vulnerable to "24 years
and +0.25 °C/decade is too little signal." This analysis answers that: mass does
not track interannual thermal or hydric anomalies either — where signal-to-noise
is far better than a 24-year trend — and where a relationship is detectable, it
runs opposite to the Amazonian one. That converts a weak null into a strong one.

---

## G2 — Species-level mass slopes

**Verdict: add. Needed to make the Amazon contrast like-for-like.**

| | Atlantic Forest (M3, this run) | Amazon (Jirinec) |
|---|---|---|
| species | 73 | 77 |
| negative point estimates | 45 (62%) | 77 (100%) |
| 95% CI entirely negative | **6 (8.2%)** | **36 (47%)** |
| 95% CI entirely positive | 2 | 0 |
| CI spans zero | 65 | 41 |

Median species slope −0.95%/decade, range −9.9 to +10.1.

**Internal validation.** The same code run on wing gives 48 negative, 9 credible
declines, 2 credible increases, 61 spanning zero of 72 — reproducing the
manuscript's published 50-tree phylogenetic numbers (49 / 9 / 2 / 61) almost
exactly. The lme4 tier tracks the phylo tier closely.

Currently the manuscript compares its *pooled* mass estimate against Jirinec's
*species tally*. This makes the comparison direct.

---

## G3 — mass:wing ratio (their headline metric)

**Verdict: add. Sharpest contrast in the whole comparison.**

7,577 shared records, 72 species. Response = ln(mass/wing), so slopes are %/decade.

| Tier | %/decade | 95% CI |
|---|---|---|
| M0 baseline | **+1.93** | +0.39 to +3.51 |
| M1 (+ observer, municipality) | +0.93 | −0.69 to +2.58 |
| M3 (full controls) | +0.38 | −1.27 to +2.06 |

Species tally at M3: **0 of 72** with CI entirely negative (1 positive, 71 span
zero); only 20 of 72 even have a negative point estimate.
Jirinec: **53 of 77 (69%)** entirely negative, all 77 negative means.

Two things worth noting. First, this is a flat contradiction of their most
prevalent result, on their own chosen metric. Second, the M0 → M3 collapse is
itself a provenance-control story: the unadjusted ratio trend is significant,
positive, and entirely an artifact of observer and site turnover.

(The ratio remains statistically inferior to the manuscript's bivariate
delta-method isometry contrast — it conflates numerator, denominator and their
covariance. It is worth reporting *because* it is the comparator's metric, not
because it is the better construction.)

---

## G4 — Phylogenetic signal in the slopes

**Verdict: clean supporting result; include as supplementary.**

| Slopes | Moran's I (E = −0.014) | median p | Pagel's λ | median p |
|---|---|---|---|---|
| mass, M0 | +0.026 | 0.011 | 0.088 | 0.143 |
| mass, M3 | +0.000 | 0.368 | 0.034 | 0.520 |
| wing, M0 | −0.011 | 0.839 | 0.006 | 0.889 |
| wing, M3 | −0.019 | 0.770 | 0.000 | 1.000 |
| ln(mass:wing), M3 | −0.027 | 0.446 | 0.000 | 1.000 |

Pooled across the same 50 trees used by the manuscript. Essentially **no
phylogenetic structure in the temporal responses** — the weak mass M0 signal is
inconsistent between the two statistics and disappears under controls.

Jirinec reported "a high degree of phylogenetic correlation" with clades differing
in rate. Our interspecific heterogeneity is unstructured, which is what sampling
noise around a null looks like. This supports the manuscript's "idiosyncratic and
species-specific" framing while also cautioning against reading biology into the
species-slope spread.

Note the comparison is not symmetric in our favour by accident: their signal test
was run on HMSC Beta coefficients from models that *already contained phylogeny*
in the hierarchical layer, so it is partly circular. Neither M0 nor M3 here
contains a phylogenetic term.

---

## G5 — Vertical foraging stratum as a moderator

**Verdict: does not survive multiplicity. Report as a hypothesis at most.**

66 of 73 species matched to EltonTraits: understory 23, ground 21, midhigh 17,
canopy 5. Twelve year × stratum interaction tests (3 traits × 2 tiers × 2
parameterisations). Eleven are null.

The exception: ln(mass:wing), M3, continuous above-understorey index,
β = −0.0051 (95% CI −0.0099 to −0.0003, t = −2.06). Higher-foraging species
show marginally steeper wing-loading decline — the same metric and direction as
Jirinec's single trait result (mass:wing declining for midstory species). But it
is one of twelve tests at nominal p ≈ 0.04, and the categorical version of the
same test is null (t = −1.37).

Worth a supplementary line as an intriguing parallel. Not a finding.

---

## G6 — Design-matched restricted tier

**Verdict: reassuring. Supplementary.**

Approximating their single-site, single-team design: 8 contributors with ≥8-year
span and ≥150 records → 4,114 records, 73 species, 1995–2017.

| Trait | M0 | M1 | M3 |
|---|---|---|---|
| mass (%/dec) | −1.48 | −1.46 | **−1.30** [−3.01, +0.44] |
| wing (%/dec) | −1.73 | −0.49 | **−0.29** [−1.40, +0.81] |
| ln(mass:wing) (%/dec) | +0.97 | +0.04 | **−0.11** [−1.90, +1.71] |

All null. The mass point estimate is somewhat more negative than the full sample
(−1.30 vs −0.55) with a wider interval, as expected on a third of the data. Note
that even *within* long-span contributors the wing M0 → M1 collapse is large
(−1.73 → −0.49): the raw trend is between-contributor structure, not within.

The null is not an artifact of pooling heterogeneous short-span contributors.

---

## G7 — Their symmetric coverage filter

**Verdict: no material change. One line in the supplement.**

Median year 2010. Applying Jirinec's ≥5-records-before-and-after criterion:
67 of 73 species pass for mass, 62 of 73 for wing. So 6 (mass) and 11 (wing)
species in our sample would not have entered their analysis.

Refit on the passing species: mass M3 −0.77%/decade [−2.54, +1.03]; wing M3
−0.89 [−1.83, +0.06]; ln(mass:wing) M3 +0.34 [−1.31, +2.02]. Indistinguishable
from the full sample. The filter difference drives nothing.

---

## If any of this goes into the manuscript

1. Refit the headline models with `--trees 50` (glmmTMB `propto` phylogenetic
   tier) — everything above is the lme4 tier.
2. G1 needs its own short Methods paragraph on the within/between climate
   decomposition and the annual-lag substitution. Both are defensible but both
   are departures from the comparator and must be stated.
3. G3 should be reported *alongside*, not instead of, the existing isometry
   contrast, with a sentence on why the ratio is the weaker construction.
4. G5's single hit must be presented with its multiplicity context or not at all.
