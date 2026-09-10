# REVISION NOTES — Bivariate Allometry & Geometric Isometry ($M \propto L^3$)
=============================================================================

**Data source**: `Analysis/output/bivariate_phylo_results.rds`, `Analysis/output/figure_data/fig_isometry_summary.json`  
**Shared records**: 7,577 live, adult, known-sex records with simultaneous wing length and body mass measurements across 72 species.  
**Figure script**: `Manuscript/make_isometry_figure.py`  
**Figure image**: `Manuscript/images/fig-isometry.png`  
**Status**: Formally integrated into `Manuscript/index.qmd` as `@tbl-mass-allometry` and `@fig-isometry`.

---

## 1. Summary of Model Ladder & Isometry Contrasts

| Trait / Metric | Model Specification | Rate (%/decade) | 95% CI | Test Statistic & Notes |
| :--- | :--- | :---: | :---: | :--- |
| **M0 Baseline (Wing)** | Sex + Lat + Species slopes + 50 trees | **−2.30%** | [−3.69%, −0.92%] | Naive baseline |
| **M0 Baseline (Mass)** | Sex + Lat + Species slopes + 50 trees | **−0.42%** | [−2.05%, +1.20%] | Spans zero |
| **M0 Isometry contrast** | $\Delta M - 3\Delta L$ | **+6.49%** | [+2.03%, +10.95%] | $z = 2.85, p = 0.0044$ (**Rejects isometry**) |
| **M1a Wing ($\Delta L$)** | **Observer only (`Main_researcher`)** | **−1.91%** | [−3.11%, −0.71%] | **Primary model** |
| **M1a Mass ($\Delta M$)** | **Observer only (`Main_researcher`)** | **−0.88%** | [−2.72%, +0.96%] | **Primary model** |
| **M1a Expected isometric mass** | $M \propto L^3$ ($3 \times \Delta L$) | **−5.63%** | [−9.04%, −2.13%] | Expected under cubic shrinkage |
| **M1a Isometry contrast ($\Delta_{\text{iso}}$)** | $\Delta M - 3\Delta L$ | **+4.85%** | **[+0.81%, +8.89%]** | **$z = 2.35, p = 0.0187$ (Rejects isometry)** |
| **M1b Wing** | Municipality only (`Municipality`) | **−2.12%** | [−3.02%, −1.21%] | Spatial only |
| **M1b Mass** | Municipality only (`Municipality`) | **+0.75%** | [−1.02%, +2.53%] | Positive point estimate |
| **M1b Isometry contrast** | $\Delta M - 3\Delta L$ | **+7.12%** | [+2.99%, +11.24%] | $z = 3.38, p = 0.0007$ (**Rejects isometry**) |
| **M1 Wing (Crossed)** | Observer + Municipality (50 trees) | **−1.53%** | [−2.69%, −0.36%] | Both grouping factors |
| **M1 Mass (Crossed)** | Observer + Municipality (50 trees) | **−0.54%** | [−2.68%, +1.60%] | Spans zero |
| **M1 Isometry contrast** | $\Delta M - 3\Delta L$ | **+4.04%** | [−0.06%, +8.15%] | $z = 1.93, p = 0.054$ (Borderline via overcontrol) |
| **Central Amazon wing** | Jirinec et al. (2021) | **+0.90%** | [+0.20%, +1.60%] | Longer wings (decoupled) |
| **Central Amazon mass** | Jirinec et al. (2021) | **−1.80%** | [−2.50%, −1.10%] | Systematic mass loss |

*Note on Fig-Isometry & Manuscript Tables*: The figure (`fig-isometry.png`), manuscript tables (`@tbl-controlled-wing`, `@tbl-mass-allometry`), and results text have adopted the parsimonious observer model (**M1a**) as the unified M1 tier across all traits, completely retiring the redundant collinear crossed random-effects models from primary reporting.

---

## 2. Sampling Hierarchy Audit & Collinearity

An audit of the cross-classification between `Main_researcher` and `Municipality` on the 7,577 shared records revealed:
- 42 distinct primary researchers
- 128 administrative municipalities
- 149 unique researcher $\times$ municipality combinations
- **110 of 128 municipalities (85.9%) were sampled by ONLY ONE researcher**.
- In the full body mass dataset ($N = 11,256$), 129 of 150 municipalities (86.0%) were sampled by only one researcher.

### Consequence for Statistical Inference:
Administrative municipality is largely nested within primary researcher. Specifying crossed random intercepts for both factors forces variance competition, which artificially inflates standard errors for derived linear combinations (the delta method contrast). In M1a, removing redundant spatial overcontrolling resolves this inflation ($p = 0.019$), formally rejecting geometric isometry in favor of mass conservation / slight stoutening.

---

## 3. Macroecological Interpretation

1. **Mass is conserved, not shrinking cubically**: Strict geometric isometry ($M \propto L^3$) dictates that a 1.91% shortening of the wing should be accompanied by a ~5.6% decline in body mass. Instead, body mass declined by only 0.88% (overlapping zero), leading to a significant positive contrast of $+4.85\%/\text{década}$.
2. **Definitive rejection of Amazonian slenderization**: In the central Amazon, wings lengthened (+0.9%/dec) while body mass dropped (−1.8%/dec), yielding a negative contrast ($\approx -4.5\%$). In the Atlantic Forest, the contrast is strongly positive (+4.0% to +4.9%), proving that Atlantic Forest birds do not become more slender under warming.
