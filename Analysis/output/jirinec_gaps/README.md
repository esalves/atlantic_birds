# Jirinec-gap analysis — exploratory side-report

**Status: exploratory. Nothing here feeds the manuscript.** Every file in this
directory is produced by `Analysis/scripts/jirinec_gap_analysis.R`, which reads
the manuscript's analytical sample but writes only to
`Analysis/output/jirinec_gaps/` and `Analysis/figures/jirinec_gaps/`. No file
that `Manuscript/index.qmd` reads is touched.

## Reproduce

```
Rscript Analysis/scripts/jirinec_gap_analysis.R            # lme4 tier (~14 min)
Rscript Analysis/scripts/jirinec_gap_analysis.R --trees 50 # + glmmTMB propto phylo tier
Rscript Analysis/scripts/jirinec_gaps_report.R             # -> jirinec_gaps_tables.md
```

## Files

| File | Contents |
|---|---|
| `jirinec_gaps_results.rds` | all model tables, species slopes, signal tests, coverage |
| `jirinec_gaps_tables.md`   | flat markdown dump of every table |
| `REPORT.md`                | the interpretive write-up |

## What is tested

Seven analyses Jirinec et al. (2021, *Sci Adv* 7:eabk1743) performed on the
BDFFP Amazonian data that the Atlantic Forest manuscript does not:

- **G1** lagged climate covariates in the *mean* models
- **G2** species-level mass slopes
- **G3** mass:wing ratio as a response
- **G4** phylogenetic signal in the *slopes*
- **G5** vertical foraging stratum as a moderator
- **G6** design-matched restricted (single-team) tier
- **G7** their symmetric ≥5-before/≥5-after coverage filter

## Caveats that apply to every number here

1. **lme4 tier only.** The phylogenetic tier (`--trees 50`) has not been run.
   The G2 wing tally reproduces the manuscript's published 50-tree numbers
   almost exactly (9 credible declines, 2 increases, 61 spanning zero of 72),
   so the lme4 tier tracks the phylo tier closely — but anything quoted in a
   manuscript must be refit with `--trees 50`.
2. **Annual, not seasonal, climate lags.** Jirinec's dry/wet-season lags are not
   transferable: our captures run year-round across a biome with no common dry
   season. Lags here are calendar-year 0/1/2.
3. **Climate is noisier than theirs.** They had one ERA5 cell at one site; we
   have downscaled CRU-TS at 357 localities across ~15° of latitude.
4. **Species slopes are shrunk BLUPs**, so signal tests on them (G4) are biased
   by shrinkage. This applies equally to Jirinec's HMSC Beta medians — and their
   test is additionally circular, because HMSC already had phylogeny in the
   hierarchical layer that generated those coefficients. Ours is not: neither M0
   nor M3 contains a phylogenetic term.
5. **G5 runs 12 interaction tests.** Its one nominal hit does not survive any
   multiplicity correction.
