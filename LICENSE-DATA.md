# Licence for the data, results and text

The derived data sets, fitted result objects, tables, figures and manuscript
text in this repository are licensed under the
**Creative Commons Attribution 4.0 International licence (CC BY 4.0)**.

> <https://creativecommons.org/licenses/by/4.0/>
> Legal code: <https://creativecommons.org/licenses/by/4.0/legalcode>

You are free to share and adapt this material for any purpose, including
commercially, provided you give appropriate credit, link to the licence, and
indicate if changes were made.

Covered by this licence:

- `Analysis/data/derived/` — the analytical data sets
- `Analysis/output/` — fitted result objects, tables and `figure_data/`
- `Analysis/figures/` and `Manuscript/images/`
- the manuscript and supplementary text

The code is licensed separately under the MIT licence; see `LICENSE`.

## Please cite

Santos, E. S. A., A. Mizuno, S. Ortega and G. Burin. Morphological stability
and subtle phenotypic shifts in Atlantic Forest birds across two decades of
climate warming.
<!-- TODO(submission): replace with the published/preprint citation and DOI. -->

## Third-party sources

We can license the derived material as above because the primary data sources
place no restrictions on reuse, or place only an attribution requirement that
CC BY 4.0 carries through. Each source must still be cited in its own right.

| Source | Used for | Terms |
| :--- | :--- | :--- |
| **ATLANTIC BIRD TRAITS** — Rodrigues et al. 2019, *Ecology* 100(6):e02647, [10.1002/ecy.2647](https://doi.org/10.1002/ecy.2647) | all morphological records | "No copyright or proprietary restrictions are associated with the use of this data set. Please cite this data paper when the data are used in publications or teaching and educational activities." |
| **EltonTraits 1.0** — Wilman et al. 2014, *Ecology* 95(7):2027, [10.1890/13-1917.1](https://doi.org/10.1890/13-1917.1) | dietary reliance on invertebrates | Copyright restrictions: none. Proprietary restrictions: none. Cite the data paper. |
| **AvesData / a complete and dynamic tree of birds** — McTavish et al. 2025, *PNAS* 122(18):e2409658122, [10.1073/pnas.2409658122](https://doi.org/10.1073/pnas.2409658122) | the 50 phylogenetic topologies | **CC BY 4.0** — attribution required, and carried through by the licence above. |
| **WorldClim 2.1 / CRU-TS 4.09** — Fick & Hijmans 2017; Harris et al. 2020 | temperature, precipitation and SPEI at capture localities | "The data are freely available for academic use and other non-commercial use… **Redistribution or commercial use is not allowed without prior permission.**" See the carve-out below. |

### WorldClim carve-out

WorldClim's terms forbid redistribution without prior permission, so **no
file containing per-record or per-grid-cell WorldClim values is distributed
here**. Specifically, `Analysis/data/derived/passer90_climate.rds` is excluded
from this repository and from the archived deposit.

To rebuild it, download WorldClim 2.1 historical monthly weather data yourself
(<https://www.worldclim.org/data/monthlywth.html>, approximately 13 GB) and run
`Analysis/scripts/climate_extraction.R`. This is only needed to re-fit the
temperature-anomaly models; every number in the main text reproduces without
it.

`Analysis/output/climate_trends.rds` **is** included. It holds fitted
regression coefficients summarising the climate trend across 357 localities —
a research result rather than a redistribution of the underlying grids, which
WorldClim's terms expressly contemplate publishing.
