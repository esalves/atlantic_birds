# ATLANTIC ANTS (raw, immutable)

Downloaded 2026-09-09 from the LEEClab GitHub repository
<https://github.com/LEEClab/Atlantic_Ants> (branch `main`):

| File | Source URL | Size |
|---|---|---|
| `Atlantic_Ants_dataset_v1_2021-05-26d.zip` | `https://raw.githubusercontent.com/LEEClab/Atlantic_Ants/main/DATASET/Atlantic_Ants_dataset_v1_2021-05-26d.zip` | 5.3 MB (unzips to `ATLANTIC_ANTS_dataset.txt`, 117 MB, 178,976 records x 57 fields, plus `ATLANTIC_ANTS_references.docx`) |
| `aants_metadata.pdf` | `https://raw.githubusercontent.com/LEEClab/Atlantic_Ants/main/METADATA/aants_metadata.pdf` | 2.6 MB (Ecology data-paper metadata: field definitions) |

The zip is the file the analysis reads. `atlantic_ants_index.R` unzips the
dataset into a temporary directory at run time, so nothing else needs to live
here; if an `extracted/` folder is present it is only a convenience copy of the
zip contents and can be deleted.

Note: the repository `.gitignore` excludes `*.zip`, so the data file is NOT
tracked by git. On a fresh clone re-download it with

```bash
curl -L -o Analysis/data/raw/atlantic_ants/Atlantic_Ants_dataset_v1_2021-05-26d.zip \
  https://raw.githubusercontent.com/LEEClab/Atlantic_Ants/main/DATASET/Atlantic_Ants_dataset_v1_2021-05-26d.zip
```

(5,270,850 bytes; 178,976 records after unzipping).

## Citation (required, CC-BY 4.0)

Silva, R. R., Martello, F., Feitosa, R. M., Silva, O. G. M., Prado, L. P.,
Brandão, C. R. F., ... Ribeiro, M. C. (2022). ATLANTIC ANTS: a data set of ants
in Atlantic Forests of South America. *Ecology*, 103(2), e3580.
<https://doi.org/10.1002/ecy.3580>

Licence: Creative Commons Attribution 4.0 (CC-BY), as stated in the data paper
and the GitHub repository. Any use of derived numbers in the manuscript must cite
Silva et al. (2022).

## File format facts (established 2026-09-09, REVISION_PLAN.md section 7)

- Tab-delimited, Latin-1 encoded, CRLF line endings, unquoted free text
  (`Method.Description`, `Habitat.Description` contain embedded quotes), so read
  with `data.table::fread(sep = "\t", quote = "", encoding = "Latin-1")`.
- `Pitfall.Number`, `Winkler.Number` and `Total.Ant.Abundance` are read as
  character (they contain non-numeric tokens for some records) and must be
  coerced with `as.integer()` / `as.numeric()`.
- `ma.lmt` = 1 flags records inside the Atlantic Forest limit as delimited by
  the data authors; `Exclude` = 1 flags records the authors ask users to drop;
  `fl.nm` is the contributor file each record came from (the unit we use as
  "contributor programme").
- `Measurement.Type` is `abundance` for a minority of records (about 14 % of the
  Atlantic Forest subset); most records are occurrence only.

Nothing in this folder is edited by the analysis. Derived tables go to
`Analysis/data/derived/ants_*.rds`, results to `Analysis/output/ants_*.rds`.
