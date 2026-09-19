# Biology Open submission checklist

Against *Manuscript preparation* (BiO, last updated 22 July 2026). Article type
assumed: **Research Article**.

## Done

| Requirement | Limit | Status |
|---|---|---|
| Title | ≤150 characters | 115 ✓ |
| Running title | ≤32 characters | 31 ✓ ("Atlantic Forest bird morphology") |
| Keywords | 3–6 | 6 ✓ (was 7) |
| Abstract | ≤200 words, no subheadings, no references | 200 ✓ |
| Summary statement | 15–30 words | 23 ✓ |
| Manuscript length | ≤8000 words incl. figure legends, excl. references and Materials & Methods | **6,439** ✓ |
| Display items | ≤8 | 6 figures, 0 main-text tables ✓ |
| Introduction without subheadings | — | ✓ |
| Results broken up by subheadings | — | ✓ (9) |
| Materials and Methods (named, sectioned) | — | ✓ renamed from "Methods" |
| **Use of AI tools** — dedicated M&M section | required | ✓ drafted — **needs author verification** |
| Competing interests | 'No competing interests declared' if none | ✓ exact wording used |
| Data and resource availability | named exactly | ✓ renamed from "Data availability" |
| Acknowledgements | — | ✓ |
| Author contributions | — | ✓ drafted from CRediT roles |
| Species names italicised | — | ✓ |

## Outstanding — author action required

1. **Funding.** BiO requires the official funder name from the Crossref Open Funder
   Registry plus all grant numbers, from *all* authors. If there is none, the
   required wording is exactly: "This research received no specific grant from any
   funding agency in the public, commercial or not-for-profit sectors."
   Currently reads "*To be completed*".
2. **Repository DOI.** The Materials and Methods still cite
   `10.5281/zenodo.atlantic_birds`, which is a placeholder, not a resolvable DOI.
   BiO requires the repository name and a persistent identifier, and asks that the
   dataset also appear in the reference list. Mint the Zenodo DOI and replace both.
3. **Verify the AI-tools statement.** BiO requires this and has a separate AI
   policies page that should be read before submitting. The drafted text describes
   the scope recorded in the repository (31 of 128 commits between 2026-06-12 and
   2026-09-19 are AI co-authored) and discloses the two errors found and corrected.
   The corresponding author must confirm it is complete and accurate.
4. **Confirm competing interests and author contributions with all co-authors.**
5. **Decide on the preprint.** The YAML still carries
   `citation: container-title: bioRxiv`. Remove it if not preprinting.

## Deferred to revision (BiO first submission is format-free)

BiO states that at first submission a manuscript may be in any format, and that a
single PDF containing all text and figures is acceptable. The following are
therefore not blocking now but will be required if the paper is invited back:

- Figures not embedded in the text; **figure legends collected at the end** of the
  manuscript file, first sentence of each in bold.
- 1.5 line spacing, continuous line numbering, page numbers.
- Figures as EPS/PDF/SVG vector art, Arial or Helvetica, RGB, ≤180 × 210 mm,
  12 pt bold uppercase panel letters, 8 pt body labelling.
- Supplementary information collated into a **single PDF** with its own reference
  list; currently `supplementary.docx` / `.html`.
- Reference list alphabetical by first author surname, Harvard (name, date) in
  text — verify the Quarto CSL produces exactly this.
- Submission checklist form (required with the manuscript).

## Notes

- Materials and Methods do not count toward the 8000-word limit, which is why the
  countable total is 6,439 despite a 3,600-word M&M.
- A plain-language summary is present in the YAML. BiO does not ask for one; it is
  harmless but could be removed.
