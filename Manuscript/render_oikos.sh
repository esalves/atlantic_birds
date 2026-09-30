#!/usr/bin/env bash
# Render the three Oikos initial-submission files into submission_oikos/.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p submission_oikos
quarto render index.qmd --profile oikos-anon --to docx
cp _manuscript/index.docx submission_oikos/manuscript_anonymised.docx
quarto render index.qmd --profile oikos-title --to docx
cp _manuscript/index.docx submission_oikos/title_page.docx
# supplementary is outside the manuscript project: strip_authors.lua drops its author block
quarto render supplementary.qmd --to docx -M anonymous=true --output supplementary_anonymised.docx
mv supplementary_anonymised.docx submission_oikos/
python3 strip_docx_author.py submission_oikos/manuscript_anonymised.docx submission_oikos/supplementary_anonymised.docx
echo "Wrote submission_oikos/{manuscript_anonymised,title_page,supplementary_anonymised}.docx"
