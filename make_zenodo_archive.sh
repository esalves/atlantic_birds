#!/usr/bin/env bash
#
# Assemble the curated Zenodo deposit.
#
# Builds a staging tree from files tracked in git HEAD, so whatever you get is
# exactly what is committed — never a stray local file. Writes a MANIFEST.txt
# with a SHA-256 for every file, then a single .zip.
#
#   ./make_zenodo_archive.sh              # -> dist/atlantic_birds_<short-sha>.zip
#   ./make_zenodo_archive.sh v1.0.0       # -> dist/atlantic_birds_v1.0.0.zip
#
# WHAT IS EXCLUDED, AND WHY
#   Analysis/data/raw/worldclim, terraclimate, AvesDataLite*
#                              third-party bulk (~14 GB); download scripts ship
#   Analysis/data/derived/passer90_climate.rds
#                              holds per-record WorldClim values; WorldClim
#                              forbids redistribution without prior permission.
#                              Rebuild with Analysis/scripts/climate_extraction.R.
#   Analysis/output/models/    ~2 GB of brms objects superseded by the glmmTMB
#                              engine; no reported number depends on them
#   archive/                   superseded and retired analyses; preserved in the
#                              git history of the public repository instead
#   jirinec_gap_analysis.R, jirinec_gaps_report.R, Analysis/output/jirinec_gaps/
#                              an exploratory side-report, never part of the
#                              manuscript; shipping it would invite readers to
#                              treat it as a result of this paper
#   Analysis/scripts/atlantic_birds_ms*            legacy Rmd + 440 MB knitr cache
#   Analysis/scripts/body_mass_descriptive.html    legacy 4.2 MB render
#   REVISION_*.md, REVIEW_RESPONSE_PLAN.md, CODE_REVIEW.md
#                              internal planning and peer-review correspondence
#
# Everything needed to reproduce every number in the main text is included.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"

VERSION="${1:-$(git rev-parse --short HEAD)}"
STAGE_NAME="atlantic_birds_${VERSION}"
DIST="$REPO/dist"
STAGE="$DIST/$STAGE_NAME"

if ! git diff-index --quiet HEAD -- 2>/dev/null; then
  echo "WARNING: the working tree has uncommitted changes." >&2
  echo "         The deposit is built from git HEAD, so those will NOT be included." >&2
fi

rm -rf "$STAGE"
mkdir -p "$STAGE"

echo "Assembling $STAGE_NAME from $(git rev-parse --short HEAD) ..."

# Pull the curated paths straight out of HEAD.
INCLUDE=(
  'Analysis/scripts'
  'Analysis/data/derived'
  'Analysis/data/raw/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv'
  'Analysis/data/raw/BirdFuncDat.txt'
  'Analysis/output'
  'Analysis/figures'
  'Manuscript'
  'README.md'
  'REPO_STRUCTURE.md'
  'LICENSE'
  'LICENSE-DATA.md'
  'CITATION.cff'
  '.zenodo.json'
  'make_zenodo_archive.sh'
)

git archive --format=tar HEAD -- "${INCLUDE[@]}" | tar -x -C "$STAGE"

# Drop what git tracks but the deposit should not carry.
rm -rf \
  "$STAGE/Analysis/output/models" \
  "$STAGE/Analysis/output/jirinec_gaps" \
  "$STAGE/Analysis/scripts/atlantic_birds_ms_cache" \
  "$STAGE/Analysis/scripts/atlantic_birds_ms_files" \
  "$STAGE/Manuscript/_freeze" \
  "$STAGE/Manuscript/.quarto"
rm -f \
  "$STAGE/Analysis/data/derived/passer90_climate.rds" \
  "$STAGE/Analysis/scripts/jirinec_gap_analysis.R" \
  "$STAGE/Analysis/scripts/jirinec_gaps_report.R" \
  "$STAGE/Analysis/scripts/atlantic_birds_ms.Rmd" \
  "$STAGE/Analysis/scripts/body_mass_descriptive.html" \
  "$STAGE/Analysis/scripts/body_mass_descriptive.qmd"
find "$STAGE" -name '.DS_Store' -delete
find "$STAGE" -name 'REVISION_*.md' -delete

# Refuse to ship anything WorldClim-derived per-record.
if find "$STAGE" -name 'passer90_climate.rds' | grep -q .; then
  echo "ERROR: a WorldClim-derived per-record file is still staged. Aborting." >&2
  exit 1
fi

# Refuse to ship the exploratory Jirinec side-report.
if find "$STAGE" -iname '*jirinec*' | grep -q .; then
  echo "ERROR: Jirinec gap-analysis material is still staged. Aborting." >&2
  find "$STAGE" -iname '*jirinec*' >&2
  exit 1
fi

# Every input the manuscript reads must be present, or the deposit cannot render.
MISSING=0
while read -r need; do
  [ -e "$STAGE/$need" ] || { echo "ERROR: missing required input: $need" >&2; MISSING=1; }
done <<'EOF'
Analysis/data/derived/passer90.rda
Analysis/output/descriptive_summary.rds
Analysis/output/controlled_wing_phylo_results.rds
Analysis/output/mass_results.rds
Analysis/output/multitrait_phylo_results.rds
Analysis/output/multitrait_results.rds
Analysis/output/bivariate_phylo_results.rds
Analysis/output/variance_phylo_results.rds
Analysis/output/diet_interaction_phylo.rds
Analysis/output/climate_trends.rds
Analysis/output/effect_scale.rds
Analysis/output/referee_reruns/referee_reruns.rds
Analysis/output/figure_data/fig_isometry_summary.json
Manuscript/index.qmd
Manuscript/supplementary.qmd
Manuscript/references.bib
EOF
[ "$MISSING" -eq 0 ] || { echo "Aborting: the deposit would not render." >&2; exit 1; }

# Manifest with checksums.
{
  echo "# $STAGE_NAME"
  echo "# git commit: $(git rev-parse HEAD)"
  echo "# built:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "#"
  echo "# sha256  size(bytes)  path"
  cd "$STAGE"
  find . -type f ! -name MANIFEST.txt | LC_ALL=C sort | while read -r f; do
    printf '%s  %s  %s\n' \
      "$(shasum -a 256 "$f" | cut -d' ' -f1)" \
      "$(wc -c < "$f" | tr -d ' ')" \
      "${f#./}"
  done
} > "$STAGE/MANIFEST.txt"

( cd "$DIST" && rm -f "$STAGE_NAME.zip" && zip -qr "$STAGE_NAME.zip" "$STAGE_NAME" )

echo
echo "  files:   $(grep -cv '^#' "$STAGE/MANIFEST.txt")"
echo "  staged:  $(du -sh "$STAGE" | cut -f1)"
echo "  zip:     $DIST/$STAGE_NAME.zip ($(du -h "$DIST/$STAGE_NAME.zip" | cut -f1))"
echo
echo "Next: reserve the Zenodo DOI, write it into Manuscript/index.qmd"
echo "      (Data and resource availability + Reproducible software environment),"
echo "      re-render, rebuild this archive, then publish."
