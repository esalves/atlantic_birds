#!/usr/bin/env bash
# run_bayes_totoro.sh
# Parallel job x tree runner for the Stan tier on Totoro (after Santiago Ortega's
# run_totoro.sh). Resumable: bayes_fit_job.R exits at once if its tree file exists.
#
# Usage: run_bayes_totoro.sh <job-list.txt> [first_tree] [last_tree] [extra flags...]
#   job-list.txt: one job .rds path per line (relative to the repo root)
#   extra flags go to bayes_fit_job.R, e.g. --phylo-slope
# Env: CONC (concurrent fits, default 40; each uses 4 cores -> 160 cores)
#
# Example:
#   nohup Analysis/scripts/run_bayes_totoro.sh Analysis/output/bayes/jobs_priority.txt 1 50 \
#     > ~/logs/bayes_priority.log 2>&1 &
set -u
cd ~/atlantic_birds_rev
LIST=$1; T0=${2:-1}; T1=${3:-50}; shift $(( $# < 3 ? $# : 3 )); EXTRA="$*"
CONC=${CONC:-40}
mkdir -p ~/logs/bayes
run_one() {
  local job=$1 t=$2
  local tag; tag=$(basename "$job" .rds)
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /usr/bin/Rscript Analysis/scripts/bayes_fit_job.R \
    --job "$job" --tree "$t" $EXTRA > ~/logs/bayes/"${tag}${EXTRA// /}_tree$(printf %02d "$t").log" 2>&1
  echo "$tag tree $t exit $?"
}
export -f run_one
export EXTRA
# compile once, before 40 processes race to do it
/usr/bin/Rscript -e 'invisible(cmdstanr::cmdstan_model("Analysis/scripts/stan/phylo_lmm_marginal.stan"))' || exit 1
echo "[start] $(date) list=$LIST trees $T0..$T1 conc=$CONC extra='$EXTRA'"
# trees in the outer loop so every job gets its first trees early
for t in $(seq "$T0" "$T1"); do grep -v '^\s*$' "$LIST" | sed "s/$/ $t/"; done |
  xargs -P "$CONC" -L 1 bash -c 'run_one "$0" "$1"'
echo "[done] $(date)"
