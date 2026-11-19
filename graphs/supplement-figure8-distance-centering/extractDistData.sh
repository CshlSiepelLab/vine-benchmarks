#!/usr/bin/bash
# Extract per-pair distance data for the centering supplement figure.
#
# For each (ntaxa, replicate, method) this reruns `evalTrees` on the sampled tree file
# (tree.<rep>.<method>.nwk) to recompute the posterior-mean pairwise distances, then
# joins them to the simulation ground-truth (tree.<rep>.true.dist.txt) via
# buildCalibrationTable.py.  Result: data/calibration_data.csv .
#
# The .nwk tree samples themselves are produced upstream by the simulation Makefile
# (see $SIMS_ROOT/Makefile and ../common.mk); this script only re-derives distances.
set -euo pipefail

# ---- lab-specific absolute paths (edit here only) -------------------------------
MAIN_DIR=/local/storage/no-backup/vine-benchmarks
SIMS_ROOT=$MAIN_DIR/dna_sims/hky_300sites
EVALTREES=$MAIN_DIR/bin/vine/bin/evalTrees
# ---------------------------------------------------------------------------------
SIZES="25,50,100"          # taxon counts shown in the figure (panels A,B,C)
NSAMP=10                    # replicates per size
METHODS="var beast mrbayes"
NPROC=24

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW="$SCRIPT_DIR/data/raw"
mkdir -p "$RAW"

# 1. Recompute per-pair posterior-mean distances with evalTrees (parallel).
run_one() {
  local d=$1 rep=$2 m=$3
  local nwk="$SIMS_ROOT/${d}taxa/tree.${rep}.${m}.nwk"
  [ -f "$nwk" ] || return 0
  "$EVALTREES" "$nwk" 2>/dev/null | grep -v '^#' | awk 'NF>=3{print $1,$2,$3}' \
      > "$RAW/${d}taxa.${rep}.${m}.dist"
}
export -f run_one
export SIMS_ROOT EVALTREES RAW

: > "$RAW/.jobs"
IFS=',' read -ra SZ <<< "$SIZES"
for d in "${SZ[@]}"; do
  for rep in $(seq 1 "$NSAMP"); do
    for m in $METHODS; do echo "$d $rep $m" >> "$RAW/.jobs"; done
  done
done
xargs -P "$NPROC" -n3 bash -c 'run_one "$@"' _ < "$RAW/.jobs"

# 2. Join to ground-truth and emit the tidy table.
python3 "$SCRIPT_DIR/buildCalibrationTable.py" \
  --sims-root "$SIMS_ROOT" \
  --raw-dir "$RAW" \
  --sizes "$SIZES" \
  --nsamp "$NSAMP" \
  --methods "$(echo "$METHODS" | tr ' ' ',')" \
  --out "$SCRIPT_DIR/data/calibration_data.csv"

echo "done -> $SCRIPT_DIR/data/calibration_data.csv"
