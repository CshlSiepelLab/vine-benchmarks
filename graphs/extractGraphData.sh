#!/bin/bash

# Unified data extraction for HKY (300, 10000) and JC69 and GTR graphs.
# Usage: bash extractGraphData.sh [hky300|hky10k|jc69|gtr]

stats="/local/storage/no-backup/vine-benchmarks/bin/stats"

set -e

if [[ $# -lt 1 ]]; then
  echo "Usage: bash extractGraphData.sh [hky300|hky10k|jc69|gtr|treeprior]" >&2
  exit 1
fi

raw_model="$1"

model="$raw_model"

case "$model" in
  hky300|hky10k|jc69|gtr|treeprior) ;;
  *)
    echo "Model must be 'hky300', 'hky10k', 'jc69', 'gtr', or 'treeprior'." >&2
    exit 1
    ;;
esac

# Resolve absolute paths
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
graphs_dir="$script_dir"

if [[ "$model" == "hky300" ]]; then
  ROOT="/local/storage/no-backup/vine-benchmarks/dna_sims/hky_300sites"
  OUT_DIR="$graphs_dir/hky300-data"
elif [[ "$model" == "hky10k" ]]; then
  ROOT="/local/storage/no-backup/vine-benchmarks/dna_sims/hky_10000sites"
  OUT_DIR="$graphs_dir/hky10k-data"
elif [[ "$model" == "gtr" ]]; then
  ROOT="/local/storage/no-backup/vine-benchmarks/dna_sims/gtr_gamma_300sites"
  OUT_DIR="$graphs_dir/gtr-data"
elif [[ "$model" == "treeprior" ]]; then
  ROOT="/local/storage/no-backup/vine-benchmarks/dna_sims/hky_300sites_treeprior"
  OUT_DIR="$graphs_dir/hky300-treeprior-data"
else
  ROOT="/local/storage/no-backup/vine-benchmarks/dna_sims/jc69_300sites"
  OUT_DIR="$graphs_dir/jc69-data"
fi

mkdir -p "$OUT_DIR"

dist_sizes=( "10" "25" "50" "100" "250")
# Size sets
if [[ "$model" == "hky300" ]]; then
  mf_sizes=( "10" "25" "50" "100" "250" "500" "1000")
  lnl_sizes=( "10" "25" "50" "100" "250" "500" "1000")
  time_sizes=( "10" "25" "50" "100" "250" "500" "1000")
  sizes=( "10" "25" "50" "100" "250" "500" "1000")
elif [[ "$model" == "gtr" ]]; then
  lnl_sizes=( "10" "25" "50" "100")
  time_sizes=( "10" "25" "50" "100")
  mf_sizes=( "10" "25" "50" "100")
  sizes=( "10" "25" "50" "100")
elif [[ "$model" == "treeprior" ]]; then
  lnl_sizes=( "10" "25" "50" "100" "250" "500" "1000")
  time_sizes=( "10" "25" "50" "100" "250" "500" "1000")
  mf_sizes=( "10" "25" "50" "100" "250" "500" "1000")
  sizes=( "10" "25" "50" "100" "250" "500" "1000")
else
  sizes=( "10" "25" "50" "100")
  if [[ "$model" == "jc69" ]]; then
    lnl_sizes=( "10" "15" "20" "25" "50" "100" )
    time_sizes=( "10" "15" "20" "25" "50" "100" )
    mf_sizes=( "10" "15" "20" "25" "50" "100" )
  else
    # hky10k
    lnl_sizes=( "10" "25" "50" "100" )
    time_sizes=( "10" "25" "50" "100" )
    mf_sizes=( "10" "25" "50" "100" )
  fi
fi



dims=(
    "archive.D2"
    "archive.D3"
    "archive.D4"
    "archive.D5"
    "archive.D6"
    "archive.D7"
    "archive.D8"    
)

path_lnl() {
  local s="$1"
  echo "$ROOT/${s}taxa/eval.all.lnl.txt"
}

path_lnl_updated() {
  local s="$1"
  echo "$ROOT/${s}taxa/updated.all.lnl.txt"
}

path_time() {
  local s="$1"
  echo "$ROOT/${s}taxa/eval.all.time.txt"
}

path_mf() {
  local s="$1"
  echo "$ROOT/${s}taxa/eval.all.mf.txt"
}

path_mf_alt() {
  local s="$1"
  echo "$ROOT/${s}taxa/eval.all.mf.txt"
}

path_dist() {
  local s="$1"
  echo "$ROOT/${s}taxa/eval.all.dist.txt"
}

path_entropy() {
  local s="$1"
  echo "$ROOT/${s}taxa/eval.all.ent.txt"
}
path_entropy_flows() {
  local s="$1"
  echo "$ROOT/${s}taxa/vineAlt/eval.all.ent.txt"
}

path_dist_flows() {
  local s="$1"
  echo "$ROOT/${s}taxa/vineAlt/eval.all.dist.txt"
}

append_lnl_stats() {
  local file="$1"
  local method="$2"
  local values
  values="$(
    awk -v method="$method" '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          name = tolower($i)
          if (name == "true") true_col = i
          if (name == tolower(method)) method_col = i
        }
        next
      }
      NR <= 11 && method_col && method_col <= NF {
        print $method_col - $true_col
      }
    ' "$file"
  )"
  if [[ -n "$values" ]]; then
    printf "%s\n" "$values" | "$stats" | sed 's/,//g' | \
      awk '{printf "%f\t%f", $2, $4}'
  else
    printf "0\t0"
  fi
}

# ---------------- lnlSummary.txt ----------------
if [[ "$model" == "jc69" ]]; then
  printf "%s\n" \
"ntaxa	ave	NJ	std	vine	std	beast	std	beast-beagle	std	\
mrbayes	std	mrbayes-beagle	std	\
dodonaphy	std	geophy	std	vaiphy	std" \
  | expand -t 1 > "$OUT_DIR/lnlSummary.txt"
elif [[ "$model" == "hky300" || "$model" == "hky10k" || "$model" == "treeprior" ]]; then
  printf "%s\n" \
"ntaxa	ave	NJ	std	vine	std	beast	std	beast-beagle	std	\
mrbayes	std	mrbayes-beagle	std" \
  | expand -t 1 > "$OUT_DIR/lnlSummary.txt"
else
  printf "%s\n" \
"ntaxa	ave	NJ	std	vine	std	beast	std	mrbayes	std" \
  > "$OUT_DIR/lnlSummary.txt"
fi

for s in "${lnl_sizes[@]}"; do
  printf "%s\t" "$s" >> "$OUT_DIR/lnlSummary.txt"
  head -11 "$(path_lnl "$s")" | tail -10 | awk '{print $1}' | \
    $stats | sed 's/,//g' | awk '{printf "%f\t", $2}' \
    >> "$OUT_DIR/lnlSummary.txt"

  if [[ "$model" == "jc69" ]]; then
    methods=(nj vine beast beast-beagle mrbayes mrbayes-beagle dodonaphy geophy vaiphy)
  elif [[ "$model" == "hky300" || "$model" == "hky10k" || "$model" == "treeprior" ]]; then
    methods=(nj vine beast beast-beagle mrbayes mrbayes-beagle)
  else
    methods=(nj vine beast mrbayes)
  fi
  for i in "${!methods[@]}"; do
    append_lnl_stats "$(path_lnl "$s")" "${methods[$i]}" \
      >> "$OUT_DIR/lnlSummary.txt"
    if [[ "$i" -lt $((${#methods[@]} - 1)) ]]; then
      printf "\t" >> "$OUT_DIR/lnlSummary.txt"
    else
      printf "\n" >> "$OUT_DIR/lnlSummary.txt"
    fi
  done
done

# ---------------- timeSummary.txt ----------------
if [[ "$model" == "jc69" ]]; then
  printf "%s\n" \
"ntaxa	vine	std	beast	std	beast-beagle	std	mrbayes	std	\
mrbayes-beagle	std	dodonaphy	std	geophy	std	vaiphy	std" \
  | expand -t 1 > "$OUT_DIR/timeSummary.txt"
elif [[ "$model" == "hky300" || "$model" == "hky10k" || "$model" == "treeprior" ]]; then
  printf "%s\n" \
"ntaxa	vine	std	beast	std	beast-beagle	std	mrbayes	std	\
mrbayes-beagle	std" \
  | expand -t 1 > "$OUT_DIR/timeSummary.txt"
else
  printf "%s\n" \
"ntaxa	vine	std	beast	std	mrbayes	std	raxml	std" \
  > "$OUT_DIR/timeSummary.txt"
fi

for s in "${time_sizes[@]}"; do
  printf "%s\t" "$s" >> "$OUT_DIR/timeSummary.txt"
  if [[ "$model" == "jc69" ]]; then
    if [[ "$s" == "25" || "$s" == "50" || "$s" == "100" ]]; then
      for col in 4 2; do
        head -11 "$(path_time "$s")" | tail -10 | \
          awk -v c="$col" '{print $c}' | $stats | sed 's/,//g' | \
          awk '{printf "%f\t%f\t", $2, $4}' \
          >> "$OUT_DIR/timeSummary.txt"
      done
      printf "0\t0\t" >> "$OUT_DIR/timeSummary.txt"
      head -11 "$(path_time "$s")" | tail -10 | \
        awk '{print $3}' | $stats | sed 's/,//g' | \
        awk '{printf "%f\t%f\t", $2, $4}' \
        >> "$OUT_DIR/timeSummary.txt"
      printf "0\t0\t0\t0\t0\t0\t0\t0\n" \
        >> "$OUT_DIR/timeSummary.txt"
    else
      for col in 6 2 3 4 5 7 8 9; do
        head -11 "$(path_time "$s")" | tail -10 | \
          awk -v c="$col" '{print $c}' | $stats | sed 's/,//g' | \
          awk '{printf "%f\t%f", $2, $4}' \
          >> "$OUT_DIR/timeSummary.txt"
        if [[ "$col" != "9" ]]; then
          printf "\t" >> "$OUT_DIR/timeSummary.txt"
        else
          printf "\n" >> "$OUT_DIR/timeSummary.txt"
        fi
      done
    fi
  elif [[ "$model" == "hky300" || "$model" == "hky10k" || "$model" == "treeprior" ]]; then
    for col in 6 2 3 4 5; do
      head -11 "$(path_time "$s")" | tail -10 | \
        awk -v c="$col" '{print $c}' | $stats | sed 's/,//g' | \
        awk '{printf "%f\t%f", $2, $4}' \
        >> "$OUT_DIR/timeSummary.txt"
      if [[ "$col" != "5" ]]; then
        printf "\t" >> "$OUT_DIR/timeSummary.txt"
      else
        printf "\n" >> "$OUT_DIR/timeSummary.txt"
      fi
    done
  else
    for col in 4 2 3 5; do
      head -11 "$(path_time "$s")" | tail -10 | \
        awk -v c="$col" '{print $c}' | $stats | sed 's/,//g' | \
        awk '{printf "%f\t%f", $2, $4}' \
        >> "$OUT_DIR/timeSummary.txt"
      if [[ "$col" -lt 5 ]]; then
        printf "\t" >> "$OUT_DIR/timeSummary.txt"
      else
        printf "\n" >> "$OUT_DIR/timeSummary.txt"
      fi
    done
  fi
done

# ---------------- speedSummary.txt ----------------
if [[ "$model" == "jc69" || "$model" == "hky300" ||
      "$model" == "hky10k" || "$model" == "treeprior" ]]; then
  printf "ntaxa\tbeast\tstd\tbeast-beagle\tstd\tmrbayes\tstd\tmrbayes-beagle\tstd\n" \
    > "$OUT_DIR/speedSummary.txt"
  for s in "${sizes[@]}"; do
    printf "%s\t" "$s" >> "$OUT_DIR/speedSummary.txt"
    if [[ "$model" == "jc69" &&
          ( "$s" == "25" || "$s" == "50" || "$s" == "100" ) ]]; then
      head -11 "$(path_time "$s")" | tail -10 | \
        awk '{print $2/$4}' | $stats | sed 's/,//g' | \
        awk '{printf "%f\t%f\t", $2, $4}' \
        >> "$OUT_DIR/speedSummary.txt"
      printf "0\t0\t" >> "$OUT_DIR/speedSummary.txt"
      head -11 "$(path_time "$s")" | tail -10 | \
        awk '{print $3/$4}' | $stats | sed 's/,//g' | \
        awk '{printf "%f\t%f\t", $2, $4}' \
        >> "$OUT_DIR/speedSummary.txt"
      printf "0\t0\n" >> "$OUT_DIR/speedSummary.txt"
    else
      for col in 2 3 4 5; do
        head -11 "$(path_time "$s")" | tail -10 | \
          awk -v c="$col" '{print $c/$6}' | \
          $stats | sed 's/,//g' | \
          awk '{printf "%f\t%f", $2, $4}' \
          >> "$OUT_DIR/speedSummary.txt"
        if [[ "$col" -lt 5 ]]; then
          printf "\t" >> "$OUT_DIR/speedSummary.txt"
        else
          printf "\n" >> "$OUT_DIR/speedSummary.txt"
        fi
      done
    fi
  done
else
  printf "%s\n" \
"ntaxa	beast	std	mrbayes	std	raxml	std" \
  > "$OUT_DIR/speedSummary.txt"
  for s in "${sizes[@]}"; do
    printf "%s\t" "$s" >> "$OUT_DIR/speedSummary.txt"
    head -11 "$(path_time "$s")" | tail -10 | \
      awk '{print $2/$4}' | $stats | sed 's/,//g' | \
      awk '{printf "%f\t%f\t", $2, $4}' \
      >> "$OUT_DIR/speedSummary.txt"
    head -11 "$(path_time "$s")" | tail -10 | \
      awk '{print $3/$4}' | $stats | sed 's/,//g' | \
      awk '{printf "%f\t%f\t", $2, $4}' \
      >> "$OUT_DIR/speedSummary.txt"
    head -11 "$(path_time "$s")" | tail -10 | \
      awk '{print $5/$4}' | $stats | sed 's/,//g' | \
      awk '{printf "%f\t%f\n", $2, $4}' \
      >> "$OUT_DIR/speedSummary.txt"
  done
fi

# ---------------- mfSummary.txt ----------------
if [[ "$model" == "hky300" || "$model" == "hky10k" || "$model" == "treeprior" ]]; then
  printf "%s\n" \
"ntaxa	true	std	NJ	std	vine	std	beast	std	beast-beagle	std	\
mrbayes	std	mrbayes-beagle	std" \
    | expand -t 1 > "$OUT_DIR/mfSummary.txt"
else
  printf "%s\n" \
"ntaxa	true	std	NJ	std	vine	std	beast	std	mrbayes	std" \
    > "$OUT_DIR/mfSummary.txt"
fi
for s in "${mf_sizes[@]}"; do
  printf "%s\t" "$s" >> "$OUT_DIR/mfSummary.txt"
  if [[ "$model" == "hky300" || "$model" == "hky10k" || "$model" == "treeprior" ]]; then
    tail -1 "$(path_mf "$s")" | \
      awk '{printf "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", \
$1,$2,$3,$4,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' \
      >> "$OUT_DIR/mfSummary.txt"
  else
    tail -1 "$(path_mf "$s")" | \
      awk '{printf "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", \
$1,$2,$3,$4,$7,$8,$9,$10,$11,$12}' >> "$OUT_DIR/mfSummary.txt"
  fi
done

# HKY-10k has no vineAlt/250-taxon distance and entropy inputs. Its figure
# panels use the four summaries above, so extraction is complete here.
if [[ "$model" == "hky10k" ]]; then
  exit 0
fi

# ---------------- distSummary.txt ----------------
printf "ntaxa\tvine\tbeast\tvineflows\tvinedev\tbeastdev\tvineflowsdev\n" > "$OUT_DIR/distSummary.txt"
for s in "${dist_sizes[@]}"; do
  vine_mean="$(tail -1 "$(path_dist "$s")" | awk '{printf "%f", $5}')"
  beast_mean="$(tail -1 "$(path_dist "$s")" | awk '{printf "%f", $9}')"
  vineflows_mean="$(
    tail -1 "$(path_dist_flows "$s")" | awk '{printf "%f", $5}'
  )"

  read -r vine_sd beast_sd <<< "$(
    head -n 11 "$(path_dist "$s")" | awk '
      NR > 1 {
        sum5 += $5
        sumsq5 += $5 * $5
        sum9 += $9
        sumsq9 += $9 * $9
        n++
      }
      END {
        if (n == 0) {
          printf "0 0"
        } else {
          mean5 = sum5 / n
          mean9 = sum9 / n
          sd5 = sqrt((sumsq5 / n) - (mean5 * mean5))
          sd9 = sqrt((sumsq9 / n) - (mean9 * mean9))
          printf "%.6f %.6f", sd5, sd9
        }
      }'
  )"

  vineflows_sd="$(
    head -n 11 "$(path_dist_flows "$s")" | awk '
      NR > 1 {
        sum5 += $5
        sumsq5 += $5 * $5
        n++
      }
      END {
        if (n == 0) {
          printf "0"
        } else {
          mean5 = sum5 / n
          sd5 = sqrt((sumsq5 / n) - (mean5 * mean5))
          printf "%.6f", sd5
        }
      }'
  )"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$s" "$vine_mean" "$beast_mean" "$vineflows_mean" \
    "$vine_sd" "$beast_sd" "$vineflows_sd" \
    >> "$OUT_DIR/distSummary.txt"
done

printf "ntaxa\tvine\tbeast\tvineflows\tvinedev\tbeastdev\tvineflowsdev\n" > "$OUT_DIR/entropySummary.txt"
for s in "${dist_sizes[@]}"; do
  vine_mean="$(tail -1 "$(path_entropy "$s")" | awk '{printf "%f", $2}')"
  beast_mean="$(tail -1 "$(path_entropy "$s")" | awk '{printf "%f", $5}')"
  vineflows_mean="$(
    tail -1 "$(path_entropy_flows "$s")" | awk '{printf "%f", $2}'
  )"

  read -r vine_sd beast_sd <<< "$(
    head -n 11 "$(path_entropy "$s")" | awk '
      NR > 1 {
        sum2 += $2
        sumsq2 += $2 * $2
        sum5 += $5
        sumsq5 += $5 * $5
        n++
      }
      END {
        if (n == 0) {
          printf "0 0"
        } else {
          mean2 = sum2 / n
          mean5 = sum5 / n
          sd2 = sqrt((sumsq2 / n) - (mean2 * mean2))
          sd5 = sqrt((sumsq5 / n) - (mean5 * mean5))
          printf "%.6f %.6f", sd2, sd5
        }
      }'
  )"

  vineflows_sd="$(
    head -n 11 "$(path_entropy_flows "$s")" | awk '
      NR > 1 {
        sum2 += $2
        sumsq2 += $2 * $2
        n++
      }
      END {
        if (n == 0) {
          printf "0"
        } else {
          mean2 = sum2 / n
          sd2 = sqrt((sumsq2 / n) - (mean2 * mean2))
          printf "%.6f", sd2
        }
      }'
  )"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$s" "$vine_mean" "$beast_mean" "$vineflows_mean" \
    "$vine_sd" "$beast_sd" "$vineflows_sd" \
    >> "$OUT_DIR/entropySummary.txt"
done

if [[ "$model" == "hky300" || "$model" == "treeprior" ]]; then
  printf "ntaxa\tNJ\tstd\tvine\tstd\tbeast\tstd\tbeast-beagle\tstd\tmrbayes\tstd\tmrbayes-beagle\tstd\n" > "$OUT_DIR/rfSummary.txt"
  for i in "${!sizes[@]}" ; do
          s="${sizes[$i]}"
          printf "%s\t" $s >> "$OUT_DIR/rfSummary.txt"
          python "$script_dir/rf-mean-and-sd.py" "$ROOT/${s}taxa/eval.all.rf.txt" | awk '{printf "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' >> "$OUT_DIR/rfSummary.txt"
      done
fi

# Branch-score distance (point / posterior-mean-tree BSD, normalized by true-tree
# length).  eval.all.bsd.txt already carries the grand row (per-method mean, SD
# across replicates) after the "-----" line, so just extract that row.  Columns
# selected mirror rfSummary: NJ, vine, beast, beast-beagle, mrbayes, mrbayes-beagle.
if [[ "$model" == "hky300" || "$model" == "treeprior" ]]; then
  printf "ntaxa\tNJ\tstd\tvine\tstd\tbeast\tstd\tbeast-beagle\tstd\tmrbayes\tstd\tmrbayes-beagle\tstd\n" > "$OUT_DIR/bsdSummary.txt"
  for i in "${!sizes[@]}" ; do
          s="${sizes[$i]}"
          printf "%s\t" $s >> "$OUT_DIR/bsdSummary.txt"
          grep -A1 '^-----' $ROOT/${s}taxa/eval.all.bsd.txt | tail -1 | awk '{printf "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' >> "$OUT_DIR/bsdSummary.txt"
      done
fi

if [[ "$model" == "hky300" ]]; then
  printf "d\tvine\tstd\tbeast\ttime\tstd\tvineH\tstd\ttimeH\tstd\n" > "$OUT_DIR/dimSummary25.txt"
  for i in "${!dims[@]}" ; do
      SUB="${dims[$i]}"
      printf "%s\t" $SUB >> "$OUT_DIR/dimSummary25.txt"
      head -11 $ROOT/testdata.dim25/$SUB/eval.all.lnl.txt | tail -10 | awk '{print $4-$5}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> "$OUT_DIR/dimSummary25.txt"
      head -11 $ROOT/testdata.dim25/$SUB/eval.all.lnl.txt | tail -10 | awk '{print "0"}' | $stats | sed 's/,//g' | awk '{printf "%f\t", $2}' >> "$OUT_DIR/dimSummary25.txt"    
      head -11 $ROOT/testdata.dim25/$SUB/eval.all.time.txt | tail -10 | awk '{print $4}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> "$OUT_DIR/dimSummary25.txt"
      head -11 $ROOT/testdata.dim25H/$SUB/eval.all.lnl.txt | tail -10 | awk '{print $4-$5}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> "$OUT_DIR/dimSummary25.txt"
      head -11 $ROOT/testdata.dim25H/$SUB/eval.all.time.txt | tail -10 | awk '{print $4}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\n", $2, $4}' >> "$OUT_DIR/dimSummary25.txt"
  done

  printf "d\tvine\tstd\tbeast\ttime\tstd\tvineH\tstd\ttimeH\tstd\n" > "$OUT_DIR/dimSummary50.txt"
  for i in "${!dims[@]}" ; do
      SUB="${dims[$i]}"
      printf "%s\t" $SUB >> "$OUT_DIR/dimSummary50.txt"
      head -11 $ROOT/testdata.dim50/$SUB/eval.all.lnl.txt | tail -10 | awk '{print $4-$5}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> "$OUT_DIR/dimSummary50.txt"
      head -11 $ROOT/testdata.dim50/$SUB/eval.all.lnl.txt | tail -10 | awk '{print "0"}' | $stats | sed 's/,//g' | awk '{printf "%f\t", $2}' >> "$OUT_DIR/dimSummary50.txt"    
      head -11 $ROOT/testdata.dim50/$SUB/eval.all.time.txt | tail -10 | awk '{print $4}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> "$OUT_DIR/dimSummary50.txt"
      head -11 $ROOT/testdata.dim50H/$SUB/eval.all.lnl.txt | tail -10 | awk '{print $4-$5}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> "$OUT_DIR/dimSummary50.txt"
      head -11 $ROOT/testdata.dim50H/$SUB/eval.all.time.txt | tail -10 | awk '{print $4}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\n", $2, $4}' >> "$OUT_DIR/dimSummary50.txt"
  done

fi

echo "Wrote summaries to: $OUT_DIR"
