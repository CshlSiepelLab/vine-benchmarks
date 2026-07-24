#!/bin/bash

ROOT=/local/storage/no-backup/vine-benchmarks/crispr_sims

sizes=(
    "10"
    "25"
    "50"
    "100"
    "250"
    "500"
    "1000"
)

for metric in rf bsd; do
    output="${metric}Summary.txt"
    printf "ntaxa\tvine\tstd\tlaml\tstd\tbeam\tstd\n" > "$output"
    for s in "${sizes[@]}"; do
        printf "%s\t" "$s" >> "$output"
        # Aggregate row layout: true mean/sd, vine mean/sd, laml mean/sd,
        # beam mean/sd. Select the three inferred methods.
        awk '/^-----/{getline; printf "%f\t%f\t%f\t%f\t%f\t%f\n", $3, $4, $5, $6, $7, $8}' \
            "$ROOT/${s}taxa/eval.all.${metric}.txt" >> "$output"
    done
done
