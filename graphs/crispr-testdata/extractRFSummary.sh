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

printf "ntaxa\tvine\tstd\tlaml\tstd\tbeam\tstd\n" > rfSummary.txt
for i in "${!sizes[@]}" ; do
    s="${sizes[$i]}"
    printf "%s\t" $s >> rfSummary.txt
    python rf-mean-and-sd.py $ROOT/${s}taxa/eval.all.rf.txt | awk '{printf "%f\t%f\t%f\t%f\t%f\t%f\n", $3, $4, $5, $6, $7, $8}' >> rfSummary.txt
done

