#!/bin/bash

ROOT=/local/storage/no-backup/vine-benchmarks/crispr_sims

stats=/local/storage/no-backup/vine-benchmarks/bin/stats

sizes=(
    "10"
    "25"
    "50"
    "100"
    "250"
    "500"
    "1000"
)

printf "ntaxa\tlaml\tstd\tvine\tstd\tbeam\tstd\n" > lnlSummary.txt
for i in "${!sizes[@]}" ; do
    s="${sizes[$i]}"
    printf "%s\t" $s >> lnlSummary.txt
    head -11 $ROOT/${s}taxa/summary.lnl.txt | tail -10 | awk '{print $2}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> lnlSummary.txt
    head -11 $ROOT/${s}taxa/summary.lnl.txt | tail -10 | awk '{print $3}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> lnlSummary.txt
    head -11 $ROOT/${s}taxa/summary.lnl.txt | tail -10 | awk '{print $4}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\n", $2, $4}' >> lnlSummary.txt
done

printf "ntaxa\tlaml\tstd\tvine\tstd\tbeam\tstd\n" > timeSummary.txt
for i in "${!sizes[@]}" ; do
    s="${sizes[$i]}"
    printf "%s\t" $s >> timeSummary.txt
    head -11 $ROOT/${s}taxa/summary.time.txt | tail -10 | awk '{print $2}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> timeSummary.txt
    head -11 $ROOT/${s}taxa/summary.time.txt | tail -10 | awk '{print $3}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> timeSummary.txt
    head -11 $ROOT/${s}taxa/summary.time.txt | tail -10 | awk '{print $4}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\n", $2, $4}' >> timeSummary.txt
done

printf "ntaxa\tlaml\tstd\tbeam\tstd\n" > speedSummary.txt
for i in "${!sizes[@]}" ; do
    s="${sizes[$i]}"
    printf "%s\t" $s >> speedSummary.txt
    head -11 $ROOT/${s}taxa/summary.time.txt | tail -10 | awk '{print $2/$3}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\t", $2, $4}' >> speedSummary.txt
    head -11 $ROOT/${s}taxa/summary.time.txt | tail -10 | awk '{print $4/$3}' | $stats | sed 's/,//g' | awk '{printf "%f\t%f\n", $2, $4}' >> speedSummary.txt
done

