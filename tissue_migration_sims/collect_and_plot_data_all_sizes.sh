#!/bin/bash

R_SRC="/local/storage/no-backup/vine-benchmarks/r/src"
RPLOTTING_SIF="/local/storage/no-backup/vine-benchmarks/containers/rplotting/rplotting.sif"

# Collect runtime data
outfileTime=eval.all.time.allSizes.txt
rm -f $outfileTime
for file in $(find ./ -maxdepth 2 -name "eval.all.time.txt"); do
    # Get header
    if [ ! -f $outfileTime ]; then
        head -n 1 $file | sed "s/^/size\\t/" > $outfileTime
    fi
    # Get taxa size
    taxaSize=$(basename $(dirname $file))
    # Get content, skippng the last two lines of summary stats
    tail -n +2 $file | head -n -2 | sed "s/^/$taxaSize\\t/" >> $outfileTime
done

# Make the runtime plot
singularity exec --bind $R_SRC:/mnt/scripts --bind ./:/mnt/files $RPLOTTING_SIF Rscript /mnt/scripts/plot_runtime_bars_allSizes.R /mnt/files/eval.all.time.allSizes.txt /mnt/files/eval.all.time.allSizes.pdf



# # Collect lnl data
# outfileLnl=eval.all.lnl.allSizes.txt
# rm -f $outfileLnl
# for file in $(find ./ -maxdepth 2 -name "eval.all.lnl.txt"); do
#     # Get header
#     if [ ! -f $outfileLnl ]; then
#         head -n 1 $file | sed "s/^/size\\t/" > $outfileLnl
#     fi
#     # Get taxa size
#     taxaSize=$(basename $(dirname $file))
#     # Get content, skippng the last two lines of summary stats
#     tail -n +2 $file | head -n -2 | sed "s/^/$taxaSize\\t/" >> $outfileLnl
# done

# # Make the lnl plot
# singularity exec --bind $R_SRC:/mnt/scripts --bind ./:/mnt/files $RPLOTTING_SIF Rscript /mnt/scripts/plot_lnl_allSizes.R /mnt/files/eval.all.lnl.allSizes.txt /mnt/files/eval.all.lnl.allSizes.pdf

# # Make the lnl plot with only vine vs. beam (no laml)
# singularity exec --bind $R_SRC:/mnt/scripts --bind ./:/mnt/files $RPLOTTING_SIF Rscript /mnt/scripts/plot_lnl_allSizes_vineBeamOnly.R /mnt/files/eval.all.lnl.allSizes.txt /mnt/files/eval.all.lnl.allSizes.vineBeamOnly.pdf

