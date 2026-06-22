#!/bin/bash

VINE_NSAMP=1000
VAROPT_ROOT="-M 200 -c 50 --covar CONST -v 0 -s $VINE_NSAMP --hky85 -H "

bestargs=(
    "$VAR_OPTROOT -D 2 "
    "$VAR_OPTROOT -D 3 "
    "$VAR_OPTROOT -D 4 "
    "$VAR_OPTROOT -D 5 "
    "$VAR_OPTROOT -D 6 "
    "$VAR_OPTROOT -D 7 "
    "$VAR_OPTROOT -D 8 "
)

dirs=(
    "archive.D2"
    "archive.D3"
    "archive.D4"
    "archive.D5"
    "archive.D6"
    "archive.D7"
    "archive.D8"
)

printf "run\tvinelnl\tbeastlnl\tvinemf\tbeastmf\tviner2\tvine95\tbeastr2\tbeast95\tvinerf\tbeastrf\ttime\n" > summary.txt

for i in "${!bestargs[@]}"; do
    args="${bestargs[$i]}"
    dir="${dirs[$i]}"
    rm -f tree.*.var.nwk tree.*.var.nwk.log eval.all.*.txt tree.*.var-time tree.*.varlnl tree.*.time tree.*.lnl
    export VAROPT="$VAROPT_ROOT $args"
			
    make -j 10
    
    # extract key statistics for summary table
    printf "$dir\t" >> summary.txt
    tail -1 eval.all.lnl.txt | awk '{printf "%f\t%f\t", $4, $5}' >> summary.txt
    tail -1 eval.all.mf.txt | awk '{printf "%f\t%f\t", $7, $9}' >> summary.txt
    tail -1 eval.all.dist.txt | awk '{printf "%f\t%f\t%f\t%f\t", $3, $5, $7, $9}' >> summary.txt
    tail -1 eval.all.rf.txt | awk '{printf "%f\t%f\t", $7, $9}' >> summary.txt
    tail -1 eval.all.time.txt | awk '{printf "%f\n", $4}' >> summary.txt
			
    # archive files
    mkdir -p $dir
    mv eval.all.*.txt $dir
done
