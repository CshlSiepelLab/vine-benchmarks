#!/usr/bin/bash
# Build data/elbo_max.csv for the Taylor-vs-MonteCarlo ELBO supplement figure.
#
# Source: the Taylor variational runs that were re-run with --batchsize 100 so
# that the final Monte Carlo pass (LNL_mc, reported on the "# Reverting ..." line)
# is a 100-sample estimate of E_q[lnL] at the SAME converged (mu, Sigma) as the
# Taylor/hybrid estimate (LNL). Both share the identical analytic KLD, so
#     ELBO_taylor = LNL      - KLD        (rev "ELB")
#     ELBO_mc     = LNL_mc   - KLD
# i.e. ELBO_taylor - ELBO_mc = LNL - LNL_mc, a same-(mu,Sigma) estimator comparison.
#
# Run from this directory.
set -euo pipefail
SIMDIR=/DATA/local/storage/no-backup/vine-benchmarks/dna_sims/hky_300sites
OUT=data/elbo_max.csv
mkdir -p data
echo "ntaxa,rep,elbo_taylor,elbo_mc,kld" > "$OUT"
for sz in 25 50; do
  for n in $(seq 1 10); do
    log="$SIMDIR/taylor-${sz}taxa/tree.${n}.var.nwk.log"
    line=$(grep '# Reverting' "$log")
    # parse "ELB: X, LNL: Y, ... KLD: K, ... LNL_mc: M"
    elb=$(echo "$line"    | sed -E 's/.*ELB: (-?[0-9.]+).*/\1/')
    kld=$(echo "$line"    | sed -E 's/.*KLD: (-?[0-9.]+).*/\1/')
    lnlmc=$(echo "$line"  | sed -E 's/.*LNL_mc: (-?[0-9.]+).*/\1/')
    elbo_mc=$(python3 -c "print(f'{$lnlmc - $kld:.2f}')")
    echo "${sz},${n},${elb},${elbo_mc},${kld}" >> "$OUT"
  done
done
echo "wrote $OUT"
