#!/usr/bin/bash
# Reruns VINE on a single alignment with a sweep of the variance-regularization strength
# (--var-reg / -v) to show that the CENTERING of the posterior mean is essentially
# invariant to it, while interval WIDTH / coverage changes a lot.  (The benchmark itself
# runs with -v 0, which disables the anti-collapse regularizer; coverage is discussed
# separately in the main paper.)
#
# Prints, per setting: cov95 cov50 medSD meanW95 meanErr RMSE slope  (join vs ground-truth).
set -euo pipefail

MAIN_DIR=/local/storage/no-backup/vine-benchmarks
SIMS_ROOT=$MAIN_DIR/dna_sims/hky_300sites
VINE=$MAIN_DIR/bin/vine/bin/vine
EVALTREES=$MAIN_DIR/bin/vine/bin/evalTrees

SIZE=${1:-10}          # taxon count
REP=${2:-1}            # replicate
NSAMP=1000             # posterior samples
SETTINGS=(0 1 2 5)     # -v multipliers to test

FA=$SIMS_ROOT/${SIZE}taxa/tree.${REP}.fa
TRUE=$SIMS_ROOT/${SIZE}taxa/tree.${REP}.true.dist.txt
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

for v in "${SETTINGS[@]}"; do
  "$VINE" "$FA" --hky85 -s "$NSAMP" -v "$v" > "$WORK/v$v.nwk" 2>/dev/null
done

python3 - "$WORK" "$EVALTREES" "$TRUE" "${SETTINGS[*]}" <<'PY'
import sys, subprocess, math, numpy as np
work, evt, true, settings = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split()
tm={}
for ln in open(true):
    s=ln.split()
    if not s or ln.lstrip().startswith('#') or s[0].lower().startswith('leaf'): continue
    tm[(s[0],s[1])]=float(s[2])
print(f"{'-v':>4s} | {'cov95':>6s} {'cov50':>6s} | {'medSD':>7s} {'meanW95':>8s} | {'meanErr':>8s} {'RMSE':>6s} {'slope':>6s}")
for v in settings:
    out=subprocess.run([evt, f"{work}/v{v}.nwk"], capture_output=True, text=True).stdout
    c95=c50=n=0; sd=[]; w=[]; e=[]; tv=[]; ev=[]
    for ln in out.splitlines():
        s=ln.split()
        if len(s)<11 or ln.lstrip().startswith('#') or s[0].lower().startswith('leaf'): continue
        k=(s[0],s[1])
        if k not in tm: continue
        mean,ssd=float(s[2]),float(s[3]); lo95,hi95=float(s[7]),float(s[8]); lo50,hi50=float(s[9]),float(s[10])
        t=tm[k]; n+=1; c95+=lo95<=t<=hi95; c50+=lo50<=t<=hi50
        sd.append(ssd); w.append(hi95-lo95); e.append(mean-t); tv.append(t); ev.append(mean)
    slope=np.polyfit(tv,ev,1)[0]; rmse=math.sqrt(np.mean(np.array(e)**2))
    print(f"{v:>4s} | {c95/n:6.3f} {c50/n:6.3f} | {np.median(sd):7.4f} {np.mean(w):8.4f} | {np.mean(e):+8.4f} {rmse:6.3f} {slope:6.3f}")
PY
