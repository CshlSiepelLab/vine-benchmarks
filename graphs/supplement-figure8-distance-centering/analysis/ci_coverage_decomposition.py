#!/usr/bin/env python3
"""Decompose 95% CI under-coverage into CENTERING (biased mean) vs WIDTH (too-narrow)
components.  z=(mean-true)/sd : mean(z)~0 => centered ; sd(z)>1 => intervals too narrow.
Reruns evalTrees to obtain the SD/CI columns (not stored in the calibration table).

NOTE: coverage/interval-width is discussed separately in the main paper (with -v settings
other than 0).  This script is provided for completeness; the supplement figure itself is
about centering only.
"""
import os, sys, math, argparse, subprocess
import numpy as np, pandas as pd

METHODS=[("var","vine"),("beast","beast"),("mrbayes","mrbayes")]

def true_map(path):
    m={}
    for ln in open(path):
        s=ln.split()
        if not s or ln.lstrip().startswith("#") or s[0].lower().startswith("leaf"): continue
        m[(s[0],s[1])]=float(s[2])
    return m

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sims-root", default="/local/storage/no-backup/vine-benchmarks/dna_sims/hky_300sites")
    ap.add_argument("--evaltrees", default="/local/storage/no-backup/vine-benchmarks/bin/vine/bin/evalTrees")
    ap.add_argument("--sizes", default="10,25,50")
    ap.add_argument("--nsamp", type=int, default=10)
    a=ap.parse_args()
    sizes=[int(s) for s in a.sizes.split(",")]

    rows=[]
    for nt in sizes:
        d=f"{nt}taxa"
        for rep in range(1,a.nsamp+1):
            tf=os.path.join(a.sims_root,d,f"tree.{rep}.true.dist.txt")
            if not os.path.isfile(tf): continue
            tm=true_map(tf)
            for mf,mn in METHODS:
                nwk=os.path.join(a.sims_root,d,f"tree.{rep}.{mf}.nwk")
                if not os.path.isfile(nwk): continue
                out=subprocess.run([a.evaltrees,nwk],capture_output=True,text=True).stdout
                for ln in out.splitlines():
                    s=ln.split()
                    if len(s)<11 or ln.lstrip().startswith("#") or s[0].lower().startswith("leaf"): continue
                    k=(s[0],s[1])
                    if k not in tm: continue
                    mean,sd=float(s[2]),float(s[3]); lo95,hi95=float(s[7]),float(s[8])
                    rows.append((nt,mn,tm[k],mean,sd,lo95,hi95))
            print(f"  done {d} rep{rep}", file=sys.stderr)
    df=pd.DataFrame(rows,columns=["ntaxa","method","true","mean","sd","lo95","hi95"])
    df["z"]=np.where(df.sd>0,(df["mean"]-df["true"])/df.sd,np.nan)
    df["in95"]=(df["true"]>=df.lo95)&(df["true"]<=df.hi95)
    df["miss_high"]=df["true"]>df.hi95

    print("="*86)
    print("95% CI coverage decomposition")
    print("  mean(z)~0 => centered ; sd(z)>1 => too narrow ; miss_hi% => truth above upper CI")
    print("="*86)
    h=f"{'ntaxa':>6s} {'method':>8s} | {'cov95':>6s} | {'mean(z)':>8s} {'sd(z)':>7s} | {'miss_hi%':>8s} {'medSD':>7s}"
    print(h); print("-"*len(h))
    for nt in sizes:
        for m in ["vine","beast","mrbayes"]:
            s=df[(df.ntaxa==nt)&(df.method==m)]
            if s.empty: continue
            nmiss=(~s.in95).sum()
            mh=100*s.miss_high.sum()/nmiss if nmiss else float("nan")
            print(f"{nt:>6d} {m:>8s} | {s.in95.mean():6.3f} | {s.z.mean():+8.2f} {s.z.std():7.2f} | {mh:7.1f}% {s.sd.median():7.4f}")
        print()

if __name__ == "__main__":
    main()
