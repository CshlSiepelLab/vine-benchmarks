#!/usr/bin/env python3
"""Branch-level centering: match branches by bipartition between each method's tree
sample and the true tree, and compare posterior-mean branch length to the truth,
split terminal vs internal.  Shows VINE's shrinkage lives in internal (deep) edges.

Requires: dendropy.  Reads tree.<rep>.{true,var,beast,mrbayes}.nwk from <sims-root>/<N>taxa/.
"""
import os, sys, math, argparse
import numpy as np, pandas as pd
import dendropy

METHODS = [("var","vine"),("beast","beast"),("mrbayes","mrbayes")]

def true_branches(path, tns):
    t = dendropy.Tree.get(path=path, schema="newick", taxon_namespace=tns)
    t.is_rooted = False
    t.encode_bipartitions(suppress_unifurcations=True)
    out = {}
    for e in t.edges():
        if e.length is None or e.bipartition is None: continue
        bm = e.bipartition.split_bitmask
        side = bin(bm).count("1")
        out[e.bipartition] = (e.length, min(side, len(tns)-side))
    return out

def est_means(path, tns, maxt):
    trees = dendropy.TreeList.get(path=path, schema="newick", taxon_namespace=tns)
    n = min(len(trees), maxt); acc = {}
    for t in trees[:n]:
        t.is_rooted = False
        t.encode_bipartitions(suppress_unifurcations=True)
        for e in t.edges():
            if e.length is None or e.bipartition is None: continue
            r = acc.get(e.bipartition)
            if r is None: acc[e.bipartition] = [e.length, 1]
            else: r[0]+=e.length; r[1]+=1
    return acc, n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims-root", default="/local/storage/no-backup/vine-benchmarks/dna_sims/hky_300sites")
    ap.add_argument("--sizes", default="10,25,50")
    ap.add_argument("--nsamp", type=int, default=10)
    ap.add_argument("--maxtrees", type=int, default=400)
    ap.add_argument("--out-csv", default=None)
    a = ap.parse_args()
    sizes = [int(s) for s in a.sizes.split(",")]

    rows=[]
    for nt in sizes:
        d=f"{nt}taxa"
        for rep in range(1, a.nsamp+1):
            tf=os.path.join(a.sims_root,d,f"tree.{rep}.true.nwk")
            if not os.path.isfile(tf): continue
            tns=dendropy.TaxonNamespace(); tb=true_branches(tf,tns)
            for mf,mn in METHODS:
                ef=os.path.join(a.sims_root,d,f"tree.{rep}.{mf}.nwk")
                if not os.path.isfile(ef): continue
                em,ntr=est_means(ef,tns,a.maxtrees)
                for bip,(tl,size) in tb.items():
                    r=em.get(bip)
                    est=(r[0]/r[1]) if r else np.nan
                    freq=(r[1]/ntr) if r else 0.0
                    rows.append((nt,rep,mn,size,size==1,tl,est,freq))
            print(f"  done {d} rep{rep}", file=sys.stderr)
    df=pd.DataFrame(rows,columns=["ntaxa","rep","method","clade_size","terminal","true_bl","est_bl","freq"])
    if a.out_csv: df.to_csv(a.out_csv,index=False)
    p=df.dropna(subset=["est_bl"]).copy(); p["err"]=p["est_bl"]-p["true_bl"]

    print("\n"+"="*84)
    print("Aggregate relative branch-length bias  Sigma(est-true)/Sigma(true), terminal vs internal")
    print("="*84)
    for nt in sizes:
        print(f"\n-- {nt} taxa --")
        for kind,mask in [("terminal",p.terminal),("internal",~p.terminal)]:
            for m in ["vine","beast","mrbayes"]:
                s=p[(p.ntaxa==nt)&(p.method==m)&mask]
                if s.empty: continue
                rel=s["err"].sum()/s["true_bl"].sum()
                print(f"  {kind:9s} {m:8s} n={len(s):5d}  relbias={rel:+.3f}")
    print("\n"+"="*84)
    print("Internal-clade recovery frequency (fraction of sampled trees containing each true split)")
    print("="*84)
    for nt in sizes:
        for m in ["vine","beast","mrbayes"]:
            s=df[(df.ntaxa==nt)&(df.method==m)&(~df.terminal)]
            if s.empty: continue
            print(f"  {nt:>4d} {m:8s} mean recovery freq={s.freq.mean():.3f}  fully-absent={100*(s.freq==0).mean():.1f}%")

if __name__ == "__main__":
    main()
