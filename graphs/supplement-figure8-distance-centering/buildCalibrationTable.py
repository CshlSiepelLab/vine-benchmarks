#!/usr/bin/env python3
"""Build the tidy calibration table for the centering supplement figure.

For every (ntaxa, replicate, method) it reads the per-pair distance summary produced
by `evalTrees` (raw/<ntaxa>.<rep>.<method>.dist, columns: leaf1 leaf2 mean ...) and
joins it, by leaf pair, to the true pairwise distances (tree.<rep>.true.dist.txt in the
simulation directory).  Output is a long CSV:

    ntaxa,rep,method,leaf1,leaf2,true,est

where `est` is the posterior-/variational-mean pairwise distance and `true` is the
distance in the simulation ground-truth tree.  Methods use the internal keys
var/beast/mrbayes and are relabelled to vine/beast/mrbayes to match the Fig. 2 palette.
"""
import os, sys, glob, argparse, csv

MKEY = {"var": "vine", "beast": "beast", "mrbayes": "mrbayes"}

def load_true(path):
    t = {}
    with open(path) as fh:
        for ln in fh:
            s = ln.split()
            if not s or ln.lstrip().startswith("#") or s[0].lower().startswith("leaf"):
                continue
            t[(s[0], s[1])] = float(s[2])
    return t

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims-root", required=True,
                    help="dir holding <N>taxa/ simulation subdirs with tree.*.true.dist.txt")
    ap.add_argument("--raw-dir", required=True,
                    help="dir holding evalTrees output raw/<N>taxa.<rep>.<method>.dist")
    ap.add_argument("--sizes", required=True, help="comma-separated taxon counts, e.g. 25,50,100")
    ap.add_argument("--nsamp", type=int, default=10, help="replicates per size")
    ap.add_argument("--methods", default="var,beast,mrbayes")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sizes = [s.strip() for s in a.sizes.split(",")]
    methods = [m.strip() for m in a.methods.split(",")]
    nrows = 0
    with open(a.out, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["ntaxa", "rep", "method", "leaf1", "leaf2", "true", "est"])
        for sz in sizes:
            for rep in range(1, a.nsamp + 1):
                tf = os.path.join(a.sims_root, f"{sz}taxa", f"tree.{rep}.true.dist.txt")
                if not os.path.isfile(tf):
                    continue
                truth = load_true(tf)
                for m in methods:
                    rf = os.path.join(a.raw_dir, f"{sz}taxa.{rep}.{m}.dist")
                    if not os.path.isfile(rf) or os.path.getsize(rf) == 0:
                        continue
                    with open(rf) as fh:
                        for ln in fh:
                            s = ln.split()
                            if len(s) < 3:
                                continue
                            k = (s[0], s[1])
                            if k in truth:
                                w.writerow([sz, rep, MKEY.get(m, m), s[0], s[1],
                                            f"{truth[k]:.6f}", f"{float(s[2]):.6f}"])
                                nrows += 1
    print(f"wrote {a.out}: {nrows} rows", file=sys.stderr)

if __name__ == "__main__":
    main()
