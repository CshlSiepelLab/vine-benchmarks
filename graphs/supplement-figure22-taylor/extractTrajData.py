#!/usr/bin/env python3
"""Build data/elbo_traj.csv: per-iteration ELBO trajectories for a few
representative replicates, for both the Taylor (batchsize-100) and Monte Carlo
(batchsize-100) runs, at 25 and 50 taxa. Used for the convergence panels of the
supplement figure.

Each run starts from the same initialization; the two curves per replicate show
that optimizing the Taylor-approximated ELBO vs. a 100-sample Monte Carlo ELBO
converges at a similar rate to essentially the same optimum. Several replicates
are shown (trees 1-3) to demonstrate this is not cherry-picked. (This is a
comparison of the two optimization runs, not a same-(mu,Sigma) per-iteration
estimator comparison.)
"""
import re, csv, os

SIMDIR = "/DATA/local/storage/no-backup/vine-benchmarks/dna_sims/hky_300sites"
REPS = [1, 2, 3]  # representative replicates (trees) shown in the trajectory panels
CONDS = [(25, "Taylor",      f"taylor-25taxa"),
         (25, "Monte Carlo", f"montecarlo-25taxa"),
         (50, "Taylor",      f"taylor-50taxa"),
         (50, "Monte Carlo", f"montecarlo-50taxa")]

def read_traj(path):
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("state\t"):
                header = [h for h in line.split("\t") if h != ""]
                continue
            if line.startswith("#") or not line.strip() or header is None:
                continue
            parts = [p for p in line.split("\t") if p != ""]
            rec = dict(zip(header, parts))
            try:
                rows.append((int(rec["state"]), float(rec["elbo"])))
            except (KeyError, ValueError):
                pass
    return rows

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "elbo_traj.csv")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ntaxa", "method", "rep", "iter", "elbo"])
    for sz, method, d in CONDS:
        for rep in REPS:
            log = os.path.join(SIMDIR, d, f"tree.{rep}.var.nwk.log")
            for it, elbo in read_traj(log):
                w.writerow([sz, method, rep, it, f"{elbo:.4f}"])
print("wrote", out)
