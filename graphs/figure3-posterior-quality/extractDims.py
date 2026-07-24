#!/usr/bin/env python3
"""Build dimSummary25.txt / dimSummary50.txt for supplement figure S4.

For each embedding dimension and each geometry (Euclidean, hyperbolic):
  - Delta lnl (VINE - BEAST) is summarized as mean +/- SD over the 10
    replicates (population SD, matching the `stats` tool);
  - VINE runtime is summarized as median / Q1 / Q3 over the 10 replicates
    (robust to the heavy-tailed hyperbolic convergence time).

Output columns (tab-separated):
  d  dlnlE dlnlE_sd  timeE_med timeE_q1 timeE_q3  dlnlH dlnlH_sd  timeH_med timeH_q1 timeH_q3
"""
import os, numpy as np

ROOT = "/DATA/local/storage/no-backup/vine-benchmarks/dna_sims/hky_300sites"
OUT_DIR = "/local/storage/no-backup/vine-benchmarks/graphs/hky300-data"
DIMS = [f"archive.D{i}" for i in range(2, 9)]

def rows(path):
    """the 10 replicate data rows (header + 10 reps; skip any grand row)."""
    with open(path) as fh:
        ls = [l for l in fh if l.strip()]
    return [l.split() for l in ls[1:11]]

def dlnl_stats(lnl_path):
    r = rows(lnl_path)
    d = np.array([float(x[3]) - float(x[4]) for x in r])  # vine - beast
    return d.mean(), d.std()  # population SD (ddof=0), matches `stats`

def time_stats(time_path):
    r = rows(time_path)
    t = np.array([float(x[3]) for x in r])  # vine time
    return np.median(t), np.percentile(t, 25), np.percentile(t, 75)

def build(size):
    euc = os.path.join(ROOT, f"testdata.dim{size}")
    hyp = os.path.join(ROOT, f"testdata.dim{size}H")
    out = os.path.join(OUT_DIR, f"dimSummary{size}.txt")
    with open(out, "w") as fh:
        fh.write("d\tdlnlE\tdlnlE_sd\ttimeE_med\ttimeE_q1\ttimeE_q3"
                 "\tdlnlH\tdlnlH_sd\ttimeH_med\ttimeH_q1\ttimeH_q3\n")
        for sub in DIMS:
            le_m, le_s = dlnl_stats(os.path.join(euc, sub, "eval.all.lnl.txt"))
            te = time_stats(os.path.join(euc, sub, "eval.all.time.txt"))
            lh_m, lh_s = dlnl_stats(os.path.join(hyp, sub, "eval.all.lnl.txt"))
            th = time_stats(os.path.join(hyp, sub, "eval.all.time.txt"))
            fh.write(f"{sub}\t{le_m:.4f}\t{le_s:.4f}\t{te[0]:.4f}\t{te[1]:.4f}\t{te[2]:.4f}"
                     f"\t{lh_m:.4f}\t{lh_s:.4f}\t{th[0]:.4f}\t{th[1]:.4f}\t{th[2]:.4f}\n")
    print("wrote", out)

for sz in (25, 50):
    build(sz)
