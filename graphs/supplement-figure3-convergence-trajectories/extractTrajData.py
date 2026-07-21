#!/usr/bin/env python3
"""Build data/convergence_traj.csv: log-likelihood vs. wall-clock time for
vine (Taylor-approximated ELBO optimizer) vs. BEAST2,
across the hky_300sites taxa-size sweep (10-500 taxa). Used for the
vine-vs-MCMC convergence supplement figure (makeLoglikGraphs.R).

Per-iteration wall-clock time is not logged directly by any method, so we
linearly interpolate each trajectory's iteration index against the method's
*total* run time for that replicate (from eval.all.time.txt).

"""
import csv
import os
import re

SIMDIR = "../../dna_sims/hky_300sites"
TAXA_SIZES = [10, 25, 50, 100, 250, 500]
REPS = [1, 2, 3]


def read_vine_traj(path):
    """Return list of (state, ll) from a vine .var.nwk.log file."""
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("state\t"):
                header = line.split("\t")
                continue
            if line.startswith("#") or not line.strip() or header is None:
                continue
            parts = line.split("\t")
            rec = dict(zip(header, parts))
            rows.append((int(rec["state"]), float(rec["ll"])))
    return rows


def read_beast_traj(path):
    """Return list of (sample, likelihood) from a BEAST tracer .log file."""
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("Sample\t"):
                header = line.split("\t")
                continue
            if line.startswith("#") or not line.strip() or header is None:
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                continue
            rec = dict(zip(header, parts))
            rows.append((int(rec["Sample"]), float(rec["likelihood"])))
    return rows


def read_time_table(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    header = lines[0].split()
    out = {}
    for line in lines[1:]:
        toks = [t for t in re.split(r"\t+", line) if t != ""]
        if len(toks) != len(header):
            continue 
        rec = dict(zip(header, toks))
        try:
            out[int(rec["samp"])] = rec
        except (KeyError, ValueError):
            pass
    return out


def main():
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "data", "convergence_traj.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    rows_out = []
    for ntaxa in TAXA_SIZES:
        d = os.path.join(SIMDIR, f"{ntaxa}taxa")
        times = read_time_table(os.path.join(d, "eval.all.time.txt"))

        for rep in REPS:
            vine_total = float(times[rep]["vine"])
            beast_total = float(times[rep]["beast"])

            vine_traj = read_vine_traj(os.path.join(d, f"tree.{rep}.var.nwk.log"))
            vine_max = max(s for s, _ in vine_traj)
            for state, ll in vine_traj:
                t_sec = (state / vine_max) * vine_total
                rows_out.append((ntaxa, "vine", rep, f"{t_sec:.4f}", f"{ll:.4f}"))

            beast_traj = read_beast_traj(os.path.join(d, f"tree.{rep}.beast.log"))
            beast_max = max(s for s, _ in beast_traj)
            for sample, ll in beast_traj:
                t_sec = (sample / beast_max) * beast_total
                rows_out.append((ntaxa, "BEAST2", rep, f"{t_sec:.4f}", f"{ll:.4f}"))

        print(f"{ntaxa} taxa: reps {REPS}")

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ntaxa", "method", "rep", "time_sec", "loglik"])
        w.writerows(rows_out)
    print("wrote", out_path, f"({len(rows_out)} rows)")


if __name__ == "__main__":
    main()
