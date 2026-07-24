#!/usr/bin/env python3
"""Build data/convergence_traj.csv: log-likelihood vs. wall-clock time for
vine (Taylor-approximated ELBO optimizer) vs. BEAST2,
across the hky_300sites taxa-size sweep (10-500 taxa). Used for the
vine-vs-MCMC convergence supplement figure (makeLoglikGraphs.R).

Per-iteration wall-clock time is not logged directly by any method, so we
linearly interpolate each trajectory's iteration index against the method's
*total* run time for that taxa size, taken from the mean-time-per-method
summary in ../hky300-data/timeSummary.txt. That file only has the mean (and
std) over replicates, not each replicate's own individually measured time,
so every replicate at a given taxa size is normalized against the same mean
total time.
"""
import csv
import os

SIMDIR = "../../dna_sims/hky_300sites"
TIME_SUMMARY_PATH = "../hky300-data/timeSummary.txt"
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


def read_time_summary(path):
    """Return {ntaxa: {"vine": mean_sec, "beast": mean_sec}} from
    timeSummary.txt. Columns are positional (ntaxa, vine, std, beast, std,
    beast-beagle, std, mrbayes, std, mrbayes-beagle, std) -- "std" repeats,
    so this can't be read by header name."""
    out = {}
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            toks = line.split()
            if not toks:
                continue
            out[int(toks[0])] = {"vine": float(toks[1]), "beast": float(toks[3])}
    return out


def main():
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "data", "convergence_traj.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    time_summary = read_time_summary(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), TIME_SUMMARY_PATH))

    rows_out = []
    for ntaxa in TAXA_SIZES:
        d = os.path.join(SIMDIR, f"{ntaxa}taxa")
        vine_total = time_summary[ntaxa]["vine"]
        beast_total = time_summary[ntaxa]["beast"]

        for rep in REPS:
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
