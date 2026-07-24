#!/usr/bin/env python3
"""
Plot vine mixture-sweep comparisons across taxa counts.

Usage:
  python plot_posterior_quality.py [options]

  --taxa       Space-separated list of taxa counts  (default: 10 25 50 100 250)
  --root       Root directory containing {n}taxa/    (default: script directory)
  --out        Output directory for plots            (default: <root>/plots)
  --methods    Method summary prefixes               (default: vine vine_flows vine_mixture vine_mixture_flows)
  --labels     Optional plot labels for vine methods (default: inferred from --methods)
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


DEFAULT_METHODS = ["vine", "vine_flows", "vine_mixture", "vine_mixture_flows"]
METHOD_LABELS = {
    "bootstrap_nj": "Bootstrap NJ",
    "vine": "Vine",
    "vine_flows": "Vine + flows",
    "vine_mixture": "Vine + mixture",
    "vine_mixture_flows": "Vine + mixture + flows",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--taxa", nargs="+", type=int, default=[10, 25, 50, 100, 250])
    p.add_argument(
        "--methods",
        nargs="+",
        default=DEFAULT_METHODS,
        metavar="METHOD",
    )
    p.add_argument(
        "--labels",
        nargs="+",
        default=None,
        metavar="LABEL",
    )
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--baseline-root", type=Path, default=None,
                   help="Directory containing the original hky_300sites summaries")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()
    if args.labels is None:
        args.labels = [METHOD_LABELS.get(method, method.replace("_", " ").title()) for method in args.methods]
    if args.baseline_root is None:
        args.baseline_root = args.root.parent / "hky_300sites"
    if len(args.methods) != len(args.labels):
        p.error("--methods and --labels must have the same length")
    return args


def read_table(path):
    if not path.exists():
        return []
    rows, header = [], None
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if header is None:
                header = line.split()
                continue
            if set(line) == {"-"}:
                break
            try:
                rows.append({k: float(v) for k, v in zip(header, line.split())})
            except ValueError:
                continue
    return rows


def mean_sd(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan, np.nan
    if arr.size == 1:
        return float(arr[0]), 0.0
    return float(np.mean(arr)), float(np.std(arr, ddof=1))


def metric_value(row, metric, key):
    if metric not in ("lnl", "mf"):
        return row.get(key)
    true_lnl = row.get("true")
    method_lnl = row.get(key)
    if true_lnl is None or method_lnl is None or true_lnl == 0:
        return None
    return 100.0 * (method_lnl - true_lnl) / abs(true_lnl)


TABLE_NAME = {
    "ci50":    "eval.all.dist.txt",
    "ci95":    "eval.all.dist.txt",
    "rf":      "eval.all.rf.txt",
    "entropy": "eval.all.ent.txt",
    "lnl":     "eval.all.lnl.txt",
    "time":    "eval.all.time.txt",
    "mf":      "eval.all.mf.txt",
    "r2":      "eval.all.dist.txt",
    "rmse":    "eval.all.dist.txt",
    "bsd":     "eval.all.bsd.txt",
}

VINE_KEY = {
    "ci50": "vine_50CI", "ci95": "vine_95CI", "rf": "vine", "entropy": "vine_top", "lnl": "vine",
    "time": "vine", "mf": "vine", "r2": "vine_r2", "rmse": "vine_RMSE",
    "bsd": "vine",
}

BASELINE_KEYS = {
    "ci50": [
        ("BEAST2", "beast_50CI"),
        ("BEAST2 + BEAGLE", "beast-beagle_50CI"),
        ("MrBayes", "mrbayes_50CI"),
        ("MrBayes + BEAGLE", "mrbayes-beagle_50CI"),
    ],
    "ci95": [
        ("BEAST2", "beast_95CI"),
        ("BEAST2 + BEAGLE", "beast-beagle_95CI"),
        ("MrBayes", "mrbayes_95CI"),
        ("MrBayes + BEAGLE", "mrbayes-beagle_95CI"),
    ],
    "rf": [
        ("BEAST2", "beast"),
        ("BEAST2 + BEAGLE", "beast-beagle"),
        ("MrBayes", "mrbayes"),
        ("MrBayes + BEAGLE", "mrbayes-beagle"),
    ],
    "entropy": [
        ("BEAST2", "beast_top"),
        ("BEAST2 + BEAGLE", "beast-beagle_top"),
        ("MrBayes", "mrbayes_top"),
        ("MrBayes + BEAGLE", "mrbayes-beagle_top"),
    ],
    "lnl": [
        ("BEAST2", "beast"),
        ("BEAST2 + BEAGLE", "beast-beagle"),
        ("MrBayes", "mrbayes"),
        ("MrBayes + BEAGLE", "mrbayes-beagle"),
    ],
    "time": [
        ("BEAST2", "beast"),
        ("BEAST2 + BEAGLE", "beast-beagle"),
        ("MrBayes", "mrbayes"),
        ("MrBayes + BEAGLE", "mrbayes-beagle"),
    ],
    "mf": [
        ("BEAST2", "beast"),
        ("BEAST2 + BEAGLE", "beast-beagle"),
        ("MrBayes", "mrbayes"),
        ("MrBayes + BEAGLE", "mrbayes-beagle"),
    ],
    "r2": [
        ("BEAST2", "beast_r2"),
        ("BEAST2 + BEAGLE", "beast-beagle_r2"),
        ("MrBayes", "mrbayes_r2"),
        ("MrBayes + BEAGLE", "mrbayes-beagle_r2"),
    ],
    "rmse": [
        ("BEAST2", "beast_RMSE"),
        ("BEAST2 + BEAGLE", "beast-beagle_RMSE"),
        ("MrBayes", "mrbayes_RMSE"),
        ("MrBayes + BEAGLE", "mrbayes-beagle_RMSE"),
    ],
    "bsd": [
        ("BEAST2", "beast"),
        ("BEAST2 + BEAGLE", "beast-beagle"),
        ("MrBayes", "mrbayes"),
        ("MrBayes + BEAGLE", "mrbayes-beagle"),
    ],
}


def collect(args, methods, metric):
    baseline_keys = BASELINE_KEYS[metric]
    labels = [*args.labels, *[label for label, _ in baseline_keys]]
    data = {label: {"mean": [], "sd": []} for label in labels}
    for n in args.taxa:
        baseline_vals = {label: None for label, _ in baseline_keys}
        for method_dir, label in zip(methods, args.labels):
            path = args.root / f"{n}taxa" / f"{method_dir}.{TABLE_NAME[metric]}"
            rows = read_table(path)
            vk = VINE_KEY[metric]
            vals = [v for r in rows for v in [metric_value(r, metric, vk)] if v is not None]
            m, s = mean_sd(vals)
            data[label]["mean"].append(m)
            data[label]["sd"].append(s)
            for baseline_label, key in baseline_keys:
                bvals = [v for r in rows for v in [metric_value(r, metric, key)] if v is not None]
                if baseline_vals[baseline_label] is None and bvals:
                    baseline_vals[baseline_label] = bvals
        for baseline_label, _ in baseline_keys:
            m, s = mean_sd(baseline_vals[baseline_label] or [])
            data[baseline_label]["mean"].append(m)
            data[baseline_label]["sd"].append(s)
    return labels, data


METHOD_COLORS = {
    "Bootstrap NJ": "#A0CBE8",
    "Vine": "#F28E2B",
    "Vine + flows": "#76B7B2",
    "Vine + mixture": "#EDC948",
    "Vine + mixture + flows": "#A0CBE8",
    "BEAST2": "#59A14F",
    "BEAST2 + BEAGLE": "#8BC184",
    "MrBayes": "#E15759",
    "MrBayes + BEAGLE": "#E98A8C",
}

FALLBACK_COLORS = ["#76B7B2", "#EDC948", "#BAB0AC"]


def grouped_bar(ax, taxa, labels, colors, data, ylabel, ylim=None, hline=None, logy=False):
    x = np.arange(len(taxa))
    width = min(0.16, 0.8 / len(labels))
    offsets = (np.arange(len(labels)) - (len(labels) - 1) / 2) * width
    for label, offset in zip(labels, offsets):
        means = np.asarray(data[label]["mean"], dtype=float)
        sds = np.asarray(data[label]["sd"], dtype=float)
        if logy:
            means = np.where(means > 0, means, np.nan)
            sds = np.where(np.isfinite(means), sds, np.nan)
        ax.bar(
            x + offset, means, width, label=label, color=colors[label],
            zorder=2,
        )
        ax.errorbar(
            x + offset, means, yerr=sds, fmt="none", ecolor="black",
            elinewidth=1, capsize=0, zorder=4,
        )
    if hline is not None:
        ax.axhline(hline, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in taxa])
    ax.set_xlabel("Number of Taxa (n)")
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    if logy:
        ax.set_yscale("log")
    ax.grid(axis="y", color="#bdbdbd", linewidth=0.8)
    ax.grid(axis="x", color="#bdbdbd", linewidth=0.8)
    ax.set_axisbelow(True)


def collect_kl(args, methods, table_name, column, extra_baselines):
    """KL(vine-method || reference) plus fixed extra baseline-vs-reference bars.

    column: name of the KL column in the summary tables (e.g. "kl_pdist" or
    "kl_topo"). extra_baselines: list of (label, filename) pairs, each
    filename relative to the {n}taxa/ directory, containing a precomputed
    KL summary table.
    """
    labels = [*args.labels, *[label for label, _ in extra_baselines]]
    data = {label: {"mean": [], "sd": []} for label in labels}
    for n in args.taxa:
        for method_dir, label in zip(methods, args.labels):
            path = args.root / f"{n}taxa" / f"{method_dir}.{table_name}"
            rows = read_table(path)
            vals = [v for r in rows for v in [r.get(column)] if v is not None]
            m, s = mean_sd(vals)
            data[label]["mean"].append(m)
            data[label]["sd"].append(s)
        for label, filename in extra_baselines:
            rows = read_table(args.root / f"{n}taxa" / filename)
            vals = [v for r in rows for v in [r.get(column)] if v is not None]
            m, s = mean_sd(vals)
            data[label]["mean"].append(m)
            data[label]["sd"].append(s)
    return labels, data


def main():
    args = parse_args()
    outdir = args.out or (args.root / "plots")
    outdir.mkdir(exist_ok=True)

    methods = args.methods

    fig, axes = plt.subplots(3, 5, figsize=(24.0, 10.8))

    for metric, ax, ylabel, kwargs in [
        ("time",    axes[0, 0], "Runtime (seconds)",    {"logy": True}),
        ("lnl",     axes[0, 1], "% change from true\nlog likelihood", {}),
        ("mf",      axes[0, 2], "% change from true\nheld-out log likelihood", {}),
        ("rf",      axes[0, 3], "Normalized Robinson-Foulds\ndistance to true tree", {}),
        ("bsd",     axes[0, 4], "Normalized branch-score\ndistance to true tree", {}),
        ("ci50",    axes[1, 0], "50% CI coverage of\npairwise distances to truth", {"ylim": (0, 1.05)}),
        ("ci95",    axes[1, 1], "95% CI coverage of\npairwise distances to truth", {"ylim": (0, 1.05)}),
        ("entropy", axes[1, 2], "Posterior entropy of\ntree topology", {}),
        ("r2",      axes[1, 3], "R² of pairwise\ndistances to truth", {"ylim": (0, 1.05)}),
        ("rmse",    axes[1, 4], "RMSE of pairwise\ndistances to truth", {}),
    ]:
        labels, data = collect(args, methods, metric)
        colors = {
            label: METHOD_COLORS.get(label, FALLBACK_COLORS[i % len(FALLBACK_COLORS)])
            for i, label in enumerate(labels)
        }
        grouped_bar(ax, args.taxa, labels, colors, data, ylabel, **kwargs)

    for ax, ylabel, table_name, column, extra_baselines in [
        (axes[2, 0], "Mean marginal pairwise-distance\nKL divergence from BEAST2 posterior",
         "eval.all.kl_pdist_beast.txt", "kl_pdist", [
             ("BEAST2 + BEAGLE", "beast-beagle.eval.all.kl_pdist_beast.txt"),
             ("MrBayes", "mrbayes.eval.all.kl_pdist_beast.txt"),
             ("MrBayes + BEAGLE", "mrbayes-beagle.eval.all.kl_pdist_beast.txt"),
         ]),
        (axes[2, 1], "Mean marginal pairwise-distance\nKL divergence from MrBayes posterior",
         "eval.all.kl_pdist_mrbayes.txt", "kl_pdist", [
             ("MrBayes + BEAGLE", "mrbayes-beagle.eval.all.kl_pdist_mrbayes.txt"),
             ("BEAST2", "beast.eval.all.kl_pdist_mrbayes.txt"),
             ("BEAST2 + BEAGLE", "beast-beagle.eval.all.kl_pdist_mrbayes.txt"),
         ]),
        (axes[2, 2], "Mean per-split topology KL\ndivergence from BEAST2 posterior",
         "eval.all.kl_topo_beast.txt", "kl_topo", [
             ("BEAST2 + BEAGLE", "beast-beagle.eval.all.kl_topo_beast.txt"),
             ("MrBayes", "mrbayes.eval.all.kl_topo_beast.txt"),
             ("MrBayes + BEAGLE", "mrbayes-beagle.eval.all.kl_topo_beast.txt"),
         ]),
        (axes[2, 3], "Mean per-split topology KL\ndivergence from MrBayes posterior",
         "eval.all.kl_topo_mrbayes.txt", "kl_topo", [
             ("MrBayes + BEAGLE", "mrbayes-beagle.eval.all.kl_topo_mrbayes.txt"),
             ("BEAST2", "beast.eval.all.kl_topo_mrbayes.txt"),
             ("BEAST2 + BEAGLE", "beast-beagle.eval.all.kl_topo_mrbayes.txt"),
         ]),
    ]:
        labels, data = collect_kl(args, methods, table_name, column, extra_baselines)
        colors = {
            label: METHOD_COLORS.get(label, FALLBACK_COLORS[i % len(FALLBACK_COLORS)])
            for i, label in enumerate(labels)
        }
        grouped_bar(ax, args.taxa, labels, colors, data, ylabel, ylim=(0, None))

    axes[2, 4].axis("off")

    axes[0, 1].ticklabel_format(axis="y", style="plain", useOffset=False)

    handles, lbls = axes[0, 0].get_legend_handles_labels()
    kl_handles, kl_lbls = axes[2, 0].get_legend_handles_labels()
    for h, l in zip(kl_handles, kl_lbls):
        if l not in lbls:
            handles.append(h)
            lbls.append(l)
    # Keep Bootstrap NJ first, but give it an empty lower legend slot so the
    # remaining (method, method+BEAGLE) pairs stay vertically aligned.
    if "Bootstrap NJ" in lbls:
        i = lbls.index("Bootstrap NJ")
        bootstrap_label = lbls.pop(i)
        bootstrap_handle = handles.pop(i)
        lbls[0:0] = [bootstrap_label, ""]
        handles[0:0] = [bootstrap_handle, Line2D([], [], linestyle="none")]
    # Matplotlib fills legends column-first, so each adjacent pair is stacked.
    ncol = len(lbls) // 2
    fig.legend(handles, lbls, loc="lower center", ncol=ncol, frameon=False,
               bbox_to_anchor=(0.5, 0.01))
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(outdir / "hky300_summary.pdf")
    plt.close(fig)
    print(f"Wrote {outdir / 'hky300_summary.pdf'}")


if __name__ == "__main__":
    main()
