#!/usr/bin/env python3
"""Align and plot multiple evalTrees RF-MDS coordinate tables."""

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = [
    ("beast", "beast"),
    ("beast-beagle", "beast+beagle"),
    ("mrbayes", "mrbayes"),
    ("mrbayes-beagle", "mrbayes + beagle"),
    ("vine", "vine"),
    ("vine_flows", "vine + flows"),
    ("bootstrap_nj", "bootstrap nj"),
]
METHOD_ORDER = {method: index for index, (method, _) in enumerate(METHODS)}
METHOD_LABEL = dict(METHODS)
METHOD_COLOR = {
    "beast": "#59A14F",
    "beast-beagle": "#8BC184",
    "mrbayes": "#E15759",
    "mrbayes-beagle": "#E98A8C",
    "vine": "#F28E2B",
    "vine_flows": "#76B7B2",
    "bootstrap_nj": "#A0CBE8",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", nargs="+", required=True, type=Path,
                        help="One or more evalTrees *.mds.tsv files")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--width-per-panel", type=float, default=3.0)
    parser.add_argument("--height", type=float, default=2.5,
                        help="Height in inches per replicate row")
    return parser.parse_args()


def read_coordinates(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"tree", "x", "y"}.issubset(rows[0]):
        raise ValueError("{} is not an evalTrees RF-MDS table".format(path))
    return np.asarray([[float(row["x"]), float(row["y"])] for row in rows])


def oriented_principal_basis(points):
    """Return principal axes with a deterministic sign for each axis."""
    centered = points - points.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    basis = vt.T
    projected = centered.dot(basis)
    for axis in range(2):
        score = np.sum(projected[:, axis] ** 3)
        if abs(score) < 1e-12:
            score = projected[np.argmax(np.abs(projected[:, axis])), axis]
        if score < 0:
            basis[:, axis] *= -1
    return basis


def align_to_reference(points, reference_basis):
    centered = points - points.mean(axis=0)
    basis = oriented_principal_basis(points)
    return centered.dot(basis).dot(reference_basis.T)


def file_method(path):
    suffix = ".pairwise_rf.mds.tsv"
    stem = path.name[:-len(suffix)] if path.name.endswith(suffix) else path.stem
    for method, _ in METHODS:
        if stem.endswith("." + method):
            return method
    raise ValueError("cannot determine method from {}".format(path))


def file_replicate(path):
    match = re.search(r"(?:^|\.)tree\.(\d+)\.", path.name)
    if not match:
        raise ValueError("cannot determine replicate from {}".format(path))
    return int(match.group(1))


def main():
    args = parse_args()
    by_replicate = {}
    for path in args.coordinates:
        replicate = file_replicate(path)
        method = file_method(path)
        if method in by_replicate.setdefault(replicate, {}):
            raise ValueError("duplicate {} input for replicate {}".format(
                method, replicate))
        by_replicate[replicate][method] = path

    replicates = sorted(by_replicate)
    methods = [method for method, _ in METHODS
               if method in by_replicate[replicates[0]]]
    expected = set(methods)
    for replicate in replicates:
        if set(by_replicate[replicate]) != expected:
            raise ValueError("replicate {} does not have the same methods as "
                             "replicate {}".format(replicate, replicates[0]))

    aligned_rows = []
    for replicate in replicates:
        clouds = [read_coordinates(by_replicate[replicate][method])
                  for method in methods]
        reference_basis = oriented_principal_basis(clouds[0])
        aligned_rows.append([
            align_to_reference(points, reference_basis) for points in clouds
        ])

    nrows = len(replicates)
    ncols = len(methods)
    fig, axes = plt.subplots(
        nrows, ncols, squeeze=False,
        figsize=(args.width_per_panel * ncols, args.height * nrows),
        sharex="row", sharey="row",
    )

    for row_index, (replicate, row) in enumerate(zip(replicates, aligned_rows)):
        combined = np.vstack(row)
        xcenter = 0.5 * (combined[:, 0].min() + combined[:, 0].max())
        ycenter = 0.5 * (combined[:, 1].min() + combined[:, 1].max())
        half_span = max(np.ptp(combined[:, 0]), np.ptp(combined[:, 1])) * 0.55
        half_span = max(half_span, 0.05)
        xlim = (xcenter - half_span, xcenter + half_span)
        ylim = (ycenter - half_span, ycenter + half_span)
        for column_index, (method, points) in enumerate(zip(methods, row)):
            ax = axes[row_index, column_index]
            ax.scatter(points[:, 0], points[:, 1], s=32, alpha=0.65,
                       color=METHOD_COLOR[method], linewidths=0)
            if row_index == 0:
                ax.set_title(METHOD_LABEL[method])
            if row_index == nrows - 1:
                ax.set_xlabel("MDS1")
            if column_index == 0:
                ax.set_ylabel("replicate {}\n\nMDS2".format(replicate))
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_aspect("equal", adjustable="box")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    fig.tight_layout(pad=0.5, w_pad=0.25, h_pad=0.35)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out)
    plt.close(fig)


if __name__ == "__main__":
    main()
