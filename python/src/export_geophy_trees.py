#!/usr/bin/env python3
"""Draw Newick trees from a fitted GeoPhy state (no training)."""

from __future__ import annotations

import argparse

import geophy as gp
from geophy.utils import init_seed, parse_config
import torch


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-prefix", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--count", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    init_seed(args.seed)
    torch.set_num_threads(1)
    config = parse_config(f"{args.state_prefix}.latest.yaml")
    reader = gp.StateReader(config, f"{args.state_prefix}.latest.pt")
    states = reader.sample_states(mc_samples=args.count)
    trees = states["samples"]["utree_samples"]
    branch_lens = states["samples"]["branch_lengths"]
    names = reader._seq_data.names
    if len(trees) != args.count:
        raise SystemExit(f"expected {args.count} GeoPhy trees, got {len(trees)}")
    with open(args.output, "w", encoding="utf-8") as handle:
        for tree, blens in zip(trees, branch_lens):
            # tree.branch_lengths (as_ete_tree's default) is a stale embedding-
            # geometry estimate GeoPhyModel.sample_states() bakes into each
            # TreeMetric before the branch-length model runs; blens (the same
            # tensor mean_tlen/the ELBO are computed from) is the actual fitted
            # posterior draw, so it must be passed explicitly.
            newick = tree.as_ete_tree(taxon_names=names, branch_lengths=blens).write(format=1)
            handle.write(newick + "\n")


if __name__ == "__main__":
    main()
