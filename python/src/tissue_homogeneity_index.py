#!/usr/bin/env python3

import sys
import math
import numpy as np
import pandas as pd
import dendropy


def same_tissue_probability(tissues):
    """
    Calculate the probability that a randomly chosen pair of cells from the 
    given list of tissues belongs to the same tissue type (sampling without replacement)
    across all unique tissues and their empirical frequencies.
    """
    n = len(tissues)
    freqs = pd.Series(tissues).value_counts()
    return sum((c / n) * ((c - 1) / (n - 1)) for c in freqs)


# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.cass.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.var.mcc.nwk"
# newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.laml_trees.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.beam.mcc.nwk"
# tissues_csv = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.tissues.expanded.csv"

newick_file = sys.argv[1]
tissues_csv = sys.argv[2]

# read tree
tree = dendropy.Tree.get(path=newick_file, schema="newick", rooting="default-rooted", preserve_underscores=True)

# If there are no branch lengths (cassiopeia), then make all lengths 1
branch_lengths_present = any(edge.length is not None for edge in tree.postorder_edge_iter())
if not branch_lengths_present:
    for edge in tree.postorder_edge_iter():
        edge.length = 1.0

# Read tissues
tissue_df = pd.read_csv(tissues_csv, header=None, names=["cell", "tissue"])
tissue_map = tissue_df.set_index("cell")["tissue"]

# Get tree tip labels
tip_labels = [leaf.taxon.label for leaf in tree.leaf_node_iter()]

# Global expected same-tissue probability
global_tissues = tissue_map.loc[tip_labels].to_numpy()
expected = same_tissue_probability(global_tissues)

# For each node, calculate the observed same-tissue pair fraction of that node's tips and compute the ratio to expected
clade_scores = []
for node in tree.postorder_node_iter():
    if node.is_leaf():
        continue    # Skip tips
    
    leaves = node.leaf_nodes()
    if len(leaves) < 2:
        continue    # Skip final nodes before tips, since we need at least 2 cells to form a pair for the statistic
    
    clade_tip_labels = [leaf.taxon.label for leaf in leaves]
    clade_tissues = tissue_map.loc[clade_tip_labels].to_numpy()
    
    observed = same_tissue_probability(clade_tissues)
    clade_scores.append(observed / expected)

thi = float(np.mean(clade_scores))
print(f"THI\t{thi:.4f}\tinternal_clades={len(clade_scores)}\texpected_same_tissue_pair_prob={expected:.6f}")
