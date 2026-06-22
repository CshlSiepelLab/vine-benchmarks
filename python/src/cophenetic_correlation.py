#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import dendropy
from itertools import combinations


def mutation_distance(vec1, vec2, missing_value="-1", ignore_missing=False):
    if ignore_missing:
        mask = (vec1 != missing_value) & (vec2 != missing_value)
        if mask.sum() == 0:
            return None
        return np.mean(vec1[mask] != vec2[mask])
    else:
        return np.mean(vec1 != vec2)    # Simple normalized Hamming distance here. Using sum (actual Hamming distance) instead of mean gives the same result since the correlation is scale invariant
    

# newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.cass.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.var.mcc.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.laml_trees.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.beam.mcc.nwk"
# mutation_tsv = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.tsv"
# missing_value = "-1"


newick_file = sys.argv[1]
mutation_tsv = sys.argv[2]
missing_value = str(sys.argv[3])


# Read in the tree
tree = dendropy.Tree.get(path=newick_file, schema="newick", rooting="default-rooted", preserve_underscores=True)

# If there are no branch lengths (cassiopeia), then make all lengths 1
branch_lengths_present = any(edge.length is not None for edge in tree.postorder_edge_iter())
if not branch_lengths_present:
    for edge in tree.postorder_edge_iter():
        edge.length = 1.0

# Compute the distance matrix implied by the tree
pdm = tree.phylogenetic_distance_matrix()   # Patristic distance matrix between tip pairs

# Read in the mutation matrix
df = pd.read_csv(mutation_tsv, sep="\t", index_col=0, dtype=str)

# Keep only shared samples
tree_samples = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
shared = [sample for sample in df.index if sample in tree_samples]
df = df.loc[shared]

# Map sample name -> taxon object in the tree
taxon_namespace = {taxon.label: taxon for taxon in tree.taxon_namespace}

# Compute the mutatin and tree distances
tree_dists = []
mut_dists = []
for a, b in combinations(shared, 2):
    d_tree = pdm.distance(taxon_namespace[a], taxon_namespace[b])
    
    vec1 = df.loc[a].to_numpy()
    vec2 = df.loc[b].to_numpy()
    d_mut = mutation_distance(vec1, vec2, missing_value=missing_value, ignore_missing=False)
    
    if d_mut is not None:
        tree_dists.append(d_tree)
        mut_dists.append(d_mut)

tree_dists = np.array(tree_dists, dtype=float)
mut_dists = np.array(mut_dists, dtype=float)

r = np.corrcoef(tree_dists, mut_dists)[0, 1] # Off diagonal entry is the correlation coefficient for this 2x2 matrix
print(f"CopheneticCorrelation\t{r:.4f}\tn_pairs={len(tree_dists)}")