#!/usr/bin/env python3

import sys
import pandas as pd
import dendropy


def fitch_hartigan_parsimony_count(tree, tissue_map, include_origin_branch=True, origin_state="0"):
    score = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            label = node.taxon.label
            node.state_set = {tissue_map[label]}
            continue
        
        child_sets = [child.state_set for child in node.child_node_iter()]
        if len(child_sets) == 2:
            intersection = set.intersection(*child_sets)
        elif include_origin_branch:
            origin_set = {origin_state}   # Assumed state at the origin of the tree
            child_set = child_sets[0]
            intersection = origin_set.intersection(child_set)
        else:
            RuntimeError("Unexpected number of child sets. Tree might not be binary.")
        
        if intersection:
            node.state_set = intersection
        else:
            node.state_set = set.union(*child_sets)
            score += 1
            
    return score


# newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.cass.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.var.mcc.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.laml_trees.nwk"
# # newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.beam.mcc.nwk"
# mutation_tsv = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.tsv"
# origin_state = "0"

newick_file = sys.argv[1]
mutation_tsv = sys.argv[2]
origin_state = str(sys.argv[3])

tree = dendropy.Tree.get(path=newick_file, schema="newick", rooting="default-rooted", preserve_underscores=True)

mut_df = pd.read_csv(mutation_tsv, sep="\t", index_col=0, dtype=str)

total_parsimony = 0
for site in mut_df.columns:
    state_map = mut_df[site].to_dict()
    total_parsimony += fitch_hartigan_parsimony_count(tree, state_map, include_origin_branch=True, origin_state=origin_state)

print(f"MinParsimonyAcrossSites\t{total_parsimony}")
