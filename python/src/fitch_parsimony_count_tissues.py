#!/usr/bin/env python3

import sys
import pandas as pd
import dendropy


def fitch_hartigan_parsimony_count(tree, tissue_map, include_origin_branch=True, primary_tissue="LL"):
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
            origin_set = {primary_tissue}   # Assumed state at the origin of the tree
            child_set = child_sets[0]
            intersection = origin_set.intersection(child_set)
        else:
            RuntimeError("Unexpected number of child sets. Tree might not be binary.")
        
        
        if intersection:
            node.state_set = intersection
        else:
            node.state_set = set.union(*child_sets)
            score += 1  # Increment score for each transition needed when there is no intersection possible
            
    return score


# newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.cass.nwk"
# newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.var.mcc.nwk"
# newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.laml_trees.nwk"
# newick_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.beam.mcc.nwk"
# tissues_csv = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.70.tissues.expanded.csv"
# primary_tissue = "LL"

newick_file = sys.argv[1]
tissues_csv = sys.argv[2]
primary_tissue = str(sys.argv[3])

# Read tree
tree = dendropy.Tree.get(path=newick_file, schema="newick", rooting="default-rooted", preserve_underscores=True)

# Read tissues
tissue_df = pd.read_csv(tissues_csv, header=None, names=["cell", "tissue"])
tissue_map = tissue_df.set_index("cell")["tissue"]

migration_count = fitch_hartigan_parsimony_count(tree, tissue_map, include_origin_branch=True, primary_tissue=primary_tissue)
print(f"MinParsimonyMigrations\t{migration_count}")
