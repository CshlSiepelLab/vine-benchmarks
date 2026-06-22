
import argparse
import pandas as pd
import cassiopeia as cas
from ete3 import Tree
        

parser = argparse.ArgumentParser()
parser.add_argument("character_matrix_csv", help="Path to input character matrix CSV file")
parser.add_argument("outfile", help="Output nwk file for saving the inferred tree")
args = parser.parse_args()

character_matrix_csv = args.character_matrix_csv
outfile = args.outfile

# Read in matrix
final_matrix = pd.read_csv(character_matrix_csv, sep=",", index_col=0)
final_matrix.index = final_matrix.index.astype(str)

# Solve cassiopeia greedy
reconstructed_tree = cas.data.CassiopeiaTree(character_matrix=final_matrix, missing_state_indicator=-1)
greedy_solver = cas.solver.VanillaGreedySolver()
greedy_solver.solve(reconstructed_tree)

# Make ete3 tree to write newick with internal node labels
connections = reconstructed_tree.edges
tree = Tree.from_parent_child_table(connections)

# Rename internal nodes from cassiopeia defaults
i = 0
for node in tree.traverse():
    if not node.is_leaf():
        node.name = f"node{i}"
        i += 1

with open(outfile, "w") as it:
    # it.write(tree.write(format=8))    # Write node and leaf names, but no branch lengths
    it.write(tree.write(format=9))  # Write leaf names only and no branch lengths
