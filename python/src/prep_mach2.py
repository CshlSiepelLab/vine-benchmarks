
import sys
from ete3 import Tree
import pandas as pd

leaf_labeled_tree = sys.argv[1]
leaf_labels_csv = str(sys.argv[2])
primary_tissue = str(sys.argv[3])
output_file_labeling = sys.argv[4]
output_file_tree = sys.argv[5]
output_file_colors = sys.argv[6]


# use tree to get edge list and branch lengths
try:
    tree = Tree(leaf_labeled_tree, format=5) # Tree with tip labels and branch lengths
except:
    tree = Tree(leaf_labeled_tree, format=1) # Fallback to general format

# Name internal nodes if not already named and remove tissue labels for internal node names, if they exist
i=1
for node in tree.traverse():
    if node.is_root():
        node.name = "root"
    elif not node.is_leaf():
        current_name = node.name
        new_name = current_name
        # Special renaming for laml outputs with polytomies resolved and named "node_dummy"
        if new_name == "node" and current_name.split("_")[1] == "dummy":
            new_name = f"nodedummy{i}"
            i += 1
        elif new_name == "":
            new_name = f"node{i}"
            i += 1
        node.name = new_name

# Get edge list
edges = []
for node in tree.traverse():
    if not node.is_leaf():
        for child in node.children:
            edges.append((node.name, child.name))

# Read in tissues to list
leaf_label = pd.read_csv(leaf_labels_csv, names=["leaf", "tissue"])
tissues = leaf_label["tissue"].unique().tolist()

# Fix when primary tissue is not a leaf label, but required in coloring scheme for MACHINA to run
if primary_tissue not in tissues:
    tissues.append(primary_tissue)

# Sort tissues to have primary tissue first and others in alphanumeric order, for consistent coloring scheme across runs
tissues = sorted(tissues, key=lambda x: (x != primary_tissue, x))

i = 1
color_map = {}
for tissue in tissues:
    color_map[tissue] = i
    i += 1

# output files
with open(output_file_tree, "w") as file:
    for edge in edges:
        file.write(f"{edge[0]}\t{edge[1]}\n")

leaf_label.to_csv(output_file_labeling, sep="\t", index=False, header=False)

with open(output_file_colors, "w") as file:
    for key, value in color_map.items():
        file.write(f"{key}\t{value}\n")
