
import sys
import pandas as pd
import numpy as np
from ete3 import Tree


def get_sep(file):
    with open(file, "r") as f:
        line = f.readline()
        if "\t" in line:
            sep = "\t"
        elif "," in line:
            sep = ","
        else:
            sep = " "
    return sep


def get_site_category(label):
    site_category = ""
    if label == primary_tissue:
        site_category = "primary"
    else:
        site_category = "metastasis"
    return site_category


treefile = sys.argv[1]
tissues = sys.argv[2]
primary_tissue = sys.argv[3]
edges_outputfile = sys.argv[4]
metadata_outputfile = sys.argv[5]

# use tree to get edge list and branch lengths
try:
    tree = Tree(treefile, format=5) # Tree with tip labels and branch lengths
except:
    tree = Tree(treefile, format=1) # Fallback to general format

# Name internal nodes if not already named
i=1
for node in tree.traverse():
    if node.is_leaf() == False:
        if node.name == "":
            node.name = f"node{i}"
            i = i + 1

# Name the root node
tree.get_tree_root().name = "root"

edges = pd.DataFrame(columns=["node1", "node2"])
branch_lengths = {}
cluster_labels = []
cluster_index_label = {}

# Set origin info
cluster_index_label["origin"] = 0
branch_lengths[cluster_index_label["origin"]] = 0.0
cluster_labels.append("origin")

# Set root info
cluster_index_label[tree.get_tree_root().name] = 1
branch_lengths[cluster_index_label[tree.get_tree_root().name]] = (tree.get_tree_root().dist)
cluster_labels.append(tree.get_tree_root().name)

# Add origin to root edge
edges.loc[len(edges)] = [cluster_index_label["origin"], cluster_index_label[tree.get_tree_root().name]]

# set all other node info by children traversal
i = 2
for node in tree.traverse():
    if node.is_leaf() == True:
        name = node.name
    else:
        node_name = node.name
        children = node.children
        for child in children:
            child_name = child.name
            cluster_labels.append(child_name)
            cluster_index_label[child_name] = i
            branch_lengths[i] = child.dist
            edges.loc[len(edges)] = [cluster_index_label[node_name], i]
            i = i + 1

### output edges
edges.to_csv(edges_outputfile, sep=" ", index=False, header=False)


# read in tissues and format metadata tsv
tissues_df = pd.read_csv(tissues, sep=get_sep(tissues), names=["id", "tissue"])
tissues_df = pd.concat([tissues_df, pd.DataFrame({"id": ["origin"], "tissue": [primary_tissue]})], ignore_index=True)   # Add origin as known primary
tissues_dict = dict(zip(tissues_df["id"].astype(str), tissues_df["tissue"].astype(str)))
unique_tissues = set(tissues_dict.values())
if primary_tissue not in unique_tissues:
    unique_tissues.add(primary_tissue)
tissue_to_int = {tissue: i for i, tissue in enumerate(unique_tissues)}

### output metadata tsv
with open(metadata_outputfile, "w") as file:
    file.write(
        "\t".join(
            [
                "anatomical_site_index",
                "anatomical_site_label",
                "cluster_index",
                "cluster_label",
                "present",
                "site_category",
                "num_mutations",
            ]
        )
    )
    for cluster_label in cluster_labels:
        cluster_index = cluster_index_label[cluster_label]
        num_mutations = branch_lengths[cluster_index]
        for anatomical_site_label in unique_tissues:
            anatomical_site_index = tissue_to_int[anatomical_site_label]
            if anatomical_site_label != primary_tissue:
                site_category = "metastasis"
            else:
                site_category = "primary"
            if (
                cluster_label in tissues_dict
                and tissues_dict[cluster_label] == anatomical_site_label
            ):
                present = 1
            else:
                present = 0
            file.write(
                f"\n{anatomical_site_index}\t{anatomical_site_label}\t{cluster_index}\t{cluster_label}\t{present}\t{site_category}\t{num_mutations}"
            )
