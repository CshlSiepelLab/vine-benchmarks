
import sys
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from ete3 import Tree
import dendropy
from copy import deepcopy
import networkx as nx

DEFAULT_COLORS = [
    "#006400",
    "#FF0000",
    "#0000CD",
    "#FFA500",
    "#800080",
    "#808080",
    "#FFC0CB",
    "#ADD8E6",
    "#A52A2A",
    "#FFFF00",
] * 3


def get_migration_counts(tree):
    tissues = []
    migration_counts = {}
    for node in tree.traverse():
        if node.is_root():
            continue
        else:
            node_tissue = node.name.split("_")[-1]
            parent_tissue = node.up.name.split("_")[-1]
            if node_tissue not in tissues:
                tissues.append(node_tissue)
            if parent_tissue not in tissues:
                tissues.append(parent_tissue)
            if node_tissue == parent_tissue:
                continue
            migration = f"{parent_tissue}_{node_tissue}"
            if migration not in migration_counts:
                migration_counts[migration] = 1
            else:
                migration_counts[migration] += 1
    return migration_counts, sorted(tissues)


def process_tree(filepath):
    # read in tree files to ete3 tree
    tree = Tree(filepath, format=8)
    # set tree root to primary
    tree.get_tree_root().name = f"0_{primary_tissue}"
    # get counts of migration events in a dict with source_recipient tissue key and count integer value
    counts, tissues = get_migration_counts(tree)
    return counts, tissues


newick = sys.argv[1]
primary_tissue = sys.argv[2]
outfile = sys.argv[3]


# read in the true tree to ete3 and get the migration counts
migration_graph, all_tissues = process_tree(newick)

# find all tissues to set the node colors
all_tissues = sorted(list(set(all_tissues) - {primary_tissue}))
custom_colors = {
    node: color
    for node, color in zip(all_tissues, DEFAULT_COLORS[0 : len(all_tissues)])
    if node != primary_tissue
}
all_tissues = [primary_tissue] + all_tissues
custom_colors[primary_tissue] = "black"

# plot the probability graph with edge thicknesses proportional to the probability
G = nx.MultiDiGraph()
for node in all_tissues:
    G.add_node(
        node,
        color=custom_colors[node],
        shape="box",
        fillcolor="white",
        penwidth=3.0,
        fontsize=32,
    )
for edge, num in migration_graph.items():
    source, target = edge.split("_")
    label = ""
    if num > 1:
        label = f"{num}"
    G.add_edge(
        source,
        target,
        color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
        penwidth=3,
        label=label,
        fontsize=24,
    )
dot = nx.nx_pydot.to_pydot(G)
dot.write_pdf(outfile)