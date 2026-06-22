#!/usr/bin/env python3

import argparse
from collections import Counter
import dendropy


def migration_graph_signature(tree, use_multiplicity=False):
    """
    Build a canonical representation of the migration graph for one tree.

    A migration is counted for each edge parent -> child where
    parent.state != child.state.

    If use_multiplicity is False:
        graph = set of unique directed edges, e.g.
        (("LL","M1"), ("M1","M2"))

    If use_multiplicity is True:
        graph = sorted tuple of (src, dst, count), e.g.
        (("LL","M1",2), ("M1","M2",5))
    """
    edges = []
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        parent = node.parent_node
        parent_state = str(parent.annotations.get_value("state"))
        child_state = str(node.annotations.get_value("state"))
        if parent_state is None or child_state is None:
            continue
        if parent_state != child_state:
            edges.append((parent_state, child_state))
    if use_multiplicity:
        edge_counts = Counter(edges)
        return tuple(sorted((src, dst, count) for (src, dst), count in edge_counts.items()))
    else:
        return tuple(sorted(set(edges)))


nexus_file = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.4.var.nex"

trees = dendropy.TreeList.get(path=nexus_file, schema="nexus", preserve_underscores=True, extract_comment_metadata=True)


graph_counts = Counter()

for tree in trees:
    sig = migration_graph_signature(tree, use_multiplicity=False)
    print(sig)
    graph_counts[sig] += 1

total_trees = sum(graph_counts.values())

print(f"Read {total_trees} trees")
print(f"Found {len(graph_counts)} distinct migration graphs\n")
