import os
import sys
import pandas as pd
import numpy as np


def process_migrations_file(filepath):
    '''
    Get true migration graph counts from simulation output .migrations file with source -> recipient
    migration event format without any file header.
    '''
    counts = {}
    with open(filepath, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if "->" not in line:
                raise ValueError(f"Bad line (missing '->'): {raw.rstrip()}")
            source, recipient = [x.strip() for x in line.split("->", 1)]
            migration = f"{source}_{recipient}"
            counts[migration] = counts.get(migration, 0) + 1
    return counts


def calculate_metrics(true_counts, inferred_counts):
    TP = 0
    FP = 0
    FN = 0
    all_keys = set(true_counts.keys()).union(set(inferred_counts.keys()))
    for key in all_keys:
        if key in inferred_counts:
            inferred_count = inferred_counts[key]
            if key in true_counts:
                true_count = true_counts[key]
                if inferred_count >= true_count:
                    TP += true_count
                    FP += inferred_count - true_count
                else:
                    TP += inferred_count
                    FN += true_count - inferred_count
            else:
                FP += inferred_count
        else:
            FN += true_counts[key]
    # compute precision as TP/(TP + FP) and recall as TP/(TP + FN)
    if (TP + FP) != 0:
        precision = TP / (TP + FP)
    else:
        precision = 0
    if (TP + FN) != 0:
        recall = TP / (TP + FN)
    else:
        recall = 0
    # calculate F1 score (2((precision * recall)/(precision + recall)))
    if precision + recall == 0:
        f1 = 0
    else:
        f1 = 2 * ((precision * recall) / (precision + recall))
    return f1, recall, precision


def posterior_threshold_metrics(posterior_prob_graph, true_counts, t=0.50):
    thresholds = [j for j in np.arange(0, 1.0, 0.01)]
    rows = []
    for thresh in thresholds:
        thresh_counts = {key: value for key, value in posterior_prob_graph.items() if value > thresh}
        if len(thresh_counts) == 0:
            f1 = 0.0
            recall = 0.0
            precision = 1.0
        else:
            edges = ["_".join(edge.split("_")[:-1]) for edge in thresh_counts.keys()]
            thresh_counts = {}
            for edge in edges:
                if edge not in thresh_counts:
                    thresh_counts[edge] = 1
                else:
                    thresh_counts[edge] += 1
            f1, recall, precision = calculate_metrics(true_counts, thresh_counts)
        rows.append({"threshold": thresh, "precision": precision, "recall": recall, "f1": f1,})
    thresh_df = pd.DataFrame(rows)
    return thresh_df


true_graph_file = sys.argv[1]
prob_graph_file = sys.argv[2]
outputfile = sys.argv[3]


# Read in true graph
true_counts = process_migrations_file(true_graph_file)

# Read in inferred graph
prob_graph = {}
with open(prob_graph_file, "r") as file:
    for line in file.readlines():
        line = line.strip().split(",")
        prob_graph[line[0]] = float(line[1])

# Calculate performance metrics
thresh_df = posterior_threshold_metrics(prob_graph, true_counts)

# Write outputs
thresh_df.to_csv(outputfile, index=False)

