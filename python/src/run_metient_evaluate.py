
import os, sys
import gzip
import pickle
from metient import metient as met
import pandas as pd
import numpy as np
import math
import torch


tree = sys.argv[1]
metadata = sys.argv[2]
patient = sys.argv[3]
primary = sys.argv[4]
output_dir = sys.argv[5]
outfile = sys.argv[6]
threads = int(sys.argv[7])
solve_polytomies = sys.argv[8].lower() == 'true' # We need to set solve_polytomies to run oon nly trees with < 150 nodes in metient v0.1.3.5.3

torch.set_num_threads(threads)
torch.set_num_interop_threads(threads)


print_config = met.PrintConfig(visualize=True, verbose=True, k_best_trees=5)
weights = met.Weights.pancancer_genetic_uniform_weighting()  # Use default weights based on Metient author's suggestion for non-human data
met.evaluate_label_clone_tree(tree, metadata, weights, print_config, output_dir, patient, solve_polytomies=solve_polytomies)

# Sort through results to obtain samples from the solution space
with gzip.open(os.path.join(output_dir, f"{patient}_{primary}.pkl.gz"), "rb") as f:
    pckl = pickle.load(f)

# Obtain samples from the solution space
num_samples = len(pckl["node_labels"])
tissues = pckl["anatomical_sites"]
migration_graphs = np.empty(num_samples, dtype=dict)
losses = np.empty(num_samples)
for i in range(num_samples):
    V = pckl["node_labels"][i]
    parents = pckl['parents'][i]
    A = met.adjacency_matrix_from_parents(parents)
    G = met.migration_graph(V, A)
    migration_graphs[i] = pd.DataFrame(G.numpy(), columns=tissues, index=tissues).transpose().to_dict()
    losses[i] = float(pckl["losses"][i])

# # Output individual migration graph solutions
# with open(f"{output_dir}/{patient}_{primary}_migration_graphs.txt", "w") as file:
#     file.write(f"loss\tmigration_graph\n")
#     for l, g in zip(losses, migration_graphs):
#         file.write(f"{l}\t{g}\n")

# # Output top solution as if it has 1.0 probabilities
# # This was done to overcome the poor calibration of a probabilistic summary graph for Quinn CP4 with only 3 divergent metient outputs
# outputfile="/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.4.cass.metient_LL.topgraph.csv"
# with open(outputfile, "w") as file:
#     for source, targets_dict in migration_graphs[0].items():
#         for target, edge_count in targets_dict.items():
#             if edge_count > 0:
#                 for n in range(1, int(edge_count) + 1):
#                     migration = f"{source}_{target}_{n}"
#                     file.write(f"{migration},1.0\n")

# Get edge-wise probability graph from all migration graph solutions
all_graphs = [(loss, graph) for loss, graph in zip(losses, migration_graphs)]

# Convert solution losses to probabilities
max_loss = max(solution[0] for solution in all_graphs)
min_loss = min(solution[0] for solution in all_graphs)
min_max_denominator = max_loss - min_loss if max_loss != min_loss else 1.0  # Prevent division by zero when all losses are equal
all_graphs = [((loss - min_loss)/min_max_denominator, counts) for loss, counts in all_graphs]  # Min-max scale losses
temp = 0.5  # Temperature parameter for softmax (fixed to 0.5 from Metient authors suggestion)
prob_denominator = sum([math.exp(-loss/temp) for loss, counts in all_graphs])
all_graphs = [(math.exp(-loss/temp)/prob_denominator, counts) for loss, counts in all_graphs]   # Convert to probabilities by temperature-scaled softmax

# Build the consensus graph
metient_prob_graph = {}
prob = None
for loss, solution in all_graphs:
    prob = loss
    for source_tissue, targets_dict in solution.items():
        for target_tissue, edge_count in targets_dict.items():
            if edge_count > 0:
                for n in range(1, int(edge_count) + 1):
                    migration = f"{source_tissue}_{target_tissue}_{n}"
                    if migration not in metient_prob_graph:
                        metient_prob_graph[migration] = prob
                    else:
                        metient_prob_graph[migration] += prob

# Sort descending by probability
metient_prob_graph = dict(sorted(metient_prob_graph.items(), key=lambda item: item[1], reverse=True))

# Output the consensus graph
with open(outfile, "w") as f:
    for migration, probability in metient_prob_graph.items():
        f.write(f"{migration},{probability}\n")


