
import argparse
import cassiopeia as cas
# import numpy as np


parser = argparse.ArgumentParser(description="Simulate CRISPR barcodes and tree.")
parser.add_argument("--in_tree", type=str, required=True, help="Input path for the ground truth tree in newick format")
parser.add_argument("--out_matrix", type=str, required=True, help="Output path for character matrix")
# parser.add_argument("--out_tree", type=str, required=True, help="Output path for ground truth tree")
# parser.add_argument("--birth_rate", type=float, default=0.4, help="Birth rate for tree process")
# parser.add_argument("--death_rate", type=float, default=0.05, help="Death rate for tree process")
# parser.add_argument("--num_tips", type=int, default=50, help="Number of tips in the tree")
parser.add_argument("--num_sites", type=int, default=30, help="Number of barcode sites")
parser.add_argument("--mut_rate", type=float, default=0.1, help="Mutation rate for barcode matrix process")
parser.add_argument("--heritable_silencing_rate", type=float, default=0.0001, help="Heritable silencing rate")
parser.add_argument("--stochastic_silencing_rate", type=float, default=0.01, help="Stochastic silencing rate")

args = parser.parse_args()

in_tree = args.in_tree
out_matrix = args.out_matrix
# out_tree = args.out_tree
# birth_rate = args.birth_rate
# death_rate = args.death_rate
# num_tips = args.num_tips
num_sites = args.num_sites
mut_rate = args.mut_rate
heritable_silencing_rate = args.heritable_silencing_rate
stochastic_silencing_rate = args.stochastic_silencing_rate

# # Commenting this out in favor of using a fixed input tree instead of simulating the tree here
# # Simulate tree using birth-death process
# simulator = cas.sim.BirthDeathFitnessSimulator(
#     birth_waiting_distribution=lambda scale: np.random.exponential(scale),
#     initial_birth_scale=1 / birth_rate,
#     death_waiting_distribution=lambda: np.random.exponential(1 / death_rate),
#     num_extant=num_tips,
# )
# ground_truth_tree = simulator.simulate_tree()

# Read in fixed input tree
ground_truth_tree = cas.data.CassiopeiaTree()
ground_truth_tree.populate_tree(tree=in_tree)

# Indel outcome priors per site
# Indel outcomes are the same across all sites
# We specify fixed values here so we can use the same values as fixed edit rates during inference
# Alternatively, it is possible to sample these within the Cas9LineageTracingDataSimulator from scratch each time,
# but it is more challenging to extract those values to ensure that the same values are used during inference
state_priors = {
    1: 0.001891385794602738,
    2: 0.0024762237066947427,
    3: 0.006434777250119938,
    4: 0.022626187451291142,
    5: 0.0004519611853084519,
    6: 0.004252307867693467,
    7: 0.008695199487743388,
    8: 0.005287228198826541,
    9: 0.00852515654186136,
    10: 0.0028375427438811205,
    11: 0.01630715871347483,
    12: 0.004283226603781587,
    13: 0.0067464292466787685,
    14: 0.016033678821811725,
    15: 0.01917491749069289,
    16: 0.0036585141765295664,
    17: 0.004891872021313965,
    18: 0.00867199597266662,
    19: 0.011140999534685158,
    20: 0.010581282853573013,
    21: 0.004541688131273875,
    22: 0.0007617862007893749,
    23: 0.036165552581349764,
    24: 0.007900631443893674,
    25: 0.00019973184290800428,
    26: 0.009877241632563307,
    27: 0.005210678254720317,
    28: 0.004181565880363428,
    29: 0.017589109966820247,
    30: 0.0017611954710861511,
    31: 0.00022959343560415528,
    32: 0.002994537359426538,
    33: 0.000801068502385592,
    34: 0.014061233100234731,
    35: 0.015163857946017672,
    36: 0.04639414188685866,
    37: 0.001504774488150319,
    38: 0.014571829241040163,
    39: 0.00011550294775905987,
    40: 0.0017412028611578365,
    41: 0.0019912235849203347,
    42: 0.015435544626458327,
    43: 0.0026174926141904324,
    44: 0.03034937733094538,
    45: 0.005414023081912718,
    46: 0.015567352766629855,
    47: 0.0018764294412664411,
    48: 0.03822133912883424,
    49: 0.005016017821502824,
    50: 0.017923582049993315,
    51: 0.00707709038325807,
    52: 0.01894418099476499,
    53: 0.006552348025361078,
    54: 0.011483251612581855,
    55: 0.00213045546278859,
    56: 0.00044723689376124275,
    57: 0.028216265577151605,
    58: 0.007500191290250902,
    59: 9.64495624736785e-06,
    60: 0.014456941979911852,
    61: 0.0015131231024375432,
    62: 0.04034276177940931,
    63: 0.007743384492752599,
    64: 0.0012231748469220054,
    65: 0.018702920184101884,
    66: 0.0008980883166958691,
    67: 0.005709122901667588,
    68: 0.014609172657832641,
    69: 0.008227991029964507,
    70: 0.00027099882320938255,
    71: 0.006800810951957274,
    72: 0.0016024909826768872,
    73: 0.002702309545993545,
    74: 0.013068054552918141,
    75: 0.021222580867047328,
    76: 0.006778296749981686,
    77: 0.005021657046780748,
    78: 0.03394700190346978,
    79: 0.01064798652834624,
    80: 0.0015709747044370804,
    81: 0.006723215374919598,
    82: 0.00918494570764104,
    83: 0.0007148514517145434,
    84: 0.033386853032898484,
    85: 0.004952562215118869,
    86: 0.0016703080977210618,
    87: 0.008780194621411974,
    88: 0.0021936132072342554,
    89: 0.0037054891861055513,
    90: 0.013720737264473232,
    91: 0.02129789157325353,
    92: 0.013693745961487449,
    93: 0.018491991117775546,
    94: 0.0011056635477868833,
    95: 0.0008658270062053632,
    96: 0.008272143321129521,
    97: 0.001439189221306065,
    98: 0.011779293293448219,
    99: 0.032353702218426296,
    100: 0.021099922150975396,
}

# Simulate Cas9 lineage tracing data on the ground truth tree
# Docs at https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.sim.Cas9LineageTracingDataSimulator.html
lt_sim = cas.sim.Cas9LineageTracingDataSimulator(
    number_of_cassettes=num_sites,
    size_of_cassette=1,
    mutation_rate=mut_rate,
    state_priors=state_priors,
    heritable_silencing_rate=heritable_silencing_rate,
    stochastic_silencing_rate=stochastic_silencing_rate,
)
lt_sim.overlay_data(ground_truth_tree)

# Write the crispr barcode cells x sites character matrix to file
ground_truth_tree.character_matrix.to_csv(out_matrix, sep=",")

# # Write ground truth tree to file
# nwk = ground_truth_tree.get_newick(record_branch_lengths=True, record_node_names=False)
# with open(out_tree, "w") as f:
#     f.write(nwk)
