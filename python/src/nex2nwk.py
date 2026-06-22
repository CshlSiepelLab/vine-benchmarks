from Bio import Phylo
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("nexus_file")
parser.add_argument("newick_file")
args = parser.parse_args()

trees = Phylo.parse(args.nexus_file, "nexus")
Phylo.write(trees, args.newick_file, "newick")
