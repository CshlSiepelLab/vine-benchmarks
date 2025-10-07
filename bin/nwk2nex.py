from Bio import Phylo
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("newick_file")
parser.add_argument("nexus_file")
args = parser.parse_args()

Phylo.write(Phylo.read(args.newick_file, "newick"), args.nexus_file, "nexus")
