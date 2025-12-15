#!/usr/bin/env python

from dendropy.simulate.treesim import birth_death_tree
import argparse

# defaults
brate_def = 1
drate_def = 0.5
ntips_def = 100

desc = "Generate a Newick-formatted tree using a birth-death process"

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-b", "--brate", default=brate_def, help = "birthrate (default %f)" % brate_def)
parser.add_argument("-d", "--drate", default=drate_def, help = "deathrate (default %f)" % drate_def)
parser.add_argument("-n", "--ntips", default=ntips_def, help = "number of tips (default %d)" % ntips_def)
args = parser.parse_args()
    
# Simulate a starting tree
tree = birth_death_tree(birth_rate=float(args.brate), death_rate=float(args.drate), num_extant_tips=int(args.ntips))

# force it to be unrooted
tree.is_rooted = False
tree.update_bipartitions()

# Convert the dendropy tree to a Newick string and print to stdout
newick_str = tree.as_string(schema="newick")
print(newick_str)

# remove the [%R]
# add more parens?
# root the tree?
