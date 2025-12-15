#!/usr/bin/env python

import dendropy
import argparse
import sys
from random import sample

desc = "Convert format of multi-tree file"

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("file", help="File containing one or more trees in Newick format, one per line")
parser.add_argument("-i", "--informat", default="newick", help="input file format, 'newick' or 'nexus' (default 'newick')")
parser.add_argument("-o", "--outformat", default="newick", help="output file format, 'newick' or 'nexus' (default 'newick')")
parser.add_argument("-s", "--sample", default=-1, help="sample specified number of trees from full set")
args = parser.parse_args()

#trees = dendropy.TaxonNamespace()
trees = dendropy.TreeList()

try:
    trees.read_from_path(args.file,schema=args.informat)
except Exception as e:
    print("Error reading %s" % args.file)
    print(e)
    exit(1)

if int(args.sample) > 0:
    sampidx = sample(range(0, len(trees) - 1), int(args.sample))
    treesamp = dendropy.TreeList([trees[i] for i in sampidx])
else:
    treesamp = trees.clone(depth=0)

treesamp.write_to_stream(sys.stdout, schema=args.outformat)

