#!/usr/bin/env python

# FIXME: should have a check to make sure all trees have same set of taxa

import dendropy
import argparse
import sys
import math
import scipy.stats

from dendropy.calculate import treecompare

dmeas_def = "WRF"
desc = "Compute distances between trees. If sets are provided, all pairs are compared.  Mean and standard deviation of distance measure are reported."

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("file1", help="File containing one or more trees")
parser.add_argument("file2", help="File containing one or more trees")
parser.add_argument("-d", "--dmeas", default=dmeas_def,
                    help = "distance measure: WRF (weighted Robinson Foulds), URF (unweighted Robinson Foulds), EUC (Euclidean), or PWD (correlation of pairwise distances) (default %s)" % dmeas_def)
parser.add_argument("-f1", "--format1", default="newick", help = "format for file1, 'newick' or 'nexus' (default 'newick')")
parser.add_argument("-f2", "--format2", default="newick", help = "format for file2, 'newick' or 'nexus' (default 'newick')")
parser.add_argument("-q", "--quiet", action='store_true', help = "output mean distance only (useful in scripting)")
args = parser.parse_args()

if args.dmeas not in ["WRF", "URF", "EUC", "PWD"]:
    parser.print_usage()
    exit(1)

# for use in pairwise distance comparison; compute pearson correlation from distance matrices
def dmatrix_pearson(tree1, tree2):
    dmat1 = tree1.phylogenetic_distance_matrix()
    dmat2 = tree2.phylogenetic_distance_matrix()
    d1 = []
    d2 = []
    for idx1, taxon1 in enumerate(tree1.taxon_namespace):
        for taxon2 in tree1.taxon_namespace:
            d1.append(dmat1.patristic_distance(taxon1, taxon2))
            d2.append(dmat2.patristic_distance(taxon1, taxon2))

    return scipy.stats.pearsonr(d1, d2).statistic
     
# establish common namespace
tns = dendropy.TaxonNamespace()

trees1 = dendropy.TreeList(taxon_namespace=tns)
try:
    trees1.read_from_path(args.file1,schema=args.format1,rooting="force-rooted")
except Exception as e:
    print("Error reading %s" % args.file1)
    print(e)
    exit(1)

try:
    trees2 = dendropy.TreeList(taxon_namespace=tns)
    trees2.read_from_path(args.file2,schema=args.format2,rooting="force-rooted")
except Exception as e:
    print("Error reading %s" % args.file2)
    print(e)
    exit(1)

totd = 0
totd2 = 0
for tree1 in trees1:
    for tree2 in trees2:
        if args.dmeas == "WRF":
            d = treecompare.weighted_robinson_foulds_distance(tree1, tree2)
        elif args.dmeas == "URF":
            d = treecompare.symmetric_difference(tree1, tree2)
        elif args.dmeas == "EUC":
            d = treecompare.euclidean_distance(tree1, tree2)
        else:     # PWD; in this case need full distance matrices
            d = dmatrix_pearson(tree1, tree2)
            
        totd += d
        totd2 += d*d

N = len(trees1)*len(trees2)

if (args.quiet == False):
    print("%d tree(s) from %s, %d tree(s) from %s" %
        (len(trees1), args.file1, len(trees2), args.file2), end=": ")
    print("mean: %.2f, sd: %.2f" % (totd/N, math.sqrt(totd2/N - (totd/N)*(totd/N))))
else:
    print("%.2f" % (totd/N), end="\t")
