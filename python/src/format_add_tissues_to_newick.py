
import os, sys

from graphposterior.tree_utils import annotate_tree_with_tissues


nwk_file = sys.argv[1]
tissues_file = sys.argv[2]
out_file = sys.argv[3]

annotate_tree_with_tissues(nwk_file, tissues_file, out_file)