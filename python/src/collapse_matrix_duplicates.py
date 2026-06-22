
import argparse
import pandas as pd
from graphposterior.matrix_utils import collapse_character_matrix, expand_clones_with_multiple_tissues

parser = argparse.ArgumentParser(description="Collapse duplicate rows in a character matrix.")
parser.add_argument("-m", help="Path to the input character matrix CSV.")
parser.add_argument("-om", help="Path to write the collapsed character matrix CSV.")
parser.add_argument("-on", help="Path to write the name mapping TSV.")
parser.add_argument("-t", default=None, help="Optional path to the tissues CSV.")
parser.add_argument("-ot", default=None, help="Optional path to write the collapsed tissues TSV.")

args = parser.parse_args()
matrix_csv_file = args.m
tissues_csv_file = args.t
output_matrix_csv_file = args.om
output_name_mapping_tsv_file = args.on
output_tissues_tsv_file = args.ot

# Read in matrix
char_matrix_df = pd.read_csv(matrix_csv_file, index_col=0)

# Make index names string type in case they are numeric to avoid issues with mapping later
char_matrix_df.index = char_matrix_df.index.astype(str)

# Remove NaN columns
char_matrix_df = char_matrix_df.dropna(axis=1, how='all')

# Optionally read in tissues
if tissues_csv_file:
    with open (tissues_csv_file, 'r') as f:
        tissues_dict = {}
        for line in f:
            parts = line.strip().split(',')
            clone_id = str(parts[0])
            tissue = str(parts[1])
            tissues_dict[clone_id] = tissue
else:
    tissues_dict = None

collapsed_matrix, maps = collapse_character_matrix(char_matrix_df, tissues_dict)

# Write out matrix and other files optionally
collapsed_matrix.to_csv(output_matrix_csv_file)

if output_name_mapping_tsv_file:
    map_original_barcodes = maps['group_to_originals']
    with open(output_name_mapping_tsv_file, 'w') as f:
        for group_id, original_barcodes in map_original_barcodes.items():
            f.write(f"{group_id}\t{original_barcodes}\n")
            
if output_tissues_tsv_file and tissues_dict:
    group_to_tissues = maps['group_to_tissues']
    with open(output_tissues_tsv_file, 'w') as f:
        f.write("group_name\ttissues\n")
        for group_id, tissues in group_to_tissues.items():
            f.write(f"{group_id}\t{tissues}\n")
            
# Expand clones with multiple tissues optionally
if tissues_csv_file and output_tissues_tsv_file:
    tissues_df = pd.DataFrame(maps['group_to_tissues'].items(), columns=['group_name', 'tissues'])
    expanded_matrix_df, expanded_tissues_df = expand_clones_with_multiple_tissues(collapsed_matrix, tissues_df)
    expanded_matrix_df.to_csv(output_matrix_csv_file.replace('.csv', '.expanded.csv'))
    expanded_tissues_df.to_csv(output_tissues_tsv_file.replace('.tsv', '.expanded.csv'), index=False, header=False)
