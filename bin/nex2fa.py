
import argparse
from Bio import AlignIO
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Convert NEXUS to FASTA format.")
parser.add_argument("nexus", help="Input NEXUS file")
parser.add_argument("fasta", help="Output FASTA file")
args = parser.parse_args()

nexus_file = args.nexus
fasta_file = args.fasta

# Read the NEXUS alignment
alignment = AlignIO.read(nexus_file, "nexus")

# Write each sequence to FASTA
SeqIO.write(alignment, fasta_file, "fasta")
