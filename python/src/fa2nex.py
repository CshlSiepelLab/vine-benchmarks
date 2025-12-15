
import argparse
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

parser = argparse.ArgumentParser(description="Convert FASTA to NEXUS format.")
parser.add_argument("fasta", help="Input FASTA file")
parser.add_argument("nexus", help="Output NEXUS file")
args = parser.parse_args()

fasta = args.fasta
nexus = args.nexus

# Read and annotate records
records = list(SeqIO.parse(fasta, "fasta"))
for record in records:
    record.annotations["molecule_type"] = "DNA"

# Build alignment
alignment = MultipleSeqAlignment(records)

# Write to NEXUS
AlignIO.write(alignment, nexus, "nexus")
