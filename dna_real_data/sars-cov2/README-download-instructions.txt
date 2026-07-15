# steps for downloading nextstrain data and setting up subsampled versions
# raw files are quite large so these commands are kept separate from makefile

# install augur and seqkit
# e.g., pip install nextstrain-augur, brew install seqkit

# download data from nextstrain 

# these were the original commands used
#curl -L -O https://data.nextstrain.org/files/ncov/open/100k/sequences.fasta.xz
#curl -L -O https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz
#curl -L -O https://data.nextstrain.org/files/ncov/open/100k/metadata.tsv.xz

# however, the archived versions from Feb 25, 2026 appear to be available here
curl -L -o sequences.fasta.xz https://nextstrain-data.s3.amazonaws.com/files/ncov/open/100k/sequences.fasta.xz?versionId=VE0xefHcGH31NBnekZFT8Mw1ILlWwCY
curl -L -o aligned.fasta.xz https://nextstrain-data.s3.amazonaws.com/files/ncov/open/aligned.fasta.xz?versionId=10lRec1ySVWc2mLoPd1b7LOb5OLbas_i
curl -L -o metadata.tsv.xz https://nextstrain-data.s3.amazonaws.com/files/ncov/open/100k/metadata.tsv.xz?versionId=qewluWuTeTNFkoULKvma8MY0.2AlstLe

# decompress
xz -d sequences.fasta.xz
xz -d aligned.fasta.xz
xz -d metadata.tsv.xz

# filter to obtain large and small data sets
augur filter   --metadata metadata.tsv   --sequences sequences.fasta --output-sequences filtered.fasta   --output-metadata filtered.tsv --min-length 29000   --exclude-ambiguous-dates-by any
augur filter   --metadata filtered.tsv   --sequences aligned.fasta --group-by region year month   --sequences-per-group 3   --subsample-seed 1 --output-sequences large-subset.fasta   --output-metadata large-subset.tsv
augur filter   --metadata large-subset.tsv   --sequences large-subset.fasta   --group-by region year month   --sequences-per-group 1   --subsample-seed 2   --output-sequences small-subset.fasta   --output-metadata small-subset.tsv

# deleted original sequences and metadata, and compressed filtered.fasta and filtered.tsv to save space

