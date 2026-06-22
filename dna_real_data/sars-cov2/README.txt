steps for setting up subsampled data set

# install augur and seqkit
# pip install nextstrain-augur
# brew install seqkit
curl -L -O https://data.nextstrain.org/files/ncov/open/100k/sequences.fasta.xz
curl -L -O https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz
curl -L -O https://data.nextstrain.org/files/ncov/open/100k/metadata.tsv.xz
xz -d sequences.fasta.xz
xz -d aligned.fasta.xz
xz -d metadata.tsv.xz
augur filter   --metadata metadata.tsv   --sequences sequences.fasta --output-sequences filtered.fasta   --output-metadata filtered.tsv --min-length 29000   --exclude-ambiguous-dates-by any
# deleted original sequences and metadata to save space
augur filter   --metadata filtered.tsv   --sequences aligned.fasta --group-by region year month   --sequences-per-group 3   --subsample-seed 1 --output-sequences large-subset.fasta   --output-metadata large-subset.tsv
augur filter   --metadata large-subset.tsv   --sequences large-subset.fasta   --group-by region year month   --sequences-per-group 1   --subsample-seed 2   --output-sequences small-subset.fasta   --output-metadata small-subset.tsv
# compressed filtered.fasta and filtered.tsv to save space

# pilot run
vine --hky vine_250.fasta --logf vine-small.log --mean vine-small.mean.nwk -s 1000 -j 8 > vine-small.nwk

# also try --dgamma 4 --gtr

