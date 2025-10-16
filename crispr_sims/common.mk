export SHELL=/usr/bin/bash

# Absolute paths
MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT_SUFFIX := crispr_sims

# Relative paths
ROOT := $(MAIN_DIR)/$(ROOT_SUFFIX)
PHAST_BIN := $(MAIN_DIR)/phast/bin
OTHER_BIN := $(MAIN_DIR)/bin

all: tree.1.true.nwk

# Leave the DNA tree birth/death simulation parameters alone for now
tree.%.true.nwk: 
	$(OTHER_BIN)/bdTree3 -b 1 -d 0.5 --oversample-k 3 --height 5 --min-edge 0.02 --expected-height $(EXPHEIGHT) --no-stem --ucln-sd 0.6 --target-stat median -n $(NTAXA) | sed 's/\[\&[UR]\] //' > $@

# Run crispr barcode simulation to get the mutation overlay for the tree
# Need to change the python command to a binary executable with cassiopeia installed to remove the local dependency
tree.%.indels.csv: tree.%.true.nwk
	python $(OTHER_BIN)/simulateCrisprBarcodes.py \
		--in_tree $< \
		--out_matrix $@ \
		--num_sites $(NSITES) \
		--mut_rate $(MUTRATE) \
		--heritable_silencing_rate $(HERITABLESILENCINGRATE) \
		--stochastic_silencing_rate $(STOCHASICSILENCINGRATE)

# Run cassiopeia-greedy to get the starting tree for laml and vine
tree.%.cass.nwk: tree.%.indels.csv
	python $(OTHER_BIN)/cassiopeiaGreedy.py $< $@

# Run laml
tree.%.laml_trees.nwk: tree.%.indels.csv tree.%.cass.nwk
	run_laml -c tree.$*.indels.csv -t tree.$*.cass.nwk -o tree.$*.laml --topology_search --noDropout

# Reformat the indel matrix to a tsv
tree.%.indels.tsv: tree.%.indels.csv
	$(OTHER_BIN)/csv-to-tsv.sh $< > $@

# Run vine
tree.%.var.nwk tree.%.var.log tree.%.var-time: tree.%.indels.tsv tree.%.cass.nwk 
	/usr/bin/time -o tree.$*.var-time $(PHAST_BIN)/vine \
		$(VAROPT) \
		-i CRISPR tree.$*.indels.tsv \
		--log tree.$*.var.log \
		--tree tree.$*.cass.nwk \
		--mean tree.$*.mean.nwk \
		> $@ || { rm -f tree.$*.var.* tree.$*.var-time tree.$*.mean.nwk ; exit 0; }

# Extract lnls for all methods per tree
tree.%.lnl: tree.%.laml_trees.nwk tree.%.var.log
	grep '^Negative' tree.$*.laml_params.txt | awk '{printf "%d\t-%f\t", $*, $$2}' > $@
	tail -1 tree.$*.var.log | awk '{printf "%f\n", $$11}' >> $@

# Extract times for all methods per tree
tree.%.time: tree.%.laml_trees.nwk tree.%.var-time
	grep '^Runtime' tree.$*.laml.log | tail -1 | awk '{printf "%d\t%f\t", $*, $$3}' > $@
	head -1 tree.$*.var-time | awk '{printf "%s\n", $$1}' | sed 's/user//' >> $@

clean:
	rm -rf tree.*.true.nwk tree.*.laml* tree.*.var* tree.*.indels.* tree.*.cass.nwk tree.*.lnl tree.*.time
