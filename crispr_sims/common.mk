export SHELL=/usr/bin/bash

# Prevent Make from removing files it thinks it should clean up
.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# Absolute paths
MAIN_DIR := /local/storage/no-backup/vine-benchmarks

# Relative paths
BIN := $(MAIN_DIR)/bin
PHAST_BIN := $(BIN)/phast/bin
VINE_BIN := $(BIN)/vine/bin
CONTAINERS := $(MAIN_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
LAML_SIF := $(CONTAINERS)/laml/laml.sif

TREES := $(shell seq -f tree.%.0f.true.nwk 1 $(NSAMP))
INDELS := $(patsubst tree.%.true.nwk,tree.%.indels.csv,$(TREES))
LNLS := $(patsubst tree.%.true.nwk,tree.%.lnl,$(TREES))
TIMES := $(patsubst tree.%.true.nwk,tree.%.time,$(TREES))
EVALRF := $(patsubst tree.%.true.nwk,tree.%.rf,$(TREES))
LAML_TREES := $(patsubst tree.%.true.nwk,tree.%.laml_trees.nwk,$(TREES))

all: summary.time.txt summary.lnl.txt

simulate: $(TREES) $(INDELS)
laml: $(LAML_TREES)

tree.%.true.nwk:
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/python/src/simulateCellTree.py \
		--out_tree /mnt/crispr_sims/$(NTAXA)taxa/tree.$*.true.nwk \
		--num_tips $(NTAXA) \
		--birth_rate 0.075 \
		--death_rate 0.005 \
		--desired_time 54

# Run crispr barcode simulation to get the mutation overlay for the tree
tree.%.indels.csv: tree.%.true.nwk
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/python/src/simulateCrisprBarcodes.py \
		--in_tree /mnt/crispr_sims/$(NTAXA)taxa/tree.$*.true.nwk \
		--out_matrix /mnt/crispr_sims/$(NTAXA)taxa/tree.$*.indels.csv \
		--num_cassettes $(NCASETTES) \
		--cassette_size $(CASETTESIZE) \
		--mut_rate 0.01 \
		--heritable_silencing_rate 0.0001 \
		--stochastic_silencing_rate 0.0

# Reformat the indel matrix to a tsv
tree.%.indels.tsv: tree.%.indels.csv
	$(BIN)/csv-to-tsv.sh $< > $@

# Run cassiopeia-greedy to get the starting tree for laml and vine
tree.%.cass.nwk: tree.%.indels.csv
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/python/src/cassiopeiaGreedy.py /mnt/crispr_sims/$(NTAXA)taxa/tree.$*.indels.csv $@

# Run laml
tree.%.laml_trees.nwk: tree.%.indels.csv tree.%.cass.nwk
	singularity exec --bind $(CURDIR):/mnt $(LAML_SIF) \
	run_laml -c /mnt/tree.$*.indels.csv -t /mnt/tree.$*.cass.nwk -o /mnt/tree.$*.laml --topology_search --noDropout --noultrametric

# Run vine
tree.%.var.nwk tree.%.var.log tree.%.var-time: tree.%.indels.tsv tree.%.cass.nwk
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine \
		$(VAROPT) \
		-i CRISPR tree.$*.indels.tsv \
		--log tree.$*.var.log \
		--tree tree.$*.cass.nwk \
		--mean tree.$*.mean.nwk > tree.$*.var.nwk

# Extract lnls for all methods per tree
tree.%.lnl: tree.%.laml_trees.nwk tree.%.var.log
	grep '^Negative' tree.$*.laml_params.txt | awk '{printf "%d\t-%f\t", $*, $$2}' > $@
	tail -1 tree.$*.var.log | awk '{printf "%f\n", $$11}' >> $@

# Extract times for all methods per tree
tree.%.time: tree.%.laml_trees.nwk tree.%.var-time tree.%.var.nwk
	grep '^Runtime' tree.$*.laml.log | tail -1 | awk '{printf "%d\t%f\t", $*, $$3}' > $@
	head -1 tree.$*.var-time | awk '{printf "%s\n", $$1}' | sed 's/user//' >> $@

# Summary of lnls and times for all methods
summary.lnl.txt: $(LNLS)
	cat $(LNLS) | awk 'BEGIN{printf "cp\tlaml\tvine\tdiff\n"} {x1 += $$2; x2 += $$3; printf "%d\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$3-$$2} END {printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\n", x1/NR, x2/NR, (x2-x1)/NR}' > $@

summary.time.txt: $(TIMES)
	cat $(TIMES) | awk 'BEGIN{printf "cp\tlaml\tvine\tspeedup\n"} NF==3 {x1 += $$2; x2 += $$3; printf "%d\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$2/$$3} END {printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\n", x1/NR, x2/NR, x1/x2}' > $@

# evalTrees
tree.%.var.rf.txt: tree.%.var.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.var.nwk -t tree.$*.true.nwk > $@

tree.%.true.rf.txt: tree.%.true.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.true.nwk -t tree.$*.true.nwk > $@

tree.%.laml.rf.txt: tree.%.laml_trees.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.laml_trees.nwk -t tree.$*.true.nwk > $@

tree.%.rf: tree.%.true.rf.txt tree.%.var.rf.txt tree.%.laml.rf.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
	done

eval.all.rf.txt: $(EVALRF)
	echo "true (sd) vine (sd) laml (sd)" > tmp
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{x1 += $$1; x1s += ($$2 * $$2); x2 += $$3; x2s += ($$4 * $$4); x3 += $$5; x3s += ($$6 * $$6); print $$0} END {printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), x3/(NR-1), sqrt(x3s/(NR-1))}' tmp > $@
	rm -f tmp

clean:
	rm -rf tree.*.* summary.*.txt eval.all.*.txt

clean_laml:
	rm -rf tree.*.laml* eval.all.*.txt summary.*.txt tree.*.lnl tree.*.time tree.*.rf
