export SHELL=/usr/bin/bash

# Prevent Make from removing files it thinks it should clean up
#.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# Absolute paths
MAIN_DIR := /local/storage/no-backup/vine-benchmarks

# Relative paths
BIN := $(MAIN_DIR)/bin
PHAST_BIN := $(BIN)/phast/bin
VINE_BIN := $(BIN)/vine/bin
#VINE_BIN := /home/staklins/projects/vine_project/vine/bin
BEAST := $(BIN)/beast/bin/beast
CONTAINERS := $(MAIN_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
LAML_SIF := $(CONTAINERS)/laml/laml.sif
BEAM_TEMPLATE := $(MAIN_DIR)/crispr_sims/beam_template.xml

MCMC_ESS_THRESHOLD := 400

# for collapsing
PYTHON_SRC := $(MAIN_DIR)/python/src
BEAM_SUP := $(CONTAINERS)/beam_sup/beam_sup.sif

TREES := $(shell seq -f tree.%.0f.true.nwk 1 $(NSAMP))
INDELS := $(patsubst tree.%.true.nwk,tree.%.indels.csv,$(TREES))
LNLS := $(patsubst tree.%.true.nwk,tree.%.lnl,$(TREES))
TIMES := $(patsubst tree.%.true.nwk,tree.%.time,$(TREES))
EVALRF := $(patsubst tree.%.true.nwk,tree.%.rf,$(TREES))
LAML_TREES := $(patsubst tree.%.true.nwk,tree.%.laml_trees.nwk,$(TREES))
VINE_TREES := $(patsubst tree.%.true.nwk,tree.%.var.nwk,$(TREES))
BEAM_TERM := $(patsubst tree.%.true.nwk,tree.%.beam.term,$(TREES))
CHECK_LNLS := $(patsubst tree.%.true.nwk,tree.%.checklnl,$(TREES))

all: summary.time.txt summary.lnl.txt eval.all.rf.txt

simulate: $(TREES) $(INDELS)
laml: $(LAML_TREES)
vine: $(VINE_TREES)
beam: $(BEAM_TERM)
check: $(CHECK_LNLS) summary.checklnl.txt

vineclean:
	rm -f $(VINE_TREES) $(LNLS) $(TIMES) summary.time.txt summary.lnl.txt

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

# Collapse dups
#tree.%.collapsed.csv tree.%.collapsing_map.tsv: tree.%.indels.csv
#	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(BEAM_SUP) python3 /mnt/scripts/collapse_matrix_duplicates.py \
#	 -m /mnt/files/tree.$*.indels.csv \
#	 -om /mnt/files/tree.$*.collapsed.csv \
#	 -on /mnt/files/tree.$*.collapsing_map.tsv	

# Reformat the indel matrix to a tsv
tree.%.indels.tsv: tree.%.indels.csv
	$(BIN)/csv-to-tsv.sh $< > $@

# Same for collapsed version
#tree.%.collapsed.tsv: tree.%.collapsed.csv
#	$(BIN)/csv-to-tsv.sh $< > $@

# Run cassiopeia-greedy to get the starting tree for laml
tree.%.cass.nwk: tree.%.indels.csv
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/python/src/cassiopeiaGreedy.py /mnt/crispr_sims/$(NTAXA)taxa/tree.$*.indels.csv /mnt/crispr_sims/$(NTAXA)taxa/tree.$*.cass.nwk

# Run LAML
# Added --noultrametric to the LAML v1.0.5 run since it gives ultrametric trees. Seems like a bug they will need to fix.
tree.%.laml_trees.nwk tree.%.laml-time: tree.%.indels.csv tree.%.cass.nwk
	singularity exec --bind $(CURDIR):/mnt $(LAML_SIF) /usr/bin/time -o /mnt/tree.$*.laml-time \
	run_laml -c /mnt/tree.$*.indels.csv -t /mnt/tree.$*.cass.nwk -o /mnt/tree.$*.laml --topology_search --noDropout --noultrametric

# Run vine
tree.%.var.nwk tree.%.var.log tree.%.var-time: tree.%.indels.tsv
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine \
		$(VAROPT) \
		-i CRISPR tree.$*.indels.tsv \
		--logf tree.$*.var.log \
		--mean tree.$*.mean.nwk > tree.$*.var.nwk

# use the collapsed version instead
#tree.%.var.nwk tree.%.var.log tree.%.var-time: tree.%.collapsed.tsv 
#	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine \
#		$(VAROPT) \
#		-i CRISPR tree.$*.collapsed.tsv \
#		--log tree.$*.var.log \
#		--mean tree.$*.mean.nwk > tree.$*.var.nwk

# removing this for now
#		--tree tree.$*.cass.nwk \


# Run BEAM
tree.%.beam.trees tree.%.beam.log tree.%.beam.term: tree.%.indels.tsv tree.%.cass.nwk $(BEAM_TEMPLATE)
	cp $(BEAM_TEMPLATE) tree.$*.beam.xml; \
	$(BEAST) \
		-working \
		$$( [ -f tree.$*.beam.xml.state ] && echo -resume || echo -overwrite ) \
		-threads 1 \
		-seed 1 \
		-D outfileLog=tree.$*.beam.log \
		-D outfileTrees=tree.$*.beam.trees \
		-D barcodeMatrix=tree.$*.indels.tsv \
		-D startingNewickFile=tree.$*.cass.nwk \
		-D mcmcLength=$(MCMCLENGTH) \
		-D sampleFreq=$(MCMCSAMPLEFREQ) \
		tree.$*.beam.xml > tree.$*.beam.term

# Convert beam nexus output to newick
BURNIN := 0.1
tree.%.beam.nwk tree.%.beam.burninRemoved.nwk: tree.%.beam.trees
	$(BIN)/nex2nwk tree.$*.beam.trees tree.$*.beam.nwk
	numTrees=$$(wc -l < tree.$*.beam.nwk) ; \
	burnin=$$(awk -v n=$$numTrees -v b=$(BURNIN) 'BEGIN{printf "%d\n", n*b}') ; \
	tail -n +$$((burnin+1)) tree.$*.beam.nwk > tree.$*.beam.burninRemoved.nwk

# Extract lnls for all methods per tree
tree.%.lnl: tree.%.laml_params.txt tree.%.var.log tree.%.beam.log
	grep '^Negative' tree.$*.laml_params.txt | awk '{printf "%d\t-%f\t", $*, $$2}' > $@
	tail -1 tree.$*.var.log | awk '{printf "%f\t", $$19}' >> $@
	grep -v '^#' tree.$*.beam.log | grep -v '^Sample' | awk '{printf "%f\n", $$3}' | sort -nr | head -1 >> $@

# Extract times for all methods per tree
tree.%.time: tree.%.laml_trees.nwk tree.%.var-time tree.%.var.nwk tree.%.beam.term
	grep '^Runtime' tree.$*.laml.log | tail -1 | awk '{printf "%d\t%f\t", $*, $$3}' > $@
	head -1 tree.$*.var-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@
	beam=$$(awk '/\/Msamples/{ \
			n_last=$$1; x_last=$$NF; sub(/\/Msamples.*/,"",x_last); \
			if (($$4+0)>$(MCMC_ESS_THRESHOLD)) { n=n_last; x=x_last; found=1; exit } \
		} \
		END{ \
			if(!found){n=n_last; x=x_last} \
			h=m=s=0; \
			if(match(x,/([0-9]+)h/,a)) h=a[1]; \
			if(match(x,/([0-9]+)m/,a)) m=a[1]; \
			if(match(x,/([0-9]+)s/,a)) s=a[1]; \
			printf "%.3f", (n/1e6) * (3600*h + 60*m + s) \
		}' tree.$*.beam.term) ; \
	echo -e "\t$$beam" >> $@

# Summary of lnls and times for all methods
summary.lnl.txt: $(LNLS)
	cat $(LNLS) | awk 'BEGIN{printf "cp\tlaml\tvine\t\tbeam\tlaml_diff\tbeam_diff\n"} {x1 += $$2; x2 += $$3; x3 += $$4; printf "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$4, $$3-$$2, $$3-$$4} END {printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", x1/NR, x2/NR, x3/NR, (x2-x1)/NR, (x2-x3)/NR}' > $@

summary.time.txt: $(TIMES)
	cat $(TIMES) | awk 'BEGIN{printf "cp\tlaml\tvine\t\tbeam\tlaml_speedup\tbeam_speedup\n"} {x1 += $$2; x2 += $$3; x3 += $$4; printf "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$4, $$2/$$3, $$4/$$3} END {printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", x1/NR, x2/NR, x3/NR, x1/x2, x3/x2}' > $@

# evalTrees
tree.%.var.rf.txt: tree.%.var.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.var.nwk -t tree.$*.true.nwk > $@

tree.%.true.rf.txt: tree.%.true.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.true.nwk -t tree.$*.true.nwk > $@

tree.%.laml.rf.txt: tree.%.laml_trees.nwk tree.%.true.nwk
	sed 's/^\[&R\]//g' tree.$*.laml_trees.nwk > tmp.nwk
	$(VINE_BIN)/evalTrees tmp.nwk -t tree.$*.true.nwk > $@
	rm tmp.nwk

tree.%.beam.rf.txt: tree.%.beam.burninRemoved.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.beam.burninRemoved.nwk -t tree.$*.true.nwk > $@

tree.%.rf: tree.%.true.rf.txt tree.%.var.rf.txt tree.%.laml.rf.txt tree.%.beam.rf.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
	done

# Normalized RF: mean and std divided by (n-3), n = taxa from true tree
# (Newick: n_leaves = commas + 1)
eval.all.rf.txt: $(EVALRF)
	nwk="$(patsubst %.rf,%.true.nwk,$(firstword $(EVALRF)))"; \
	n=$$(grep -o ',' "$$nwk" | wc -l); n=$$((n + 1)); \
	div=$$((n - 3)); [ $$div -lt 1 ] && div=1; \
	echo "true (sd) vine (sd) laml (sd) beam (sd)" > tmp; \
	for file in $^ ; do \
		awk -v d="$$div" '{printf "%f\t%f\t", $$2/d, $$3/d}' \
		  "$$file" >> tmp; \
		echo >> tmp; \
	done; \
	awk '{x1+=$$1; x1s+=$$2*$$2; x2+=$$3; x2s+=$$4*$$4; \
	  x3+=$$5; x3s+=$$6*$$6; x4+=$$7; x4s+=$$8*$$8; print $$0} END { \
	  printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", \
	  x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), \
	  x3/(NR-1), sqrt(x3s/(NR-1)), x4/(NR-1), sqrt(x4s/(NR-1))}' \
	  tmp > $@; \
	rm -f tmp

tree.%.checklnl: tree.%.indels.tsv tree.%.laml_trees.nwk
	$(VINE_BIN)/crisprLnl tree.$*.indels.tsv tree.$*.laml_trees.nwk tree.$*.laml_params.txt > $@

summary.checklnl.txt: $(CHECK_LNLS)
	printf "%17s %15s %15s %15s %15s\n" "file" "VINElnl" "LAMLlnl" "diff" "%diff" > $@
	for file in $^ ; do \
		printf "%17s " $${file} >> $@ ;\
		awk '{printf "%15.4f ", $$NF} END {printf "\n"}' $${file} >> $@ ;\
	done

clean:
	rm -rf tree.*.* summary.*.txt eval.all.*.txt

clean_laml:
	rm tree.*.laml* eval.all.*.txt summary.*.txt tree.*.lnl tree.*.time tree.*.rf

clean_vine:
	rm tree.*.var* tree.*.mean* eval.all.*.txt summary.*.txt tree.*.lnl tree.*.time tree.*.rf

clean_beam:
	rm tree.*.beam* eval.all.*.txt summary.*.txt tree.*.lnl tree.*.time tree.*.rf

clean_not_converged_beam:
	for num in $$(seq 1 $(NSAMP)); do \
		if [ -f tree.$$num.beam.term ]; then \
			ess=$$(grep "Msamples" tree.$$num.beam.term | tail -n 1 | awk '{print $$4}'); \
			if [ -z "$$ess" ] || awk -v ess="$$ess" -v threshold="$(MCMC_ESS_THRESHOLD)" 'BEGIN{exit !(ess < threshold)}'; then \
				echo "Deleting tree.$$num (ESS=$$ess)"; \
				rm -f tree.$$num.beam.term; \
			else \
				echo "Keeping tree.$$num (ESS=$$ess)"; \
			fi; \
		fi; \
	done

archive_summary:
	archive_dir=archive_summary.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv summary.* $$archive_dir; \
	mv tree.*.time $$archive_dir; \
	mv tree.*.lnl $$archive_dir; \
	mv tree.*.rf $$archive_dir

archive_beam:
	archive_dir=archive_beam.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv summary.* $$archive_dir; \
	mv tree.*.time $$archive_dir; \
	mv tree.*.lnl $$archive_dir; \
	mv tree.*.rf $$archive_dir; \
	mv tree.*.beam* $$archive_dir

archive_laml:
	archive_dir=archive_laml.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv summary.* $$archive_dir; \
	mv tree.*.time $$archive_dir; \
	mv tree.*.lnl $$archive_dir; \
	mv tree.*.rf $$archive_dir; \
	mv tree.*.laml* $$archive_dir

archive_vine:
	archive_dir=archive_vine.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv eval.all.*.txt $$archive_dir; \
	mv summary.* $$archive_dir; \
	mv tree.*.time $$archive_dir; \
	mv tree.*.lnl $$archive_dir; \
	mv tree.*.rf $$archive_dir; \
	mv tree.*.var* $$archive_dir
