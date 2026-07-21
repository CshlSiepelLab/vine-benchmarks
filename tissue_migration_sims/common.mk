export SHELL=/usr/bin/bash

# Give beast2 access to beagle
export LD_LIBRARY_PATH := /local/storage/no-backup/vine-benchmarks/lib:$(LD_LIBRARY_PATH)

# Prevent Make from removing files it thinks it should clean up
#.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS: tree.%.beam.trees tree.%.beam.log tree.%.xml tree.%.xml.state

# Path to BEAM template
BEAM_TEMPLATE := /local/storage/no-backup/vine-benchmarks/tissue_migration_sims/beam_template.xml

BEAM_ESS_TARGET := 400

# Absolute paths
MAIN_DIR := /local/storage/no-backup/vine-benchmarks

# Relative paths
PYTHON_SRC := $(MAIN_DIR)/python/src
R_SRC := $(MAIN_DIR)/r/src
BIN := $(MAIN_DIR)/bin
# VINE_BIN := $(BIN)/vine/bin
VINE_BIN := /home/staklins/projects/vine_project/vine/bin
BEAST := $(BIN)/beast/bin/beast
TREEANNOTATOR := $(BIN)/beast/bin/treeannotator
CONTAINERS := $(MAIN_DIR)/containers
GRAPHPOSTERIOR_SIF := $(CONTAINERS)/graphposterior/graphposterior.sif
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
METIENT_SIF := $(CONTAINERS)/metient/metient.sif
MACH2_SIF := $(CONTAINERS)/mach2/mach2.sif
LAML_SIF := $(CONTAINERS)/laml/laml.sif
RPLOTTING_SIF := $(CONTAINERS)/rplotting/rplotting.sif

# VAROPT?=-s 1000 -v 3 --cov DIST --planar-flow --radial-flow
VAROPT?=-s 10000 -v 0

TREE_IDS := $(shell seq 1 $(NSAMP))
TREES := $(addsuffix _cell_tree.nwk,$(addprefix tree.,$(TREE_IDS)))
MATRIXCSVS := $(addsuffix _indel_character_matrix.csv,$(addprefix tree.,$(TREE_IDS)))
TISSUECSVS := $(addsuffix _cell_tree.labeling.csv,$(addprefix tree.,$(TREE_IDS)))
TRUEGRAPHPLOTS := $(addsuffix _cell_tree.migrations.pdf,$(addprefix tree.,$(TREE_IDS)))
CASTREES := $(addsuffix .cass.nwk,$(addprefix tree.,$(TREE_IDS)))
LAMLTREES := $(addsuffix .laml_trees.nwk,$(addprefix tree.,$(TREE_IDS)))
VINELOGS := $(addsuffix .var.log,$(addprefix tree.,$(TREE_IDS)))
VINEPROBGRAPHS := $(addsuffix .var_probability_graph.csv,$(addprefix tree.,$(TREE_IDS)))
BEAMXMLS := $(addsuffix .beam.xml,$(addprefix tree.,$(TREE_IDS)))
BEAMTERMS := $(addsuffix .beam.term,$(addprefix tree.,$(TREE_IDS)))
BEAMPROBGRAPHS := $(addsuffix .beam_probability_graph.csv,$(addprefix tree.,$(TREE_IDS)))
VINEPRFS := $(addsuffix .var.precision_recall.csv,$(addprefix tree.,$(TREE_IDS)))
BEAMPRFS := $(addsuffix .beam.precision_recall.csv,$(addprefix tree.,$(TREE_IDS)))
METIENTMETA := $(addsuffix .metient.metadata.txt,$(addprefix tree.,$(TREE_IDS)))
METIENTPRFS := $(addsuffix .metient.precision_recall.csv,$(addprefix tree.,$(TREE_IDS)))
MACH2PRFS := $(addsuffix .mach2.precision_recall.csv,$(addprefix tree.,$(TREE_IDS)))
TIMES := $(addsuffix .time,$(addprefix tree.,$(TREE_IDS)))
50PRFS := $(addsuffix .50_prf.txt,$(addprefix tree.,$(TREE_IDS)))
LNLS := $(addsuffix .lnl,$(addprefix tree.,$(TREE_IDS)))
LAMLTANGLEGRAMS := $(addsuffix .laml.tanglegram.pdf,$(addprefix tree.,$(TREE_IDS)))
VINETANGLEGRAMS := $(addsuffix .var.tanglegram.pdf,$(addprefix tree.,$(TREE_IDS)))
BEAMTANGLEGRAMS := $(addsuffix .beam.tanglegram.pdf,$(addprefix tree.,$(TREE_IDS)))

# LAML, Metient, MACH2, and BEAM are intentionally omitted for the largest datasets.
# Keep the summary schemas stable by writing NA for their fixed columns below.
AUX_METHODS_ENABLED := $(if $(filter 500 1000,$(NTAXA)),,1)
BEAM_ENABLED := $(if $(filter 500 1000,$(NTAXA)),,1)
COMMON_MK := $(lastword $(MAKEFILE_LIST))
LNL_AUX_DEPS := $(if $(AUX_METHODS_ENABLED),tree.%.laml_params.txt)
LNL_BEAM_DEPS := $(if $(BEAM_ENABLED),tree.%.beam.log)
TIME_AUX_DEPS := $(if $(AUX_METHODS_ENABLED),tree.%.laml-time tree.%.metient-time tree.%.mach2-time)
TIME_BEAM_DEPS := $(if $(BEAM_ENABLED),tree.%.beam.term)
PRF_AUX_DEPS := $(if $(AUX_METHODS_ENABLED),tree.%.metient.precision_recall.csv tree.%.mach2.precision_recall.csv)
PRF_BEAM_DEPS := $(if $(BEAM_ENABLED),tree.%.beam.precision_recall.csv)
THRESH_AUX_PRFS := $(if $(AUX_METHODS_ENABLED),$(METIENTPRFS) $(MACH2PRFS))
THRESH_BEAM_PRFS := $(if $(BEAM_ENABLED),$(BEAMPRFS))

all: $(TRUEGRAPHPLOTS) eval.all.time.pdf eval.all.50_prf.pdf eval.all.thresh_prf.pdf eval.all.lnl.pdf
simulate: $(TREES) $(MATRIXCSVS) $(TISSUECSVS) $(TRUEGRAPHPLOTS)
cass: $(CASTREES)
ifeq ($(AUX_METHODS_ENABLED),1)
laml: $(LAMLTREES)
prep_metient: $(METIENTMETA)
metient: $(METIENTPRFS)
mach2: $(MACH2PRFS)
else
laml prep_metient metient mach2:
	@echo "$@ is disabled for NTAXA=$(NTAXA)"
endif
vine: $(VINEPROBGRAPHS)
ifeq ($(BEAM_ENABLED),1)
prep_beam: $(BEAMXMLS)
beam: $(BEAMTERMS) $(BEAMPROBGRAPHS)
else
prep_beam beam:
	@echo "$@ is disabled for NTAXA=$(NTAXA)"
endif

# Simulate data
tree.%_cell_tree.nwk tree.%_cell_tree.labeling tree.%_cell_tree.vertex.labeling tree.%_cell_tree.migrations tree.%_indel_character_matrix.tsv:
	singularity exec --bind ./:/mnt $(GRAPHPOSTERIOR_SIF) \
	run_met_cancer_barcode_simulation \
		--outdir /mnt/ \
		--outprefix tree.$* \
		--num_generations 250 \
		--migration_rate 1e-6 \
		--num_cells_downsample $(NTAXA) \
		--num_possible_tissues 10 \
		--num_sites 3 \
		--num_barcodes 15 \
		--mutationrate 0.0025 \
		--heritable_silencing_rate 0.0001 \
		--stochastic_silencing_rate 0.01

tree.%_cell_tree.tissue_labels.nwk: tree.%_cell_tree.nwk tree.%_cell_tree.vertex.labeling
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) \
		python3 /mnt/scripts/format_add_tissues_to_newick.py /mnt/files/tree.$*_cell_tree.nwk /mnt/files/tree.$*_cell_tree.vertex.labeling /mnt/files/$@

tree.%_cell_tree.migrations.pdf: tree.%_cell_tree.tissue_labels.nwk
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) \
		python3 /mnt/scripts/plot_migration_graph_from_newick.py /mnt/files/$< P /mnt/files/$@

tree.%_indel_character_matrix.csv: tree.%_indel_character_matrix.tsv
	sed -e 's/\t/,/g' $< > $@

tree.%_cell_tree.labeling.csv: tree.%_cell_tree.labeling
	sed -e 's/ /,/g' $< > $@

tree.%.cass.nwk: tree.%_indel_character_matrix.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(CASSIOPEIA_SIF) \
	python /mnt/scripts/cassiopeiaGreedy.py /mnt/files/$< /mnt/files/$@

# Run VINE in tissue migration mode
tree.%.var.log tree.%.var-time tree.%.var.nwk tree.%.var.dot tree.%.var.nex: tree.%_indel_character_matrix.tsv tree.%_cell_tree.labeling.csv
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine $(VAROPT) -i CRISPR tree.$*_indel_character_matrix.tsv --logf tree.$*.var.log --mean tree.$*.var.mean.nwk \
	--migration tree.$*_cell_tree.labeling.csv --migration-rate-prior 0 --primary P --sample-graphs tree.$*.var.dot --labeled-trees tree.$*.var.nex > tree.$*.var.nwk

# Get MCC tree from VINE output
tree.%.var.mcc.nex tree.%.var.mcc.nwk: tree.%.var.nex
	$(TREEANNOTATOR) -lowMem true -height keep -burnin 0 $^ | \
    sed 's/\[[^]]*\]//g' > tree.$*.var.mcc.nex
	$(BIN)/nex2nwk tree.$*.var.mcc.nex tree.$*.var.mcc.nwk

# Plot VINE vs true tree tanglegram
# Seemingly need to remove the origin node here to prevent the tree read in step hanging for a long time
tree.%.var.tanglegram.pdf: tree.%_cell_tree.nwk tree.%_cell_tree.labeling.csv tree.%.var.mcc.nwk
	sed 's/)1:0;/);/' tree.$*_cell_tree.nwk > tree.$*_cell_tree.nwk.tmp
	sed 's/):0.00000;/);/' tree.$*.var.mcc.nwk > tree.$*.var.mcc.nwk.tmp
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_tanglegram.R \
		/mnt/files/tree.$*_cell_tree.nwk.tmp "True" /mnt/files/tree.$*.var.mcc.nwk.tmp "VINE" /mnt/files/tree.$*_cell_tree.labeling.csv /mnt/files/tree.$*.var.tanglegram.pdf "P"
	rm tree.$*_cell_tree.nwk.tmp tree.$*.var.mcc.nwk.tmp

# Process VINE results to get an edge-wise consensus probability graph
tree.%.var_probability_graph.csv: tree.%.var.nex
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/beam_postprocessing.py \
		/mnt/files/tree.$*.var.nex state P 1 0.0 /mnt/files/tree.$*.var

# Get VINE precision, recall, and F1 scores from edge-wise probability graph
tree.%.var.precision_recall.csv: tree.%.var_probability_graph.csv tree.%_cell_tree.migrations
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/get_performance_metrics.py \
		/mnt/files/tree.$*_cell_tree.migrations \
		/mnt/files/tree.$*.var_probability_graph.csv \
		/mnt/files/tree.$*.var.precision_recall.csv

# Format BEAM XML input file
tree.%.beam.xml: tree.%_cell_tree.labeling.csv $(BEAM_TEMPLATE)
	cp $(BEAM_TEMPLATE) $@
	allTissues=$$(cut -d, -f2 tree.$*_cell_tree.labeling.csv | sort | uniq | tr '\n' ' ' | sed 's/ *$$//'); \
	numAllTissues=$$(echo $${allTissues} | wc -w); \
	if ! echo $${allTissues} | grep -q 'P'; then \
		numAllTissues=$$((numAllTissues + 1)); \
	fi; \
	numTissueRates=$$(( ((numAllTissues * numAllTissues) - numAllTissues) / 2 )); \
	sed -i "s|numTissueRates|$${numTissueRates}|g" $@; \
	numOtherTissues=$$(cut -d, -f2 tree.$*_cell_tree.labeling.csv | sort | uniq | sed 's/P//' | grep -v '^$$' | wc -l); \
	rootTissueFreqs="1"; \
	for i in $$(seq 1 $$numOtherTissues); do \
		rootTissueFreqs="$${rootTissueFreqs} 0"; \
	done; \
	sed -i "s|rootTissueFreqs|$${rootTissueFreqs}|g" $@

# Run BEAM - requires some tuning here of the mcmc length and then sampling frequency
# The seed is fixed to permit resuming runs that don't converge
tree.%.beam.trees tree.%.beam.log tree.%.beam.term: tree.%_indel_character_matrix.tsv tree.%_cell_tree.labeling.csv tree.%.beam.xml tree.%.cass.nwk
	$(BEAST) \
		-working \
		$$( [ -f tree.$*.beam.xml.state ] && echo -resume || echo -overwrite ) \
		-threads 1 \
		-seed 1 \
		-D outfileLog=tree.$*.beam.log \
		-D outfileTrees=tree.$*.beam.trees \
		-D barcodeMatrix=tree.$*_indel_character_matrix.tsv \
		-D tissueCsv=tree.$*_cell_tree.labeling.csv \
		-D startingNewickFile=tree.$*.cass.nwk \
		-D mcmcLength=$(MCMCLENGTH) \
		-D sampleFreq=$(MCMCSAMPLEFREQ) \
		tree.$*.beam.xml > tree.$*.beam.term

# Get MCC tree from BEAM output
tree.%.beam.mcc.nex tree.%.beam.mcc.nwk: tree.%.beam.trees
	$(TREEANNOTATOR) -lowMem true -height keep -burnin 10 $^ | \
    sed 's/\[[^]]*\]//g' > tree.$*.beam.mcc.nex
	$(BIN)/nex2nwk tree.$*.beam.mcc.nex tree.$*.beam.mcc.nwk

# Plot BEAM vs true tree tanglegram
# Seemingly need to remove the origin node here to prevent the tree read in step hanging for a long time
tree.%.beam.tanglegram.pdf: tree.%_cell_tree.nwk tree.%_cell_tree.labeling.csv tree.%.beam.mcc.nwk
	sed 's/)1:0;/);/' tree.$*_cell_tree.nwk > tree.$*_cell_tree.nwk.tmp2
	sed 's/):0.00000;/);/' tree.$*.beam.mcc.nwk > tree.$*.beam.mcc.nwk.tmp
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_tanglegram.R \
		/mnt/files/tree.$*_cell_tree.nwk.tmp2 "True" /mnt/files/tree.$*.beam.mcc.nwk.tmp "BEAM" /mnt/files/tree.$*_cell_tree.labeling.csv /mnt/files/tree.$*.beam.tanglegram.pdf "P"
	rm tree.$*_cell_tree.nwk.tmp2 tree.$*.beam.mcc.nwk.tmp
	
tree.%.beam_probability_graph.csv: tree.%.beam.trees tree.%.beam.term
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/beam_postprocessing.py \
		/mnt/files/tree.$*.beam.trees location P 250 0.10 /mnt/files/tree.$*.beam

# Get BEAM precision, recall, and F1 scores from edge-wise probability graph
tree.%.beam.precision_recall.csv: tree.%.beam_probability_graph.csv tree.%_cell_tree.migrations
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/get_performance_metrics.py \
		/mnt/files/tree.$*_cell_tree.migrations \
		/mnt/files/tree.$*.beam_probability_graph.csv \
		/mnt/files/tree.$*.beam.precision_recall.csv

# Run LAML to get an input tree for Metient and MACH2
# Added --noultrametric to the LAML v1.0.5 run since it gives ultrametric trees. Seems like a bug they will need to fix.
# Remove the origin node from the LAML tree to be compatible with Metient and MACH2 input prep requirements
tree.%.laml_trees.nwk tree.%.laml-time tree.%.laml_trees_no_origin.nwk: tree.%_indel_character_matrix.csv tree.%.cass.nwk
	singularity exec --bind $(CURDIR):/mnt $(LAML_SIF) /usr/bin/time -o /mnt/tree.$*.laml-time \
	run_laml -c /mnt/tree.$*_indel_character_matrix.csv -t /mnt/tree.$*.cass.nwk -o /mnt/tree.$*.laml --topology_search --keep_polytomies --noultrametric; \
	fixed_newick=$$(cat tree.$*.laml_trees.nwk | sed 's/\[&R\]//g' | sed "s/);//g"); \
	fixed_newick="$${fixed_newick:1});"; \
	last_branch=$$(echo "$$fixed_newick" | grep -oP ':[0-9]+\.[0-9]+(?=\)\;)' | tail -1 | tr -d ':'); \
	second_last_branch=$$(echo "$$fixed_newick" | grep -oP ':[0-9]+\.[0-9]+(?=\))' | tail -2 | head -1 | tr -d ':'); \
	new_branch=$$(echo "$$last_branch + $$second_last_branch" | bc); \
	modified_newick=$$(echo "$$fixed_newick" | sed -r "s/$$second_last_branch\):$$last_branch/0$$new_branch/g"); \
	echo "$$modified_newick" > tree.$*.laml_trees_no_origin.nwk

# Plot LAML vs true tree tanglegram
# Some adjustments needed to handle the origin branch here
tree.%.laml.tanglegram.pdf: tree.%_cell_tree.nwk tree.%_cell_tree.labeling.csv tree.%.laml_trees.nwk
	sed 's/)1:0;/);/' tree.$*_cell_tree.nwk > tree.$*_cell_tree.nwk.tmp3
	cat tree.$*.laml_trees.nwk | sed 's/\[&R\]//g' | sed 's/):[0-9.]*)\?;/);/' > tree.$*.laml_trees.nwk.tmp
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_tanglegram.R \
		/mnt/files/tree.$*_cell_tree.nwk.tmp3 "True" /mnt/files/tree.$*.laml_trees.nwk.tmp "LAML" /mnt/files/tree.$*_cell_tree.labeling.csv /mnt/files/tree.$*.laml.tanglegram.pdf "P"
	rm tree.$*_cell_tree.nwk.tmp3 tree.$*.laml_trees.nwk.tmp

# Prepare metient input files
tree.%.metient.tree.txt tree.%.metient.metadata.txt: tree.%.laml_trees_no_origin.nwk tree.%_cell_tree.labeling.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/format_metient.py \
		/mnt/files/tree.$*.laml_trees_no_origin.nwk \
		/mnt/files/tree.$*_cell_tree.labeling.csv \
		P \
		/mnt/files/tree.$*.metient.tree.txt \
		/mnt/files/tree.$*.metient.metadata.txt

# Run Metient-evaluate
# Resolve polytomies only if the tree has fewer than 150 taxa, which is a requirement by Metient
tree.%.metient_probability_graph.csv tree.%.metient-time tree.%.metient.term: tree.%.metient.tree.txt tree.%.metient.metadata.txt
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(METIENT_SIF) /usr/bin/time -o /mnt/files/tree.$*.metient-time \
	python3 /mnt/scripts/run_metient_evaluate.py \
		/mnt/files/tree.$*.metient.tree.txt \
		/mnt/files/tree.$*.metient.metadata.txt \
		tree.$*.metient \
		P \
		/mnt/files/ \
		/mnt/files/tree.$*.metient_probability_graph.csv \
		1 \
		$$(if [ $$(expr 2 \* $(NTAXA) - 1) -lt 150 ]; then echo "True"; else echo "False"; fi) > tree.$*.metient.term

# Get metient precision, recall, and F1 scores from edge-wise probability graph
tree.%.metient.precision_recall.csv: tree.%.metient_probability_graph.csv tree.%_cell_tree.migrations
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/get_performance_metrics.py \
		/mnt/files/tree.$*_cell_tree.migrations \
		/mnt/files/tree.$*.metient_probability_graph.csv \
		/mnt/files/tree.$*.metient.precision_recall.csv

# Prepare mach2 input files
tree.%.mach2.labeling tree.%.mach2.tree tree.%.mach2.colors: tree.%.laml_trees_no_origin.nwk tree.%_cell_tree.labeling.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/prep_mach2.py \
		/mnt/files/tree.$*.laml_trees_no_origin.nwk \
		/mnt/files/tree.$*_cell_tree.labeling.csv \
		P \
		/mnt/files/tree.$*.mach2.labeling \
		/mnt/files/tree.$*.mach2.tree \
		/mnt/files/tree.$*.mach2.colors

# Run mach2
# NOTE: An active gurobi license must be available locally with its path set in the variable GRB_LICENSE_FILE 
tree.%.mach2.P-G-0.graph tree.%.mach2-time: tree.%.mach2.labeling tree.%.mach2.tree tree.%.mach2.colors
	singularity exec --bind $(GRB_LICENSE_FILE):/mnt/gurobi.lic --env GRB_LICENSE_FILE=/mnt/gurobi.lic --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files \
	$(MACH2_SIF) /usr/bin/time -o /mnt/files/tree.$*.mach2-time mach2 \
	/mnt/files/tree.$*.mach2.tree \
	/mnt/files/tree.$*.mach2.labeling \
	--colormap /mnt/files/tree.$*.mach2.colors -p P --log -o /mnt/files/tree.$*.mach2 -t 1;
	for file in tree.$*.mach2/*; do \
		mv $$file ./tree.$*.mach2.$$(basename $$file); \
	done; \
	rm -rf tree.$*.mach2

# Summarize mach2 results
tree.%.mach2_probability_graph.csv tree.%.mach2.graph_results_combined.txt: tree.%.mach2.P-G-0.graph
	@results=$$(find ./ -maxdepth 1 -name "tree.$*.mach2.*.graph"); \
	echo "result_num,source,target,count" > tree.$*.mach2.graph_results_combined.txt; \
	for file in $$results; do \
		result_num=$$(basename "$$file" | cut -d'-' -f3 | cut -d'.' -f1); \
		while IFS= read -r line; do \
			echo "$$result_num,$$(echo "$$line" | tr ' ' ',' | tr '\t' ',')" >> tree.$*.mach2.graph_results_combined.txt; \
		done < "$$file"; \
	done; \
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) \
		python3 /mnt/scripts/convert_mach2_results_to_consensus_graph.py \
		/mnt/files/tree.$*.mach2.graph_results_combined.txt \
		/mnt/files/tree.$*.mach2_probability_graph.csv

# Get mach2 precision, recall, and F1 scores from edge-wise probability graph
tree.%.mach2.precision_recall.csv: tree.%.mach2_probability_graph.csv tree.%_cell_tree.migrations
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/get_performance_metrics.py \
		/mnt/files/tree.$*_cell_tree.migrations \
		/mnt/files/tree.$*.mach2_probability_graph.csv \
		/mnt/files/tree.$*.mach2.precision_recall.csv


# Extract log likelihoods, using NA when methods are disabled.
tree.%.lnl: tree.%.var.log $(LNL_BEAM_DEPS) $(LNL_AUX_DEPS) $(COMMON_MK)
	@printf "%s\t" "$*" > $@; \
	if [ "$(AUX_METHODS_ENABLED)" = "1" ]; then \
		grep '^Negative' tree.$*.laml_params.txt | awk '{printf "-%f\t", $$2}' >> $@; \
	else \
		printf "NA\t" >> $@; \
	fi
	tail -1 tree.$*.var.log | awk '{printf "%f\t%f\t", $$11, $$21}' >> $@
	@if [ "$(BEAM_ENABLED)" = "1" ]; then \
		grep -v '^#' tree.$*.beam.log | grep -v '^Sample' | \
			awk '{a=$$(NF-4); b=$$(NF-3); s=a+b; if (NR==1 || s>max) {max=s; A=a; B=b}} END {printf "%f\t%f\n", A, B}' >> $@; \
	else \
		printf "NA\tNA\n" >> $@; \
	fi

eval.all.lnl.txt: $(LNLS)
	awk 'function avg(sum, n) { return n ? sprintf("%.2f", sum/n) : "NA" } \
	     BEGIN { print "cp\tlaml_tree\tvine_tree\tvine_mig\tbeam_tree\tbeam_mig" } \
	     { print; if ($$2 != "NA") { x1 += $$2; n1++ } \
	       x2 += $$3; x3 += $$4; n++; \
	       if ($$5 != "NA") { x4 += $$5; n4++ } \
	       if ($$6 != "NA") { x5 += $$6; n5++ } } \
	     END { print "-----------------------------------------"; \
	       printf "all\t%s\t%.2f\t%.2f\t%s\t%s\n", \
	              avg(x1, n1), x2/n, x3/n, avg(x4, n4), avg(x5, n5) }' $(LNLS) > $@

eval.all.lnl.pdf: eval.all.lnl.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_lnl.R /mnt/files/eval.all.lnl.txt /mnt/files/eval.all.lnl.pdf

# Get times per simulation
# For BEAM, we take the time at which ESS > $(BEAM_ESS_TARGET) ("convergence" point), not the total arbitrary chain length time
# For other methods we use the total runtime as reported by the time command
tree.%.time: tree.%.var-time $(TIME_BEAM_DEPS) $(TIME_AUX_DEPS) $(COMMON_MK)
	@{ \
	printf "%s\t" "$*" > $@; \
	vine=$$(head -n 1 tree.$*.var-time | awk '{printf "%.3f", $$1}'); \
	laml="NA"; metient="NA"; mach2="NA"; beam="NA"; beam_improvement="NA"; \
	if [ "$(AUX_METHODS_ENABLED)" = "1" ]; then \
		laml=$$(head -n 1 tree.$*.laml-time | awk '{printf "%.3f", $$1}'); \
		metient=$$(head -n 1 tree.$*.metient-time | awk '{printf "%.3f", $$1}'); \
		mach2=$$(head -n 1 tree.$*.mach2-time | awk '{printf "%.3f", $$1}'); \
	fi; \
	if [ "$(BEAM_ENABLED)" = "1" ]; then \
		beam=$$(awk '/\/Msamples/{ \
		n_last=$$1; x_last=$$NF; sub(/\/Msamples.*/,"",x_last); \
		if (($$3+0)>$(BEAM_ESS_TARGET)) { n=n_last; x=x_last; found=1; exit } \
	} \
	END{ \
		if(!found){n=n_last; x=x_last} \
		h=m=s=0; \
		if(match(x,/([0-9]+)h/,a)) h=a[1]; \
		if(match(x,/([0-9]+)m/,a)) m=a[1]; \
		if(match(x,/([0-9]+)s/,a)) s=a[1]; \
		printf "%.3f", (n/1e6) * (3600*h + 60*m + s) \
		}' tree.$*.beam.term); \
		beam_improvement=$$(echo "scale=6; $$beam / $$vine" | bc -l | awk '{printf "%.3f", $$1}'); \
	fi; \
	echo -e "$$vine\t$$laml\t$$metient\t$$mach2\t$$beam\t$$beam_improvement" >> $@; \
	}

# Get overall summary files with times from all simulations
eval.all.time.txt: $(TIMES)
	awk 'function avg(sum, n) { return n ? sprintf("%.2f", sum/n) : "NA" } \
	     BEGIN{print "sim\tvine\tlaml\tmetient\tmach2\tbeam\tbeam_improvement"} \
	     FNR==1 { \
	       print; \
	       sumV+=$$2; if ($$3 != "NA") { sumL+=$$3; nL++ } \
	       if ($$4 != "NA") { sumM+=$$4; nM++ } \
	       if ($$5 != "NA") { sum2+=$$5; n2++ } \
	       if ($$6 != "NA") { sumB+=$$6; nB++ } \
	       if ($$7 != "NA") { sumIB+=$$7; nIB++ } n++; \
	     } \
	     END{ \
	       if(n>0){ \
	         print "-----------------------------------------"; \
	         printf "avg\t%.2f\t%s\t%s\t%s\t%s\t%s\n", \
	                sumV/n, avg(sumL,nL), avg(sumM,nM), avg(sum2,n2), \
	                avg(sumB,nB), avg(sumIB,nIB); \
	       } \
	     }' $(TIMES) > $@



eval.all.time.pdf: eval.all.time.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_runtime_bars.R /mnt/files/eval.all.time.txt /mnt/files/eval.all.time.pdf

# Get precision, recall, and f1 scores for 0.50 probability threshold graphs for quick comparisons
THRESH := 0.5
tree.%.50_prf.txt: tree.%.var.precision_recall.csv $(PRF_BEAM_DEPS) $(PRF_AUX_DEPS) $(COMMON_MK)
	@SIM="$*"; \
	awk -F, -v sim="$$SIM" 'NR==1{ \
		print "sim\tmethod\tprecision\trecall\tf1"; \
	} NR>1{ \
		if($$1=="$(THRESH)"){ \
			print sim"\tVINE\t"$$2"\t"$$3"\t"$$4; \
		} \
	}' tree.$*.var.precision_recall.csv > $@; \
	if [ "$(BEAM_ENABLED)" = "1" ]; then \
		awk -F, -v sim="$$SIM" 'NR==1{next} NR>1{ \
			if($$1=="$(THRESH)"){ \
				print sim"\tBEAM\t"$$2"\t"$$3"\t"$$4; \
			} \
		}' tree.$*.beam.precision_recall.csv >> $@; \
	fi; \
	if [ "$(AUX_METHODS_ENABLED)" = "1" ]; then \
		awk -F, -v sim="$$SIM" 'NR==1{next} NR>1{ \
		if($$1=="$(THRESH)"){ \
			print sim"\tMETIENT\t"$$2"\t"$$3"\t"$$4; \
		} \
		}' tree.$*.metient.precision_recall.csv >> $@; \
		awk -F, -v sim="$$SIM" 'NR==1{next} NR>1{ \
		if($$1=="$(THRESH)"){ \
			print sim"\tMACH2\t"$$2"\t"$$3"\t"$$4; \
		} \
		}' tree.$*.mach2.precision_recall.csv >> $@; \
	fi

eval.all.50_prf.txt: $(50PRFS)
	awk 'FNR==1 && NR==1 {print "sim\tmethod\tprecision\trecall\tf1"; next} FNR==1 {next} {print}' tree.*.50_prf.txt > $@

eval.all.50_prf.pdf: eval.all.50_prf.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_f1_box.R /mnt/files/eval.all.50_prf.txt /mnt/files/eval.all.50_prf.pdf

eval.all.thresh_prf.txt: $(VINEPRFS) $(THRESH_BEAM_PRFS) $(THRESH_AUX_PRFS) $(COMMON_MK)
	echo -e "sim\tmethod\tthreshold\tprecision\trecall\tf1" > $@; \
	for file in $(filter %.precision_recall.csv,$^); do \
		method=$$(basename "$$file" | cut -d'.' -f3); \
		if [ "$$method" == "var" ]; then method="VINE"; \
		elif [ "$$method" == "beam" ]; then method="BEAM"; \
		elif [ "$$method" == "metient" ]; then method="METIENT"; \
		elif [ "$$method" == "mach2" ]; then method="MACH2"; \
		else method="UNKNOWN"; fi; \
		sim=$$(basename "$$file" | cut -d'.' -f2); \
		awk -F, -v sim="$$sim" -v method="$$method" 'NR>1{ print sim"\t"method"\t"$$1"\t"$$2"\t"$$3"\t"$$4; }' $$file >> $@; \
	done

eval.all.thresh_prf.pdf: eval.all.thresh_prf.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_precision_recall.R /mnt/files/eval.all.thresh_prf.txt /mnt/files/eval.all.thresh_prf.pdf

clean:
	rm -f tree.*

clean_beam:
	rm -f tree.*.beam* tree.*.time tree.*.50_prf.txt eval.*

clean_not_converged_beam:
	for cp in $(TREE_IDS); do \
		if [ -f tree.$$cp.beam.term ]; then \
			ess=$$(grep "Msamples" tree.$$cp.beam.term | tail -n 1 | awk '{print $$6}'); \
			if [ -z "$$ess" ] || awk -v ess="$$ess" -v target=$(BEAM_ESS_TARGET) 'BEGIN{exit !(ess < target)}'; then \
				echo "Deleting tree.$$cp (ESS=$$ess)"; \
				rm -f tree.$$cp.beam.term \
				      tree.$$cp.beam_probability_graph.csv \
				      tree.$$cp.beam*.pdf; \
			fi; \
		fi; \
	done

clean_vine:
	rm -f tree.*.var* eval.all.*

archive_all:
	archive_dir=archive_all.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv tree.* $$archive_dir; \
	mv eval.* $$archive_dir;

archive_vine:
	archive_dir=archive_vine.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv tree.*.var_probability_graph.csv $$archive_dir; \
	mv tree.*.var.precision_recall.csv $$archive_dir; \
	mv tree.*.var-time $$archive_dir; \
	mv tree.*.var.dot $$archive_dir; \
	mv tree.*.var.nwk $$archive_dir; \
	mv tree.*.var.nex $$archive_dir; \
	mv tree.*.var.log $$archive_dir; \
	mv tree.*.time $$archive_dir; \
	mv tree.*.50_prf.txt $$archive_dir; \
	mv eval.* $$archive_dir;

archive_beam:
	archive_dir=archive_beam.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv tree.*.beam* $$archive_dir; \
	mv tree.*.time $$archive_dir; \
	mv tree.*.50_prf.txt $$archive_dir; \
	mv eval.* $$archive_dir;

archive_metient:
	archive_dir=archive_metient.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv tree.*.metient* $$archive_dir; \
	mv tree.*.time $$archive_dir; \
	mv tree.*.50_prf.txt $$archive_dir; \
	mv eval.* $$archive_dir;
