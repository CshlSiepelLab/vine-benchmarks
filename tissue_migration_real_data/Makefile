export SHELL=/usr/bin/bash

# Give beast2 access to beagle
export LD_LIBRARY_PATH := /local/storage/no-backup/vine-benchmarks/lib:$(LD_LIBRARY_PATH)

THREADS := 1

# Prevent Make from removing files it thinks it should clean up
.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS: quinn.%.beam.term quinn.%.beam.log quinn.%.beam.trees

BEAM_ESS_TARGET := 400

# Path to raw data
QUINNROOT := ../crispr_real_data/quinn_cps/5k

# Path to BEAM template
BEAM_TEMPLATE := /local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/beam_template.xml

# Absolute paths
MAIN_DIR := /local/storage/no-backup/vine-benchmarks

# Relative paths
PYTHON_SRC := $(MAIN_DIR)/python/src
BIN := $(MAIN_DIR)/bin
R_SRC := $(MAIN_DIR)/r/src
VINE := $(BIN)/vine/bin/vine
BEAST := $(BIN)/beast/bin/beast
TREEANNOTATOR := $(BIN)/beast/bin/treeannotator
CONTAINERS := $(MAIN_DIR)/containers
GRAPHPOSTERIOR_SIF := $(CONTAINERS)/graphposterior/graphposterior.sif
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
MACH2_SIF := $(CONTAINERS)/mach2/mach2.sif
LAML_SIF := $(CONTAINERS)/laml/laml.sif
RPLOTTING_SIF := $(CONTAINERS)/rplotting/rplotting.sif
METIENT_SIF := $(CONTAINERS)/metient/metient.sif


# Subset of CPs for migration analysis in the VINE paper
# Excluded in the original paper barcode quality control: 5 16 18 25 33 38 39 41 50 53 65 69 75 81 87 88 93
# Exclude > 250 taxa: 4 6 7 8 9 10 11 12 13 14 15 17 19 20 21
# Exclude because only one tissue observed: 22 29 46 78 85 94 96
QUINNCPS := 23 24 26 27 28 30 31 32 34 35 36 37 40 42 43 44 45 47 48 49 51 52 54 55 56 57 58 59 60 61 62 63 64 66 67 68 70 71 72 73 74 76 77 79 80 82 83 84 86 89 90 91 92 95 97 98 99 100

# Need to exclude only from MACH2 runs (due to a bug in their code?): 31

MATRICES := $(QUINNCPS:%=quinn.%.tsv)
MATRICES_CSV := $(QUINNCPS:%=quinn.%.csv)
FASTA := $(QUINNCPS:%=quinn.%.fasta)
TISSUES := $(QUINNCPS:%=quinn.%.tissues)
CASSTREES := $(QUINNCPS:%=quinn.%.cass.nwk)
VARTREES := $(QUINNCPS:%=quinn.%.var.nwk)
VARLOG := $(QUINNCPS:%=quinn.%.var.log)
VARPROBGRAPH := $(QUINNCPS:%=quinn.%.var_probability_graph.csv)
VARTERM := $(QUINNCPS:%=quinn.%.var.term)
BEAMPROBGRAPH := $(QUINNCPS:%=quinn.%.beam_probability_graph.csv)
BEAMTERM := $(QUINNCPS:%=quinn.%.beam.term)
LAMLTREES := $(QUINNCPS:%=quinn.%.laml_trees_no_origin.nwk)
MACH2COLORS := $(QUINNCPS:%=quinn.%.mach2.colors)
MACH2PROBGRAPHS := $(QUINNCPS:%=quinn.%.mach2_probability_graph.csv)
CASSMACH2PROBGRAPHS := $(QUINNCPS:%=quinn.%.cass.mach2_probability_graph.csv)
METIENTMETA := $(QUINNCPS:%=quinn.%.metient.metadata.txt)
METIENTPROBGRAPHS := $(QUINNCPS:%=quinn.%.metient_probability_graph.csv)
CASSMETIENTPROBGRAPHS := $(QUINNCPS:%=quinn.%.cass.metient_probability_graph.csv)
LAMLTANGLEGRAM := $(QUINNCPS:%=quinn.%.var.laml.tanglegram.pdf)
CASSTANGLEGRAM := $(QUINNCPS:%=quinn.%.var.cass.tanglegram.pdf)
BEAMTANGLEGRAM := $(QUINNCPS:%=quinn.%.var.beam.tanglegram.pdf)
TIMES := $(QUINNCPS:%=quinn.%.time)
LNLS := $(QUINNCPS:%=quinn.%.lnl)
STATS := $(QUINNCPS:%=quinn.%.stats)
GRAPHSTATS := $(QUINNCPS:%=quinn.%.graph.stats)

VAROPT?=-M 400 -c 100 -v 0 

all: $(VARPROBGRAPH) $(BEAMPROBGRAPH) $(METIENTPROBGRAPHS) # $(MACH2PROBGRAPHS)
tanglegrams: $(LAMLTANGLEGRAM) $(CASSTANGLEGRAM)
cass: $(CASSTREES)
laml: $(LAMLTREES)
prepmach2: $(MACH2COLORS)
mach2: $(MACH2PROBGRAPHS)
cassmach2: $(CASSMACH2PROBGRAPHS)
prepmetient: $(METIENTMETA)
metient: $(METIENTPROBGRAPHS)
cassmetient: $(CASSMETIENTPROBGRAPHS)
vine: $(VARPROBGRAPH)
beam: $(BEAMPROBGRAPH)
vine_infer: $(VARLOG)
beam_infer: $(BEAMTERM)
source_data: $(MATRICES) $(CASSTREES)

# Copy source data and reformat
quinn.%.original.csv:
	sed 's/,r[0-9]*$$//' $(QUINNROOT)/`echo $@ | awk -F\. '{print $$2}'`/matrix.csv > $@

# Extract the tissue labels and simplify them to remove multiple tissue samples (e.g. RE and RW -> R; M1, M2, M3 -> M)
quinn.%.original.tissues: quinn.%.original.csv
	tail -n +2 $< | cut -d, -f1 | awk -F. '{print $$0 "," $$1}' > $@
# 	awk -F, '{ gsub(/[0-9]/,"",$$2); if ($$2 ~ /^R/) $$2="R"; print $$1 "," $$2 }' $@ > $@.tmp && mv $@.tmp $@

quinn.%.expanded.csv quinn.%.tissues.expanded.csv quinn.%.csv quinn.%.tissues.tsv quinn.%.collapsing_map.tsv: quinn.%.original.csv quinn.%.original.tissues
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/collapse_matrix_duplicates.py \
	 -m /mnt/files/quinn.$*.original.csv \
	 -om /mnt/files/quinn.$*.csv \
	 -on /mnt/files/quinn.$*.collapsing_map.tsv	\
	 -t /mnt/files/quinn.$*.original.tissues \
	 -ot /mnt/files/quinn.$*.tissues.tsv

quinn.%.tsv: quinn.%.expanded.csv
	$(BIN)/csv-to-tsv.sh $< > $@

quinn.%.cass.nwk: quinn.%.expanded.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(CASSIOPEIA_SIF) \
	python /mnt/scripts/cassiopeiaGreedy.py /mnt/files/quinn.$*.expanded.csv /mnt/files/quinn.$*.cass.nwk

# Plot VINE MCC circular tree
quinn.%.cass.circular.pdf: quinn.%.tissues.expanded.csv quinn.%.cass.nwk
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_circular_tree.R \
		/mnt/files/quinn.$*.cass.nwk "CASSIOPEIA-GREEDY" /mnt/files/quinn.$*.tissues.expanded.csv /mnt/files/quinn.$*.cass.circular.pdf "LL"

# Added --noultrametric to the LAML v1.0.5 run since it gives ultrametric trees. Seems like a big they will need to fix.
# Need to replace _ characters that LAML cannot accept
quinn.%.laml_trees.nwk quinn.%.laml-time: quinn.%.expanded.csv quinn.%.cass.nwk
	singularity exec --bind $(CURDIR):/mnt $(LAML_SIF) /usr/bin/time -o /mnt/quinn.$*.laml-time \
	run_laml -c /mnt/quinn.$*.expanded.csv -t /mnt/quinn.$*.cass.nwk -o /mnt/quinn.$*.laml --topology_search --keep_polytomies --noultrametric

# Reformat laml tree to remove origin node so it can be used to prepare MACH2 input files
quinn.%.laml_trees_no_origin.nwk: quinn.%.laml_trees.nwk
	fixed_newick=$$(cat quinn.$*.laml_trees.nwk | sed 's/\[&R\]//g' | sed "s/);//g"); \
	fixed_newick="$${fixed_newick:1});"; \
	last_branch=$$(echo "$$fixed_newick" | grep -oP ':[0-9]+\.[0-9]+(?=\)\;)' | tail -1 | tr -d ':'); \
	second_last_branch=$$(echo "$$fixed_newick" | grep -oP ':[0-9]+\.[0-9]+(?=\))' | tail -2 | head -1 | tr -d ':'); \
	new_branch=$$(echo "$$last_branch + $$second_last_branch" | bc); \
	modified_newick=$$(echo "$$fixed_newick" | sed -r "s/$$second_last_branch\):$$last_branch/0$$new_branch/g"); \
	echo "$$modified_newick" > quinn.$*.laml_trees_no_origin.nwk

# Prepare mach2 input files from LAML
quinn.%.mach2.labeling quinn.%.mach2.tree quinn.%.mach2.colors: quinn.%.laml_trees_no_origin.nwk quinn.%.tissues.expanded.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/prep_mach2.py \
		/mnt/files/quinn.$*.laml_trees_no_origin.nwk \
		/mnt/files/quinn.$*.tissues.expanded.csv \
		LL \
		/mnt/files/quinn.$*.mach2.labeling \
		/mnt/files/quinn.$*.mach2.tree \
		/mnt/files/quinn.$*.mach2.colors

# Run mach2 from LAML
# NOTE: An active gurobi license must be available locally with its path set in the variable GRB_LICENSE_FILE 
quinn.%.mach2.LL-G-0.graph quinn.%.mach2-time: quinn.%.mach2.labeling quinn.%.mach2.tree quinn.%.mach2.colors
	singularity exec --bind $(GRB_LICENSE_FILE):/mnt/gurobi.lic --env GRB_LICENSE_FILE=/mnt/gurobi.lic --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files \
	$(MACH2_SIF) /usr/bin/time -o /mnt/files/quinn.$*.mach2-time mach2 \
	/mnt/files/quinn.$*.mach2.tree \
	/mnt/files/quinn.$*.mach2.labeling \
	--colormap /mnt/files/quinn.$*.mach2.colors -p LL --log -o /mnt/files/quinn.$*.mach2 -t $(THREADS);
	for file in quinn.$*.mach2/*; do \
		mv $$file ./quinn.$*.mach2.$$(basename $$file); \
	done; \
	rm -rf quinn.$*.mach2

# Summarize mach2 results from LAML
quinn.%.mach2_probability_graph.csv quinn.%.mach2.graph_results_combined.txt: quinn.%.mach2.LL-G-0.graph
	@results=$$(find ./ -maxdepth 1 -name "quinn.$*.mach2.*.graph"); \
	echo "result_num,source,target,count" > quinn.$*.mach2.graph_results_combined.txt; \
	for file in $$results; do \
		result_num=$$(basename "$$file" | cut -d'-' -f3 | cut -d'.' -f1); \
		while IFS= read -r line; do \
			echo "$$result_num,$$(echo "$$line" | tr ' ' ',' | tr '\t' ',')" >> quinn.$*.mach2.graph_results_combined.txt; \
		done < "$$file"; \
	done; \
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) \
		python3 /mnt/scripts/convert_mach2_results_to_consensus_graph.py \
		/mnt/files/quinn.$*.mach2.graph_results_combined.txt \
		/mnt/files/quinn.$*.mach2_probability_graph.csv

# Prepare mach2 input files from Cassiopeia-greedy tree
quinn.%.cass.mach2.labeling quinn.%.cass.mach2.tree quinn.%.cass.mach2.colors: quinn.%.cass.nwk quinn.%.tissues.expanded.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/prep_mach2.py \
		/mnt/files/quinn.$*.cass.nwk \
		/mnt/files/quinn.$*.tissues.expanded.csv \
		LL \
		/mnt/files/quinn.$*.cass.mach2.labeling \
		/mnt/files/quinn.$*.cass.mach2.tree \
		/mnt/files/quinn.$*.cass.mach2.colors

# Run mach2 from Cassiopeia-greedy tree
# NOTE: An active gurobi license must be available locally with its path set in the variable GRB_LICENSE_FILE 
# Since MACH2 requires the specified primary tissue to be in the observed labeling, we use the most common tissue when the LL primary is not observed in the clone
quinn.%.cass.mach2.LL-G-0.graph quinn.%.cass.mach2-time: quinn.%.cass.mach2.labeling quinn.%.cass.mach2.tree quinn.%.cass.mach2.colors
	if awk '{print $$2}' quinn.$*.cass.mach2.labeling | grep -q "^LL$$"; then \
		primary_tissue="LL"; \
	else \
		primary_tissue=$$(awk '{print $$2}' quinn.$*.cass.mach2.labeling | sort | uniq -c | sort -rn | awk 'NR==1 {print $$2}'); \
	fi; \
	singularity exec --bind $(GRB_LICENSE_FILE):/mnt/gurobi.lic --env GRB_LICENSE_FILE=/mnt/gurobi.lic --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files \
	$(MACH2_SIF) /usr/bin/time -o /mnt/files/quinn.$*.cass.mach2-time mach2 \
	/mnt/files/quinn.$*.cass.mach2.tree \
	/mnt/files/quinn.$*.cass.mach2.labeling \
	--colormap /mnt/files/quinn.$*.cass.mach2.colors -p $$primary_tissue --log -o /mnt/files/quinn.$*.cass.mach2 -t $(THREADS);
	for file in quinn.$*.cass.mach2/*; do \
		mv $$file ./quinn.$*.cass.mach2.$$(basename $$file); \
	done; \
	rm -rf quinn.$*.cass.mach2
	touch quinn.$*.cass.mach2.LL-G-0.graph # Create empty graph file for cases where mach2 outputs more than 10 solutions so the 0 is either 00 or 00 in the output name

# Summarize mach2 results from Cassiopea-greedy tree
quinn.%.cass.mach2_probability_graph.csv quinn.%.cass.mach2.graph_results_combined.txt: quinn.%.cass.mach2.LL-G-0.graph
	@results=$$(find ./ -maxdepth 1 -name "quinn.$*.cass.mach2.*.graph"); \
	echo "result_num,source,target,count" > quinn.$*.cass.mach2.graph_results_combined.txt; \
	for file in $$results; do \
		result_num=$$(basename "$$file" | cut -d'-' -f3 | cut -d'.' -f1); \
		while IFS= read -r line; do \
			echo "$$result_num,$$(echo "$$line" | tr ' ' ',' | tr '\t' ',')" >> quinn.$*.cass.mach2.graph_results_combined.txt; \
		done < "$$file"; \
	done; \
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) \
		python3 /mnt/scripts/convert_mach2_results_to_consensus_graph.py \
		/mnt/files/quinn.$*.cass.mach2.graph_results_combined.txt \
		/mnt/files/quinn.$*.cass.mach2_probability_graph.csv

quinn.%.cass.mach2_90_graph_matrix.pdf: quinn.%.cass.mach2_probability_graph.csv
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_graph_matrix_realData.R \
		/mnt/files/quinn.$*.cass.mach2_probability_graph.csv 0.90 /mnt/files/quinn.$*.cass.mach2_90_graph_matrix.pdf "LL"

quinn.4.cass.mach2.LL-G-0.graph_matrix.pdf: quinn.4.cass.mach2.LL-G-0.graph.csv
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_graph_matrix_realData.R \
		/mnt/files/quinn.4.cass.mach2.LL-G-0.graph.csv 0.90 /mnt/files/quinn.4.cass.mach2.LL-G-0.graph_matrix.pdf "LL"


# Prepare metient input files from LAML
quinn.%.metient.tree.txt quinn.%.metient.metadata.txt: quinn.%.laml_trees_no_origin.nwk quinn.%.tissues.expanded.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/format_metient.py \
		/mnt/files/quinn.$*.laml_trees_no_origin.nwk \
		/mnt/files/quinn.$*.tissues.expanded.csv \
		LL \
		/mnt/files/quinn.$*.metient.tree.txt \
		/mnt/files/quinn.$*.metient.metadata.txt

# Run Metient-evaluate from LAML
quinn.%.metient_probability_graph.csv quinn.%.metient-time quinn.%.metient.term: quinn.%.metient.tree.txt quinn.%.metient.metadata.txt quinn.%.expanded.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(METIENT_SIF) /usr/bin/time -o /mnt/files/quinn.$*.metient-time \
	python3 /mnt/scripts/run_metient_evaluate.py \
		/mnt/files/quinn.$*.metient.tree.txt \
		/mnt/files/quinn.$*.metient.metadata.txt \
		quinn.$*.metient \
		LL \
		/mnt/files/ \
		/mnt/files/quinn.$*.metient_probability_graph.csv \
		$(THREADS) \
		$$(if [ $$(expr 2 \* $$(tail -n +2 quinn.$*.expanded.csv | wc -l) - 1) -lt 150 ]; then echo "True"; else echo "False"; fi) > quinn.$*.metient.term

# Prepare metient input files from Cassiopeia-greedy
quinn.%.cass.metient.tree.txt quinn.%.cass.metient.metadata.txt: quinn.%.cass.nwk quinn.%.tissues.expanded.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/format_metient.py \
		/mnt/files/quinn.$*.cass.nwk \
		/mnt/files/quinn.$*.tissues.expanded.csv \
		LL \
		/mnt/files/quinn.$*.cass.metient.tree.txt \
		/mnt/files/quinn.$*.cass.metient.metadata.txt

# Run Metient-evaluate from Cassiopeia-greedy tree
quinn.%.cass.metient_probability_graph.csv quinn.%.cass.metient-time quinn.%.cass.metient.term: quinn.%.cass.metient.tree.txt quinn.%.cass.metient.metadata.txt quinn.%.expanded.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(METIENT_SIF) /usr/bin/time -o /mnt/files/quinn.$*.cass.metient-time \
	python3 /mnt/scripts/run_metient_evaluate.py \
		/mnt/files/quinn.$*.cass.metient.tree.txt \
		/mnt/files/quinn.$*.cass.metient.metadata.txt \
		quinn.$*.cass.metient \
		LL \
		/mnt/files/ \
		/mnt/files/quinn.$*.cass.metient_probability_graph.csv \
		$(THREADS) \
		$$(if [ $$(expr 2 \* $$(tail -n +2 quinn.$*.expanded.csv | wc -l) - 1) -lt 150 ]; then echo "True"; else echo "False"; fi) > quinn.$*.cass.metient.term

quinn.%.cass.metient_90_graph_matrix.pdf: quinn.%.cass.metient_probability_graph.csv
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_graph_matrix_realData.R \
		/mnt/files/quinn.$*.cass.metient_probability_graph.csv 0.90 /mnt/files/quinn.$*.cass.metient_90_graph_matrix.pdf "LL"

# quinn.4.cass.metient_LL.topgraph_matrix.pdf: quinn.4.cass.metient_LL.topgraph.csv
# 	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_graph_matrix_realData.R \
# 		/mnt/files/quinn.4.cass.metient_LL.topgraph.csv 0.90 /mnt/files/quinn.4.cass.metient_LL.topgraph_matrix.pdf "LL"

# Run vine in tissue migration mode
quinn.%.var.log quinn.%.var-time quinn.%.var.nwk quinn.%.var.dot quinn.%.var.nex: quinn.%.tsv quinn.%.tissues.expanded.csv
	/usr/bin/time -o quinn.$*.var-time $(VINE) $(VAROPT) -j $(THREADS) -i CRISPR quinn.$*.tsv --logfile quinn.$*.var.log --mean quinn.$*.mean.nwk \
	--migration quinn.$*.tissues.expanded.csv --migration-rate-prior 0 --primary LL --sample-graphs quinn.$*.var.dot --labeled-trees quinn.$*.var.nex > quinn.$*.var.nwk

# Get MCC tree from VINE output
quinn.%.var.mcc.nex quinn.%.var.mcc.nwk: quinn.%.var.nex
	$(TREEANNOTATOR) -lowMem true -height keep -burnin 0 $^ | \
    sed 's/\[[^]]*\]//g' > quinn.$*.var.mcc.nex
	$(BIN)/nex2nwk quinn.$*.var.mcc.nex quinn.$*.var.mcc.nwk

# Plot VINE MCC circular tree
quinn.%.var.mcc.circular.pdf: quinn.%.tissues.expanded.csv quinn.%.var.mcc.nwk
	sed 's/):0.00000;/);/' quinn.$*.var.mcc.nwk > quinn.$*.var.mcc.nwk.tmp
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_circular_tree.R \
		/mnt/files/quinn.$*.var.mcc.nwk.tmp "VINE" /mnt/files/quinn.$*.tissues.expanded.csv /mnt/files/quinn.$*.var.mcc.circular.pdf "LL"
	rm quinn.$*.var.mcc.nwk.tmp

# Plot VINE vs CASS tree tanglegram
quinn.%.var.cass.tanglegram.pdf: quinn.%.cass.nwk quinn.%.tissues.expanded.csv quinn.%.var.mcc.nwk
	sed 's/):0.00000;/);/' quinn.$*.var.mcc.nwk > quinn.$*.var.mcc.nwk.tmp
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_tanglegram.R \
		/mnt/files/quinn.$*.var.mcc.nwk.tmp "VINE" /mnt/files/quinn.$*.cass.nwk "CASSIOPEIA-GREEDY" /mnt/files/quinn.$*.tissues.expanded.csv /mnt/files/quinn.$*.var.cass.tanglegram.pdf "LL"
	rm quinn.$*.var.mcc.nwk.tmp


# Plot VINE vs LAML tree tanglegram
quinn.%.var.laml.tanglegram.pdf: quinn.%.laml_trees.nwk quinn.%.tissues.expanded.csv quinn.%.var.mcc.nwk
	cat quinn.$*.laml_trees.nwk | sed 's/\[&R\]//g' | sed 's/):[0-9.]*)\?;/);/' > quinn.$*.laml_trees.nwk.tmp
	sed 's/):0.00000;/);/' quinn.$*.var.mcc.nwk > quinn.$*.var.mcc.nwk.tmp
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_tanglegram.R \
		/mnt/files/quinn.$*.var.mcc.nwk.tmp "VINE" /mnt/files/quinn.$*.laml_trees.nwk.tmp "LAML" /mnt/files/quinn.$*.tissues.expanded.csv /mnt/files/quinn.$*.var.laml.tanglegram.pdf "LL"
	rm quinn.$*.laml_trees.nwk.tmp quinn.$*.var.mcc.nwk.tmp

quinn.%.var_probability_graph.csv: quinn.%.var.nex
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/beam_postprocessing.py \
		/mnt/files/quinn.$*.var.nex state LL 1 0.0 /mnt/files/quinn.$*.var

quinn.%.var_90_graph_matrix.pdf: quinn.%.var_probability_graph.csv
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_graph_matrix_realData.R \
		/mnt/files/quinn.$*.var_probability_graph.csv 0.90 /mnt/files/quinn.$*.var_90_graph_matrix.pdf "LL"


quinn.%.beam.xml: quinn.%.tissues.expanded.csv quinn.%.cass.nwk $(BEAM_TEMPLATE)
	cp $(BEAM_TEMPLATE) $@
	allTissues=$$(cut -d, -f2 quinn.$*.tissues.expanded.csv | sort | uniq | tr '\n' ' ' | sed 's/ *$$//'); \
	numAllTissues=$$(echo $${allTissues} | wc -w); \
	if ! echo $${allTissues} | grep -q 'LL'; then \
		numAllTissues=$$((numAllTissues + 1)); \
	fi; \
	numTissueRates=$$(( ((numAllTissues * numAllTissues) - numAllTissues) / 2 )); \
	sed -i "s|numTissueRates|$${numTissueRates}|g" $@; \
	numOtherTissues=$$(cut -d, -f2 quinn.$*.tissues.expanded.csv | sort | uniq | sed 's/LL//' | grep -v '^$$' | wc -l); \
	rootTissueFreqs="1"; \
	for i in $$(seq 1 $$numOtherTissues); do \
		rootTissueFreqs="$${rootTissueFreqs} 0"; \
	done; \
	sed -i "s|rootTissueFreqs|$${rootTissueFreqs}|g" $@

# Run BEAM - requires some tuning here of the mcmc length and then sampling frequency
# The seed is fixed to permit resuming runs that don't converge
quinn.%.beam.trees quinn.%.beam.log quinn.%.beam.term: quinn.%.tsv quinn.%.tissues.expanded.csv quinn.%.beam.xml quinn.%.cass.nwk
	$(BEAST) \
		-working \
		$$( [ -f quinn.$*.beam.xml.state ] && echo -resume || echo -overwrite ) \
		-threads 5 \
		-seed 1 \
		-D outfileLog=quinn.$*.beam.log \
		-D outfileTrees=quinn.$*.beam.trees \
		-D barcodeMatrix=quinn.$*.tsv \
		-D tissueCsv=quinn.$*.tissues.expanded.csv \
		-D startingNewickFile=quinn.$*.cass.nwk \
		-D mcmcLength=25000000 \
		-D sampleFreq=2500 \
		quinn.$*.beam.xml > quinn.$*.beam.term

quinn.%.beam_probability_graph.csv: quinn.%.beam.trees quinn.%.beam.term
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/beam_postprocessing.py \
		/mnt/files/quinn.$*.beam.trees location LL 54 0.10 /mnt/files/quinn.$*.beam

# Get MCC tree from BEAM output
quinn.%.beam.mcc.nex quinn.%.beam.mcc.nwk: quinn.%.beam.trees
	$(TREEANNOTATOR) -lowMem true -height keep -burnin 10 $^ | \
    sed 's/\[[^]]*\]//g' > quinn.$*.beam.mcc.nex
	$(BIN)/nex2nwk quinn.$*.beam.mcc.nex quinn.$*.beam.mcc.nwk

# Plot VINE vs CASS tree tanglegram
quinn.%.var.beam.tanglegram.pdf: quinn.%.beam.mcc.nwk quinn.%.tissues.expanded.csv quinn.%.var.mcc.nwk
	sed 's/):0.00000;/);/' quinn.$*.var.mcc.nwk > quinn.$*.var.mcc.nwk.tmp
	sed 's/):0.00000;/);/' quinn.$*.beam.mcc.nwk > quinn.$*.beam.mcc.nwk.tmp
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_tanglegram.R \
		/mnt/files/quinn.$*.var.mcc.nwk.tmp "VINE" /mnt/files/quinn.$*.beam.mcc.nwk.tmp "BEAM" /mnt/files/quinn.$*.tissues.expanded.csv /mnt/files/quinn.$*.var.beam.tanglegram.pdf "LL"
	rm quinn.$*.var.mcc.nwk.tmp quinn.$*.beam.mcc.nwk.tmp

# Get times per simulation
# For BEAM, we take the time at which ESS > $(BEAM_ESS_TARGET) ("convergence" point), not the total arbitrary chain length time
# For other methods we use the total runtime as reported by the time command
quinn.%.time: quinn.%.var-time quinn.%.laml-time quinn.%.metient-time quinn.%.beam.term	# TODO: Need to change metient to laml-metient after I run them tonight (3/20/26)
	@{ \
	printf "%s\t" "$*" > $@; \
	vine=$$(head -n 1 quinn.$*.var-time | awk '{printf "%.3f", $$1}'); \
	laml=$$(head -n 1 quinn.$*.laml-time | awk '{printf "%.3f", $$1}'); \
	metient=$$(head -n 1 quinn.$*.metient-time | awk '{printf "%.3f", $$1}'); \
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
	}' quinn.$*.beam.term); \
	echo -e "$$vine\t$$laml\t$$metient\t$$beam" >> $@; \
	}

# Get overall summary files with times from all CPs
eval.all.time.txt: $(TIMES)
	awk 'BEGIN{print "sim\tvine\tlaml\tmetient\tbeam"} \
	     FNR==1 { \
	       print; \
	       sumV+=$$2; sumL+=$$3; sumM+=$$4; sumB+=$$5; n++; \
	     } \
	     END{ \
	       if(n>0){ \
	         print "-----------------------------------------"; \
	         printf "avg\t%.2f\t%.2f\t%.2f\t%.2f\n", \
	                sumV/n, sumL/n, sumM/n, sumB/n; \
	       } \
	     }' $(TIMES) > $@

eval.all.time.pdf: eval.all.time.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_runtime_bars_realData.R /mnt/files/eval.all.time.txt /mnt/files/eval.all.time.pdf

# Extract lnls for BEAM and VINE
quinn.%.lnl: quinn.%.laml_params.txt quinn.%.var.log quinn.%.beam.log
	grep '^Negative' quinn.$*.laml_params.txt | awk '{printf "%d\t-%f\t", $*, $$2}' > $@
	tail -1 quinn.$*.var.log | awk '{printf "%f\t%f\t", $$11, $$21}' >> $@
	grep -v '^#' quinn.$*.beam.log | grep -v '^Sample' | \
		awk '{a=$$(NF-4); b=$$(NF-3); s=a+b; if (NR==1 || s>max) {max=s; A=a; B=b}} END {printf "%f\t%f\n", A, B}' >> $@

eval.all.lnl.txt: $(LNLS)
	cat $(LNLS) | awk 'BEGIN{printf "cp\tlaml_tree\tvine_tree\tvine_mig\tbeam_tree\tbeam_mig\n"} {x1 += $$2; x2 += $$3; x3 += $$4; x4 += $$5; x5 += $$6; printf "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$4, $$5, $$6} END {printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", x1/NR, x2/NR, x3/NR, x4/NR, x5/NR}' > $@

eval.all.lnl.pdf: eval.all.lnl.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_lnl_allSizes_vineBeamOnly_realData.R /mnt/files/eval.all.lnl.txt /mnt/files/eval.all.lnl.pdf

# Get quantitative tree comparisons
quinn.%.stats: quinn.%.cass.nwk quinn.%.laml_trees.nwk quinn.%.beam.mcc.nwk quinn.%.var.mcc.nwk quinn.%.tissues.expanded.csv quinn.%.tsv
	fields=("$*"); \
	for tree in quinn.$*.cass.nwk quinn.$*.laml_trees.nwk quinn.$*.beam.mcc.nwk quinn.$*.var.mcc.nwk; do \
		echo "Calculating stats for $$tree"; \
		cc=$$(singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/cophenetic_correlation.py /mnt/files/$$tree /mnt/files/quinn.$*.tsv "-1" | awk '{print $$2}'); \
		pm=$$(singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/fitch_parsimony_count_mutations.py /mnt/files/$$tree /mnt/files/quinn.$*.tsv "0" | awk '{print $$2}'); \
		thi=$$(singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/tissue_homogeneity_index.py /mnt/files/$$tree /mnt/files/quinn.$*.tissues.expanded.csv | awk '{print $$2}'); \
		pt=$$(singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/fitch_parsimony_count_tissues.py /mnt/files/$$tree /mnt/files/quinn.$*.tissues.expanded.csv "LL" | awk '{print $$2}'); \
		fields+=($$cc $$pm $$thi $$pt); \
	done; \
	printf '%s\t' "$${fields[@]}" | sed 's/\t$$/\n/' > $@

eval.all.stats.txt: $(STATS)
	cat $(STATS) | awk 'BEGIN{printf "cp\tcass_cc\tcass_pm\tcass_thi\tcass_pt\tlaml_cc\tlaml_pm\tlaml_thi\tlaml_pt\tbeam_cc\tbeam_pm\tbeam_thi\tbeam_pt\tvine_cc\tvine_pm\tvine_thi\tvine_pt\n"} {x1+=$$2;x2+=$$3;x3+=$$4;x4+=$$5;x5+=$$6;x6+=$$7;x7+=$$8;x8+=$$9;x9+=$$10;x10+=$$11;x11+=$$12;x12+=$$13;x13+=$$14;x14+=$$15;x15+=$$16;x16+=$$17;printf "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",$$1,$$2,$$3,$$4,$$5,$$6,$$7,$$8,$$9,$$10,$$11,$$12,$$13,$$14,$$15,$$16,$$17} END{printf "-----------------------------------------\navg\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",x1/NR,x2/NR,x3/NR,x4/NR,x5/NR,x6/NR,x7/NR,x8/NR,x9/NR,x10/NR,x11/NR,x12/NR,x13/NR,x14/NR,x15/NR,x16/NR}' > $@

eval.all.stats.pdf: eval.all.stats.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_tree_stats_realData.R /mnt/files/eval.all.stats.txt /mnt/files/eval.all.stats.pdf --hide-cass-cc

# Get quantitative migration graph comparisons
quinn.%.graph.stats: quinn.%.metient_probability_graph.csv quinn.%.beam_probability_graph.csv quinn.%.var_probability_graph.csv
	threshold=0.5; \
	fields=("$*"); \
	singularity_cmd="singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF)"; \
	for graph in $^; do \
		echo "Calculating stats for $$graph"; \
		read nm ncm nss <<< $$($$singularity_cmd python3 /mnt/scripts/count_migration_comigrations_seedingsites.py /mnt/files/$$graph $$threshold | awk '{print $$2, $$4, $$6}'); \
		fields+=($$nm $$ncm $$nss); \
	done; \
	printf '%s\t' "$${fields[@]}" | sed 's/\t$$/\n/' > $@

eval.all.graph.stats.txt: $(GRAPHSTATS)
	cat $(GRAPHSTATS) | awk 'BEGIN{printf "cp\tmetient_mig\tmetient_comig\tmetient_seedsites\tbeam_mig\tbeam_comig\tbeam_seedsites\tvar_mig\tvar_comig\tvar_seedsites\n"} {x1+=$$2;x2+=$$3;x3+=$$4;x4+=$$5;x5+=$$6;x6+=$$7;x7+=$$8;x8+=$$9;x9+=$$10;printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$$1,$$2,$$3,$$4,$$5,$$6,$$7,$$8,$$9,$$10} END{printf "-----------------------------------------\navg\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",x1/NR,x2/NR,x3/NR,x4/NR,x5/NR,x6/NR,x7/NR,x8/NR,x9/NR}' > $@

eval.all.graph.stats.pdf: eval.all.graph.stats.txt
	singularity exec --bind $(R_SRC):/mnt/scripts --bind ./:/mnt/files $(RPLOTTING_SIF) Rscript /mnt/scripts/plot_graph_stats_realData.R /mnt/files/eval.all.graph.stats.txt /mnt/files/eval.all.graph.stats.pdf


clean:
	rm -f quinn.*

clean_not_converged_beam:
	for cp in $(QUINNCPS); do \
		if [ -f quinn.$$cp.beam.term ]; then \
			ess=$$(grep "Msamples" quinn.$$cp.beam.term | tail -n 1 | awk '{print $$3}'); \
			if [ -z "$$ess" ] || awk -v ess="$$ess" -v target=$(BEAM_ESS_TARGET) 'BEGIN{exit !(ess < target)}'; then \
				echo "Deleting quinn.$$cp (ESS=$$ess)"; \
				rm -f quinn.$$cp.beam.term \
				      quinn.$$cp.beam_probability_graph.csv \
				      quinn.$$cp.beam*.pdf; \
			fi; \
		fi; \
	done

# Target to keep resuming BEAM runs until they converge, rather than having to know the mcmc chain length ahead of time. This is a little clunky but works.
.PHONY: repeat_beam_runs_until_done
J ?= 10
repeat_beam_runs_until_done:
	@while ! $(MAKE) -q; do $(MAKE) beam_infer -j $(J); $(MAKE) clean_not_converged_beam; done; $(MAKE) beam -j $(J); echo "Nothing to be done. Exiting."

archive_all:
	archive_dir=archive_all.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv quinn.* $$archive_dir;

archive_vine:
	archive_dir=archive_vine.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv quinn.*.var* $$archive_dir; \
	mv eval.all* $$archive_dir; \
	mv quinn.*.time $$archive_dir; \
	mv quinn.*.stats $$archive_dir; \
	mv quinn.*.graph.stats $$archive_dir; \
	mv quinn.*.lnl $$archive_dir;

