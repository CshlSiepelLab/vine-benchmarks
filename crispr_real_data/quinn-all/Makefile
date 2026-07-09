export SHELL=/usr/bin/bash

QUINNROOT := ../quinn_cps/5k

# Prevent Make from removing files it thinks it should clean up
.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# Absolute paths
MAIN_DIR := /local/storage/no-backup/vine-benchmarks

# Relative paths
PYTHON_SRC := $(MAIN_DIR)/python/src
BIN := $(MAIN_DIR)/bin
PHAST_BIN := $(BIN)/phast/bin
VINE_BIN := $(BIN)/vine/bin
#VINE_BIN := /home/staklins/projects/vine_project/vine/bin
CONTAINERS := $(MAIN_DIR)/containers
LAML_SIF := $(CONTAINERS)/laml/laml.sif
GRAPHPOSTERIOR_SIF := $(CONTAINERS)/graphposterior/graphposterior.sif
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif

# Excluded in the original paper barcode quality control: 5 16 18 25 33 38 39 41 50 53 65 69 75 81 87 88 93
# Originally excluded because they were too large due to VINE memory requirements, but now can be solved: 1 2 3
QUINNCPS := 4 6 7 8 9 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 28 29 30 31 32 34 35 36 37 40 42 43 44 45 46 47 48 49 51 52 54 55 56 57 58 59 60 61 62 63 64 66 67 68 70 71 72 73 74 76 77 78 79 80 82 83 84 85 86 89 90 91 92 94 95 96 97 98 99 100


MATRICES := $(QUINNCPS:%=quinn.%.tsv)
MATRICES_CSV := $(QUINNCPS:%=quinn.%.csv)
CASSTREES := $(QUINNCPS:%=quinn.%.cass.nwk)
LAMLTREES := $(QUINNCPS:%=quinn.%.laml_trees.nwk)
VARTREES := $(QUINNCPS:%=quinn.%.var.nwk)
VARLOG := $(QUINNCPS:%=quinn.%.var.log)
TIMES := $(QUINNCPS:%=quinn.%.time)
LNLS := $(QUINNCPS:%=quinn.%.lnl)
CHECKLNLS := $(QUINNCPS:%=quinn.%.checklnl)

VAROPT?=-v 0 

all: $(VARLOG) $(LAMLTREES) $(VARTREES) $(TIMES) $(LNLS) summary.lnl.txt summary.time.txt

laml: $(LAMLTREES)
vine: $(VARLOG) $(VARTREES)
infer: $(VARLOG) $(LAMLTREES)
check: $(CHECKLNLS) summary.checklnl.txt

vineclean:
	rm -rf $(VARLOG) $(VARTREES) $(TIMES) $(LNLS) summary.lnl.txt summary.time.txt

# these targets just copy source data and reformat
quinn.%.original.csv:
	sed 's/,r[0-9]*$$//' $(QUINNROOT)/`echo $@ | awk -F\. '{print $$2}'`/matrix.csv > $@

quinn.%.csv quinn.%.collapsing_map.tsv: quinn.%.original.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(GRAPHPOSTERIOR_SIF) python3 /mnt/scripts/collapse_matrix_duplicates.py \
	 -m /mnt/files/quinn.$*.original.csv \
	 -om /mnt/files/quinn.$*.csv \
	 -on /mnt/files/quinn.$*.collapsing_map.tsv	

quinn.%.tsv: quinn.%.csv
	$(BIN)/csv-to-tsv.sh $< > $@

quinn.%.cass.nwk: quinn.%.csv
	singularity exec --bind $(PYTHON_SRC):/mnt/scripts --bind ./:/mnt/files $(CASSIOPEIA_SIF) \
	python /mnt/scripts/cassiopeiaGreedy.py /mnt/files/quinn.$*.csv /mnt/files/quinn.$*.cass.nwk

source-data: $(MATRICES) $(CASSTREES)

# Added --noultrametric to the LAML v1.0.5 run since it gives ultrametric trees. Seems like a big they will need to fix.
quinn.%.laml_trees.nwk quinn.%.laml-time: quinn.%.csv quinn.%.cass.nwk
	singularity exec --bind $(CURDIR):/mnt $(LAML_SIF) /usr/bin/time -o /mnt/quinn.$*.laml-time \
	run_laml -c /mnt/quinn.$*.csv -t /mnt/quinn.$*.cass.nwk -o /mnt/quinn.$*.laml --topology_search --noDropout --noultrametric

laml: $(LAMLTREES)

quinn.%.var.log quinn.%.var.nwk quinn.%.var-time quinn.%.mean.nwk: quinn.%.tsv
	/usr/bin/time -o quinn.$*.var-time $(VINE_BIN)/vine $(VAROPT) -i CRISPR quinn.$*.tsv --logf quinn.$*.var.log --mean quinn.$*.mean.nwk > quinn.$*.var.nwk

quinn.%.lnl: quinn.%.laml_trees.nwk quinn.%.var.log
	grep '^Negative' quinn.$*.laml_params.txt | awk '{printf "%d\t-%f\t", $*, $$2}' > $@
	tail -1 quinn.$*.var.log | awk '{printf "%f\n", $$19}' >> $@

quinn.%.time: quinn.%.laml_trees.nwk quinn.%.var-time
	grep '^Runtime' quinn.$*.laml.log | tail -1 | awk '{printf "%d\t%f\t", $*, $$3}' > $@
	awk '{printf "%f\n", $$1}' quinn.$*.var-time >> $@

summary.lnl.txt: $(LNLS)
	cat $(LNLS) | awk 'BEGIN{printf "cp\tlaml\tvine\tdiff\n"} {x1 += $$2; x2 += $$3; n++; printf "%s\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$3-$$2} END {printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\n", x1/n, x2/n, (x2-x1)/n}' > $@

summary.time.txt: $(TIMES)
	cat $(TIMES) | awk 'BEGIN{printf "cp\tlaml\tvine\tspeedup\n"} NF==3 {x1 += $$2; x2 += $$3; n++; printf "%s\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$2/$$3} END {printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\n", x1/n, x2/n, x1/x2}' > $@

quinn.%.tr: quinn.%.var.log
	grep -v '^#' $^ > $@

var: $(VARLOG)

quinn.%.checklnl: quinn.%.indels.tsv quinn.%.laml_trees.nwk
	$(VINE_BIN)/crisprLnl quinn.$*.tsv quinn.$*.laml_trees.nwk quinn.$*.laml_params.txt > $@

summary.checklnl.txt: $(CHECKLNLS)
	printf "%17s %15s %15s %15s %15s\n" "file" "VINElnl" "LAMLlnl" "diff" "%diff" > $@
	for file in $^ ; do \
		printf "%17s " $${file} >> $@ ;\
		awk '{printf "%15.4f ", $$NF} END {printf "\n"}' $${file} >> $@ ;\
	done

clean:
	rm -f $(MATRICES) $(CASSTREES) $(LAMLTREES) $(VARLOG) $(TIMES) $(LNLS) quinn.*.var.log quinn.*.mean.nwk quinn.*.csv quinn.*.laml_* summary.lnl.txt summary.time.txt

.PRECIOUS: $(MATRICES) $(CASSTREES) $(LAMLTREES) $(VARLOG) $(MATRICES_CSV)


archive_all:
	archive_dir=archive_all.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv quinn.* $$archive_dir; \
	mv summary.lnl.txt $$archive_dir; \
	mv summary.time.txt $$archive_dir;

archive_vine:
	archive_dir=archive_vine.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv quinn.*.var.log $$archive_dir; \
	mv quinn.*.var.nwk $$archive_dir; \
	mv quinn.*.var-time $$archive_dir; \
	mv summary.lnl.txt $$archive_dir; \
	mv summary.time.txt $$archive_dir;


