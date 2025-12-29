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
CONTAINERS := $(MAIN_DIR)/containers
LAML_SIF := $(CONTAINERS)/laml/laml.sif

# Copy of below new cutoff CPs just for testing since some in the full set fail with VINE for now due to a bug
#QUINNCPS := 7 10 14 15 17 19 20 21 22 23 24 26 27 28 29 30 31 32 34 35 36 37 40 42 43 44 45 46 47 48 49 51 52 54 55 56 57 58 59 60 61 62 63 64 66 67 68 70 71 72 73 74 76 77 78 79 80 82 83 84 85 86 89 90 91 92 94 95 96 97 98 99 100

# New cutoff for CPs with around ~650 tips or less.
# QUINNCPS := 6 7 8 9 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 28 29 30 31 32 34 35 36 37 40 42 43 44 45 46 47 48 49 51 52 54 55 56 57 58 59 60 61 62 63 64 66 67 68 70 71 72 73 74 76 77 78 79 80 82 83 84 85 86 89 90 91 92 94 95 96 97 98 99 100

# Old cutoff for only CPs with number >= 30
QUINNCPS := 30 31 32 34 35 36 37 40 42 43 44 45 46 47 48 49 51 52 54 55 56 57 58 59 60 61 62 63 64 66 67 68 70 71 72 73 74 76 77 78 79 80 82 83 84 85 86 89 90 91 92 94 95 96 97 98 99 100


MATRICES := $(QUINNCPS:%=quinn.%.tsv)
MATRICES_CSV := $(QUINNCPS:%=quinn.%.csv)
CASSTREES := $(QUINNCPS:%=quinn.%.cass.nwk)
LAMLTREES := $(QUINNCPS:%=quinn.%.laml_trees.nwk)
VARLOG := $(QUINNCPS:%=quinn.%.var.log)
TIMES := $(QUINNCPS:%=quinn.%.time)
LNLS := $(QUINNCPS:%=quinn.%.lnl)

VAROPT?=-M 400 -b 10 -c 100 -D 6 -v 0 --taylor

all: $(VARLOG) $(LAMLTREES) $(TIMES) $(LNLS) summary.lnl.txt summary.time.txt

laml: $(LAMLTREES)
vine: $(VARLOG)
infer: $(VARLOG) $(LAMLTREES)

# these targets just copy source data and reformat
quinn.%.csv:
	sed 's/,r[0-9]*$$//' $(QUINNROOT)/`echo $@ | awk -F\. '{print $$2}'`/matrix.csv > $@

quinn.%.tsv: quinn.%.csv
	$(BIN)/csv-to-tsv.sh $< > $@

quinn.%.cass.nwk:
	cp $(QUINNROOT)/`echo $@ | awk -F\. '{print $$2}'`/cassiopeia_greedy_inferred.nwk $@

source-data: $(MATRICES) $(CASSTREES)

# Added --noultrametric to the LAML v1.0.5 run since it gives ultrametric trees. Seems like a big they will need to fix.
quinn.%.laml_trees.nwk: quinn.%.csv quinn.%.cass.nwk
	singularity exec --bind $(CURDIR):/mnt $(LAML_SIF) \
	run_laml -c /mnt/quinn.$*.csv -t /mnt/quinn.$*.cass.nwk -o /mnt/quinn.$*.laml --topology_search --noDropout --noultrametric

laml: $(LAMLTREES)

quinn.%.var.log quinn.%.var-time: quinn.%.tsv quinn.%.cass.nwk
	/usr/bin/time -o quinn.$*.var-time $(VINE_BIN)/vine $(VAROPT) -i CRISPR quinn.$*.tsv --log quinn.$*.var.log --mean quinn.$*.mean.nwk > $@
#	/usr/bin/time -o quinn.$*.var-time $(VINE_BIN)/vine $(VAROPT) -i CRISPR quinn.$*.tsv --log quinn.$*.var.log --tree quinn.$*.cass.nwk --mean quinn.$*.mean.nwk > $@

quinn.%.lnl: quinn.%.laml_trees.nwk quinn.%.var.log
	grep '^Negative' quinn.$*.laml_params.txt | awk '{printf "%d\t-%f\t", $*, $$2}' > $@
	tail -1 quinn.$*.var.log | awk '{printf "%f\n", $$11}' >> $@

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

clean:
	rm -f $(MATRICES) $(CASSTREES) $(LAMLTREES) $(VARLOG) $(TIMES) $(LNLS) quinn.*.var.log quinn.*.mean.nwk quinn.*.csv quinn.*.laml_* summary.lnl.txt summary.time.txt

.PRECIOUS: $(MATRICES) $(CASSTREES) $(LAMLTREES) $(VARLOG) $(MATRICES_CSV)


archive_laml:
	archive_dir=archive_laml.$$(date +%Y-%m-%d_%H:%M:%S); \
	mkdir -p $$archive_dir; \
	mv quinn.*.laml_trees.nwk $$archive_dir; \
	mv quinn.[0-9].laml.log $$archive_dir; \
	mv quinn.*.laml_ckpt.txt $$archive_dir; \
	mv quinn.*.laml_params.txt $$archive_dir; \
	mv quinn.*.laml_annotations.txt $$archive_dir; \
	mv quinn.*.time $$archive_dir; \
	mv quinn.*.lnl $$archive_dir; \
	mv summary.lnl.txt $$archive_dir; \
	mv summary.time.txt $$archive_dir;


