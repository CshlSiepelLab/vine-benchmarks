export SHELL=/usr/bin/bash

# Prevent Make from removing files it thinks it should clean up
# .NOTINTERMEDIATE:   # conflicts with .SECONDARY on some GNU make versions
.SECONDARY:
.PRECIOUS:

# edit for local structure; this is the only place absolute paths are used
MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT_SUFFIX := dna_sims/jc69_300sites

ROOT := $(MAIN_DIR)/$(ROOT_SUFFIX)
BIN := $(MAIN_DIR)/bin
PYTHON_SRC := $(MAIN_DIR)/python/src
VINE_BIN := $(BIN)/vine/bin
PHAST_BIN := $(MAIN_DIR)/phast/bin
OTHER_BIN := $(MAIN_DIR)/bin
BEAST_BIN := $(BIN)/beast/bin
BEAST := $(BEAST_BIN)/beast
MRBAYES := $(OTHER_BIN)/mb
BEAST_TEMPLATE := $(ROOT)/beast_template.xml
DODONAPHY_SIF := $(MAIN_DIR)/dodonaphy/dodonaphy.sif
GEOPHY_SIF := $(MAIN_DIR)/geophy/geophy.sif
GEOPHY_CONFIG := $(MAIN_DIR)/geophy/default.yaml
VAIPHY_SIF := $(MAIN_DIR)/vaiphy/vaiphy.sif
BURNIN_PCT := 10
METRIC_TREE_COUNT := 1000

# Log approximately 10,000 MCMC samples for each program.
BEAST_SAMPLEFREQ := $(shell awk 'BEGIN {f=int($(BEAST_MCMCLEN)/10000); print (f > 0 ? f : 1)}')
MRBAYES_SAMPLEFREQ := $(shell awk 'BEGIN {f=int($(MRBAYES_MCMCLEN)/10000); print (f > 0 ? f : 1)}')

BEAST_TEMPLATE := $(ROOT)/beast_template.xml
TREES := $(shell seq -f tree.%.0f.true.nwk 1 $(NSAMP))
FA := $(patsubst %.true.nwk,%.fa,$(TREES))
MOD := $(patsubst %.true.nwk,%.true.mod,$(TREES))
MLMOD := $(patsubst %.true.nwk,%.ml.mod,$(TREES))
NJMOD := $(patsubst %.true.nwk,%.nj.mod,$(TREES))
NJ := $(patsubst %.true.nwk,%.nj.nwk,$(TREES))
ML := $(patsubst %.true.nwk,%.ml.nwk,$(TREES))
VAR := $(patsubst %.true.nwk,%.var.nwk,$(TREES))
EVALRF := $(patsubst tree.%.true.nwk,tree.%.rf,$(TREES))
EVALMF := $(patsubst tree.%.true.nwk,tree.%.mf,$(TREES))
EVALDIST := $(patsubst tree.%.true.nwk,tree.%.dist,$(TREES))
LNL := $(patsubst tree.%.true.nwk,tree.%.lnl,$(TREES))
TIME := $(patsubst tree.%.true.nwk,tree.%.time,$(TREES))
VARLOG := $(patsubst %.true.nwk,%.var.nwk.log,$(TREES))
VARTIME := $(patsubst %.true.nwk,%.var-time,$(TREES))
TRACER := $(patsubst %.true.nwk,%.tr,$(TREES))

XML := $(patsubst %.true.nwk,%.beast.xml,$(TREES))
BEASTLOG := $(patsubst %.true.nwk,%.beast.log,$(TREES))
BEASTXML := $(patsubst %.true.nwk,%.beast.xml,$(TREES))
XMLSTATE := $(patsubst %.true.nwk,%.beast.xml.state,$(TREES))
BEASTTREES := $(patsubst %.true.nwk,%.beast-tree.trees,$(TREES))
BEASTNWK := $(patsubst %.true.nwk,%.beast.nwk,$(TREES))
BEASTLNL := $(patsubst %.true.nwk,%.beastlnl,$(TREES))
BEASTTERM := $(patsubst %.true.nwk,%.beast.term,$(TREES))
BEASTBEAGLELOG := $(patsubst %.true.nwk,%.beast-beagle.log,$(TREES))
BEASTBEAGLENWK := $(patsubst %.true.nwk,%.beast-beagle.nwk,$(TREES))

MRBAYESLOG := $(patsubst %.true.nwk,%.mrbayes.nex.p,$(TREES))
MRBAYESBEAGLELOG := $(patsubst %.true.nwk,%.mrbayes-beagle.nex.p,$(TREES))
MRBAYESBEAGLENWK := $(patsubst %.true.nwk,%.mrbayes-beagle.nwk,$(TREES))

# evalTrees stuff
FAHELDOUT := $(patsubst %.true.nwk,%.heldout.fa,$(TREES))

.PHONY: beast-beagle beast-beagle-nwk mrbayes-beagle mrbayes-beagle-nwk

all: eval.all.lnl.txt eval.all.rf.txt eval.all.mf.txt eval.all.time.txt eval.all.dist.txt

mcmc: $(BEASTLOG) $(BEASTBEAGLELOG) $(MRBAYESLOG) $(MRBAYESBEAGLELOG)
beast: $(BEASTLOG)
beast-beagle: $(BEASTBEAGLELOG)
mrbayes: $(MRBAYESLOG)
mrbayes-beagle: $(MRBAYESBEAGLELOG)
vine: $(VARLOG)

tree.%.true.nwk: 
	$(OTHER_BIN)/bdTree3 -b 1 -d 0.5 --oversample-k 3 --height 5 --min-edge 0.02 --expected-height $(EXPHEIGHT) --no-stem --ucln-sd 0.6 --target-stat median -n $(NTAXA) | sed 's/\[\&[UR]\] //' > $@


tree.%.fa: tree.%.true.nwk
		cp ../base-jc69.mod tmp.mod
		echo -n "TREE: " >> tmp.mod
		cat $< >> tmp.mod
		$(PHAST_BIN)/base_evolve --nsites $(NSITES) tmp.mod > $@
		rm tmp.mod

tree.%.heldout.fa: tree.%.true.nwk
		cp ../base-jc69.mod tmp.mod
		echo -n "TREE: " >> tmp.mod
		cat $< >> tmp.mod
		$(PHAST_BIN)/base_evolve --nsites $(NSITES) tmp.mod > $@
		rm tmp.mod

tree.%.nj.nwk: tree.%.fa
		$(PHAST_BIN)/vine --nj-only $< > $@

tree.%.ml.mod: tree.%.nj.nwk tree.%.fa
	$(PHAST_BIN)/phyloFit --subst-mod JC69 --tree $^ -o tree.$*.ml

tree.%.ml.nwk: tree.%.ml.mod
	$(PHAST_BIN)/tree_doctor --tree-only $^ > $@

tree.%.var.nwk tree.%.var-time tree.%.var.nwk.log: tree.%.fa 
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine $< -l tree.$*.var.nwk.log $(VAROPT) --mean tree.$*.mean.nwk > tree.$*.var.nwk

tree.%.beast.xml:
	cp "$(BEAST_TEMPLATE)" $@

tree.%.beast.term tree.%.beast-tree.trees tree.%.beast.log: tree.%.beast.xml tree.%.fa
	rm -f tree.$*.beast-tree.trees tree.$*.beast.log
	"$(BEAST)" -java -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$(BEAST_SAMPLEFREQ) -D printfreq=$(MCMC_PRINTFREQ) $< > tree.$*.beast.term

tree.%.beast-beagle.xml:
	cp "$(BEAST_TEMPLATE)" $@

tree.%.beast-beagle.term tree.%.beast-beagle-tree.trees tree.%.beast-beagle.log: tree.%.beast-beagle.xml tree.%.fa
	rm -f tree.$*.beast-beagle-tree.trees tree.$*.beast-beagle.log
	"$(BEAST)" -beagle -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$(BEAST_SAMPLEFREQ) -D printfreq=$(MCMC_PRINTFREQ) $< > tree.$*.beast-beagle.term

tree.%.beast.nwk: tree.%.beast-tree.trees
	$(BEAST_BIN)/logcombiner -log $^ -o thinned.trees.$*.beast \
	  -burnin $(BURNIN_PCT) -resample 5000
	python3 "$(PYTHON_SRC)/time2subs.py" thinned.trees.$*.beast tmp.$*.beast.nex
	$(OTHER_BIN)/convertTrees -i nexus tmp.$*.beast.nex > $@
	rm -f tmp.$*.beast.nex

tree.%.beast-beagle.nwk: tree.%.beast-beagle-tree.trees
	THIN=$$(( ($(BEAST_MCMCLEN) / $(BEAST_SAMPLEFREQ)) * (100 - $(BURNIN_PCT)) / 100 / $(METRIC_TREE_COUNT) )); \
	if [ "$$THIN" -lt 1 ]; then THIN=1; fi; \
	RESAMPLE=$$(( $(BEAST_SAMPLEFREQ) * THIN )); \
	$(BEAST_BIN)/logcombiner -log $^ -o thinned.trees.$*.beast-beagle \
	  -burnin $(BURNIN_PCT) -resample $$RESAMPLE
	python3 "$(PYTHON_SRC)/time2subs.py" thinned.trees.$*.beast-beagle tmp.$*.beast-beagle.nex
	$(OTHER_BIN)/convertTrees -i nexus tmp.$*.beast-beagle.nex > tmp.$*.beast-beagle.nwk
	N=$$(awk 'NF{n++} END{print n}' tmp.$*.beast-beagle.nwk); \
	TARGET=$(METRIC_TREE_COUNT); \
	if [ "$$N" -lt "$$TARGET" ]; then TARGET=$$N; fi; \
	awk -v n="$$N" -v target="$$TARGET" \
	  'NF{j++; if (int(j*target/n) > int((j-1)*target/n)) print}' \
	  tmp.$*.beast-beagle.nwk > $@
	rm -f tmp.$*.beast-beagle.nex tmp.$*.beast-beagle.nwk

tree.%.nex: tree.%.fa
	$(OTHER_BIN)/fa2nex $< $@

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	$(OTHER_BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex --mcmc_length $(MRBAYES_MCMCLEN) --model JC69 \
		--sample_freq $(MRBAYES_SAMPLEFREQ) --print_freq $(MCMC_PRINTFREQ) --diagn_freq $(MCMC_PRINTFREQ)

tree.%.mrbayes-beagle.nex: tree.%.nex
	$(OTHER_BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes-beagle.nex --mcmc_length $(MRBAYES_MCMCLEN) --model JC69 \
		--sample_freq $(MRBAYES_SAMPLEFREQ) --print_freq $(MCMC_PRINTFREQ) --diagn_freq $(MCMC_PRINTFREQ) --use_beagle

# Run MrBayes
tree.%.mrbayes.term tree.%.mrbayes.nex.p tree.%.mrbayes.nex.t: tree.%.mrbayes.nex
	$(MRBAYES) tree.$*.mrbayes.nex > tree.$*.mrbayes.term

tree.%.mrbayes-beagle.term tree.%.mrbayes-beagle.nex.p tree.%.mrbayes-beagle.nex.t: tree.%.mrbayes-beagle.nex
	$(MRBAYES) tree.$*.mrbayes-beagle.nex > tree.$*.mrbayes-beagle.term

# Get mrbayes tree in nexus format
tree.%.mrbayes.nwk: tree.%.mrbayes.nex.t
	SKIP=$$(awk -v p=$(BURNIN_PCT) \
	  '/^[[:space:]]*tree[[:space:]]+/{c++} END{printf "%d", int(c*p/100)}' $<); \
	awk -v skip="$$SKIP" '\
	  /^[[:space:]]*tree[[:space:]]+/{ \
	    if (++c <= skip) next; \
	    if ((c - skip) % 200) next \
	  } \
	  1' \
	  $< > tree.$*.mrbayes.thinned.t
	$(OTHER_BIN)/convertTrees -i nexus \
	  tree.$*.mrbayes.thinned.t \
	  | sed 's/^\[&[^]]*\]\s*//' > $@
	rm -f tree.$*.mrbayes.thinned.t

tree.%.mrbayes-beagle.nwk: tree.%.mrbayes-beagle.nex.t
	SKIP=$$(awk -v p=$(BURNIN_PCT) \
	  '/^[[:space:]]*tree[[:space:]]+/{c++} END{printf "%d", int(c*p/100)}' $<); \
	N=$$(awk '/^[[:space:]]*tree[[:space:]]+/{c++} END{print c}' $<); \
	REMAINING=$$(( N - SKIP )); \
	TARGET=$(METRIC_TREE_COUNT); \
	if [ "$$REMAINING" -lt "$$TARGET" ]; then TARGET=$$REMAINING; fi; \
	awk -v skip="$$SKIP" -v remaining="$$REMAINING" -v target="$$TARGET" '\
	  /^[[:space:]]*tree[[:space:]]+/{ \
	    if (++c <= skip) next; \
	    j = c - skip; \
	    if (int(j*target/remaining) == int((j-1)*target/remaining)) next \
	  } \
	  1' \
	  $< > tree.$*.mrbayes-beagle.thinned.t
	$(OTHER_BIN)/convertTrees -i nexus \
	  tree.$*.mrbayes-beagle.thinned.t \
	  | sed 's/^\[&[^]]*\]\s*//' > $@
	rm -f tree.$*.mrbayes-beagle.thinned.t

# Get NJ tree in nexus format
tree.%.nj.nex: tree.%.nj.nwk
	$(OTHER_BIN)/nwk2nex $< $@

# Run dodonaphy
tree.%.dodonaphy.term tree.%.dodonaphy.elbo.txt tree.%.dodonaphy-time: tree.%.nex tree.%.nj.nex
	mkdir tree.$*.dodonaphy
	cp $< tree.$*.dodonaphy/tree.$*.nex
	singularity exec --bind $(CURDIR):/mnt $(DODONAPHY_SIF) \
		env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
		/usr/bin/time -o /mnt/tree.$*.dodonaphy-time \
		dodo \
		--path_root /mnt/tree.$*.dodonaphy \
		--path_dna /mnt/tree.$*.nex  \
		--start /mnt/tree.$*.nj.nex \
		--infer vi \
		--model JC69 \
		--embed up \
		--temp 0.00001 \
		--prior exponential \
		--connect nj \
		--boosts 3 \
		--importance 3 \
		--curv -100 \
		--learn 0.1 \
		--epochs 2000 \
		--draws 1000 \
		--overwrite > tree.$*.dodonaphy.term
	rm -f tree.$*.dodonaphy/tree.$*.nex
	find ./tree.$*.dodonaphy/ -type f -exec bash -c 'for file; do mv "$$file" ./tree.$*.dodonaphy.$$(basename "$$file"); done' _ {} +
	rm -rf tree.$*.dodonaphy


# Run GeoPhy
tree.%.geophy.term tree.%.geophy-time tree.%.geophy.eval.latest.txt: tree.%.nex
	cp $< tree.$*.geophy.nex
	cp $(GEOPHY_CONFIG) tree.$*.geophy.config.yaml
	/usr/bin/time -o tree.$*.geophy-time \
		singularity exec \
		--bind $(CURDIR):/mnt \
		$(GEOPHY_SIF) \
		bash /opt/app/scripts/run_geophy.sh \
		-i /mnt/tree.$*.geophy.nex \
		-o /mnt/tree.$*.geophy \
		-c /mnt/tree.$*.geophy.config.yaml > tree.$*.geophy.term
	rm -f tree.$*.geophy.nex

tree.%.vaiphy.term: tree.%.nex
	mkdir -p tree.$*
	cp tree.$*.nex tree.$*/tree.$*.nex
	/usr/bin/time -o tree.$*.vaiphy-time \
	singularity exec \
		--bind $(CURDIR):/mnt \
		--bind /local/storage/no-backup/vine-benchmarks/VaiPhy:/opt/vaiphy \
		--pwd /mnt \
		--env PYTHONWARNINGS=ignore \
		$(VAIPHY_SIF) \
		python -W ignore \
			/opt/vaiphy/src/main.py \
			--data_path /mnt/ \
			--result_path /mnt/ \
			--max_iter 200 \
			--n_particles 128 \
		--dataset tree.$* > tree.$*.vaiphy.log 2> tree.$*.vaiphy.err
	echo done > $@

# extract training likelihoods
tree.%.true.mod: tree.%.true.nwk tree.%.fa
	cp ../base-jc69.mod tmp.mod
	echo -n "TREE: " >> tmp.mod
	cat $< >> tmp.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tmp.mod -o tree.$*.true tree.$*.fa
	rm tmp.mod

tree.%.nj.mod: tree.%.nj.nwk tree.%.fa
	cp ../base-jc69.mod tmp.mod
	echo -n "TREE: " >> tmp.mod
	cat $< >> tmp.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tmp.mod -o tree.$*.nj tree.$*.fa
	rm tmp.mod

tree.%.modlnl: tree.%.true.mod tree.%.nj.mod tree.%.ml.mod 
	rm -f $@
	for file in $^ ; do \
		echo -n "$${file} " >> $@ ;\
		grep LNL $${file} | awk '{print $$2}' >> $@ ;\
	done

tree.%.varlnl: tree.%.var.nwk.log 
	echo -n "$^ " > $@
	tail -1 $^ | awk '{print $$11}' >> $@ 

tree.%.beastlnl: tree.%.beast.log
	echo -n "$^ " > $@
	grep -v '^#' $^ | grep -v '^Sample' | awk '{print $$3}' | sort -nr | head -1 >> $@

tree.%.beast-beaglelnl: tree.%.beast-beagle.log
	echo -n "$^ " > $@
	grep -v '^#' $^ | grep -v '^Sample' | awk '{print $$3}' | sort -nr | head -1 >> $@

tree.%.mrbayeslnl: tree.%.mrbayes.nex.p tree.%.mrbayes.term
	echo -n "$< " > $@
	grep -v '^\[' $< | grep -v '^Gen' | awk '{print $$2}' | sort -gr | head -1 | awk '{printf "%.6f\n", $$1}' >> $@

tree.%.mrbayes-beaglelnl: tree.%.mrbayes-beagle.nex.p tree.%.mrbayes-beagle.term
	echo -n "$< " > $@
	grep -v '^\[' $< | grep -v '^Gen' | awk '{print $$2}' | sort -gr | head -1 | awk '{printf "%.6f\n", $$1}' >> $@

tree.%.dodonaphylnl: tree.%.dodonaphy.elbo.txt
	echo -n "$^ " > $@
	awk 'NR==1 || $$1>max {max=$$1} END {printf "%.6f\n", max}' $^ >> $@

tree.%.geophylnl: tree.%.geophy.eval.latest.txt
	echo -n "$^ " > $@
	awk 'NR==2 {printf "%.6f\n", $$6}' $^ >> $@

tree.%.vaiphylnl: tree.%.vaiphy.log
	echo -n "$^ " > $@
	# Prefer maximum sampled LL (selected, bifurcated); fall back to selected
	( grep -F "Selected tree log-likelihood (bifurcated):" $^ || true ) \
		| tail -1 \
		| awk '{printf "%.6f\n", $$NF}' >> $@
 
tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.beast-beaglelnl tree.%.mrbayeslnl tree.%.mrbayes-beaglelnl tree.%.dodonaphylnl tree.%.geophylnl tree.%.vaiphylnl
	cat $^ | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@


eval.all.lnl.txt: $(LNL)
	echo -e "true\tnj\tml\tvine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle\tdodonaphy\tgeophy\tvaiphy" > tmp
	for file in $^ ; do \
		awk '{val=$$2; gsub(/,/, "", val); printf "%s\t", val}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{for(i=1;i<=11;i++) x[i]+=$$i; print $$0} END {printf "-----\n"; for(i=1;i<=11;i++) printf "%f%s", x[i]/(NR-1), (i<11 ? "\t" : "\n")}' tmp > $@
	rm tmp

# Create a version where each row's values are offset by its true value
updated.all.lnl.txt: eval.all.lnl.txt
	awk 'NR==1{print; next} \
		$$1=="-----"{print; next} \
		{t=$$1; for(i=1;i<=NF;i++) $$i=$$i-t; print}' $< > $@


# Extract timing info
tree.%.time: tree.%.beast.term tree.%.beast-beagle.term tree.%.var-time tree.%.mrbayes.term tree.%.mrbayes-beagle.term tree.%.dodonaphy-time tree.%.geophy-time tree.%.vaiphy-time
	echo -e "samp\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle\tvine\tdodonaphy\tgeophy\tvaiphy" > $@; \
	beast_time=$$(grep '^Total calculation time' tree.$*.beast.term | awk '{print $$4}'); \
	printf "$*\t%s\t" "$$beast_time" >> $@; \
	beast_beagle_time=$$(grep '^Total calculation time' tree.$*.beast-beagle.term | awk '{print $$4}'); \
	printf "%s\t" "$$beast_beagle_time" >> $@; \
	mrbayes_time=$$(grep 'Analysis used' tree.$*.mrbayes.term | awk '{printf "%s\t", $$3}'); \
	printf "%s\t" "$$mrbayes_time" >> $@; \
	mrbayes_beagle_time=$$(grep 'Analysis used' tree.$*.mrbayes-beagle.term | awk '{printf "%s\t", $$3}'); \
	printf "%s\t" "$$mrbayes_beagle_time" >> $@; \
	head -1 tree.$*.var-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@; \
	head -1 tree.$*.dodonaphy-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@; \
	head -1 tree.$*.geophy-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@; \
	head -1 tree.$*.vaiphy-time | awk '{printf "%s\n", $$1}' | sed 's/user//' >> $@; \

eval.all.time.txt: $(TIME)
	awk 'FNR==1 && NR==1 {print; next} FNR==2 {print; for(i=2;i<=NF;i++) sum[i]+=$$i; n++} END {if(n>0){printf "-----------------------------------------\nall"; for(i=2;i<=NF;i++) printf "\t%.2f", sum[i]/n; printf "\n"}}' $(TIME) > $@

# evalTrees stuff
# (1) modelFit
# extract kappa from vine log
tree.%.var.mf.txt: tree.%.var.nwk tree.%.heldout.fa tree.%.var.nwk.log
	$(VINE_BIN)/evalTrees tree.$*.var.nwk -f tree.$*.heldout.fa > $@

# use true kappa
tree.%.true.mf.txt: tree.%.true.nwk tree.%.heldout.fa
	$(VINE_BIN)/evalTrees tree.$*.true.nwk -f tree.$*.heldout.fa > $@

# use ML kappa
tree.%.nj.mf.txt: tree.%.nj.nwk tree.%.heldout.fa tree.%.ml.mod
	$(VINE_BIN)/evalTrees tree.$*.nj.nwk -f tree.$*.heldout.fa > $@

# use ML kappa
tree.%.ml.mf.txt: tree.%.ml.nwk tree.%.heldout.fa tree.%.ml.mod
	$(VINE_BIN)/evalTrees tree.$*.ml.nwk -f tree.$*.heldout.fa > $@

# use posterior mean kappa from beast log
tree.%.beast.mf.txt: tree.%.beast.nwk tree.%.heldout.fa tree.%.beast.log
	$(VINE_BIN)/evalTrees tree.$*.beast.nwk -f tree.$*.heldout.fa > $@

tree.%.beast-beagle.mf.txt: tree.%.beast-beagle.nwk tree.%.heldout.fa tree.%.beast-beagle.log
	$(VINE_BIN)/evalTrees tree.$*.beast-beagle.nwk -f tree.$*.heldout.fa > $@

tree.%.mrbayes.mf.txt: tree.%.mrbayes.nwk tree.%.heldout.fa
	$(VINE_BIN)/evalTrees tree.$*.mrbayes.nwk -f tree.$*.heldout.fa > $@

tree.%.mrbayes-beagle.mf.txt: tree.%.mrbayes-beagle.nwk tree.%.heldout.fa
	$(VINE_BIN)/evalTrees tree.$*.mrbayes-beagle.nwk -f tree.$*.heldout.fa > $@

tree.%.vaiphy.nwk: tree.%.vaiphy.log
	# Prefer bifurcated tree; fall back to non-bifurcated if needed
	( grep -F "Selected tree Newick (bifurcated):" $< || true ) \
		| tail -1 \
		| sed 's/.*Selected tree Newick (bifurcated):[[:space:]]*//' \
		> $@ ; \
	if [ ! -s "$@" ]; then \
		grep -F "Selected tree Newick:" $< \
			| tail -1 \
			| sed 's/.*Selected tree Newick:[[:space:]]*//' > $@ ; \
	fi
		
tree.%.vaiphy.renamed.nwk: tree.%.vaiphy.log tree.%.vaiphy.nwk
	python3 \
	  /local/storage/no-backup/vine-benchmarks/bin/rename_vaiphy_tree.py \
	  --log tree.$*.vaiphy.log \
	  --in tree.$*.vaiphy.nwk \
	  --out $@

tree.%.mf: tree.%.true.mf.txt tree.%.nj.mf.txt tree.%.ml.mf.txt tree.%.var.mf.txt tree.%.beast.mf.txt tree.%.beast-beagle.mf.txt tree.%.mrbayes.mf.txt tree.%.mrbayes-beagle.mf.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
	done

# (2) distances
tree.%.true.dist.txt: tree.%.true.nwk 
	$(VINE_BIN)/evalTrees tree.$*.true.nwk > $@

tree.%.var.dist.txt: tree.%.var.nwk tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.var.nwk > tmp.$*.var.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.var.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f tmp.$*.var.dist

tree.%.nj.dist.txt: tree.%.nj.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.nj.nwk > tmp.$*.nj.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.nj.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f tmp.$*.nj.dist

tree.%.ml.dist.txt: tree.%.ml.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.ml.nwk > tmp.$*.ml.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.ml.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f tmp.$*.ml.dist

tree.%.beast.dist.txt: tree.%.beast.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.beast.nwk > tmp.$*.beast.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.beast.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f tmp.$*.beast.dist

tree.%.beast-beagle.dist.txt: tree.%.beast-beagle.nwk tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.beast-beagle.nwk > tmp.$*.beast-beagle.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.beast-beagle.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f tmp.$*.beast-beagle.dist

tree.%.dist: tree.%.ml.dist.txt tree.%.var.dist.txt tree.%.beast.dist.txt tree.%.beast-beagle.dist.txt
	paste $^ > $@

eval.all.dist.txt: $(EVALDIST)
	cat $^ | awk '{print $$1, $$2, $$6, $$7, $$8, $$9, $$11, $$12, $$13, $$14, $$16, $$17, $$18, $$19}' | awk 'BEGIN {printf "ML_r2 ML_RMSE vine_r2 vine_RMSE vine_95CI vine_50CI beast_r2 beast_RMSE beast_95CI beast_50CI beast-beagle_r2 beast-beagle_RMSE beast-beagle_95CI beast-beagle_50CI\n"} {for(i=1;i<=NF;i++) x[i]+=$$i; print $$0} END{printf "-----\n"; for(i=1;i<=14;i++) printf "%f%s", x[i]/NR, (i<14 ? "\t" : "\n")}' > $@

# (3) RF dist
tree.%.var.rf.txt: tree.%.var.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.var.nwk -t tree.$*.true.nwk > $@

tree.%.true.rf.txt: tree.%.true.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.true.nwk -t tree.$*.true.nwk > $@

tree.%.nj.rf.txt: tree.%.nj.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.nj.nwk -t tree.$*.true.nwk > $@

tree.%.ml.rf.txt: tree.%.ml.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.ml.nwk -t tree.$*.true.nwk > $@

tree.%.beast.rf.txt: tree.%.beast.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.beast.nwk -t tree.$*.true.nwk > $@

tree.%.beast-beagle.rf.txt: tree.%.beast-beagle.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.beast-beagle.nwk -t tree.$*.true.nwk > $@

tree.%.mrbayes.rf.txt: tree.%.mrbayes.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.mrbayes.nwk -t tree.$*.true.nwk > $@

tree.%.mrbayes-beagle.rf.txt: tree.%.mrbayes-beagle.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.mrbayes-beagle.nwk -t tree.$*.true.nwk > $@

tree.%.rf: tree.%.true.rf.txt tree.%.nj.rf.txt tree.%.ml.rf.txt tree.%.var.rf.txt tree.%.beast.rf.txt tree.%.beast-beagle.rf.txt tree.%.mrbayes.rf.txt tree.%.mrbayes-beagle.rf.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
	done

eval.all.rf.txt: $(EVALRF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > tmp
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{for(i=1;i<=NF;i+=2){x[i]+=$$i; xs[i]+=($$(i+1)*$$(i+1))} print $$0} END {printf "-----\n"; for(i=1;i<=15;i+=2) printf "%f\t%f%s", x[i]/(NR-1), sqrt(xs[i]/(NR-1)), (i<15 ? "\t" : "\n")}' tmp > $@
	rm -f tmp

eval.all.mf.txt: $(EVALMF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > tmp
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{for(i=1;i<=NF;i+=2){if($$i!="nan") x[i]+=$$i; xs[i]+=($$(i+1)*$$(i+1))} print $$0} END {printf "-----\n"; for(i=1;i<=15;i+=2) printf "%f\t%f%s", x[i]/(NR-1), sqrt(xs[i]/(NR-1)), (i<15 ? "\t" : "\n")}' tmp > $@
	rm -f tmp

# for use in debugging
tracer: $(TRACER)

tree.%.tr: tree.%.var.nwk.log
	grep -v '^#' $^ > $@

raxml_eliminate:
	rm -rf tree.*.raxml*
	rm -rf tree.*.lnl
	rm -rf tree.*.time
clean:
	rm -rf $(TREES) $(FA) $(MOD) $(ML) $(MLMOD) $(NJMOD) $(NJ) $(VAR) $(VARNEX) $(EVALURF) $(LNL) $(VARLOG) $(VARTIME) tree.*.mean*.nwk tree.*.lnl.diffs tree.*.varlnl tree.*.modlnl eval.all.*.txt tree.*.beast* *.mf.txt *.rf.txt *.dist.txt  tree.*.time tree.*.lnl tree.*.mf tree.*.rf $(FAHELDOUT) $(TRACER)
	rm -rf tree.*.mrbayes* tree.*.nex
	rm -rf tree.*.raxml*
	rm -rf tree.*.dodonaphy*
	rm -rf tree.*.geophy*


clean_vine:
	rm -rf tree.*.var* eval.all.*.txt tree.*.mf tree.*.rf tree.*.dist tree.*.lnl updated.all.lnl.txt

clean-beast-beagle:
	rm -rf tree.*.beast-beagle* thinned.trees.*.beast-beagle tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt

clean-mrbayes-beagle:
	rm -rf tree.*.mrbayes-beagle* tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt


# Some useful rules for archiving results more efficiently
archive_mcmc:
	archive_dir="archive.beast_mrbayes_$$(date +%Y-%m-%d_%H:%M:%S)"; \
	mkdir $$archive_dir; \
	mv tree.*.beast* $$archive_dir/; \
	mv tree.*.mrbayes* $$archive_dir/; \
	mv eval.all.*.txt $$archive_dir/; \
	mv tree.*.lnl $$archive_dir/; \
	mv tree.*.time $$archive_dir/; \
	mv tree.*.rf $$archive_dir/; \
	mv tree.*.mf $$archive_dir/; \
	mv tree.*.dist $$archive_dir/; \
	mv beast_ess_runtime_scale_factor* $$archive_dir/; \
	mv mrbayes_ess_runtime_scale_factor* $$archive_dir/

archive_vine:
	archive_dir="archive.vine_$$(date +%Y-%m-%d_%H:%M:%S)"; \
	mkdir $$archive_dir; \
	mv tree.*.var* $$archive_dir/; \
	mv eval.all.*.txt $$archive_dir/; \
	mv tree.*.lnl $$archive_dir/; \
	mv tree.*.time $$archive_dir/; \
	mv tree.*.rf $$archive_dir/; \
	mv tree.*.mf $$archive_dir/; \
	mv tree.*.dist $$archive_dir/
