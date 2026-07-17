export SHELL=/usr/bin/bash

# Prevent Make from removing files it thinks it should clean up
#.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# edit for local structure; this is the only place absolute paths are used
MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT_SUFFIX := dna_sims/hky_10000sites

ROOT := $(MAIN_DIR)/$(ROOT_SUFFIX)
PYTHON_SRC := $(MAIN_DIR)/python/src
BIN := $(MAIN_DIR)/bin
PHAST_BIN := $(BIN)/phast/bin
VINE_BIN := $(BIN)/vine/bin
BEAST_BIN := $(BIN)/beast/bin
BEAST := $(BEAST_BIN)/beast
MRBAYES := $(BIN)/mb
BEAST_TEMPLATE := $(ROOT)/beast_template.xml
BURNIN_PCT := 10
METRIC_TREE_COUNT := 1000

TREES := $(shell seq -f tree.%.0f.true.nwk 1 $(NSAMP))
FA := $(patsubst %.true.nwk,%.fa,$(TREES))
MOD := $(patsubst %.true.nwk,%.true.mod,$(TREES))
MLMOD := $(patsubst %.true.nwk,%.ml.mod,$(TREES))
NJMOD := $(patsubst %.true.nwk,%.nj.mod,$(TREES))
NJ := $(patsubst %.true.nwk,%.nj.nwk,$(TREES))
ML := $(patsubst %.true.nwk,%.ml.nwk,$(TREES))
VAR := $(patsubst %.true.nwk,%.var.nwk,$(TREES))
EVALRF := $(patsubst tree.%.true.nwk,tree.%.rf,$(TREES))
EVALBSD := $(patsubst tree.%.true.nwk,tree.%.bsd,$(TREES))
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

.PHONY: beast-beagle beast-beagle-nwk mrbayes-beagle mrbayes-beagle-nwk bsd clean-mf

all: eval.all.lnl.txt eval.all.rf.txt eval.all.mf.txt eval.all.time.txt eval.all.dist.txt

bsd: eval.all.bsd.txt
vine: $(VAR)
infer: $(BEASTLOG) $(BEASTBEAGLELOG) $(MRBAYESLOG) $(MRBAYESBEAGLELOG) $(VAR)
beast: $(BEASTLOG)
beast-beagle: $(BEASTBEAGLELOG)
mrbayes: $(MRBAYESLOG)
mrbayes-beagle: $(MRBAYESBEAGLELOG)


tree.%.true.nwk: 
	$(BIN)/bdTree3 -b 1 -d 0.5 --oversample-k 3 --height 5 --min-edge 0.02 --expected-height $(EXPHEIGHT) --no-stem --ucln-sd 0.6 --target-stat median -n $(NTAXA) | sed 's/\[\&[UR]\] //' > $@

tree.%.fa: tree.%.true.nwk
	cp ../base-hky.mod tree.$*.fa.tmp.mod
	echo -n "TREE: " >> tree.$*.fa.tmp.mod
	cat $< >> tree.$*.fa.tmp.mod
	$(PHAST_BIN)/base_evolve --nsites $(NSITES) tree.$*.fa.tmp.mod > $@
	rm -f tree.$*.fa.tmp.mod

tree.%.heldout.fa: tree.%.true.nwk
	cp ../base-hky.mod tree.$*.heldout.fa.tmp.mod
	echo -n "TREE: " >> tree.$*.heldout.fa.tmp.mod
	cat $< >> tree.$*.heldout.fa.tmp.mod
	$(PHAST_BIN)/base_evolve --nsites $(NSITES) tree.$*.heldout.fa.tmp.mod > $@
	rm -f tree.$*.heldout.fa.tmp.mod

tree.%.nj.nwk: tree.%.fa
	$(VINE_BIN)/vine --nj-only $< > $@

tree.%.ml.mod: tree.%.nj.nwk tree.%.fa
	$(PHAST_BIN)/phyloFit --subst-mod HKY85 --tree $^ -o tree.$*.ml

tree.%.ml.nwk: tree.%.ml.mod
	$(PHAST_BIN)/tree_doctor --tree-only $^ > $@

tree.%.var.nwk tree.%.var-time tree.%.var.nwk.log: tree.%.fa 
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine $< -l tree.$*.var.nwk.log $(VAROPT) --mean tree.$*.mean.nwk > tree.$*.var.nwk

tree.%.beast.xml:
	cp "$(BEAST_TEMPLATE)" $@

tree.%.beast.term tree.%.beast-tree.trees tree.%.beast.log: tree.%.beast.xml tree.%.fa
	rm -f tree.$*.beast-tree.trees tree.$*.beast.log
	"$(BEAST)" -java -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$(BEAST_MCMC_SAMPLEFREQ) -D printfreq=$(MCMC_PRINTFREQ) $< > tree.$*.beast.term

tree.%.beast-beagle.xml:
	cp "$(BEAST_TEMPLATE)" $@

tree.%.beast-beagle.term tree.%.beast-beagle-tree.trees tree.%.beast-beagle.log: tree.%.beast-beagle.xml tree.%.fa
	rm -f tree.$*.beast-beagle-tree.trees tree.$*.beast-beagle.log
	"$(BEAST)" -beagle -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$(BEAST_MCMC_SAMPLEFREQ) -D printfreq=$(MCMC_PRINTFREQ) $< > tree.$*.beast-beagle.term

tree.%.beast.nwk: tree.%.beast-tree.trees
	THIN=$$(( ($(BEAST_MCMCLEN) / $(BEAST_MCMC_SAMPLEFREQ)) * (100 - $(BURNIN_PCT)) / 100 / $(METRIC_TREE_COUNT) )); \
	if [ "$$THIN" -lt 1 ]; then THIN=1; fi; \
	RESAMPLE=$$(( $(BEAST_MCMC_SAMPLEFREQ) * THIN )); \
	$(BEAST_BIN)/logcombiner -log $^ -o thinned.trees.$*.beast \
	  -burnin $(BURNIN_PCT) -resample $$RESAMPLE
	python3 "$(PYTHON_SRC)/time2subs.py" thinned.trees.$*.beast tree.$*.beast.tmp.nex
	$(BIN)/convertTrees -i nexus tree.$*.beast.tmp.nex > tree.$*.beast.tmp.nwk
	N=$$(awk 'NF{n++} END{print n}' tree.$*.beast.tmp.nwk); \
	TARGET=$(METRIC_TREE_COUNT); \
	if [ "$$N" -lt "$$TARGET" ]; then TARGET=$$N; fi; \
	awk -v n="$$N" -v target="$$TARGET" \
	  'NF{j++; if (int(j*target/n) > int((j-1)*target/n)) print}' \
	  tree.$*.beast.tmp.nwk > $@
	rm -f tree.$*.beast.tmp.nex tree.$*.beast.tmp.nwk thinned.trees.$*.beast

tree.%.beast-beagle.nwk: tree.%.beast-beagle-tree.trees
	THIN=$$(( ($(BEAST_MCMCLEN) / $(BEAST_MCMC_SAMPLEFREQ)) * (100 - $(BURNIN_PCT)) / 100 / $(METRIC_TREE_COUNT) )); \
	if [ "$$THIN" -lt 1 ]; then THIN=1; fi; \
	RESAMPLE=$$(( $(BEAST_MCMC_SAMPLEFREQ) * THIN )); \
	$(BEAST_BIN)/logcombiner -log $^ -o thinned.trees.$*.beast-beagle \
	  -burnin $(BURNIN_PCT) -resample $$RESAMPLE
	python3 "$(PYTHON_SRC)/time2subs.py" thinned.trees.$*.beast-beagle tree.$*.beast-beagle.tmp.nex
	$(BIN)/convertTrees -i nexus tree.$*.beast-beagle.tmp.nex > tree.$*.beast-beagle.tmp.nwk
	N=$$(awk 'NF{n++} END{print n}' tree.$*.beast-beagle.tmp.nwk); \
	TARGET=$(METRIC_TREE_COUNT); \
	if [ "$$N" -lt "$$TARGET" ]; then TARGET=$$N; fi; \
	awk -v n="$$N" -v target="$$TARGET" \
	  'NF{j++; if (int(j*target/n) > int((j-1)*target/n)) print}' \
	  tree.$*.beast-beagle.tmp.nwk > $@
	rm -f tree.$*.beast-beagle.tmp.nex tree.$*.beast-beagle.tmp.nwk

tree.%.nex: tree.%.fa
	$(BIN)/fa2nex $< $@

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	$(BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex --mcmc_length $(MRBAYES_MCMCLEN) --model HKY \
		--sample_freq $(MRBAYES_MCMC_SAMPLEFREQ) --print_freq $(MCMC_PRINTFREQ) --diagn_freq $(MCMC_PRINTFREQ)

tree.%.mrbayes-beagle.nex: tree.%.nex
	$(BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes-beagle.nex --mcmc_length $(MRBAYES_MCMCLEN) --model HKY \
		--sample_freq $(MRBAYES_MCMC_SAMPLEFREQ) --print_freq $(MCMC_PRINTFREQ) --diagn_freq $(MCMC_PRINTFREQ) --use_beagle

# Run MrBayes
tree.%.mrbayes.term tree.%.mrbayes.nex.p tree.%.mrbayes.nex.t: tree.%.mrbayes.nex
	$(MRBAYES) tree.$*.mrbayes.nex > tree.$*.mrbayes.term

tree.%.mrbayes-beagle.term tree.%.mrbayes-beagle.nex.p tree.%.mrbayes-beagle.nex.t: tree.%.mrbayes-beagle.nex
	$(MRBAYES) tree.$*.mrbayes-beagle.nex > tree.$*.mrbayes-beagle.term

# Get mrbayes tree in nexus format
tree.%.mrbayes.nwk: tree.%.mrbayes.nex.t
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
	  $< > tree.$*.mrbayes.thinned.t
	$(BIN)/convertTrees -i nexus tree.$*.mrbayes.thinned.t \
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
	$(BIN)/convertTrees -i nexus tree.$*.mrbayes-beagle.thinned.t \
	  | sed 's/^\[&[^]]*\]\s*//' > $@
	rm -f tree.$*.mrbayes-beagle.thinned.t

# Run raxml
tree.%.raxml.term: tree.%.fa
	rm -f $@
	sed 's/> />/g' $< > tree.$*.raxml.fa
	$(BIN)/raxml-ng --msa tree.$*.raxml.fa --model HKY+F --prefix tree.$* --search1 --threads 1 > tree.$*.raxml.term
	rm -f tree.$*.raxml.fa

# extract training likelihoods
tree.%.true.mod: tree.%.true.nwk tree.%.fa
	cp ../base-hky.mod tree.$*.true.tmp.mod
	echo -n "TREE: " >> tree.$*.true.tmp.mod
	cat $< >> tree.$*.true.tmp.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tree.$*.true.tmp.mod -o tree.$*.true tree.$*.fa
	rm -f tree.$*.true.tmp.mod

tree.%.nj.mod: tree.%.nj.nwk tree.%.fa
	cp ../base-hky.mod tree.$*.nj.tmp.mod
	echo -n "TREE: " >> tree.$*.nj.tmp.mod
	cat $< >> tree.$*.nj.tmp.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tree.$*.nj.tmp.mod -o tree.$*.nj tree.$*.fa
	rm -f tree.$*.nj.tmp.mod

tree.%.modlnl: tree.%.true.mod tree.%.nj.mod tree.%.ml.mod 
	rm -f $@
	for file in $^ ; do \
		echo -n "$${file} " >> $@ ;\
		grep LNL $${file} | awk '{print $$2}' >> $@ ;\
	done

tree.%.varlnl: tree.%.var.nwk.log 
	echo -n "$^ " > $@
	tail -1 $^ | awk '{print $$11}' | sed 's/,//g' >> $@ 

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

tree.%.raxmllnl: tree.%.raxml.term
	echo -n "$^ " > $@
	grep 'Final LogLikelihood:' $^ | awk '{printf "%.6f\n", $$3}' >> $@

tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.beast-beaglelnl tree.%.mrbayeslnl tree.%.mrbayes-beaglelnl tree.%.raxmllnl
	cat $^ | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@

eval.all.lnl.txt: $(LNL)
	echo -e "true\tnj\tml\tvine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle\traxml" > $@.tmp
	for file in $^ ; do \
		awk '{printf "%s\t", $$2}' $${file} >> $@.tmp ;\
		echo >> $@.tmp ;\
	done
	awk '{for(i=1;i<=NF;i++) x[i]+=$$i; print $$0} END {printf "-----\n"; for(i=1;i<=9;i++) printf "%f%s", x[i]/(NR-1), (i<9 ? "\t" : "\n")}' $@.tmp > $@
	rm -f $@.tmp

# Create a version where each row's values are offset by its true value
updated.all.lnl.txt: eval.all.lnl.txt
	awk 'NR==1{print; next} \
		$$1=="-----"{print; next} \
		{t=$$1; for(i=1;i<=NF;i++) $$i=$$i-t; print}' $< > $@


## extract timing info
tree.%.time: tree.%.beast.term tree.%.beast-beagle.term tree.%.var-time tree.%.mrbayes.term tree.%.mrbayes-beagle.term tree.%.raxml.term
	echo -e "samp\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle\tvine\traxml" > $@; \
	beast_time=$$(grep '^Total calculation time' tree.$*.beast.term | awk '{print $$4}'); \
	printf "$*\t%s\t" "$$beast_time" >> $@; \
	beast_beagle_time=$$(grep '^Total calculation time' tree.$*.beast-beagle.term | awk '{print $$4}'); \
	printf "%s\t" "$$beast_beagle_time" >> $@; \
	mrbayes_time=$$(grep 'Analysis used' tree.$*.mrbayes.term | awk '{printf "%s\t", $$3}'); \
	printf "%s\t" "$$mrbayes_time" >> $@; \
	mrbayes_beagle_time=$$(grep 'Analysis used' tree.$*.mrbayes-beagle.term | awk '{printf "%s\t", $$3}'); \
	printf "%s\t" "$$mrbayes_beagle_time" >> $@; \
	head -1 tree.$*.var-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@; \
	grep 'Elapsed time:' tree.$*.raxml.term | awk '{printf "%s\n", $$3}' >> $@


eval.all.time.txt: $(TIME)
	awk 'FNR==1 && NR==1 {print; next} FNR==2 {print; for(i=2;i<=NF;i++) sum[i]+=$$i; n++} END {if(n>0){printf "-----------------------------------------\nall"; for(i=2;i<=NF;i++) printf "\t%.2f", sum[i]/n; printf "\n"}}' $(TIME) > $@

# evalTrees stuff
# (1) modelFit
# extract kappa from vine log
tree.%.var.mf.txt: tree.%.var.nwk tree.%.heldout.fa tree.%.var.nwk.log
	kappa=$$(awk '{for (i=1; i<=NF; i++) if ($$i == "kappa:") kappa=$$(i+1)} END {print kappa}' tree.$*.var.nwk.log); \
	$(VINE_BIN)/evalTrees tree.$*.var.nwk -f tree.$*.heldout.fa \
	  -k "$$kappa" > $@.tmp && mv $@.tmp $@

# use true kappa
tree.%.true.mf.txt: tree.%.true.nwk tree.%.heldout.fa
	$(VINE_BIN)/evalTrees tree.$*.true.nwk -f tree.$*.heldout.fa -k 4 > $@

# use ML kappa
tree.%.nj.mf.txt: tree.%.nj.nwk tree.%.heldout.fa tree.%.ml.mod
	kappa=`awk '$$1 == "BACKGROUND:" {pi_c = $$3; pi_g = $$4 } $$1<0 {print ($$3/pi_g) / ($$2/pi_c)}' tree.$*.ml.mod` ;\
	$(VINE_BIN)/evalTrees tree.$*.nj.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# use ML kappa
tree.%.ml.mf.txt: tree.%.ml.nwk tree.%.heldout.fa tree.%.ml.mod
	kappa=`awk '$$1 == "BACKGROUND:" {pi_c = $$3; pi_g = $$4 } $$1<0 {print ($$3/pi_g) / ($$2/pi_c)}' tree.$*.ml.mod` ;\
	$(VINE_BIN)/evalTrees tree.$*.ml.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# use posterior mean kappa from beast log
tree.%.beast.mf.txt: tree.%.beast.nwk tree.%.heldout.fa tree.%.beast.log
	kappa=`awk '{if (inlog) {sum += $$10; n++} ; if ($$1 == "Sample") inlog=1} END {print sum/n}' tree.$*.beast.log` ;\
	$(VINE_BIN)/evalTrees tree.$*.beast.nwk -f tree.$*.heldout.fa -k $$kappa > $@

tree.%.beast-beagle.mf.txt: tree.%.beast-beagle.nwk tree.%.heldout.fa tree.%.beast-beagle.log
	kappa=`awk '{if (inlog) {sum += $$10; n++} ; if ($$1 == "Sample") inlog=1} END {print sum/n}' tree.$*.beast-beagle.log` ;\
	$(VINE_BIN)/evalTrees tree.$*.beast-beagle.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# Use the post-burn-in posterior mean kappa for both MrBayes modes.
tree.%.mrbayes.mf.txt: tree.%.mrbayes.nwk tree.%.heldout.fa \
                       tree.%.mrbayes.nex.p
	kappa=`awk -v p=$(BURNIN_PCT) '\
	  $$1=="Gen"{\
	    for(i=1;i<=NF;i++) if($$i=="kappa") c=i; next\
	  }\
	  $$1!="Gen" && c { vals[++n]=$$c }\
	  END { \
	    if(n){ \
	      skip=int(n*p/100); \
	      if(skip>=n) skip=n-1; \
	      sum=0; cnt=0; \
	      for(i=skip+1;i<=n;i++){ sum+=vals[i]; cnt++ } \
	      if(cnt) printf "%.6f\\n", sum/cnt; \
	    } \
	  }\
	' tree.$*.mrbayes.nex.p` ;\
	$(VINE_BIN)/evalTrees tree.$*.mrbayes.nwk -f tree.$*.heldout.fa \
	  -k $$kappa > $@

tree.%.mrbayes-beagle.mf.txt: tree.%.mrbayes-beagle.nwk tree.%.heldout.fa \
                              tree.%.mrbayes-beagle.nex.p
	kappa=`awk -v p=$(BURNIN_PCT) '\
	  $$1=="Gen"{\
	    for(i=1;i<=NF;i++) if($$i=="kappa") c=i; next\
	  }\
	  $$1!="Gen" && c { vals[++n]=$$c }\
	  END { \
	    if(n){ \
	      skip=int(n*p/100); \
	      if(skip>=n) skip=n-1; \
	      sum=0; cnt=0; \
	      for(i=skip+1;i<=n;i++){ sum+=vals[i]; cnt++ } \
	      if(cnt) printf "%.6f\\n", sum/cnt; \
	    } \
	  }\
	' tree.$*.mrbayes-beagle.nex.p` ;\
	$(VINE_BIN)/evalTrees tree.$*.mrbayes-beagle.nwk -f tree.$*.heldout.fa \
	  -k $$kappa > $@

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
	$(VINE_BIN)/evalTrees tree.$*.var.nwk > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f $@.tmp

tree.%.nj.dist.txt: tree.%.nj.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.nj.nwk > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f $@.tmp

tree.%.ml.dist.txt: tree.%.ml.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.ml.nwk > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f $@.tmp

tree.%.beast.dist.txt: tree.%.beast.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.beast.nwk > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f $@.tmp

tree.%.beast-beagle.dist.txt: tree.%.beast-beagle.nwk tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.beast-beagle.nwk > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f $@.tmp

tree.%.mrbayes.dist.txt: tree.%.mrbayes.nwk tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.mrbayes.nwk > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f $@.tmp

tree.%.mrbayes-beagle.dist.txt: tree.%.mrbayes-beagle.nwk tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.mrbayes-beagle.nwk > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f $@.tmp

tree.%.dist: tree.%.ml.dist.txt tree.%.var.dist.txt tree.%.beast.dist.txt tree.%.beast-beagle.dist.txt tree.%.mrbayes.dist.txt tree.%.mrbayes-beagle.dist.txt
	paste $^ > $@

eval.all.dist.txt: $(EVALDIST)
	cat $^ | awk '{print $$1, $$2, $$6, $$7, $$8, $$9, $$11, $$12, $$13, $$14, $$16, $$17, $$18, $$19, $$21, $$22, $$23, $$24, $$26, $$27, $$28, $$29}' | awk 'BEGIN {printf "ML_r2 ML_RMSE vine_r2 vine_RMSE vine_95CI vine_50CI beast_r2 beast_RMSE beast_95CI beast_50CI beast-beagle_r2 beast-beagle_RMSE beast-beagle_95CI beast-beagle_50CI mrbayes_r2 mrbayes_RMSE mrbayes_95CI mrbayes_50CI mrbayes-beagle_r2 mrbayes-beagle_RMSE mrbayes-beagle_95CI mrbayes-beagle_50CI\n"} {for(i=1;i<=NF;i++) x[i]+=$$i; print $$0} END{printf "-----\n"; for(i=1;i<=22;i++) printf "%f%s", x[i]/NR, (i<22 ? "\t" : "\n")}' > $@

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
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > $@.tmp
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> $@.tmp ;\
		echo >> $@.tmp ;\
	done
	awk '{for(i=1;i<=NF;i+=2){x[i]+=$$i; xs[i]+=($$(i+1)*$$(i+1))} print $$0} END {printf "-----\n"; for(i=1;i<=15;i+=2) printf "%f\t%f%s", x[i]/(NR-1), sqrt(xs[i]/(NR-1)), (i<15 ? "\t" : "\n")}' $@.tmp > $@
	rm -f $@.tmp

eval.all.mf.txt: $(EVALMF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > $@.tmp
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> $@.tmp ;\
		echo >> $@.tmp ;\
	done
	awk '{for(i=1;i<=NF;i+=2){if($$i!="nan") x[i]+=$$i; xs[i]+=($$(i+1)*$$(i+1))} print $$0} END {printf "-----\n"; for(i=1;i<=15;i+=2) printf "%f\t%f%s", x[i]/(NR-1), sqrt(xs[i]/(NR-1)), (i<15 ? "\t" : "\n")}' $@.tmp > $@
	rm -f $@.tmp

# Expected normalized Kuhner-Felsenstein branch-score distance.  evalTrees
# reports the distance from every posterior draw to the true tree; extract its
# Mean rather than the separate posterior-mean-tree point estimate.
tree.%.true.bsd.txt: tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $< > $@

tree.%.nj.bsd.txt: tree.%.nj.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

tree.%.ml.bsd.txt: tree.%.ml.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

tree.%.var.bsd.txt: tree.%.var.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

tree.%.beast.bsd.txt: tree.%.beast.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

tree.%.beast-beagle.bsd.txt: tree.%.beast-beagle.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

tree.%.mrbayes.bsd.txt: tree.%.mrbayes.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

tree.%.mrbayes-beagle.bsd.txt: tree.%.mrbayes-beagle.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

tree.%.bsd: tree.%.true.bsd.txt tree.%.nj.bsd.txt tree.%.ml.bsd.txt \
    tree.%.var.bsd.txt tree.%.beast.bsd.txt tree.%.beast-beagle.bsd.txt \
    tree.%.mrbayes.bsd.txt tree.%.mrbayes-beagle.bsd.txt
	rm -f $@
	for file in $^; do \
	  echo -n "$$file     " >> $@; \
	  awk '/^Mean:/ {mean=$$2} /Reference tree length:/ {rl=$$NF} \
	       END {printf "%f\t%f\n", mean, rl}' "$$file" >> $@; \
	done

eval.all.bsd.txt: $(EVALBSD)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > $@.tmp
	for file in $^; do \
	  awk '{printf "%f\t", ($$3>0 ? $$2/$$3 : 0)}' "$$file" >> $@.tmp; \
	  echo >> $@.tmp; \
	done
	awk '{for(i=1;i<=8;i++){x[i]+=$$i; xs[i]+=$$i*$$i} print} \
	     END {print "-----"; for(i=1;i<=8;i++){m=x[i]/(NR-1); v=xs[i]/(NR-1)-m*m; \
	       if(v<0)v=0; printf "%f\t%f%s",m,sqrt(v),(i<8?"\t":"\n")}}' $@.tmp > $@
	rm -f $@.tmp

# for use in debugging
tracer: $(TRACER)

tree.%.tr: tree.%.var.nwk.log
	grep -v '^#' $^ > $@

clean:
	rm -rf $(TREES) $(FA) $(MOD) $(ML) $(MLMOD) $(NJMOD) $(NJ) $(VAR) $(VARNEX) $(EVALURF) $(LNL) $(VARLOG) $(VARTIME) tree.*.mean*.nwk tree.*.lnl.diffs tree.*.varlnl tree.*.modlnl eval.all.*.txt tree.*.beast* *.mf.txt *.rf.txt *.dist.txt  tree.*.time tree.*.lnl tree.*.mf tree.*.rf $(FAHELDOUT) $(TRACER)
	rm -rf tree.*.mrbayes* tree.*.nex
	rm -rf tree.*.raxml*

clean_mcmc:
	rm -rf tree.*.beast* tree.*.mrbayes* beast_ess* mrbayes_ess* eval.all.*

# Remove only derived MrBayes model-fit outputs. Preserve the expensive
# .nex.p and .nex.t inference results so `make eval.all.mf.txt` can reuse them.
clean-mf:
	rm -f tree.*.mrbayes.nwk tree.*.mrbayes-beagle.nwk
	rm -f tree.*.mrbayes.mf.txt tree.*.mrbayes-beagle.mf.txt
	rm -f tree.*.mf eval.all.mf.txt

clean-beast-beagle:
	rm -rf tree.*.beast-beagle* thinned.trees.*.beast-beagle tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt

clean-mrbayes-beagle:
	rm -rf tree.*.mrbayes-beagle* tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt

clean_vine:
	rm -rf tree.*.var* eval.all.*.txt tree.*.mf tree.*.rf tree.*.dist tree.*.lnl

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
