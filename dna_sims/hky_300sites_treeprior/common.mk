export SHELL=/usr/bin/bash

# Prevent Make from removing files it thinks it should clean up
#.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# edit for local structure; this is the only place absolute paths are used
MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT_SUFFIX := dna_sims/hky_300sites

ROOT := $(MAIN_DIR)/$(ROOT_SUFFIX)
BIN := $(MAIN_DIR)/bin
PYTHON_SRC := $(MAIN_DIR)/python/src
PHAST_BIN := $(BIN)/phast/bin
VINE_BIN := $(BIN)/vine/bin
BEAST_BIN := $(MAIN_DIR)/bin/beast/bin
BEAST := $(BEAST_BIN)/beast
MRBAYES := $(BIN)/mb
BEAST_TEMPLATE := $(ROOT)/beast_template.xml
CONTAINERS := $(MAIN_DIR)/containers
BURNIN_PCT := 10
METRIC_TREE_COUNT := 1000

# Log approximately 10,000 MCMC samples for each program.
BEAST_SAMPLEFREQ := $(shell awk 'BEGIN {f=int($(BEAST_MCMCLEN)/10000); print (f > 0 ? f : 1)}')
MRBAYES_SAMPLEFREQ := $(shell awk 'BEGIN {f=int($(MRBAYES_MCMCLEN)/10000); print (f > 0 ? f : 1)}')

TREES := $(shell seq -f tree.%.0f.true.nwk 1 $(NSAMP))
FA := $(patsubst %.true.nwk,%.fa,$(TREES))
MOD := $(patsubst %.true.nwk,%.true.mod,$(TREES))
MLMOD := $(patsubst %.true.nwk,%.ml.mod,$(TREES))
NJMOD := $(patsubst %.true.nwk,%.nj.mod,$(TREES))
NJ := $(patsubst %.true.nwk,%.nj.nwk,$(TREES))
ML := $(patsubst %.true.nwk,%.ml.nwk,$(TREES))
RAXML := $(patsubst %.true.nwk,%.raxml.term,$(TREES))
VAR := $(patsubst %.true.nwk,%.var.nwk,$(TREES))
EVALRF := $(patsubst tree.%.true.nwk,tree.%.rf,$(TREES))
EVALBSD := $(patsubst tree.%.true.nwk,tree.%.bsd,$(TREES))
EVALMF := $(patsubst tree.%.true.nwk,tree.%.mf,$(TREES))
EVALDIST := $(patsubst tree.%.true.nwk,tree.%.dist,$(TREES))
EVALENT := $(patsubst tree.%.true.nwk,tree.%.ent,$(TREES))
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

.PHONY: beast-beagle beast-beagle-nwk mrbayes-beagle mrbayes-beagle-nwk clean_rf clean_bsd

all: eval.all.lnl.txt eval.all.rf.txt eval.all.mf.txt eval.all.time.txt eval.all.dist.txt eval.all.ent.txt

bsd: eval.all.bsd.txt

simulate: $(FA)
max_lik: $(ML) $(RAXML)
mcmc: $(BEASTLOG) $(BEASTBEAGLELOG) $(MRBAYESLOG) $(MRBAYESBEAGLELOG)
beast: $(BEASTLOG)
beast-beagle: $(BEASTBEAGLELOG)
mrbayes: $(MRBAYESLOG)
mrbayes-beagle: $(MRBAYESBEAGLELOG)
ml: $(ML)
beastnwk: $(BEASTNWK)
beast-beagle-nwk: $(BEASTBEAGLENWK)
mrbayes-beagle-nwk: $(MRBAYESBEAGLENWK)
vine: $(VAR) $(VARLOG)

# clean up vine files for rerun
vineclean:
	rm -f $(VAR) $(VARLOG) $(LNL) $(TIME) eval.all.lnl.txt eval.all.time.txt tree.*.varlnl

vinemfclean:
	rm -f $(EVALRF) $(EVALMF) eval.all.mf.txt eval.all.rf.txt

tree.%.true.nwk: 
	$(BIN)/bdTree3 -b 1 -d 0.5 --oversample-k 3 --height 5 --min-edge 0.02 --expected-height $(EXPHEIGHT) --no-stem --ucln-sd 0.6 --target-stat median -n $(NTAXA) | sed 's/\[\&[UR]\] //' > $@

tree.%.fa: tree.%.true.nwk
	cp ../base-hky.mod tmp.$*.fa.mod
	echo -n "TREE: " >> tmp.$*.fa.mod
	cat $< >> tmp.$*.fa.mod
	$(PHAST_BIN)/base_evolve --nsites $(NSITES) tmp.$*.fa.mod > $@
	rm tmp.$*.fa.mod

tree.%.heldout.fa: tree.%.true.nwk
	cp ../base-hky.mod tmp.$*.heldout.mod
	echo -n "TREE: " >> tmp.$*.heldout.mod
	cat $< >> tmp.$*.heldout.mod
	$(PHAST_BIN)/base_evolve --nsites $(NSITES) tmp.$*.heldout.mod > $@
	rm tmp.$*.heldout.mod

tree.%.nj.nwk: tree.%.fa
	$(VINE_BIN)/vine --nj-only $< > $@

tree.%.ml.mod: tree.%.nj.nwk tree.%.fa
	$(PHAST_BIN)/phyloFit --subst-mod HKY85 --tree $^ -o tree.$*.ml

tree.%.ml.nwk: tree.%.ml.mod
	$(PHAST_BIN)/tree_doctor --tree-only $^ > $@

tree.%.var.nwk tree.%.var-time tree.%.var.nwk.log: tree.%.fa 
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine $< -l tree.$*.var.nwk.log $(VAROPT) \
	--mean tree.$*.mean.nwk > tree.$*.var.nwk

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
	  -burnin $(BURNIN_PCT) -resample $$(( $(BEAST_MCMCLEN) / 1000 ))
	python3 "$(PYTHON_SRC)/time2subs.py" thinned.trees.$*.beast tmp.$*.beast.nex
	$(BIN)/convertTrees -i nexus tmp.$*.beast.nex > $@
	rm -f tmp.$*.beast.nex

tree.%.beast-beagle.nwk: tree.%.beast-beagle-tree.trees
	THIN=$$(( ($(BEAST_MCMCLEN) / $(BEAST_SAMPLEFREQ)) * (100 - $(BURNIN_PCT)) / 100 / $(METRIC_TREE_COUNT) )); \
	if [ "$$THIN" -lt 1 ]; then THIN=1; fi; \
	RESAMPLE=$$(( $(BEAST_SAMPLEFREQ) * THIN )); \
	$(BEAST_BIN)/logcombiner -log $^ -o thinned.trees.$*.beast-beagle \
	  -burnin $(BURNIN_PCT) -resample $$RESAMPLE
	python3 "$(PYTHON_SRC)/time2subs.py" thinned.trees.$*.beast-beagle tmp.$*.beast-beagle.nex
	$(BIN)/convertTrees -i nexus tmp.$*.beast-beagle.nex > tmp.$*.beast-beagle.nwk
	N=$$(awk 'NF{n++} END{print n}' tmp.$*.beast-beagle.nwk); \
	TARGET=$(METRIC_TREE_COUNT); \
	if [ "$$N" -lt "$$TARGET" ]; then TARGET=$$N; fi; \
	awk -v n="$$N" -v target="$$TARGET" \
	  'NF{j++; if (int(j*target/n) > int((j-1)*target/n)) print}' \
	  tmp.$*.beast-beagle.nwk > $@
	rm -f tmp.$*.beast-beagle.nex tmp.$*.beast-beagle.nwk

tree.%.nex: tree.%.fa
	$(BIN)/fa2nex $< $@

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	$(BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex --mcmc_length $(MRBAYES_MCMCLEN) --model HKY \
		--sample_freq $(MRBAYES_SAMPLEFREQ) --print_freq $(MCMC_PRINTFREQ) --diagn_freq $(MCMC_PRINTFREQ)

tree.%.mrbayes-beagle.nex: tree.%.nex
	$(BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes-beagle.nex --mcmc_length $(MRBAYES_MCMCLEN) --model HKY \
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
	THIN=$$(( $(MRBAYES_MCMCLEN) / $(MRBAYES_SAMPLEFREQ) / 1000 )); \
	awk -v skip="$$SKIP" -v thin="$$THIN" '\
	  /^[[:space:]]*tree[[:space:]]+/{ \
	    if (++c <= skip) next; \
	    if ((c - skip) % thin) next \
	  } \
	  1' \
	  $< > tree.$*.mrbayes.thinned.t
	$(BIN)/convertTrees -i nexus \
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
	$(BIN)/convertTrees -i nexus \
	  tree.$*.mrbayes-beagle.thinned.t \
	  | sed 's/^\[&[^]]*\]\s*//' > $@
	rm -f tree.$*.mrbayes-beagle.thinned.t

# Run raxml
tree.%.raxml.term: tree.%.fa
	rm -f $@
	sed 's/> />/g' $< > tree.$*.raxml.fa
	$(BIN)/raxml-ng --msa tree.$*.raxml.fa --model HKY+F --prefix tree.$* --search1 --threads 1 > tree.$*.raxml.term
	rm -f tree.$*.raxml.fa

# Get NJ tree in nexus format
tree.%.nj.nex: tree.%.nj.nwk
	$(BIN)/nwk2nex $< $@

# extract training likelihoods
tree.%.true.mod: tree.%.true.nwk tree.%.fa
	cp ../base-hky.mod tmp.$*.true.mod
	echo -n "TREE: " >> tmp.$*.true.mod
	cat $< >> tmp.$*.true.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tmp.$*.true.mod -o tree.$*.true tree.$*.fa
	rm tmp.$*.true.mod

tree.%.nj.mod: tree.%.nj.nwk tree.%.fa
	cp ../base-hky.mod tmp.$*.nj.mod
	echo -n "TREE: " >> tmp.$*.nj.mod
	cat $< >> tmp.$*.nj.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tmp.$*.nj.mod -o tree.$*.nj tree.$*.fa
	rm tmp.$*.nj.mod

tree.%.modlnl: tree.%.true.mod tree.%.nj.mod tree.%.ml.mod 
	rm -f $@
	for file in $^ ; do \
		echo -n "$${file} " >> $@ ;\
		grep LNL $${file} | awk '{print $$2}' >> $@ ;\
	done

# For an apples-to-apples comparison with the MCMC methods, vine reports the
# expected (posterior-mean) data log-likelihood E_q[lnL], so beast and mrbayes
# below report the post-burnin posterior MEAN of their per-sample likelihood
# (not the single max-likelihood draw, which is an upward-biased order
# statistic). Prior is excluded in all three (beast col 3 = likelihood, not
# posterior; mrbayes col 2 = lnLike, not lnPrior; vine LNL/LNL_mc, not LPRIOR).

# vine: final unbiased MC estimate of E_q[lnL] (LNL_mc); fall back to LNL.
tree.%.varlnl: tree.%.var.nwk.log
	echo -n "$^ " > $@
	tail -1 $^ | awk '{ll="";mc=""; for(i=1;i<=NF;i++){if($$i=="LNL:")ll=$$(i+1); if($$i=="LNL_mc:")mc=$$(i+1)} v=(mc!=""?mc:ll); gsub(/,/,"",v); print v}' >> $@

tree.%.beastlnl: tree.%.beast.log
	echo -n "$^ " > $@
	awk -v p=$(BURNIN_PCT) '$$1 ~ /^[0-9]+$$/ {v[++n]=$$3} END {skip=int(n*p/100); s=0;c=0; for(i=skip+1;i<=n;i++){s+=v[i];c++} if(c) printf "%.6f\n", s/c}' $^ >> $@

tree.%.beast-beaglelnl: tree.%.beast-beagle.log
	echo -n "$^ " > $@
	awk -v p=$(BURNIN_PCT) '$$1 ~ /^[0-9]+$$/ {v[++n]=$$3} END {skip=int(n*p/100); s=0;c=0; for(i=skip+1;i<=n;i++){s+=v[i];c++} if(c) printf "%.6f\n", s/c}' $^ >> $@

tree.%.mrbayeslnl: tree.%.mrbayes.nex.p tree.%.mrbayes.term
	echo -n "$< " > $@
	awk -v p=$(BURNIN_PCT) '$$1 ~ /^[0-9]+$$/ {v[++n]=$$2} END {skip=int(n*p/100); s=0;c=0; for(i=skip+1;i<=n;i++){s+=v[i];c++} if(c) printf "%.6f\n", s/c}' $< >> $@

tree.%.mrbayes-beaglelnl: tree.%.mrbayes-beagle.nex.p tree.%.mrbayes-beagle.term
	echo -n "$< " > $@
	awk -v p=$(BURNIN_PCT) '$$1 ~ /^[0-9]+$$/ {v[++n]=$$2} END {skip=int(n*p/100); s=0;c=0; for(i=skip+1;i<=n;i++){s+=v[i];c++} if(c) printf "%.6f\n", s/c}' $< >> $@

tree.%.raxmllnl: tree.%.raxml.term
	echo -n "$^ " > $@
	grep 'Final LogLikelihood:' $^ | awk '{printf "%.6f\n", $$3}' >> $@

tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.beast-beaglelnl tree.%.mrbayeslnl tree.%.mrbayes-beaglelnl tree.%.raxmllnl tree.%.beast-beagle.log
	cat tree.$*.modlnl tree.$*.varlnl tree.$*.beastlnl tree.$*.beast-beaglelnl tree.$*.mrbayeslnl tree.$*.mrbayes-beaglelnl tree.$*.raxmllnl | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@

# When EXCLUDE is set (e.g. EXCLUDE=3 or EXCLUDE=3,5,7), those replicates are
# omitted; when EXCLUDE is empty, all replicates are included.
eval.all.lnl.txt: $(LNL)
	echo -e "true\tnj\tml\tvine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle\traxml" > lnltmp; \
	excl=",$(EXCLUDE),"; \
	for file in $(LNL); do \
	  num=$$(basename $$file | sed 's/tree\.\([0-9]*\)\.lnl/\1/'); \
	  case "$$excl" in *",$$num,"*) continue ;; esac; \
	  awk '{printf "%s\t", $$2}' $$file >> lnltmp; \
	  echo >> lnltmp; \
	done; \
	awk '{ \
	  x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; \
	  x5 += $$5; x6 += $$6; x7 += $$7; x8 += $$8; x9 += $$9; \
	  print $$0 \
	} END { \
	  n = NR - 1; \
	  if (n > 0) \
	    printf ("-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", \
	      x1/n, x2/n, x3/n, x4/n, x5/n, x6/n, x7/n, x8/n, x9/n \
	    ) \
	}' lnltmp > $@; \
	rm -f lnltmp

# Create a version where each row's values are offset by its true value
updated.all.lnl.txt: eval.all.lnl.txt
	awk 'NR==1{print; next} \
		$$1=="-----"{print; next} \
		{t=$$1; for(i=1;i<=NF;i++) $$i=$$i-t; print}' $< > $@

# extract timing info
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

# need to get kappa from mrbayes log
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
	$(VINE_BIN)/evalTrees tree.$*.var.nwk > tmp.$*.var.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.var.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm tmp.$*.var.dist

tree.%.nj.dist.txt: tree.%.nj.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.nj.nwk > tmp.$*.nj.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.nj.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm tmp.$*.nj.dist

tree.%.ml.dist.txt: tree.%.ml.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.ml.nwk > tmp.$*.ml.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.ml.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm tmp.$*.ml.dist

tree.%.beast.dist.txt: tree.%.beast.nwk  tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.beast.nwk > tmp.$*.beast.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.beast.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm tmp.$*.beast.dist

tree.%.beast-beagle.dist.txt: tree.%.beast-beagle.nwk tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees tree.$*.beast-beagle.nwk > tmp.$*.beast-beagle.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.beast-beagle.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm tmp.$*.beast-beagle.dist

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
	nwk="$(patsubst %.rf,%.true.nwk,$(firstword $(EVALRF)))"; \
	n=$$(grep -o ',' "$$nwk" | wc -l); n=$$((n + 1)); \
	div=$$((n - 3)); [ $$div -lt 1 ] && div=1; \
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > tmprf; \
	for file in $^ ; do \
		awk -v d="$$div" '{printf "%f\t%f\t", $$2/d, $$3/d}' "$$file" >> tmprf; \
		echo >> tmprf; \
	done
	awk '{for(i=1;i<=NF;i+=2){x[i]+=$$i; xs[i]+=($$(i+1)*$$(i+1))} print $$0} END {printf "-----\n"; for(i=1;i<=15;i+=2) printf "%f\t%f%s", x[i]/(NR-1), sqrt(xs[i]/(NR-1)), (i<15 ? "\t" : "\n")}' tmprf > $@
	rm -f tmprf

eval.all.mf.txt: $(EVALMF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > tmpmf
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmpmf ;\
		echo >> tmpmf ;\
	done
	awk '{for(i=1;i<=NF;i+=2){if($$i!="nan") x[i]+=$$i; xs[i]+=($$(i+1)*$$(i+1))} print $$0} END {printf "-----\n"; for(i=1;i<=15;i+=2) printf "%f\t%f%s", x[i]/(NR-1), sqrt(xs[i]/(NR-1)), (i<15 ? "\t" : "\n")}' tmpmf > $@
	rm -f tmpmf

# topological entropy
tree.%.var.ent.txt: tree.%.var.nwk
	$(VINE_BIN)/evalTrees -e tree.$*.var.nwk | awk '{printf "%f\t", $$NF} END {printf "\n"}' > $@

tree.%.beast.ent.txt: tree.%.beast.nwk  
	$(VINE_BIN)/evalTrees -e tree.$*.beast.nwk | awk '{printf "%f\t", $$NF} END {printf "\n"}' > $@

tree.%.beast-beagle.ent.txt: tree.%.beast-beagle.nwk
	$(VINE_BIN)/evalTrees -e tree.$*.beast-beagle.nwk | awk '{printf "%f\t", $$NF} END {printf "\n"}' > $@

tree.%.ent: tree.%.var.ent.txt tree.%.beast.ent.txt tree.%.beast-beagle.ent.txt
	paste $^ > $@

eval.all.ent.txt: $(EVALENT)
	cat $^ | awk 'BEGIN {printf "vine_spl\tvine_top\tvine_br\tbeast_spl\tbeast_top\tbeast_br\tbeast-beagle_spl\tbeast-beagle_top\tbeast-beagle_br\n"} {for(i=1;i<=NF;i++) x[i]+=$$i; print $$0} END {printf "-----\n"; for(i=1;i<=9;i++) printf "%f%s", x[i]/NR, (i<9 ? "\t" : "\n")}' > $@

# for use in debugging
tracer: $(TRACER)

tree.%.tr: tree.%.var.nwk.log
	grep -v '^#' $^ > $@

clean:
	rm -rf $(TREES) $(FA) $(MOD) $(ML) $(MLMOD) $(NJMOD) $(NJ) $(VAR) $(VARNEX) $(EVALURF) $(LNL) $(VARLOG) $(VARTIME) tree.*.mean*.nwk tree.*.lnl.diffs tree.*.varlnl tree.*.modlnl eval.all.*.txt tree.*.beast* *.mf.txt *.rf.txt *.dist.txt  tree.*.time tree.*.lnl tree.*.mf tree.*.rf $(FAHELDOUT) $(TRACER)
	rm -rf tree.*.mrbayes* tree.*.nex
	rm -rf tree.*.raxml*

archive_vine:
	archive_dir="archive.vine_$$(date +%Y-%m-%d_%H:%M:%S)"; \
	mkdir $$archive_dir; \
	mv tree.*.var* $$archive_dir/; \
	mv eval.all.*.txt $$archive_dir/; \
	mv tree.*.lnl $$archive_dir/; \
	mv tree.*.time $$archive_dir/; \
	mv tree.*.rf $$archive_dir/; \
	mv tree.*.mf $$archive_dir/; \
	mv tree.*.dist $$archive_dir/; \
	mv tree.*.ent $$archive_dir/

clean_vine:
	rm -rf tree.*.var* tree.*.lnl tree.*.time tree.*.rf tree.*.mf tree.*.dist tree.*.ent eval.all.*.txt

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

clean_beast:
	rm -rf tree.*.beast*  tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt

clean-beast-beagle:
	rm -rf tree.*.beast-beagle* thinned.trees.*.beast-beagle tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf tree.*.ent eval.all.*.txt

clean_mrbayes:
	rm -rf tree.*.mrbayes* tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt

clean-mrbayes-beagle:
	rm -rf tree.*.mrbayes-beagle* tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt

clean_rf:
	rm -f eval.all.rf.txt tree.*.rf tree.*.rf.txt tmprf

clean_bsd:
	rm -f eval.all.bsd.txt tree.*.bsd tree.*.bsd.txt tmpbsd

# Branch-score distance (BSD, Kuhner-Felsenstein) to the true tree.  The
# aggregate uses point/posterior-mean-tree BSD normalized by true-tree length.
tree.%.true.bsd.txt: tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.true.nwk -b tree.$*.true.nwk > $@

tree.%.nj.bsd.txt: tree.%.nj.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.nj.nwk -b tree.$*.true.nwk > $@

tree.%.ml.bsd.txt: tree.%.ml.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.ml.nwk -b tree.$*.true.nwk > $@

tree.%.var.bsd.txt: tree.%.var.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.var.nwk -b tree.$*.true.nwk > $@

tree.%.beast.bsd.txt: tree.%.beast.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.beast.nwk -b tree.$*.true.nwk > $@

tree.%.beast-beagle.bsd.txt: tree.%.beast-beagle.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.beast-beagle.nwk -b tree.$*.true.nwk > $@

tree.%.mrbayes.bsd.txt: tree.%.mrbayes.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.mrbayes.nwk -b tree.$*.true.nwk > $@

tree.%.mrbayes-beagle.bsd.txt: tree.%.mrbayes-beagle.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.mrbayes-beagle.nwk -b tree.$*.true.nwk > $@

tree.%.bsd: tree.%.true.bsd.txt tree.%.nj.bsd.txt tree.%.ml.bsd.txt tree.%.var.bsd.txt tree.%.beast.bsd.txt tree.%.beast-beagle.bsd.txt tree.%.mrbayes.bsd.txt tree.%.mrbayes-beagle.bsd.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '/Point.*BSD:/ {pt=$$NF} /Reference tree length:/ {rl=$$NF} END {printf "%f\t%f\n", pt, rl}' $${file} >> $@ ;\
	done

eval.all.bsd.txt: $(EVALBSD)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) beast-beagle (sd) mrbayes (sd) mrbayes-beagle (sd)" > tmpbsd
	for file in $^ ; do \
		awk '{printf "%f\t", ($$3>0 ? $$2/$$3 : 0)}' "$$file" >> tmpbsd ;\
		echo >> tmpbsd ;\
	done
	awk '{for(i=1;i<=8;i++){x[i]+=$$i; xs[i]+=$$i*$$i} print $$0} END {printf "-----\n"; for(i=1;i<=8;i++){m=x[i]/(NR-1); v=xs[i]/(NR-1)-m*m; if(v<0)v=0; printf "%f\t%f%s", m, sqrt(v), (i<8 ? "\t" : "\n")}}' tmpbsd > $@
	rm -f tmpbsd
