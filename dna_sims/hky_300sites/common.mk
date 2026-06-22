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
# VINE_BIN := $(BIN)/vine/bin
VINE_BIN := /home/staklins/projects/vine_project/vine/bin
BEAST_BIN := $(MAIN_DIR)/bin/beast/bin
BEAST := $(BEAST_BIN)/beast
MRBAYES := $(BIN)/mb
BEAST_TEMPLATE := $(ROOT)/beast_template.xml
CONTAINERS := $(MAIN_DIR)/containers
DODONAPHY_SIF := $(CONTAINERS)/dodonaphy/dodonaphy.sif
GEOPHY_SIF := $(CONTAINERS)/geophy/geophy.sif
GEOPHY_CONFIG := $(CONTAINERS)/geophy/default.yaml
BURNIN_PCT := 10


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

MRBAYESLOG := $(patsubst %.true.nwk,%.mrbayes.nex.p,$(TREES))

# evalTrees stuff
FAHELDOUT := $(patsubst %.true.nwk,%.heldout.fa,$(TREES))

all: eval.all.lnl.txt eval.all.rf.txt eval.all.mf.txt eval.all.time.txt eval.all.dist.txt eval.all.ent.txt

simulate: $(FA)
max_lik: $(ML) $(RAXML)
mcmc: $(BEASTLOG) $(MRBAYESLOG)
beast: $(BEASTLOG)
mrbayes: $(MRBAYESLOG)
ml: $(ML)
beastnwk: $(BEASTNWK)
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
	"$(BEAST)" -java -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$(MCMC_SAMPLEFREQ) -D printfreq=$(MCMC_PRINTFREQ) $< > tree.$*.beast.term

tree.%.beast.nwk: tree.%.beast-tree.trees
	$(BEAST_BIN)/logcombiner -log $^ -o thinned.trees.$*.beast \
	  -burnin $(BURNIN_PCT) -resample $$(( $(BEAST_MCMCLEN) / 1000 ))
	python3 "$(PYTHON_SRC)/time2subs.py" thinned.trees.$*.beast tmp.$*.beast.nex
	$(BIN)/convertTrees -i nexus tmp.$*.beast.nex > $@
	rm -f tmp.$*.beast.nex

# Old method - now doing pilot runs to determine mcmc chain length ahead of time
# # Calculate ESS-based runtime scaling factor for when the chain converged
# beast_ess_runtime_scale_factor.txt: $(BEASTLOG)
# 	$(BIN)/ess_for_dataset_replicates \
# 		--logfiles "$$(echo $(BEASTLOG) | tr ' ' ',')"  \
# 		--parameters "Tree.Length,Tree.height" \
# 		--outputfile $@ \
# 		--ess_threshold 625 \
# 		--burnin 0.1

tree.%.nex: tree.%.fa
	$(BIN)/fa2nex $< $@

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	$(BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex --mcmc_length $(MRBAYES_MCMCLEN) --model HKY \
		--sample_freq $(MCMC_SAMPLEFREQ) --print_freq $(MCMC_PRINTFREQ) --diagn_freq $(MCMC_PRINTFREQ)

# Run MrBayes
tree.%.mrbayes.term tree.%.mrbayes.nex.p tree.%.mrbayes.nex.t: tree.%.mrbayes.nex
	$(MRBAYES) tree.$*.mrbayes.nex > tree.$*.mrbayes.term

# Get mrbayes tree in nexus format
tree.%.mrbayes.nwk: tree.%.mrbayes.nex.t
	SKIP=$$(awk -v p=$(BURNIN_PCT) \
	  '/^[[:space:]]*tree[[:space:]]+/{c++} END{printf "%d", int(c*p/100)}' $<); \
	THIN=$$(( $(MRBAYES_MCMCLEN) / $(MCMC_SAMPLEFREQ) / 1000 )); \
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

# Old method - now doing pilot runs to determine mcmc chain length ahead of time
# # Calculate ESS-based runtime scaling factor for when the chain converged
# mrbayes_ess_runtime_scale_factor.txt: $(MRBAYESLOG)
# 	$(BIN)/ess_for_dataset_replicates \
# 		--logfiles "$$(echo $(MRBAYESLOG) | tr ' ' ',')"  \
# 		--parameters "TL" \
# 		--outputfile $@ \
# 		--ess_threshold 625 \
# 		--burnin 0.1


# Run raxml
tree.%.raxml.term: tree.%.fa
	rm -f $@
	sed 's/> />/g' $< > tree.$*.raxml.fa
	$(BIN)/raxml-ng --msa tree.$*.raxml.fa --model HKY+F --prefix tree.$* --search1 --threads 1 > tree.$*.raxml.term
	rm -f tree.$*.raxml.fa

# Get NJ tree in nexus format
tree.%.nj.nex: tree.%.nj.nwk
	$(BIN)/nwk2nex $< $@

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

tree.%.varlnl: tree.%.var.nwk.log 
	echo -n "$^ " > $@
	tail -1 $^ | awk '{print $$11}' | sed 's/,//g' >> $@ 

tree.%.beastlnl: tree.%.beast.log
	echo -n "$^ " > $@
	grep -v '^#' $^ | grep -v '^Sample' | awk '{print $$3}' | sort -nr | head -1 >> $@

tree.%.mrbayeslnl: tree.%.mrbayes.nex.p tree.%.mrbayes.term
	echo -n "$< " > $@
	grep -v '^\[' $< | grep -v '^Gen' | awk '{print $$2}' | sort -gr | head -1 | awk '{printf "%.6f\n", $$1}' >> $@

tree.%.raxmllnl: tree.%.raxml.term
	echo -n "$^ " > $@
	grep 'Final LogLikelihood:' $^ | awk '{printf "%.6f\n", $$3}' >> $@

tree.%.dodonaphylnl: tree.%.dodonaphy.elbo.txt
	echo -n "$^ " > $@
	awk 'NR==1 || $$1>max {max=$$1} END {printf "%.6f\n", max}' $^ >> $@

tree.%.geophylnl: tree.%.geophy.eval.latest.txt
	echo -n "$^ " > $@
	awk 'NR==2 {printf "%.6f\n", $$6}' $^ >> $@

#tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.mrbayeslnl tree.%.raxmllnl tree.%.dodonaphylnl tree.%.geophylnl
tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.mrbayeslnl tree.%.raxmllnl 
	cat $^ | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@

# When EXCLUDE is set (e.g. EXCLUDE=3 or EXCLUDE=3,5,7), those replicates are
# omitted; when EXCLUDE is empty, all replicates are included.
eval.all.lnl.txt: $(LNL)
	echo -e "true\tnj\tml\tvine\tbeast\tmrbayes\traxml" > lnltmp; \
	excl=",$(EXCLUDE),"; \
	for file in $(LNL); do \
	  num=$$(basename $$file | sed 's/tree\.\([0-9]*\)\.lnl/\1/'); \
	  case "$$excl" in *",$$num,"*) continue ;; esac; \
	  awk '{printf "%s\t", $$2}' $$file >> lnltmp; \
	  echo >> lnltmp; \
	done; \
	awk '{ \
	  x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; \
	  x5 += $$5; x6 += $$6; x7 += $$7; \
	  print $$0 \
	} END { \
	  n = NR - 1; \
	  if (n > 0) \
	    printf ("-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\n", \
	      x1/n, x2/n, x3/n, x4/n, x5/n, x6/n, x7/n \
	    ) \
	}' lnltmp > $@; \
	rm -f lnltmp

# Create a version where each row's values are offset by its true value
updated.all.lnl.txt: eval.all.lnl.txt
	awk 'NR==1{print; next} \
		$$1=="-----"{print; next} \
		{t=$$1; for(i=1;i<=NF;i++) $$i=$$i-t; print}' $< > $@

# extract timing info
#tree.%.time: tree.%.beast.term tree.%.var-time tree.%.mrbayes.term tree.%.raxml.term tree.%.dodonaphy-time tree.%.geophy-time
tree.%.time: tree.%.beast.term tree.%.var-time tree.%.mrbayes.term tree.%.raxml.term # beast_ess_runtime_scale_factor.txt mrbayes_ess_runtime_scale_factor.txt
	echo -e "samp\tbeast\tmrbayes\tvine\traxml" > $@; \
	beast_time=$$(grep '^Total calculation time' tree.$*.beast.term | awk '{print $$4}'); \
	printf "$*\t%s\t" "$$beast_time" >> $@; \
	mrbayes_time=$$(grep 'Analysis used' tree.$*.mrbayes.term | awk '{printf "%s\t", $$(3)}'); \
	printf "%s\t" "$$mrbayes_time" >> $@; \
	head -1 tree.$*.var-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@; \
	grep 'Elapsed time:' tree.$*.raxml.term | awk '{printf "%s\n", $$3}' >> $@
#	head -1 tree.$*.dodonaphy-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@
#	head -1 tree.$*.geophy-time | awk '{printf "%s\n", $$1}' | sed 's/user//' >> $@

# Old method - now doing pilot runs to determine mcmc chain length ahead of time
# @beast_scale_factor=$$(cat beast_ess_runtime_scale_factor.txt); \
# mrbayes_scale_factor=$$(cat mrbayes_ess_runtime_scale_factor.txt); \
# echo -e "samp\tbeast\tmrbayes\tvine\traxml" > $@; \
# beast_time=$$(grep '^Total calculation time' tree.$*.beast.term | awk '{print $$4}'); \
# scaled_beast_time=$$(awk -v t="$$beast_time" -v s="$$beast_scale_factor" 'BEGIN{printf "%f", t*s}'); \
# printf "$*\t%s\t" "$$scaled_beast_time" >> $@; \
# mrbayes_time=$$(grep 'Analysis used' tree.$*.mrbayes.term | awk '{printf "%s\t", $$(3)}'); \
# scaled_mrbayes_time=$$(awk -v t="$$mrbayes_time" -v s="$$mrbayes_scale_factor" 'BEGIN{printf "%f", t*s}'); \
# printf "%s\t" "$$scaled_mrbayes_time" >> $@; \

eval.all.time.txt: $(TIME)
	awk 'FNR==1 && NR==1 {print; next} FNR==2 {print; for(i=2;i<=NF;i++) sum[i]+=$$i; n++} END {if(n>0){printf "-----------------------------------------\nall"; for(i=2;i<=NF;i++) printf "\t%.2f", sum[i]/n; printf "\n"}}' $(TIME) > $@

# evalTrees stuff
# (1) modelFit
# extract kappa from vine log
tree.%.var.mf.txt: tree.%.var.nwk tree.%.heldout.fa tree.%.var.nwk.log
	kappa=`tail -1 tree.$*.var.nwk.log | awk '{print $$19}'` ;\
	$(VINE_BIN)/evalTrees tree.$*.var.nwk -f tree.$*.heldout.fa -k $$kappa > $@

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


tree.%.mf: tree.%.true.mf.txt tree.%.nj.mf.txt tree.%.ml.mf.txt tree.%.var.mf.txt tree.%.beast.mf.txt tree.%.mrbayes.mf.txt
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

tree.%.dist: tree.%.ml.dist.txt tree.%.var.dist.txt tree.%.beast.dist.txt
	paste $^ > $@

eval.all.dist.txt: $(EVALDIST)
	cat $^ | awk '{print $$1, $$2, $$6, $$7, $$8, $$9, $$11, $$12, $$13, $$14}' | awk 'BEGIN {printf "ML_r2 ML_RMSE vine_r2 vine_RMSE vine_95CI vine_50CI beast_r2 beast_RMSE beast_95CI beast_50CI\n"} {x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; x5 += $$5; x6 += $$6; x7 += $$7; x8 += $$8; x9 += $$9; x10 += $$10; print $$0} END{printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/NR, x2/NR, x3/NR, x4/NR, x5/NR, x6/NR, x7/NR, x8/NR, x9/NR, x10/NR}' > $@

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

tree.%.mrbayes.rf.txt: tree.%.mrbayes.nwk tree.%.true.nwk
	$(VINE_BIN)/evalTrees tree.$*.mrbayes.nwk -t tree.$*.true.nwk > $@

tree.%.rf: tree.%.true.rf.txt tree.%.nj.rf.txt tree.%.ml.rf.txt tree.%.var.rf.txt tree.%.beast.rf.txt tree.%.mrbayes.rf.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
	done

eval.all.rf.txt: $(EVALRF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) mrbayes (sd)" > tmprf
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmprf ;\
		echo >> tmprf ;\
	done
	awk '{x1 += $$1; x1s += ($$2 * $$2); x2 += $$3; x2s += ($$4 * $$4); x3 += $$5; x3s += ($$6 * $$6); x4 += $$7; x4s += ($$8 * $$8); x5 += $$9; x5s += ($$10 * $$10); x6 += $$11; x6s += ($$12 * $$12); print $$0} END {printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), x3/(NR-1), sqrt(x3s/(NR-1)), x4/(NR-1), sqrt(x4s/(NR-1)), x5/(NR-1), sqrt(x5s/(NR-1)), x6/(NR-1), sqrt(x6s/(NR-1))}' tmprf > $@
	rm -f tmprf

eval.all.mf.txt: $(EVALMF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd) mrbayes (sd)" > tmpmf
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmpmf ;\
		echo >> tmpmf ;\
	done
	awk '{x1 += $$1; x1s += ($$2 * $$2); if ($$3 != "nan") x2 += $$3; x2s += ($$4 * $$4); if ($$5 != "nan") x3 += $$5; x3s += ($$6 * $$6); x4 += $$7; x4s += ($$8 * $$8); x5 += $$9; x5s += ($$10 * $$10); if ($$11 != "nan") x6 += $$11; x6s += ($$12 * $$12); print $$0} END {printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), x3/(NR-1), sqrt(x3s/(NR-1)), x4/(NR-1), sqrt(x4s/(NR-1)), x5/(NR-1), sqrt(x5s/(NR-1)), x6/(NR-1), sqrt(x6s/(NR-1))}' tmpmf > $@
	rm -f tmpmf

# topological entropy
tree.%.var.ent.txt: tree.%.var.nwk
	$(VINE_BIN)/evalTrees -e tree.$*.var.nwk | awk '{printf "%f\t", $$NF} END {printf "\n"}' > $@

tree.%.beast.ent.txt: tree.%.beast.nwk  
	$(VINE_BIN)/evalTrees -e tree.$*.beast.nwk | awk '{printf "%f\t", $$NF} END {printf "\n"}' > $@

tree.%.ent: tree.%.var.ent.txt tree.%.beast.ent.txt
	paste $^ > $@

eval.all.ent.txt: $(EVALENT)
	cat $^ | awk 'BEGIN {printf "vine_spl\tvine_top\tvine_br\t\tbeast_spl\tbeast_top\tbeast_br\n"} {x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; x5 += $$5; x6 += $$6; print $$0} END {printf "-----\n%f\t%f\t%f\t\t%f\t%f\t%f\n", x1/NR, x2/NR, x3/NR, x4/NR, x5/NR, x6/NR}' > $@

# for use in debugging
tracer: $(TRACER)

tree.%.tr: tree.%.var.nwk.log
	grep -v '^#' $^ > $@

clean:
	rm -rf $(TREES) $(FA) $(MOD) $(ML) $(MLMOD) $(NJMOD) $(NJ) $(VAR) $(VARNEX) $(EVALURF) $(LNL) $(VARLOG) $(VARTIME) tree.*.mean*.nwk tree.*.lnl.diffs tree.*.varlnl tree.*.modlnl eval.all.*.txt tree.*.beast* *.mf.txt *.rf.txt *.dist.txt  tree.*.time tree.*.lnl tree.*.mf tree.*.rf $(FAHELDOUT) $(TRACER)
	rm -rf tree.*.mrbayes* tree.*.nex
	rm -rf tree.*.raxml*
	rm -rf tree.*.dodonaphy*
	rm -rf tree.*.geophy*

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

clean_mrbayes:
	rm -rf tree.*.mrbayes* tree.*.dist tree.*.time tree.*.lnl tree.*.mf tree.*.rf eval.all.*.txt
