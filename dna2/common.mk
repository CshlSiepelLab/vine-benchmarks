export SHELL=/usr/bin/bash

# Prevent Make from removing files it thinks it should clean up
.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# edit for local structure; this is the only place absolute paths are used
MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT_SUFFIX := dna2

ROOT := $(MAIN_DIR)/$(ROOT_SUFFIX)
PHAST_BIN := $(MAIN_DIR)/phast/bin
OTHER_BIN := $(MAIN_DIR)/bin
BEAST := $(MAIN_DIR)/beast/bin/beast
MRBAYES := $(OTHER_BIN)/mb
BEAST_TEMPLATE := $(ROOT)/beast_template.xml
DODONAPHY_SIF := $(MAIN_DIR)/dodonaphy/dodonaphy.sif
GEOPHY_SIF := $(MAIN_DIR)/geophy/geophy.sif
GEOPHY_CONFIG := $(MAIN_DIR)/geophy/default.yaml

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

# evalTrees stuff
FAHELDOUT := $(patsubst %.true.nwk,%.heldout.fa,$(TREES))

all: eval.all.lnl.txt eval.all.rf.txt eval.all.mf.txt eval.all.time.txt eval.all.dist.txt

tree.%.true.nwk: 
#	bdTree.py -b 5 -d 5 -n $(NTAXA) | sed 's/\[\&[UR]\] //' > $@
	$(OTHER_BIN)/bdTree2.py -b 1 -d 0.5 --oversample-k 3 --height 5 --min-edge 0.02 --expected-height $(EXPHEIGHT) --no-stem -n $(NTAXA) | sed 's/\[\&[UR]\] //' > $@

tree.%.fa: tree.%.true.nwk
	cp ../base-hky.mod tmp.mod
	echo -n "TREE: " >> tmp.mod
	cat $< >> tmp.mod
	$(PHAST_BIN)/base_evolve --nsites $(NSITES) tmp.mod > $@
	rm tmp.mod

tree.%.heldout.fa: tree.%.true.nwk
	cp ../base-hky.mod tmp.mod
	echo -n "TREE: " >> tmp.mod
	cat $< >> tmp.mod
	$(PHAST_BIN)/base_evolve --nsites $(NSITES) tmp.mod > $@
	rm tmp.mod

tree.%.nj.nwk: tree.%.fa
	$(PHAST_BIN)/vine --nj-only $< > $@

tree.%.ml.mod: tree.%.nj.nwk tree.%.fa
	$(PHAST_BIN)/phyloFit --subst-mod HKY85 --tree $^ -o tree.$*.ml

tree.%.ml.nwk: tree.%.ml.mod
	$(PHAST_BIN)/tree_doctor --tree-only $^ > $@

tree.%.var.nwk tree.%.var-time tree.%.var.nwk.log: tree.%.fa 
	/usr/bin/time -o tree.$*.var-time $(PHAST_BIN)/vine $< -l tree.$*.var.nwk.log $(VAROPT) --mean tree.$*.mean.nwk > tree.$*.var.nwk

tree.%.beast.xml:
	cp "$(BEAST_TEMPLATE)" $@

tree.%.beast.term tree.%.beast-tree.trees tree.%.beast.log: tree.%.beast.xml tree.%.fa
	rm -f tree.$*.beast-tree.trees tree.$*.beast.log
	"$(BEAST)" -java -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) $< > tree.$*.beast.term

tree.%.beast.nwk: tree.%.beast-tree.trees
	python3 "$(ROOT)/time2subs.py" $^ tmp.nex
	$(OTHER_BIN)/convertTrees.py -i nexus tmp.nex > $@
	rm -f tmp.nex

tree.%.nex: tree.%.fa
	$(OTHER_BIN)/fa2nex $< $@

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	python3 $(OTHER_BIN)/addMrbayesModelToNex.py --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex --mcmc_length $(MRBAYES_MCMCLEN)

# Run MrBayes
tree.%.mrbayes.term tree.%.mrbayes.nex.p tree.%.mrbayes.nex.t: tree.%.mrbayes.nex
	$(MRBAYES) tree.$*.mrbayes.nex > tree.$*.mrbayes.term

# Run raxml
tree.%.raxml.term: tree.%.fa
	rm -f $@
	sed 's/> />/g' $< > tree.$*.raxml.fa
	$(OTHER_BIN)/raxml-ng --msa tree.$*.raxml.fa --model HKY+F --prefix tree.$* --search1 --threads 1 > tree.$*.raxml.term
	rm -f tree.$*.raxml.fa

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
		--temp 0.00001 \
		--prior "exponential" \
		--connect nj \
		--boosts 3 \
		--importance 3 \
		--curv -100 \
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
		  -c /mnt/tree.$*.geophy.config.yaml 
		> tree.$*.geophy.term
	rm -f tree.$*.geophy.nex

# extract training likelihoods
tree.%.true.mod: tree.%.true.nwk tree.%.fa
	cp ../base-hky.mod tmp.mod
	echo -n "TREE: " >> tmp.mod
	cat $< >> tmp.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tmp.mod -o tree.$*.true tree.$*.fa
	rm tmp.mod

tree.%.nj.mod: tree.%.nj.nwk tree.%.fa
	cp ../base-hky.mod tmp.mod
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

tree.%.mrbayeslnl: tree.%.mrbayes.nex.p tree.%.mrbayes.term
	echo -n "$< " > $@
	grep -v '^\[' $< | grep -v '^Gen' | awk '{print $$2}' | sort -nr | head -1 | awk '{printf "%.6f\n", $$1}' >> $@

tree.%.raxmllnl: tree.%.raxml.term
	echo -n "$^ " > $@
	grep 'Final LogLikelihood:' $^ | awk '{printf "%.6f\n", $$3}' >> $@

tree.%.dodonaphylnl: tree.%.dodonaphy.elbo.txt
	echo -n "$^ " > $@
	awk 'NR==1 || $$1>max {max=$$1} END {printf "%.6f\n", max}' $^ >> $@

tree.%.geophylnl: tree.%.geophy.eval.latest.txt
	echo -n "$^ " > $@
	awk 'NR==2 {printf "%.6f\n", $$4}' $^ >> $@

tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.mrbayeslnl tree.%.raxmllnl tree.%.dodonaphylnl tree.%.geophylnl
	cat $^ | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@

eval.all.lnl.txt: $(LNL)
	echo -e "true\tnj\tml\tvine\tbeast\tmrbayes\traxml\tdodonaphy\tgeophy" > tmp
	for file in $^ ; do \
		awk '{printf "%s\t", $$3}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; x5 += $$5; x6 += $$6; x7 += $$7; x8 += $$8; x9 += $$9; print $$0} END {printf ("-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), x2/(NR-1), x3/(NR-1), x4/(NR-1), x5/(NR-1), x6/(NR-1), x7/(NR-1), x8/(NR-1), x9/(NR-1)) }' tmp > $@
	rm tmp

# extract timing info
tree.%.time: tree.%.beast.term tree.%.var-time tree.%.mrbayes.term tree.%.raxml.term tree.%.dodonaphy-time tree.%.geophy-time
	echo -e "samp\tbeast\tmrbayes\tvine\traxml\tdodonaphy\tgeophy" > $@
	grep '^Total calculation time' tree.$*.beast.term | awk '{printf "$*\t%f\t", $$4}' >> $@
	grep 'Analysis used' tree.$*.mrbayes.term | awk '{printf "%s\t", $$(3)}' >> $@
	head -1 tree.$*.var-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@
	grep 'Elapsed time:' tree.$*.raxml.term | awk '{printf "%s\t", $$3}' >> $@
	head -1 tree.$*.dodonaphy-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@
	head -1 tree.$*.geophy-time | awk '{printf "%s\n", $$1}' | sed 's/user//' >> $@

eval.all.time.txt: $(TIME)
	awk 'FNR==1 && NR==1 {print; next} FNR==2 {print; for(i=2;i<=NF;i++) sum[i]+=$$i; n++} END {if(n>0){printf "-----------------------------------------\nall"; for(i=2;i<=NF;i++) printf "\t%.2f", sum[i]/n; printf "\n"}}' $(TIME) > $@

# evalTrees stuff
# (1) modelFit
# extract kappa from vine log
tree.%.var.mf.txt: tree.%.var.nwk tree.%.heldout.fa tree.%.var.nwk.log
	kappa=`tail -1 tree.$*.var.nwk.log | awk '{print $$19}'` ;\
	$(PHAST_BIN)/evalTrees tree.$*.var.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# use true kappa
tree.%.true.mf.txt: tree.%.true.nwk tree.%.heldout.fa
	$(PHAST_BIN)/evalTrees tree.$*.true.nwk -f tree.$*.heldout.fa -k 4 > $@

# use ML kappa
tree.%.nj.mf.txt: tree.%.nj.nwk tree.%.heldout.fa tree.%.ml.mod
	kappa=`awk '$$1 == "BACKGROUND:" {pi_c = $$3; pi_g = $$4 } $$1<0 {print ($$3/pi_g) / ($$2/pi_c)}' tree.$*.ml.mod` ;\
	$(PHAST_BIN)/evalTrees tree.$*.nj.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# use ML kappa
tree.%.ml.mf.txt: tree.%.ml.nwk tree.%.heldout.fa tree.%.ml.mod
	kappa=`awk '$$1 == "BACKGROUND:" {pi_c = $$3; pi_g = $$4 } $$1<0 {print ($$3/pi_g) / ($$2/pi_c)}' tree.$*.ml.mod` ;\
	$(PHAST_BIN)/evalTrees tree.$*.ml.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# use posterior mean kappa from beast log
tree.%.beast.mf.txt: tree.%.beast.nwk tree.%.heldout.fa tree.%.beast.log
	kappa=`awk '{if (inlog) {sum += $$10; n++} ; if ($$1 == "Sample") inlog=1} END {print sum/n}' tree.$*.beast.log` ;\
	$(PHAST_BIN)/evalTrees tree.$*.beast.nwk -f tree.$*.heldout.fa -k $$kappa > $@

tree.%.mf: tree.%.true.mf.txt tree.%.nj.mf.txt tree.%.ml.mf.txt tree.%.var.mf.txt tree.%.beast.mf.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
	done

# (2) distances
tree.%.true.dist.txt: tree.%.true.nwk 
	$(PHAST_BIN)/evalTrees tree.$*.true.nwk > $@

tree.%.var.dist.txt: tree.%.var.nwk tree.%.true.dist.txt
	$(PHAST_BIN)/evalTrees tree.$*.var.nwk > tmp
	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

tree.%.nj.dist.txt: tree.%.nj.nwk  tree.%.true.dist.txt
	$(PHAST_BIN)/evalTrees tree.$*.nj.nwk > tmp
	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

tree.%.ml.dist.txt: tree.%.ml.nwk  tree.%.true.dist.txt
	$(PHAST_BIN)/evalTrees tree.$*.ml.nwk > tmp
	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

tree.%.beast.dist.txt: tree.%.beast.nwk  tree.%.true.dist.txt
	$(PHAST_BIN)/evalTrees tree.$*.beast.nwk > tmp
	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

tree.%.dist: tree.%.ml.dist.txt tree.%.var.dist.txt tree.%.beast.dist.txt
	paste $^ > $@

eval.all.dist.txt: $(EVALDIST)
	cat $^ | awk '{print $$1, $$2, $$6, $$7, $$8, $$9, $$11, $$12, $$13, $$14}' | awk 'BEGIN {printf "ML_r2 ML_RMSE vine_r2 vine_RMSE vine_95CI vine_50CI beast_r2 beast_RMSE beast_95CI beast_50CI\n"} {x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; x5 += $$5; x6 += $$6; x7 += $$7; x8 += $$8; x9 += $$9; x10 += $$10; print $$0} END{printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/NR, x2/NR, x3/NR, x4/NR, x5/NR, x6/NR, x7/NR, x8/NR, x9/NR, x10/NR}' > $@

# (3) RF dist
tree.%.var.rf.txt: tree.%.var.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.var.nwk -t tree.$*.true.nwk > $@

tree.%.true.rf.txt: tree.%.true.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.true.nwk -t tree.$*.true.nwk > $@

tree.%.nj.rf.txt: tree.%.nj.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.nj.nwk -t tree.$*.true.nwk > $@

tree.%.ml.rf.txt: tree.%.ml.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.ml.nwk -t tree.$*.true.nwk > $@

tree.%.beast.rf.txt: tree.%.beast.nwk tree.%.true.nwk
	$(PHAST_BIN)/evalTrees tree.$*.beast.nwk -t tree.$*.true.nwk > $@

tree.%.rf: tree.%.true.rf.txt tree.%.nj.rf.txt tree.%.ml.rf.txt tree.%.var.rf.txt tree.%.beast.rf.txt
	rm -f $@
	for file in $^ ; do \
		echo -n "$$file     " >> $@ ;\
		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
	done

eval.all.rf.txt: $(EVALRF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd)" > tmp
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{x1 += $$1; x1s += ($$2 * $$2); x2 += $$3; x2s += ($$4 * $$4); x3 += $$5; x3s += ($$6 * $$6); x4 += $$7; x4s += ($$8 * $$8); x5 += $$9; x5s += ($$10 * $$10); print $$0} END {printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), x3/(NR-1), sqrt(x3s/(NR-1)), x4/(NR-1), sqrt(x4s/(NR-1)), x5/(NR-1), sqrt(x5s/(NR-1))}' tmp > $@
	rm -f tmp

eval.all.mf.txt: $(EVALMF)
	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd)" > tmp
	for file in $^ ; do \
		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{x1 += $$1; x1s += ($$2 * $$2); if ($$3 != "nan") x2 += $$3; x2s += ($$4 * $$4); if ($$5 != "nan") x3 += $$5; x3s += ($$6 * $$6); x4 += $$7; x4s += ($$8 * $$8); x5 += $$9; x5s += ($$10 * $$10); print $$0} END {printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), x3/(NR-1), sqrt(x3s/(NR-1)), x4/(NR-1), sqrt(x4s/(NR-1)), x5/(NR-1), sqrt(x5s/(NR-1))}' tmp > $@
	rm -f tmp

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
