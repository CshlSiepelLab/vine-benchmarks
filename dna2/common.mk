export SHELL=/usr/bin/bash

# Prevent Make from removing files it thinks it should clean up
.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# edit for local structure; this is the only place absolute paths are used
MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT_SUFFIX := dna2_becca_test_copy

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
	"$(BEAST)" -java -working -D fastapath=tree.$*.fa -D mcmclength=$(MCMCLEN) $< > tree.$*.beast.term

tree.%.beast.nwk: tree.%.beast-tree.trees
	python3 "$(ROOT)/time2subs.py" $^ tmp.nex
	$(OTHER_BIN)/convertTrees.py -i nexus tmp.nex > $@
	rm -f tmp.nex

tree.%nex: tree.%.fa
	$(OTHER_BIN)/fa2nex $< tree.$*.nex

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	python3 $(OTHER_BIN)/addMrbayesModelToNex.py --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex --mcmc_length $(MCMCLEN)

# Run MrBayes
tree.%.mrbayes.nex.term tree.%.mrbayes.nex.p tree.%.mrbayes.nex.t: tree.%.mrbayes.nex
	$(MRBAYES) tree.$*.mrbayes.nex > tree.$*.mrbayes.nex.term

# Run raxml
tree.%.raxml.term: tree.%.fa
	rm -f $@
	sed 's/> />/g' $< > tree.$*.raxml.fa
	$(OTHER_BIN)/raxml-ng --msa tree.$*.raxml.fa --model HKY+F --prefix tree.$*.raxml --search1 --threads 1 > tree.$*.raxml.term
	rm -f tree.$*.raxml.fa

# Run dodonaphy as a container to handle several pip/python dependencies
# Not yet sure what all the dodo inputs do...
# What is a good --temp parameter?
# We likely want the --connect Nj to be used, but it currently seems to have a bug?
# What is a good number of --epochs to use?
# There seems to be long post-processing after the main run, so is there a way to only time the main run?
tree.%.dodonaphy.term tree.%.dodonaphy-time: tree.%.nex
	mkdir tree.$*.dodonaphy
	cp $< tree.$*.dodonaphy/tree.$*.nex
	/usr/bin/time -o tree.$*.dodonaphy-time singularity exec --bind $(CURDIR)/tree.$*.dodonaphy:/mnt $(DODONAPHY_SIF) \
		dodo \
		--path_root /mnt/ \
		--path_dna $^ \
		--infer vi \
		--temp 0.5 \
		--prior "exponential" \
		--connect geodesics \
		--epochs 10 \
		--overwrite > tree.%.dodonaphy.term
	rm -f tree.$*.dodonaphy/tree.$*.nex

# Run GeoPhy (in Singularity) on NEXUS; captures stdout and wall time
tree.%.geophy.term tree.%.geophy-time: tree.%.nex
	mkdir -p tree.$*.geophy/out
	cp $< tree.$*.geophy/tree.$*.nex
	cp $(GEOPHY_CONFIG) tree.$*.geophy/config.yaml
	/usr/bin/time -o tree.$*.geophy-time \
		singularity exec \
		--bind $(CURDIR)/tree.$*.geophy:/mnt \
		$(GEOPHY_SIF) \
		bash /opt/app/scripts/run_geophy.sh \
		  -i /mnt/tree.$*.nex \
		  -o /mnt/out \
		  -c /mnt/config.yaml 
		> tree.$*.geophy.term
	rm -f tree.$*.geophy/tree.$*.nex

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

tree.%.mrbayeslnl: tree.%.mrbayes.nex.p tree.%.mrbayes.nex.term
	echo -n "$< " > $@
	grep -v '^\[' $< | grep -v '^Gen' | awk '{print $$2}' | sort -nr | head -1 | awk '{printf "%.6f\n", $$1}' >> $@

tree.%.raxmllnl: tree.%.raxml.term
	echo -n "$^ " > $@
	grep 'Final LogLikelihood:' $^ | awk '{printf "%.6f\n", $$3}' >> $@

tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.mrbayeslnl tree.%.raxmllnl
	cat $^ | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@

eval.all.lnl.txt: $(LNL)
	echo "true nj ml vine beast mrbayes raxml" > tmp
	for file in $^ ; do \
		awk '{printf "%s\t", $$3}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; x5 += $$5; x6 += $$6; x7 += $$7; print $$0} END {printf ("-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), x2/(NR-1), x3/(NR-1), x4/(NR-1), x5/(NR-1), x6/(NR-1), x7/(NR-1)) }' tmp > $@
	rm tmp

# extract timing info
tree.%.time: tree.%.beast.term tree.%.var-time tree.%.mrbayes.nex.term
	grep '^Total calculation time' tree.$*.beast.term | awk '{printf "$*\t%f\t", $$4}' > $@
	grep 'Analysis used' tree.$*.mrbayes.nex.term | awk '{printf "%s\t", $$(3)}' >> $@
	head -1 tree.$*.var-time | awk '{print $$1}' | sed 's/user//' >> $@

eval.all.time.txt: $(TIME)
	cat $(TIME) | awk 'BEGIN{printf "samp\tbeast\tmrbayes\tvine\tspeedup_beast\tspeedup_mrbayes\n"} {x1 += $$2; x2 += $$3; x3 += $$4; printf "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", $$1, $$2, $$3, $$4, $$2/$$4, $$3/$$4} END { printf "-----------------------------------------\nall\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", x1/NR, x2/NR, x3/NR, x1/x3, x2/x3 }' > $@

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
	rm -rf tree.*.mrbayes*
	rm -rf tree.*.raxml*
	rm -rf tree.*.dodonaphy
	rm -rf tree.*.geophy
