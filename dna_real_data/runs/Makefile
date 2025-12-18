export SHELL=/usr/bin/bash

# Method settings
# BEAST_MCMCLEN := 4000000
# MRBAYES_MCMCLEN := 4000000

tree.ds1.beast.%: BEAST_MCMCLEN := 1800000
tree.ds2.beast.%: BEAST_MCMCLEN := 1800000
tree.ds3.beast.%: BEAST_MCMCLEN := 1800000
tree.ds4.beast.%: BEAST_MCMCLEN := 1733333
tree.ds5.beast.%: BEAST_MCMCLEN := 1633333
tree.ds6.beast.%: BEAST_MCMCLEN := 1800000
tree.ds7.beast.%: BEAST_MCMCLEN := 1800000
tree.ds8.beast.%: BEAST_MCMCLEN := 1733333

tree.ds1.mrbayes.%: MRBAYES_MCMCLEN := 566667
tree.ds2.mrbayes.%: MRBAYES_MCMCLEN := 666667
tree.ds3.mrbayes.%: MRBAYES_MCMCLEN := 766667
tree.ds4.mrbayes.%: MRBAYES_MCMCLEN := 733333
tree.ds5.mrbayes.%: MRBAYES_MCMCLEN := 733333
tree.ds6.mrbayes.%: MRBAYES_MCMCLEN := 800000
tree.ds7.mrbayes.%: MRBAYES_MCMCLEN := 800000
tree.ds8.mrbayes.%: MRBAYES_MCMCLEN := 933333

NUM_MCMC_SAMPLES := 10000
MCMC_PRINTFREQ := 1000
VINE_NSAMP := 1000
VAROPT?=-M 300 -b 10 -c 100 -D 7 --hky85 -s $(VINE_NSAMP) -v 0


# Prevent Make from removing files it thinks it should clean up
#.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# edit for local structure; this is the only place absolute paths are used
MAIN_DIR := /local/storage/no-backup/vine-benchmarks

PYTHON_SRC := $(MAIN_DIR)/python/src
BIN := $(MAIN_DIR)/bin
VINE_BIN := $(BIN)/vine/bin
PHAST_BIN := $(BIN)/phast/bin
BEAST := $(BIN)/beast/bin/beast
MRBAYES := $(BIN)/mb
BEAST_TEMPLATE := $(MAIN_DIR)/dna_sims/hky_300sites/beast_template.xml
CONTAINERS := $(MAIN_DIR)/containers
DODONAPHY_SIF := $(CONTAINERS)/dodonaphy/dodonaphy.sif
GEOPHY_SIF := $(CONTAINERS)/geophy/geophy.sif
GEOPHY_CONFIG := $(CONTAINERS)/geophy/default.yaml
# VAIPHY_SIF := $(CONTAINERS)/vbpi_gnn_env/vaiphy.sif

# Use all real FASTA files in the directory with the real DS1-8 data
FA_DIR := /local/storage/no-backup/vine-benchmarks/dna_real_data/data
FA := $(sort $(wildcard $(FA_DIR)/*.fa))

# Strip paths and extensions for consistent base names
BASENAMES := $(notdir $(basename $(FA)))

NJ := $(addprefix tree.,$(addsuffix .nj.nwk,$(BASENAMES)))
ML := $(addprefix tree.,$(addsuffix .ml.nwk,$(BASENAMES)))
VAR := $(addprefix tree.,$(addsuffix .var.nwk,$(BASENAMES)))
BEASTLOG := $(addprefix tree.,$(addsuffix .beast.log,$(BASENAMES)))
MRBAYESLOG := $(addprefix tree.,$(addsuffix .mrbayes.nex.p,$(BASENAMES)))

LNL := $(addprefix tree.,$(addsuffix .lnl,$(BASENAMES)))
TIME := $(addprefix tree.,$(addsuffix .time,$(BASENAMES)))


all: eval.all.lnl.txt eval.all.time.txt

infer: $(NJ) $(ML) $(VAR) $(BEASTLOG) $(MRBAYESLOG)

tree.%.nj.nwk: $(FA_DIR)/%.fa
	$(VINE_BIN)/vine --nj-only $< > $@

tree.%.ml.mod: tree.%.nj.nwk $(FA_DIR)/%.fa
	$(PHAST_BIN)/phyloFit --subst-mod HKY85 --tree $^ -o tree.$*.ml

tree.%.ml.nwk: tree.%.ml.mod
	$(PHAST_BIN)/tree_doctor --tree-only $^ > $@

tree.%.var.nwk tree.%.var-time tree.%.var.nwk.log: $(FA_DIR)/%.fa 
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine $< -l tree.$*.var.nwk.log $(VAROPT) --mean tree.$*.mean.nwk > tree.$*.var.nwk

tree.%.beast.xml:
	cp "$(BEAST_TEMPLATE)" $@

tree.%.beast.term tree.%.beast-tree.trees tree.%.beast.log: tree.%.beast.xml $(FA_DIR)/%.fa
	rm -f tree.$*.beast-tree.trees tree.$*.beast.log
	"$(BEAST)" -java -working -D fastapath=$(FA_DIR)/$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$$(( $(BEAST_MCMCLEN) / $(NUM_MCMC_SAMPLES) )) -D printfreq=$(MCMC_PRINTFREQ) $< > tree.$*.beast.term

tree.%.beast.nwk: tree.%.beast-tree.trees
	python3 "$(PYTHON_SRC)/time2subs.py" $^ tmp.nex
	$(BIN)/convertTrees -i nexus tmp.nex > $@
	rm -f tmp.nex

tree.%.nex: $(FA_DIR)/%.fa
	$(BIN)/fa2nex $< $@

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	$(BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex --mcmc_length $(MRBAYES_MCMCLEN) --model HKY \
		--sample_freq $$(( $(MRBAYES_MCMCLEN) / $(NUM_MCMC_SAMPLES) )) --print_freq $(MCMC_PRINTFREQ) --diagn_freq $(MCMC_PRINTFREQ)

# Run MrBayes
tree.%.mrbayes.term tree.%.mrbayes.nex.p tree.%.mrbayes.nex.t: tree.%.mrbayes.nex
	$(MRBAYES) tree.$*.mrbayes.nex > tree.$*.mrbayes.term

# Get NJ tree in nexus format
tree.%.nj.nex: tree.%.nj.nwk
	$(BIN)/nwk2nex $< $@

# Run raxml
tree.%.raxml.term: $(FA_DIR)/%.fa
	rm -f $@
	sed 's/> />/g' $< > tree.$*.raxml.fa
	$(BIN)/raxml-ng --msa tree.$*.raxml.fa --model HKY+F --prefix tree.$* --search1 --threads 1 > tree.$*.raxml.term
	rm -f tree.$*.raxml.fa

# tree.%.vaiphy.term: tree.%.nex
# 	mkdir -p tree.$*
# 	cp tree.$*.nex tree.$*/tree.$*.nex
# 	/usr/bin/time -o tree.$*.vaiphy-time \
# 	singularity exec \
# 		--bind $(CURDIR):/mnt \
# 		--bind /local/storage/no-backup/vine-benchmarks/VaiPhy:/opt/vaiphy \
# 		--pwd /mnt \
# 		--env PYTHONWARNINGS=ignore \
# 		$(VAIPHY_SIF) \
# 		python -W ignore \
# 			/opt/vaiphy/src/main.py \
# 			--data_path /mnt/ \
# 			--result_path /mnt/ \
# 			--max_iter 200 \
# 			--n_particles 128 \
# 		--dataset tree.$* > tree.$*.vaiphy.log 2> tree.$*.vaiphy.err
# 	echo done > $@

tree.%.nj.mod: tree.%.nj.nwk $(FA_DIR)/%.fa
	cp base-hky.mod tmp.mod
	echo -n "TREE: " >> tmp.mod
	cat $< >> tmp.mod
	$(PHAST_BIN)/phyloFit --lnl --init-model tmp.mod -o tree.$*.nj $(FA_DIR)/$*.fa
	rm tmp.mod

tree.%.modlnl: tree.%.nj.mod tree.%.ml.mod 
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
	grep -v '^\[' $< | grep -v '^Gen' | awk '{print $$2}' | sort -gr | head -1 | awk '{printf "%.6f\n", $$1}' >> $@

tree.%.raxmllnl: tree.%.raxml.term
	echo -n "$^ " > $@
	grep 'Final LogLikelihood:' $^ | awk '{printf "%.6f\n", $$3}' >> $@

tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.mrbayeslnl tree.%.raxmllnl 
	cat $^ | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@

eval.all.lnl.txt: $(LNL)
	echo -e "ds\tnj\tml\tvine\tbeast\tmrbayes\traxml" > tmp
	for file in $^ ; do \
		printf "%s\t" `echo $${file} | awk -F\. '{print $$2}'` >> tmp ;\
		awk '{printf "%s\t", $$2}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{x1 += $$2; x2 += $$3; x3 += $$4; x4 += $$5; x5 += $$6; x6 += $$7; print $$0} END {printf ("-----\nave\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), x2/(NR-1), x3/(NR-1), x4/(NR-1), x5/(NR-1), x6/(NR-1))}' tmp > $@
	rm tmp

# extract timing info
tree.%.time: tree.%.beast.term tree.%.var-time tree.%.mrbayes.term tree.%.raxml.term
	echo -e "ds\tbeast\tvine\tmrbayes\traxml" > $@; \
	beast_time=$$(grep '^Total calculation time' tree.$*.beast.term | awk '{print $$4}'); \
	printf "$*\t%s\t" "$$beast_time" >> $@; \
	head -1 tree.$*.var-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@; \
	mrbayes_time=$$(grep 'Analysis used' tree.$*.mrbayes.term | awk '{printf "%s\t", $$(3)}'); \
	printf "%s\t" "$$mrbayes_time" >> $@; \
	grep 'Elapsed time:' tree.$*.raxml.term | awk '{printf "%s\n", $$3}' >> $@

eval.all.time.txt: $(TIME)
	awk 'FNR==1 && NR==1 {print; next} FNR==2 {print; for(i=2;i<=NF;i++) sum[i]+=$$i; n++} END {if(n>0){printf "-----------------------------------------\nall"; for(i=2;i<=NF;i++) printf "\t%.2f", sum[i]/n; printf "\n"}}' $(TIME) > $@

# # evalTrees stuff
# # (1) modelFit
# # extract kappa from vine log
# tree.%.var.mf.txt: tree.%.var.nwk tree.%.heldout.fa tree.%.var.nwk.log
# 	kappa=`tail -1 tree.$*.var.nwk.log | awk '{print $$19}'` ;\
# 	$(PHAST_BIN)/evalTrees tree.$*.var.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# # use true kappa
# tree.%.true.mf.txt: tree.%.true.nwk tree.%.heldout.fa
# 	$(PHAST_BIN)/evalTrees tree.$*.true.nwk -f tree.$*.heldout.fa -k 4 > $@

# # use ML kappa
# tree.%.nj.mf.txt: tree.%.nj.nwk tree.%.heldout.fa tree.%.ml.mod
# 	kappa=`awk '$$1 == "BACKGROUND:" {pi_c = $$3; pi_g = $$4 } $$1<0 {print ($$3/pi_g) / ($$2/pi_c)}' tree.$*.ml.mod` ;\
# 	$(PHAST_BIN)/evalTrees tree.$*.nj.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# # use ML kappa
# tree.%.ml.mf.txt: tree.%.ml.nwk tree.%.heldout.fa tree.%.ml.mod
# 	kappa=`awk '$$1 == "BACKGROUND:" {pi_c = $$3; pi_g = $$4 } $$1<0 {print ($$3/pi_g) / ($$2/pi_c)}' tree.$*.ml.mod` ;\
# 	$(PHAST_BIN)/evalTrees tree.$*.ml.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# # use posterior mean kappa from beast log
# tree.%.beast.mf.txt: tree.%.beast.nwk tree.%.heldout.fa tree.%.beast.log
# 	kappa=`awk '{if (inlog) {sum += $$10; n++} ; if ($$1 == "Sample") inlog=1} END {print sum/n}' tree.$*.beast.log` ;\
# 	$(PHAST_BIN)/evalTrees tree.$*.beast.nwk -f tree.$*.heldout.fa -k $$kappa > $@

# tree.%.mf: tree.%.true.mf.txt tree.%.nj.mf.txt tree.%.ml.mf.txt tree.%.var.mf.txt tree.%.beast.mf.txt
# 	rm -f $@
# 	for file in $^ ; do \
# 		echo -n "$$file     " >> $@ ;\
# 		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
# 	done

# # (2) distances
# tree.%.true.dist.txt: tree.%.true.nwk 
# 	$(PHAST_BIN)/evalTrees tree.$*.true.nwk > $@

# tree.%.var.dist.txt: tree.%.var.nwk tree.%.true.dist.txt
# 	$(PHAST_BIN)/evalTrees tree.$*.var.nwk > tmp
# 	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

# tree.%.nj.dist.txt: tree.%.nj.nwk  tree.%.true.dist.txt
# 	$(PHAST_BIN)/evalTrees tree.$*.nj.nwk > tmp
# 	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

# tree.%.ml.dist.txt: tree.%.ml.nwk  tree.%.true.dist.txt
# 	$(PHAST_BIN)/evalTrees tree.$*.ml.nwk > tmp
# 	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

# tree.%.beast.dist.txt: tree.%.beast.nwk  tree.%.true.dist.txt
# 	$(PHAST_BIN)/evalTrees tree.$*.beast.nwk > tmp
# 	python3 "$(ROOT)/sumDists.py" tmp tree.$*.true.dist.txt | grep -v '^#' > $@

# tree.%.dist: tree.%.ml.dist.txt tree.%.var.dist.txt tree.%.beast.dist.txt
# 	paste $^ > $@

# eval.all.dist.txt: $(EVALDIST)
# 	cat $^ | awk '{print $$1, $$2, $$6, $$7, $$8, $$9, $$11, $$12, $$13, $$14}' | awk 'BEGIN {printf "ML_r2 ML_RMSE vine_r2 vine_RMSE vine_95CI vine_50CI beast_r2 beast_RMSE beast_95CI beast_50CI\n"} {x1 += $$1; x2 += $$2; x3 += $$3; x4 += $$4; x5 += $$5; x6 += $$6; x7 += $$7; x8 += $$8; x9 += $$9; x10 += $$10; print $$0} END{printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/NR, x2/NR, x3/NR, x4/NR, x5/NR, x6/NR, x7/NR, x8/NR, x9/NR, x10/NR}' > $@

# # (3) RF dist
# tree.%.var.rf.txt: tree.%.var.nwk tree.%.true.nwk
# 	$(PHAST_BIN)/evalTrees tree.$*.var.nwk -t tree.$*.true.nwk > $@

# tree.%.true.rf.txt: tree.%.true.nwk tree.%.true.nwk
# 	$(PHAST_BIN)/evalTrees tree.$*.true.nwk -t tree.$*.true.nwk > $@

# tree.%.nj.rf.txt: tree.%.nj.nwk tree.%.true.nwk
# 	$(PHAST_BIN)/evalTrees tree.$*.nj.nwk -t tree.$*.true.nwk > $@

# tree.%.ml.rf.txt: tree.%.ml.nwk tree.%.true.nwk
# 	$(PHAST_BIN)/evalTrees tree.$*.ml.nwk -t tree.$*.true.nwk > $@

# tree.%.beast.rf.txt: tree.%.beast.nwk tree.%.true.nwk
# 	$(PHAST_BIN)/evalTrees tree.$*.beast.nwk -t tree.$*.true.nwk > $@

# tree.%.rf: tree.%.true.rf.txt tree.%.nj.rf.txt tree.%.ml.rf.txt tree.%.var.rf.txt tree.%.beast.rf.txt
# 	rm -f $@
# 	for file in $^ ; do \
# 		echo -n "$$file     " >> $@ ;\
# 		awk '$$1 == "Mean:" {printf "%f\t", $$2} $$1 == "Std:" {printf "%f\n", $$2}' $${file} >> $@ ;\
# 	done

# eval.all.rf.txt: $(EVALRF)
# 	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd)" > tmp
# 	for file in $^ ; do \
# 		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmp ;\
# 		echo >> tmp ;\
# 	done
# 	awk '{x1 += $$1; x1s += ($$2 * $$2); x2 += $$3; x2s += ($$4 * $$4); x3 += $$5; x3s += ($$6 * $$6); x4 += $$7; x4s += ($$8 * $$8); x5 += $$9; x5s += ($$10 * $$10); print $$0} END {printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), x3/(NR-1), sqrt(x3s/(NR-1)), x4/(NR-1), sqrt(x4s/(NR-1)), x5/(NR-1), sqrt(x5s/(NR-1))}' tmp > $@
# 	rm -f tmp

# eval.all.mf.txt: $(EVALMF)
# 	echo "true (sd) nj (sd) ml (sd) vine (sd) beast (sd)" > tmp
# 	for file in $^ ; do \
# 		awk '{printf "%s\t%s\t", $$2, $$3}' $${file} >> tmp ;\
# 		echo >> tmp ;\
# 	done
# 	awk '{x1 += $$1; x1s += ($$2 * $$2); if ($$3 != "nan") x2 += $$3; x2s += ($$4 * $$4); if ($$5 != "nan") x3 += $$5; x3s += ($$6 * $$6); x4 += $$7; x4s += ($$8 * $$8); x5 += $$9; x5s += ($$10 * $$10); print $$0} END {printf "-----\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), sqrt(x1s/(NR-1)), x2/(NR-1), sqrt(x2s/(NR-1)), x3/(NR-1), sqrt(x3s/(NR-1)), x4/(NR-1), sqrt(x4s/(NR-1)), x5/(NR-1), sqrt(x5s/(NR-1))}' tmp > $@
# 	rm -f tmp

# # for use in debugging
# tracer: $(TRACER)

# tree.%.tr: tree.%.var.nwk.log
# 	grep -v '^#' $^ > $@

clean:
	rm tree.*.* eval.all.*

archive_all:
	archive_dir="archive.all_$$(date +%Y-%m-%d_%H:%M:%S)"; \
	mkdir $$archive_dir; \
	find . -maxdepth 1 ! -name . ! -name "Makefil*" ! -name "archive*" !-name "base-hky.mod" -exec mv -t "$$archive_dir" {} +
