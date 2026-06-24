export SHELL=/usr/bin/bash

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
VAROPT?=-M 400 -c 100 --hky85 -s $(VINE_NSAMP) -v 0 


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

# clean up vine files for rerun
vineclean:
	rm -f $(VAR) $(LNL) $(TIME) eval.all.lnl.txt eval.all.time.txt *.varlnl *.var.nwk.log

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
	tail -1 $^ | awk '{print $$19}' >> $@ 

tree.%.beastlnl: tree.%.beast.log
	echo -n "$^ " > $@
	grep -v '^#' $^ | grep -v '^Sample' | awk '{print $$3}' | sort -nr | head -1 >> $@

tree.%.mrbayeslnl: tree.%.mrbayes.nex.p tree.%.mrbayes.term
	echo -n "$< " > $@
	grep -v '^\[' $< | grep -v '^Gen' | awk '{print $$2}' | sort -gr | head -1 | awk '{printf "%.6f\n", $$1}' >> $@

tree.%.lnl: tree.%.modlnl tree.%.varlnl tree.%.beastlnl tree.%.mrbayeslnl 
	cat $^ | awk '{if (true == 0) true = $$2; printf "%s %f\n", $$0, $$2 - true}' > $@

eval.all.lnl.txt: $(LNL)
	echo -e "ds\tnj\tml\tvine\tbeast\tmrbayes" > tmp
	for file in $^ ; do \
		printf "%s\t" `echo $${file} | awk -F\. '{print $$2}'` >> tmp ;\
		awk '{printf "%s\t", $$2}' $${file} >> tmp ;\
		echo >> tmp ;\
	done
	awk '{x1 += $$2; x2 += $$3; x3 += $$4; x4 += $$5; x5 += $$6; print $$0} END {printf ("-----\nave\t%f\t%f\t%f\t%f\t%f\n", x1/(NR-1), x2/(NR-1), x3/(NR-1), x4/(NR-1), x5/(NR-1))}' tmp > $@
	rm tmp

# extract timing info
tree.%.time: tree.%.beast.term tree.%.var-time tree.%.mrbayes.term
	echo -e "ds\tbeast\tvine\tmrbayes" > $@; \
	beast_time=$$(grep '^Total calculation time' tree.$*.beast.term | awk '{print $$4}'); \
	printf "$*\t%s\t" "$$beast_time" >> $@; \
	head -1 tree.$*.var-time | awk '{printf "%s\t", $$1}' | sed 's/user//' >> $@; \
	mrbayes_time=$$(grep 'Analysis used' tree.$*.mrbayes.term | awk '{printf "%s\t", $$(3)}'); \
	printf "%s\t" "$$mrbayes_time" >> $@; \

eval.all.time.txt: $(TIME)
	awk 'FNR==1 && NR==1 {print; next} FNR==2 {print; for(i=2;i<=NF;i++) sum[i]+=$$i; n++} END {if(n>0){printf "-----------------------------------------\nall"; for(i=2;i<=NF;i++) printf "\t%.2f", sum[i]/n; printf "\n"}}' $(TIME) > $@

clean:
	rm tree.*.* eval.all.*

archive_all:
	archive_dir="archive.all_$$(date +%Y-%m-%d_%H:%M:%S)"; \
	mkdir $$archive_dir; \
	find . -maxdepth 1 ! -name . ! -name "Makefil*" ! -name "archive*" !-name "base-hky.mod" -exec mv -t "$$archive_dir" {} +
