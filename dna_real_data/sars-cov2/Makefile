# Prevent Make from removing files it thinks it should clean up
#.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

NUM_MCMC_SAMPLES := 10000
MCMC_PRINTFREQ := 1000
BEAST_MCMCLEN := 100000000

MAIN_DIR := /local/storage/no-backup/vine-benchmarks
BIN := $(MAIN_DIR)/bin
VINE_BIN := $(BIN)/vine/bin
VAROPT?=--gtr --dgamma 4 -s 1000 -j 8
#VAROPT?=--hky -s 1000 -j 8 
PYTHON_SRC := $(MAIN_DIR)/python/src
BEAST := $(BIN)/beast/bin/beast
BEAST_TEMPLATE := $(MAIN_DIR)/dna_real_data/sars-cov2/beast_template.xml

vine: vine-small.nwk vine-large.nwk

vine-%.nwk: %-subset.fasta
	/usr/bin/time -o vine-$*.time $(VINE_BIN)/vine $(VAROPT) $< --logf vine-$*.log --mean vine-$*.mean.nwk --embedding vine-$*.emb.tsv > $@

vineclean:
	rm -f vine-*.nwk vine-*.log

beast: beast-small.nwk beast-large.nwk

beast-%.xml: 
	cp "$(BEAST_TEMPLATE)" $@

beast-%.term beast-%-tree.trees beast-%.log: beast-%.xml %-subset.fasta
	"$(BEAST)" -beagle -threads 8 -working -resume -seed 1772572215427 -D fastapath=$*-subset.fasta -D filebase=beast-$* -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$$(( $(BEAST_MCMCLEN) / $(NUM_MCMC_SAMPLES) )) -D printfreq=$(MCMC_PRINTFREQ) $< >> beast-$*.term 2>&1

beast-%.nwk: beast-%-tree.trees
	python3 "$(PYTHON_SRC)/time2subs.py" $^ tmp.nex
	$(BIN)/convertTrees -i nexus tmp.nex > $@
	rm -f tmp.nex

beastclean:
	rm -f beast-*.nwk beast-*.log beast-*.term beast-*-tree.trees beast-*.xml 
