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
# Posterior-calibration parameterization: DIST covariance, both normalizing
# flows (radial + planar), and inflated variance regularization.  Outputs are
# written alongside the default run with a ".flows" infix.
FLOWBASE?=$(VAROPT) -S DIST -F -Z
FLOWOPT?=$(FLOWBASE) -v 3
# Variance-regularization sweep.  Note -v 0 disables the regularizer, so that
# run doubles as the ablation isolating "DIST covariance + flows" from "-v".
# (vine-small.flows.nwk is the -v 3 member of this family.)
VREG_VALUES?=0 1 2
PYTHON_SRC := $(MAIN_DIR)/python/src
BEAST := $(BIN)/beast/bin/beast
BEAST_TEMPLATE := $(MAIN_DIR)/dna_real_data/sars-cov2/beast_template.xml

.PHONY: vine vine-flows vine-sweep archive-vine rerun-vine vineclean beast beastclean

vine: vine-small.nwk vine-large.nwk

vine-%.nwk: %-subset.fasta
	/usr/bin/time -o vine-$*.time $(VINE_BIN)/vine $(VAROPT) $< --logf vine-$*.log --mean vine-$*.mean.nwk --embedding vine-$*.emb.tsv > $@

# Calibrated ("flows") run.  Only the 364-taxon subset is needed for the
# posterior-uncertainty supplementary figure; `make vine-large.flows.nwk`
# builds the 1060-taxon set too.
vine-flows: vine-small.flows.nwk

vine-%.flows.nwk: %-subset.fasta
	/usr/bin/time -o vine-$*.flows.time $(VINE_BIN)/vine $(FLOWOPT) $< --logf vine-$*.flows.log --mean vine-$*.flows.mean.nwk --embedding vine-$*.flows.emb.tsv > $@

# Variance-regularization sweep: one rule per value in VREG_VALUES, producing
# vine-<size>.flows-v<N>.{nwk,log,time,mean.nwk,emb.tsv}.  Run in parallel with
#   make -j3 vine-sweep
vine-sweep: $(foreach v,$(VREG_VALUES),vine-small.flows-v$(v).nwk)

define VINE_FLOWS_RULE
vine-%.flows-v$(1).nwk: %-subset.fasta
	/usr/bin/time -o vine-$$*.flows-v$(1).time $$(VINE_BIN)/vine $$(FLOWBASE) -v $(1) $$< --logf vine-$$*.flows-v$(1).log --mean vine-$$*.flows-v$(1).mean.nwk --embedding vine-$$*.flows-v$(1).emb.tsv > $$@
endef
$(foreach v,$(VREG_VALUES),$(eval $(call VINE_FLOWS_RULE,$(v))))

VINE_RUN_FILES := $(foreach size,small large,$(foreach suffix,nwk log time mean.nwk emb.tsv,vine-$(size).$(suffix)))
VINE_ARCHIVE_DIR ?= archive.vine.$(shell date +%Y%m%d-%H%M%S)

archive-vine:
	@files="$(strip $(wildcard $(VINE_RUN_FILES)))"; \
	if [ -z "$$files" ]; then \
		echo "No Vine run outputs to archive"; \
	else \
		mkdir -p "$(VINE_ARCHIVE_DIR)"; \
		mv $$files "$(VINE_ARCHIVE_DIR)/"; \
		echo "Archived Vine run outputs in $(VINE_ARCHIVE_DIR)"; \
	fi

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
