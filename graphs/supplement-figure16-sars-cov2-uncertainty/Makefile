# SARS-CoV-2 posterior uncertainty: Vine --posterior vs BEAST 2 on the 364-taxon subset.
#
# Inputs come from the Figure 5 folder, which builds them from the analysis
# runs in ../../dna_real_data/sars-cov2 (vine_twostage-small.nwk = VINE with
# --posterior; beast-small-tree.trees converted to substitutions/site with the
# first 2000 of 10001 samples dropped):
#   beast-small.nwk, vine_twostage-small.nwk   posterior trees (panel A)
# Missing inputs are built there with `make -C ../figure5-sars-cov2 <file>`.
# Panel B's pairwise-distance summaries are made here, from the same trees with
# branch lengths x1000 (substitutions per 1000 sites): evalTrees prints 4
# decimals, which rounds many SARS-CoV-2 distances (~1e-4 per site) too
# coarsely; the relative CI width does not depend on the units.
#
# Usage: make            -> plots/sars_cov2_posterior_uncertainty.pdf
# RSCRIPT defaults to the rplotting container (needs ggplot2, patchwork, MASS).

MAIN_DIR := $(abspath ../..)
FIG5     := ../figure5-sars-cov2
FASTA    := $(MAIN_DIR)/dna_real_data/sars-cov2/small-subset.fasta
RSCRIPT  ?= apptainer exec -B $(MAIN_DIR) $(MAIN_DIR)/containers/rplotting/rplotting.sif Rscript

VINE_BIN := $(MAIN_DIR)/bin/vine/bin
INPUTS := $(FIG5)/beast-small.nwk $(FIG5)/vine_twostage-small.nwk

.PHONY: all
all: plots/sars_cov2_posterior_uncertainty.pdf

$(INPUTS):
	$(MAKE) -C $(FIG5) $(notdir $@)

# split posterior probabilities; branch lengths are per site, so supporting
# substitutions need the alignment length
splits.tsv: splits.py $(FIG5)/beast-small.nwk $(FIG5)/vine_twostage-small.nwk $(FASTA)
	sites=$$(awk '/^>/ {if (n++) exit; next} {L += length($$0)} END {print L}' $(FASTA)); \
	python3 splits.py --reference $(FIG5)/beast-small.nwk --estimate $(FIG5)/vine_twostage-small.nwk \
	  --sites $$sites > $@.tmp && mv $@.tmp $@

# pairwise-distance summaries, branch lengths x1000
%-small.x1000.dist: $(FIG5)/%-small.nwk
	perl -pe 's/:([-+0-9.eE]+)/":" . ($$1 * 1000)/ge' $< > $*-small.x1000.nwk.tmp
	$(VINE_BIN)/evalTrees $*-small.x1000.nwk.tmp > $@.tmp && mv $@.tmp $@
	rm -f $*-small.x1000.nwk.tmp

plots/sars_cov2_posterior_uncertainty.pdf: makeGraphs.R splits.tsv beast-small.x1000.dist vine_twostage-small.x1000.dist
	$(RSCRIPT) makeGraphs.R
