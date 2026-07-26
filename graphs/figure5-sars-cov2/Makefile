# Figure 5 (SARS-CoV-2) — figure generation
#
# This directory generates the five Figure 5 panels ONLY.  The data analysis
# (VINE and BEAST 2 inference, raw posteriors) lives in the analysis directory
# ../../dna_real_data/sars-cov2 (= SRCDIR); nothing large is duplicated here.
# Raw posteriors are symlinked in from SRCDIR at build time and the small
# derived summaries the plots need (MCC trees, pairwise-distance tables) are
# built locally.  Big intermediates (full .nex posteriors) are marked .SECONDARY
# so they are rebuilt on demand and never force a re-render once the panels are
# current.
#
# PARAMETERIZATION (locked in 2026-07-24): Figure 5 uses the DEFAULT VINE run
# for every panel.  Its point estimates match BEAST tightly, whereas the
# calibrated -v 1 run trades a little scatter accuracy for wider, better-
# calibrated credible intervals and is reserved for the posterior-uncertainty
# supplement (built under ../../dna_real_data/sars-cov2/graphs/uncertainty).
# See README.md.  The -v 1 panel variants are kept in archive.flows-v1-panels/.

.SECONDARY:
.PRECIOUS:

MAIN_DIR   = ../..
BIN        = $(MAIN_DIR)/bin
VINE_BIN   = $(BIN)/vine/bin
SRCDIR     = $(MAIN_DIR)/dna_real_data/sars-cov2
PYTHON_SRC = $(MAIN_DIR)/python/src

# R with the phylogenetics stack (ape, ggtree, treeio, phangorn, phytools,
# ggrastr, ggpointdensity, uwot).  FONTCONFIG_PATH points cairo at the system
# fonts so the default `sans` family resolves to Nimbus Sans (the URW Helvetica
# clone) rather than the conda env's DejaVu; see README.md.
FONTENV = FONTCONFIG_PATH=/etc/fonts FONTCONFIG_FILE=/etc/fonts/fonts.conf
RSCRIPT ?= $(FONTENV) /home/asiepel/miniforge3/envs/rfig5/bin/Rscript

TREEANN = $(BIN)/beast/bin/treeannotator -lowMem true -height keep -burnin 0

.PHONY: all graphs clean-products clean-intermediates

# ---- the five panels ------------------------------------------------------
all graphs: vine-vs-beast.densitree.pdf tanglegram.compact.pdf scatter.pdf \
            vine-large.densitree.pdf embedding-pca.pdf

vine-vs-beast.densitree.pdf: densitree-compare.R \
        vine-small.nex vine-small.mcc.nex beast-small.nex beast-small.mcc.nex \
        small-subset.tsv
	$(RSCRIPT) $<

# tanglegram.R emits the compact and the full version in one pass.
tanglegram.compact.pdf vine-vs-beast.tanglegram.pdf: tanglegram.R \
        vine-small.mcc.nex beast-small.mcc.nex small-subset.tsv
	$(RSCRIPT) $<

scatter.pdf: scatter-plot.R vine-small.dist beast-small.dist
	$(RSCRIPT) $<

vine-large.densitree.pdf: densitree-large.R \
        vine-large.nex vine-large.mcc.nex large-subset.tsv
	$(RSCRIPT) $<

embedding-pca.pdf: embedding-pca.R vine-large.emb.tsv large-subset.tsv
	$(RSCRIPT) $<

# ---- raw inputs: symlinked from the analysis dir (never copied) -----------
vine-small.nwk vine-large.nwk vine-small.flows-v1.nwk beast-small-tree.trees \
small-subset.tsv large-subset.tsv vine-large.emb.tsv:
	ln -sf $(SRCDIR)/$@ $@

# ---- derived summaries (built locally) ------------------------------------
# BEAST: raw time-tree posterior -> substitutions, burn-in removed, rates
# stripped (.nex) -> newick -> pairwise-distance table.
beast-%.nex: beast-%-tree.trees
	python3 "$(PYTHON_SRC)/time2subs.py" --skip 2000 $< tmp.$*.nex
	sed 's/\[\&rate=[.0-9]*\]//g' tmp.$*.nex > $@
	rm -f tmp.$*.nex

beast-%.nwk: beast-%.nex
	$(BIN)/convertTrees -i nexus $< > $@

# VINE: newick posterior -> nexus.
vine-%.nex: vine-%.nwk
	$(BIN)/convertTrees -o nexus $< > $@

# Maximum-clade-credibility trees for the tanglegram / densitrees.
%.mcc.nex: %.nex
	$(TREEANN) $< > $@

# Pairwise-distance tables (mean + 95% CI) for the scatter.
%.dist: %.nwk
	$(VINE_BIN)/evalTrees $< > $@

# ---- housekeeping ---------------------------------------------------------
# Remove rendered panels (force a re-render).
clean-products:
	rm -f *.pdf *.png

# Remove the big rebuildable intermediates but keep the small cached summaries
# (*.mcc.nex, *.dist) and the rendered panels.
clean-intermediates:
	rm -f beast-small.nex beast-small.nwk vine-small.nex vine-large.nex
