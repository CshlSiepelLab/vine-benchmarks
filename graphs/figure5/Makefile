.SECONDARY:
.PRECIOUS:

MAIN_DIR=../..
BIN=$(MAIN_DIR)/bin
VINE_BIN=$(BIN)/vine/bin
SRCDIR=$(MAIN_DIR)/dna_real_data/sars-cov2
PYTHON_SRC=$(MAIN_DIR)/python/src

vine-small.nwk: ${SRCDIR}/vine-small.nwk
	cp -p $< $@

vine-large.nwk: ${SRCDIR}/vine-large.nwk
	cp -p $< $@

vine-%.nex: vine-%.nwk
	$(BIN)/convertTrees -o nexus $< > $@

beast-small-tree.trees: ${SRCDIR}/beast-small-tree.trees
	cp -p $< $@

beast-large-tree.trees: ${SRCDIR}/beast-large-tree.trees
	cp -p $< $@

beast-%.nex: beast-%-tree.trees
	python3 "./python/time2subs.py" --skip 2000 $^ tmp.nex
	sed 's/\[\&rate=[.0-9]*\]//g' tmp.nex > $@
	rm -r tmp.nex

beast-%.nwk: beast-%.nex
	$(BIN)/convertTrees -i nexus $^ > $@

vine-%.mcc.nex: vine-%.nex
	$(BIN)/beast/bin/treeannotator -lowMem true -height keep -burnin 0 $^ > $@

beast-%.mcc.nex: beast-%.nex
	$(BIN)/beast/bin/treeannotator -lowMem true -height keep -burnin 0 $^ > $@

small-subset.tsv: ${SRCDIR}/small-subset.tsv
	cp -p $< $@

large-subset.tsv: ${SRCDIR}/large-subset.tsv
	cp -p $< $@

vine-large.emb.tsv: ${SRCDIR}/vine-large.emb.tsv
	cp -p $< $@

%.dist: %.nwk
	$(VINE_BIN)/evalTrees $^ > $@

# ---- figure generation ----

graphs: vine-vs-beast.densitree.pdf vine-large.densitree.pdf \
        tanglegram.compact.pdf \
        embedding-pca.pdf scatter.pdf

vine-vs-beast.densitree.pdf: densitree-compare.R \
        vine-small.nex vine-small.mcc.nex beast-small.nex beast-small.mcc.nex \
        small-subset.tsv
	Rscript $<

vine-large.densitree.pdf: densitree-large.R \
        vine-large.nex vine-large.mcc.nex \
        large-subset.tsv
	Rscript $<

tanglegram.compact.pdf: tanglegram.R \
        vine-small.mcc.nex beast-small.mcc.nex \
        small-subset.tsv
	Rscript $<

embedding-pca.pdf: embedding-pca.R vine-large.emb.tsv large-subset.tsv
	Rscript $<

scatter.pdf: scatter-plot.R vine-small.dist beast-small.dist
	Rscript $<
