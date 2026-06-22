.SECONDARY:
.PRECIOUS:

SRCDIR=~/sars-cov2
GRAPHSDIR=$(CURDIR)
PYTHON_SRC=./python

vine-%.nwk: ${SRCDIR}/vine-%.nwk
	cp -p $< $@	

vine-%.nex: vine-%.nwk
	convertTrees.py -o nexus $< > $@

beast-%.nex: ${SRCDIR}/beast-%-tree.trees
	python3 "$(PYTHON_SRC)/time2subs.py" --skip 2000 $^ tmp.nex
	sed 's/\[\&rate=[.0-9]*\]//g' tmp.nex > $@
	rm -r tmp.nex

beast-%.nwk: beast-%.nex
	convertTrees.py -i nexus $^ > $@

vine-%.mcc.nex: vine-%.nex
	treeannotator -lowMem true -height keep -burnin 0 $^ > $@

beast-%.mcc.nex: beast-%.nex
	treeannotator -lowMem true -height keep -burnin 0 $^ > $@

%.dist: %.nwk
	evalTrees $^ > $@

# ---- figure generation ----

graphs: vine-vs-beast.densitree.pdf vine-large.densitree.pdf \
        tanglegram.compact.pdf \
        embedding-pca.pdf scatter.pdf

vine-vs-beast.densitree.pdf: densitree-compare.R \
        vine-small.nex vine-small.mcc.nex beast-small.nex beast-small.mcc.nex \
        ${SRCDIR}/small-subset.tsv
	Rscript $< "$(GRAPHSDIR)" "$(SRCDIR)"

vine-large.densitree.pdf: densitree-large.R \
        vine-large.nex vine-large.mcc.nex \
        ${SRCDIR}/large-subset.tsv
	Rscript $< "$(GRAPHSDIR)" "$(SRCDIR)"

tanglegram.compact.pdf: tanglegram.R \
        vine-small.mcc.nex beast-small.mcc.nex \
        ${SRCDIR}/small-subset.tsv
	Rscript $< "$(GRAPHSDIR)" "$(SRCDIR)"

embedding-pca.pdf: embedding-pca.R ${SRCDIR}/vine-large.emb.tsv ${SRCDIR}/large-subset.tsv
	Rscript $< "$(SRCDIR)/vine-large.emb.tsv" "$(SRCDIR)/large-subset.tsv"

scatter.pdf: scatter-plot.R vine-small.dist beast-small.dist
	Rscript $< "$(GRAPHSDIR)/vine-small.dist" "$(GRAPHSDIR)/beast-small.dist"
