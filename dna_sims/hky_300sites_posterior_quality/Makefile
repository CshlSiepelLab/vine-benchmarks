SHELL := /usr/bin/bash

# Keep each Make job single-threaded; use make -jN to control job-level
# parallelism without OpenBLAS or OpenMP multiplying the CPU usage per job.
export OPENBLAS_NUM_THREADS := 1
export OMP_NUM_THREADS := 1

.DELETE_ON_ERROR:
.SECONDARY:
.PRECIOUS:
.SECONDEXPANSION:

MAIN_DIR      ?= /local/storage/no-backup/vine-benchmarks
BASE_DATA_DIR ?= $(MAIN_DIR)/dna_sims/hky_300sites
PYTHON_SRC    := $(MAIN_DIR)/python/src

# Flag here to decide if to initialize vine with ml tree or not. 1 uses the ml tree.
USE_BASE_ML_TREE ?= 0

VINE_BIN1   ?= /home/staklins/projects/vine_project/vine/bin
VINE_BIN   ?= /home/staklins/projects/vine_project/vine/bin
# VINE_BIN ?= /local/storage/no-backup/vine-benchmarks/bin/vine/bin
PYTHON     ?= /home/staklins/miniforge3/envs/pyplotting/bin/python
R_SCRIPT   ?= /home/staklins/miniforge3/envs/rplotting/bin/Rscript
MPLCONFIGDIR ?= /tmp/matplotlib-vine
MDS_PLOT_SCRIPT ?= plot_mds_comparison.py
MDS_PLOT_METHODS := beast beast-beagle mrbayes mrbayes-beagle vine vine_flows bootstrap_nj
MDS_TAXA ?= 10
MDS_SAMPLE ?= 1
MDS_FILES ?= $(foreach m,$(MDS_PLOT_METHODS),$(MDS_TAXA)taxa/tree.$(MDS_SAMPLE).$(m).pairwise_rf.mds.tsv)
MDS_PLOT ?= plots/mds/$(MDS_TAXA)taxa/tree.$(MDS_SAMPLE).mds_comparison.pdf

TAXA       ?= 10 25 50 100 250 500 1000
NSAMP      ?= 10
VINE_NSAMP ?= 1000
BOOTSTRAP_NSAMP ?= 1000
BOOTSTRAP_NSITES ?= 300

SAMPLES := $(shell seq 1 $(NSAMP))
METHODS := bootstrap_nj vine vine_flows
# vine_mixture
# vine_mixture_flows

SHARED_ARGS := --hky85 -s $(VINE_NSAMP)
FLOWS_ARGS := -v 3 --cov DIST --planar-flow --radial-flow
MIXTURE_ARGS := --n-mixture-components 3

VINE_START_TREE_ARG = $(if $(filter 1 yes true on,$(USE_BASE_ML_TREE)),--tree $(BASE_DATA_DIR)/$*.ml.nwk)

VINE_ARGS := $(SHARED_ARGS) -v 0
VINE_FLOWS_ARGS := $(SHARED_ARGS) $(FLOWS_ARGS)
VINE_MIXTURE_ARGS := $(SHARED_ARGS) $(MIXTURE_ARGS) -v 0
VINE_MIXTURE_FLOWS_ARGS := $(SHARED_ARGS) $(FLOWS_ARGS) $(MIXTURE_ARGS)

BASE_FASTA      := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).fa))
BASE_TRUE_NWK   := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).true.nwk))
BASE_TRUE_DIST  := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).true.dist.txt))
BASE_BEAST_DIST := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast.dist.txt))
BASE_BEAST_BEAGLE_DIST := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast-beagle.dist.txt))
BASE_MRBAYES_DIST := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes.dist.txt))
BASE_MRBAYES_BEAGLE_DIST := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes-beagle.dist.txt))
BASE_BEAST_RF   := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast.rf.txt))
BASE_BEAST_BEAGLE_RF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast-beagle.rf.txt))
BASE_MRBAYES_RF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes.rf.txt))
BASE_MRBAYES_BEAGLE_RF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes-beagle.rf.txt))
BASE_BEAST_ENT  := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast.ent.txt))
BASE_BEAST_BEAGLE_ENT := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast-beagle.ent.txt))
BASE_MRBAYES_ENT := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes.ent.txt))
BASE_MRBAYES_BEAGLE_ENT := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes-beagle.ent.txt))
BASE_TRUE_MOD := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).true.mod))
BASE_HELDOUT_FA := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).heldout.fa))
BASE_TRUE_MF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).true.mf.txt))
BASE_BEAST_MF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast.mf.txt))
BASE_BEAST_BEAGLE_MF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast-beagle.mf.txt))
BASE_MRBAYES_MF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes.mf.txt))
BASE_MRBAYES_BEAGLE_MF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes-beagle.mf.txt))
BASE_BEAST_LNL  := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beastlnl))
BASE_BEAST_BEAGLE_LNL := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast-beaglelnl))
BASE_MRBAYES_LNL := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayeslnl))
BASE_MRBAYES_BEAGLE_LNL := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes-beaglelnl))
BASE_BEAST_TERM := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast.term))
BASE_BEAST_BEAGLE_TERM := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast-beagle.term))
BASE_MRBAYES_TERM := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes.term))
BASE_MRBAYES_BEAGLE_TERM := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes-beagle.term))
BASE_BEAST_NWK := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast.nwk))
BASE_BEAST_BEAGLE_NWK := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).beast-beagle.nwk))
BASE_MRBAYES_NWK := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes.nwk))
BASE_MRBAYES_BEAGLE_NWK := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).mrbayes-beagle.nwk))
BASE_INPUTS := $(BASE_FASTA) $(BASE_TRUE_NWK) $(BASE_TRUE_DIST) \
	$(BASE_BEAST_DIST) $(BASE_BEAST_BEAGLE_DIST) $(BASE_MRBAYES_DIST) $(BASE_MRBAYES_BEAGLE_DIST) \
	$(BASE_BEAST_RF) $(BASE_BEAST_BEAGLE_RF) $(BASE_MRBAYES_RF) $(BASE_MRBAYES_BEAGLE_RF) \
	$(BASE_BEAST_ENT) $(BASE_BEAST_BEAGLE_ENT) $(BASE_MRBAYES_ENT) $(BASE_MRBAYES_BEAGLE_ENT) \
	$(BASE_TRUE_MOD) $(BASE_BEAST_LNL) $(BASE_BEAST_BEAGLE_LNL) $(BASE_MRBAYES_LNL) $(BASE_MRBAYES_BEAGLE_LNL) \
	$(BASE_BEAST_TERM) $(BASE_BEAST_BEAGLE_TERM) $(BASE_MRBAYES_TERM) $(BASE_MRBAYES_BEAGLE_TERM) \
	$(BASE_BEAST_NWK) $(BASE_BEAST_BEAGLE_NWK) $(BASE_MRBAYES_NWK) $(BASE_MRBAYES_BEAGLE_NWK) \
	$(BASE_HELDOUT_FA) $(BASE_TRUE_MF) $(BASE_BEAST_MF) $(BASE_BEAST_BEAGLE_MF) $(BASE_MRBAYES_MF) $(BASE_MRBAYES_BEAGLE_MF)

VINE_NWK := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(n)taxa/tree.$(i).vine.nwk))
BOOTSTRAP_NJ_NWK := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(n)taxa/tree.$(i).bootstrap_nj.nwk))
VINE_FLOWS_NWK := $(VINE_NWK:vine.nwk=vine_flows.nwk)
VINE_MIXTURE_NWK := $(VINE_NWK:vine.nwk=vine_mixture.nwk)
VINE_MIXTURE_FLOWS_NWK := $(VINE_NWK:vine.nwk=vine_mixture_flows.nwk)
METHOD_NWK_vine := $(VINE_NWK)
METHOD_NWK_bootstrap_nj := $(BOOTSTRAP_NJ_NWK)
METHOD_NWK_vine_flows := $(VINE_FLOWS_NWK)
METHOD_NWK_vine_mixture := $(VINE_MIXTURE_NWK)
METHOD_NWK_vine_mixture_flows := $(VINE_MIXTURE_FLOWS_NWK)
NWK := $(foreach m,$(METHODS),$(METHOD_NWK_$(m)))
VARIATIONAL_NWK := $(filter-out $(BOOTSTRAP_NJ_NWK),$(NWK))

# Pairwise RF/MDS evaluation uses existing posterior samples only.  Local methods
# are read from this directory; BEAST/MrBayes samples remain in BASE_DATA_DIR,
# while all result files are written here.
PAIRWISE_RF_METHODS := $(METHODS) beast beast-beagle mrbayes mrbayes-beagle
PAIRWISE_RF := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(foreach m,$(PAIRWISE_RF_METHODS),$(n)taxa/tree.$(i).$(m).pairwise_rf.tsv)))
MDS_COORDS := $(PAIRWISE_RF:.pairwise_rf.tsv=.pairwise_rf.mds.tsv)
MDS_PLOTS := $(foreach n,$(TAXA),plots/mds/$(n)taxa.mds_comparison.pdf)

DIST_SUM := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.dist.txt))
RF_SUM   := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.rf.txt))
ENT_SUM  := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.ent.txt))
LNL_SUM  := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.lnl.txt))
TIME_SUM := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.time.txt))
MF_SUM   := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.mf.txt))
BSD_SUM  := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.bsd.txt))
KL_PDIST_BEAST_SUM    := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.kl_pdist_beast.txt))
KL_PDIST_MRBAYES_SUM  := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.kl_pdist_mrbayes.txt))
KL_PDIST_BEAST_BEAGLE_SUM   := $(foreach n,$(TAXA),$(n)taxa/beast-beagle.eval.all.kl_pdist_beast.txt)
KL_PDIST_MRBAYES_BEAGLE_SUM := $(foreach n,$(TAXA),$(n)taxa/mrbayes-beagle.eval.all.kl_pdist_mrbayes.txt)
KL_PDIST_MRBAYES_VS_BEAST_SUM        := $(foreach n,$(TAXA),$(n)taxa/mrbayes.eval.all.kl_pdist_beast.txt)
KL_PDIST_MRBAYES_BEAGLE_VS_BEAST_SUM := $(foreach n,$(TAXA),$(n)taxa/mrbayes-beagle.eval.all.kl_pdist_beast.txt)
KL_PDIST_BEAST_VS_MRBAYES_SUM        := $(foreach n,$(TAXA),$(n)taxa/beast.eval.all.kl_pdist_mrbayes.txt)
KL_PDIST_BEAST_BEAGLE_VS_MRBAYES_SUM := $(foreach n,$(TAXA),$(n)taxa/beast-beagle.eval.all.kl_pdist_mrbayes.txt)
KL_TOPO_BEAST_SUM    := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.kl_topo_beast.txt))
KL_TOPO_MRBAYES_SUM  := $(foreach n,$(TAXA),$(foreach m,$(METHODS),$(n)taxa/$(m).eval.all.kl_topo_mrbayes.txt))
KL_TOPO_BEAST_BEAGLE_SUM   := $(foreach n,$(TAXA),$(n)taxa/beast-beagle.eval.all.kl_topo_beast.txt)
KL_TOPO_MRBAYES_BEAGLE_SUM := $(foreach n,$(TAXA),$(n)taxa/mrbayes-beagle.eval.all.kl_topo_mrbayes.txt)
KL_TOPO_MRBAYES_VS_BEAST_SUM        := $(foreach n,$(TAXA),$(n)taxa/mrbayes.eval.all.kl_topo_beast.txt)
KL_TOPO_MRBAYES_BEAGLE_VS_BEAST_SUM := $(foreach n,$(TAXA),$(n)taxa/mrbayes-beagle.eval.all.kl_topo_beast.txt)
KL_TOPO_BEAST_VS_MRBAYES_SUM        := $(foreach n,$(TAXA),$(n)taxa/beast.eval.all.kl_topo_mrbayes.txt)
KL_TOPO_BEAST_BEAGLE_VS_MRBAYES_SUM := $(foreach n,$(TAXA),$(n)taxa/beast-beagle.eval.all.kl_topo_mrbayes.txt)
SUMMARIES := $(DIST_SUM) $(RF_SUM) $(ENT_SUM) $(LNL_SUM) $(TIME_SUM) $(MF_SUM) $(BSD_SUM) \
	$(KL_PDIST_BEAST_SUM) $(KL_PDIST_MRBAYES_SUM) $(KL_PDIST_BEAST_BEAGLE_SUM) $(KL_PDIST_MRBAYES_BEAGLE_SUM) \
	$(KL_PDIST_MRBAYES_VS_BEAST_SUM) $(KL_PDIST_MRBAYES_BEAGLE_VS_BEAST_SUM) \
	$(KL_PDIST_BEAST_VS_MRBAYES_SUM) $(KL_PDIST_BEAST_BEAGLE_VS_MRBAYES_SUM) \
	$(KL_TOPO_BEAST_SUM) $(KL_TOPO_MRBAYES_SUM) $(KL_TOPO_BEAST_BEAGLE_SUM) $(KL_TOPO_MRBAYES_BEAGLE_SUM) \
	$(KL_TOPO_MRBAYES_VS_BEAST_SUM) $(KL_TOPO_MRBAYES_BEAGLE_VS_BEAST_SUM) \
	$(KL_TOPO_BEAST_VS_MRBAYES_SUM) $(KL_TOPO_BEAST_BEAGLE_VS_MRBAYES_SUM)

VINE_SUMMARIES := $(foreach n,$(TAXA),$(foreach s,dist rf ent lnl time mf bsd,$(n)taxa/vine.eval.all.$(s).txt))
VINE_FLOWS_SUMMARIES := $(subst /vine.eval,/vine_flows.eval,$(VINE_SUMMARIES))
VINE_MIXTURE_SUMMARIES := $(subst /vine.eval,/vine_mixture.eval,$(VINE_SUMMARIES))
VINE_MIXTURE_FLOWS_SUMMARIES := $(subst /vine.eval,/vine_mixture_flows.eval,$(VINE_SUMMARIES))
BOOTSTRAP_NJ_SUMMARIES := $(foreach s,$(SUMMARIES),$(if $(findstring /bootstrap_nj.eval.,$(s)),$(s)))

.PHONY: all check-base-data vines bootstrap_nj vine vine_flows vine_mixture vine_mixture_flows summaries pairwise_rf mds mds_plot mds_plot_one plot clean clean-vine clean-local-data clean-rf clean-mds clean-plots

all: vines summaries plot
vines: check-base-data $(NWK)
bootstrap_nj: check-base-data $(BOOTSTRAP_NJ_SUMMARIES)
vine: check-base-data $(VINE_SUMMARIES)
vine_flows: check-base-data $(VINE_FLOWS_SUMMARIES)
vine_mixture: check-base-data $(VINE_MIXTURE_SUMMARIES)
vine_mixture_flows: check-base-data $(VINE_MIXTURE_FLOWS_SUMMARIES)
summaries: check-base-data $(SUMMARIES)
pairwise_rf: $(PAIRWISE_RF)
mds: $(MDS_COORDS)
mds_plot: $(MDS_PLOTS)
mds_plot_one: $(MDS_PLOT_SCRIPT)
	@test -n "$(strip $(MDS_FILES))" || { \
	  printf 'Set MDS_FILES to one or more .mds.tsv paths\n' >&2; exit 2; }
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) $(MDS_PLOT_SCRIPT) \
	  --coordinates $(MDS_FILES) --out $(MDS_PLOT)
	@test -s $(MDS_PLOT)
plot: plots/hky300_summary.pdf \
	plots/hky300_posterior_quality_sup1.pdf plots/hky300_posterior_quality_sup2.pdf

check-base-data: $(BASE_INPUTS)
	@printf 'Using existing simulation, BEAST2, BEAST2+BEAGLE, MrBayes, and MrBayes+BEAGLE summaries from %s\n' "$(BASE_DATA_DIR)"

# Deliberately omit tree-file prerequisites here: a missing posterior should
# fail instead of allowing Make to invoke one of the inference rules above.
define PAIRWISE_RF_RULE
$(1).pairwise_rf.tsv:
	@test -s $(2) || { printf 'Missing existing posterior: %s\n' '$(2)' >&2; exit 1; }
	$(VINE_BIN)/evalTrees \
	  --rf-matrix $(1).pairwise_rf.tsv \
	  --rf-mds $(1).pairwise_rf.mds.tsv $(2)
	@test -s $(1).pairwise_rf.tsv
	@test -s $(1).pairwise_rf.mds.tsv

$(1).pairwise_rf.mds.tsv: $(1).pairwise_rf.tsv
	@if ! test -s $$@; then rm -f $$<; $(MAKE) $$<; fi
	@test -s $$@
endef

$(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(foreach m,$(METHODS),\
  $(eval $(call PAIRWISE_RF_RULE,$(n)taxa/tree.$(i).$(m),$(n)taxa/tree.$(i).$(m).nwk)))))
$(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(foreach m,beast beast-beagle mrbayes mrbayes-beagle,\
  $(eval $(call PAIRWISE_RF_RULE,$(n)taxa/tree.$(i).$(m),$(BASE_DATA_DIR)/$(n)taxa/tree.$(i).$(m).nwk)))))

define MDS_PLOT_RULE
plots/mds/$(1)taxa.mds_comparison.pdf: \
    $$(foreach i,$$(SAMPLES),$$(foreach m,$$(MDS_PLOT_METHODS),$(1)taxa/tree.$$(i).$$(m).pairwise_rf.mds.tsv)) \
    $$(MDS_PLOT_SCRIPT)
	MPLCONFIGDIR=$$(MPLCONFIGDIR) $$(PYTHON) $$(MDS_PLOT_SCRIPT) \
	  --coordinates $$(filter %.tsv,$$^) --out $$@
	@test -s $$@
endef

$(foreach n,$(TAXA),$(eval $(call MDS_PLOT_RULE,$(n))))

$(BOOTSTRAP_NJ_NWK): %.bootstrap_nj.nwk: $(BASE_DATA_DIR)/%.fa $(BASE_DATA_DIR)/%.ml.mod
	mkdir -p $(dir $@)
	tmp_fa=$$(mktemp /tmp/bootstrap_nj.XXXXXX.fa); times=$$(mktemp /tmp/bootstrap_nj.XXXXXX.times); \
	out=$@.tmp; bootstrap_time=$*.bootstrap_nj.time.tmp; nj_time=$*.nj.time.tmp; \
	trap 'rm -f "$$tmp_fa" "$$times" "$$out" "$$bootstrap_time" "$$nj_time"' EXIT; \
	rm -f "$$out"; \
	export input="$<" tmp_fa out nboot="$(BOOTSTRAP_NSAMP)" \
	  nsites="$(BOOTSTRAP_NSITES)" vine="$(VINE_BIN)/vine" times; \
	/usr/bin/bash -c ' \
	  set -e; \
	  for ((rep=1; rep<=nboot; rep++)); do \
	    seed=$$(( (RANDOM << 15) ^ RANDOM ^ rep )); \
	    awk -v seed="$$seed" -v nout="$$nsites" '\'' \
	      /^>/ {if (n) seq[n]=s; n++; header[n]=$$0; s=""; next} \
	           {s=s $$0} \
	      END { \
	        seq[n]=s; nsite=length(seq[1]); srand(seed); \
	        if (nout < 1) exit 2; \
	        for (j=1; j<=nout; j++) col[j]=int(rand()*nsite)+1; \
	        for (i=1; i<=n; i++) { \
	          print header[i]; \
	          for (j=1; j<=nout; j++) printf "%s", substr(seq[i],col[j],1); \
	          print ""; \
	        } \
	      }'\'' "$$input" > "$$tmp_fa"; \
	    start=$$(date +%s%N); \
	    "$$vine" --nj-only --hky85 "$$tmp_fa" >> "$$out" 2>/dev/null; \
	    end=$$(date +%s%N); \
	    awk -v start="$$start" -v end="$$end" '\''BEGIN {printf "%.9f\n", (end-start)/1e9}'\'' >> "$$times"; \
	  done'; \
	awk '{sum+=$$1} END {printf "%.9fuser\n", sum}' "$$times" > "$$bootstrap_time"; \
	awk 'NR==1 {printf "%.9fuser\n", $$1}' "$$times" > "$$nj_time"; \
	mv "$$out" $@; \
	mv "$$bootstrap_time" $*.bootstrap_nj.time; \
	mv "$$nj_time" $*.nj.time; \
	rm -f "$$tmp_fa" "$$times"; trap - EXIT
	kappa=$$(awk '$$1 == "BACKGROUND:" {pi_c=$$3; pi_g=$$4} $$1+0 < 0 {print ($$3/pi_g)/($$2/pi_c); exit}' $(word 2,$^)); \
	printf 'kappa: %s\n' "$$kappa" > $*.bootstrap_nj.nwk.log
	test "$$(wc -l < $@)" -eq "$(BOOTSTRAP_NSAMP)"

$(VINE_NWK): %.vine.nwk: $(BASE_DATA_DIR)/%.fa $$(if $$(filter 1 yes true on,$$(USE_BASE_ML_TREE)),$$(BASE_DATA_DIR)/$$*.ml.nwk)
	mkdir -p $(dir $@)
	rm -f $*.vine.nwk $*.vine.nwk.log $*.vine.time $*.vine.mean.nwk
	/usr/bin/time -o $*.vine.time \
	  $(VINE_BIN)/vine $< -l $*.vine.nwk.log $(VINE_ARGS) \
	  $(VINE_START_TREE_ARG) --mean $*.vine.mean.nwk > $@
	test -s $@
	test -s $*.vine.nwk.log

$(VINE_FLOWS_NWK): %.vine_flows.nwk: $(BASE_DATA_DIR)/%.fa $$(if $$(filter 1 yes true on,$$(USE_BASE_ML_TREE)),$$(BASE_DATA_DIR)/$$*.ml.nwk)
	mkdir -p $(dir $@)
	rm -f $*.vine_flows.nwk $*.vine_flows.nwk.log $*.vine_flows.time $*.vine_flows.mean.nwk
	/usr/bin/time -o $*.vine_flows.time \
	  $(VINE_BIN)/vine $< -l $*.vine_flows.nwk.log $(VINE_FLOWS_ARGS) \
	  $(VINE_START_TREE_ARG) --mean $*.vine_flows.mean.nwk > $@
	test -s $@
	test -s $*.vine_flows.nwk.log

$(VINE_MIXTURE_NWK): %.vine_mixture.nwk: $(BASE_DATA_DIR)/%.fa $$(if $$(filter 1 yes true on,$$(USE_BASE_ML_TREE)),$$(BASE_DATA_DIR)/$$*.ml.nwk)
	mkdir -p $(dir $@)
	rm -f $*.vine_mixture.nwk $*.vine_mixture.nwk.log $*.vine_mixture.time $*.vine_mixture.mean.nwk
	/usr/bin/time -o $*.vine_mixture.time \
	  $(VINE_BIN)/vine $< -l $*.vine_mixture.nwk.log $(VINE_MIXTURE_ARGS) \
	  $(VINE_START_TREE_ARG) --mean $*.vine_mixture.mean.nwk > $@
	test -s $@
	test -s $*.vine_mixture.nwk.log

$(VINE_MIXTURE_FLOWS_NWK): %.vine_mixture_flows.nwk: $(BASE_DATA_DIR)/%.fa $$(if $$(filter 1 yes true on,$$(USE_BASE_ML_TREE)),$$(BASE_DATA_DIR)/$$*.ml.nwk)
	mkdir -p $(dir $@)
	rm -f $*.vine_mixture_flows.nwk $*.vine_mixture_flows.nwk.log $*.vine_mixture_flows.time $*.vine_mixture_flows.mean.nwk
	/usr/bin/time -o $*.vine_mixture_flows.time \
	  $(VINE_BIN)/vine $< -l $*.vine_mixture_flows.nwk.log $(VINE_MIXTURE_FLOWS_ARGS) \
	  $(VINE_START_TREE_ARG) --mean $*.vine_mixture_flows.mean.nwk > $@
	test -s $@
	test -s $*.vine_mixture_flows.nwk.log

$(NWK:.nwk=.nwk.log): %.nwk.log: | %.nwk
	@test -s $@

$(NWK:.nwk=.time): %.time: | %.nwk
	@test -s $@

$(VARIATIONAL_NWK:.nwk=.varlnl): %.varlnl: %.nwk.log
	echo -n "$< " > $@
	tail -1 $< | awk '{print $$11}' | sed 's/,//g' >> $@ 
# 	printf '%s ' $< > $@
# 	tail -1 $< | awk '{ \
# 	  for (i = 1; i <= NF; i++) { \
# 	    if ($$i == "LNL:") lnl = $$(i+1); \
# 	    else if ($$i == "LNL_mc:") lnl_mc = $$(i+1); \
# 	  } \
# 	  v = (lnl_mc != "" ? lnl_mc : lnl); \
# 	  gsub(",", "", v); \
# 	  print v; \
	# 	}' >> $@

$(BOOTSTRAP_NJ_NWK:.nwk=.varlnl): %.varlnl: %.nwk %.nwk.log $$(BASE_DATA_DIR)/$$(basename $$*).fa
	kappa=$$(awk '$$1 == "kappa:" {print $$2}' $(word 2,$^)); \
	$(VINE_BIN)/evalTrees $< -f $(word 3,$^) -k "$$kappa" | \
	  awk -v file="$<" 'BEGIN {printf "%s ", file} $$1 == "Mean:" {print $$2; found=1} END {if (!found) exit 1}' > $@

$(NWK:.nwk=.lnl): %.lnl: %.varlnl \
    $$(BASE_DATA_DIR)/$$(basename $$*).true.mod \
    $$(BASE_DATA_DIR)/$$(basename $$*).beastlnl \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast-beaglelnl \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayeslnl \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes-beaglelnl
	{ echo -e "true\tvine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle"; \
	  vine_lnl=$$(awk '{print $$2}' $<); \
	  true_lnl=$$(grep LNL $(word 2,$^) | awk '{print $$2}'); \
	  beast_lnl=$$(awk '{print $$2}' $(word 3,$^)); \
	  beast_beagle_lnl=$$(awk '{print $$2}' $(word 4,$^)); \
	  mrbayes_lnl=$$(awk '{print $$2}' $(word 5,$^)); \
	  mrbayes_beagle_lnl=$$(awk '{print $$2}' $(word 6,$^)); \
	  awk -v t="$$true_lnl" -v v="$$vine_lnl" -v b="$$beast_lnl" -v bb="$$beast_beagle_lnl" -v m="$$mrbayes_lnl" -v mb="$$mrbayes_beagle_lnl" \
	      'BEGIN {printf "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", t, v, b, bb, m, mb}'; \
	} > $@

$(NWK:.nwk=.var.mf.txt): %.var.mf.txt: %.nwk %.nwk.log $$(BASE_DATA_DIR)/$$(basename $$*).heldout.fa
	kappa=$$(awk '{for (i = 1; i <= NF; i++) if ($$i == "kappa:") kappa=$$(i+1)} END {print kappa}' $*.nwk.log); \
	$(VINE_BIN)/evalTrees $*.nwk -f $(word 3,$^) -k "$$kappa" > $@.tmp && mv $@.tmp $@

$(NWK:.nwk=.mf): %.mf: %.var.mf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).true.mf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast.mf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast-beagle.mf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes.mf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes-beagle.mf.txt
	{ echo -e "true\tvine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle"; \
	  true_mf=$$(awk '$$1 == "Mean:" {print $$2}' $(word 2,$^)); \
	  vine_mf=$$(awk '$$1 == "Mean:" {print $$2}' $<); \
	  beast_mf=$$(awk '$$1 == "Mean:" {print $$2}' $(word 3,$^)); \
	  beast_beagle_mf=$$(awk '$$1 == "Mean:" {print $$2}' $(word 4,$^)); \
	  mrbayes_mf=$$(awk '$$1 == "Mean:" {print $$2}' $(word 5,$^)); \
	  mrbayes_beagle_mf=$$(awk '$$1 == "Mean:" {print $$2}' $(word 6,$^)); \
	  awk -v t="$$true_mf" -v v="$$vine_mf" -v b="$$beast_mf" -v bb="$$beast_beagle_mf" -v m="$$mrbayes_mf" -v mb="$$mrbayes_beagle_mf" \
	      'BEGIN {printf "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", t, v, b, bb, m, mb}'; \
	} > $@

$(NWK:.nwk=.time_cmp): %.time_cmp: %.time \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast.term \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast-beagle.term \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes.term \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes-beagle.term
	{ echo -e "samp\tvine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle"; \
	  vine_time=$$(head -1 $< | awk '{print $$1}' | sed 's/user//'); \
	  beast_time=$$(grep '^Total calculation time' $(word 2,$^) | awk '{print $$4}'); \
	  beast_beagle_time=$$(grep '^Total calculation time' $(word 3,$^) | awk '{print $$4}'); \
	  mrbayes_time=$$(grep 'Analysis used' $(word 4,$^) | awk '{print $$3}'); \
	  mrbayes_beagle_time=$$(grep 'Analysis used' $(word 5,$^) | awk '{print $$3}'); \
	  awk -v s="$(word 2,$(subst ., ,$(notdir $*)))" -v v="$$vine_time" -v b="$$beast_time" -v bb="$$beast_beagle_time" -v m="$$mrbayes_time" -v mb="$$mrbayes_beagle_time" \
	      'BEGIN {printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", s, v, b, bb, m, mb}'; \
	} > $@

$(NWK:.nwk=.var.dist.txt): %.var.dist.txt: %.nwk $$(BASE_DATA_DIR)/$$(basename $$*).true.dist.txt
	$(VINE_BIN)/evalTrees $< > $@.tmp
	python3 "$(PYTHON_SRC)/sumDists.py" $@.tmp $(word 2,$^) | grep -v '^#' > $@
	rm -f $@.tmp

BASENAMES := $(foreach n,$(TAXA),$(foreach i,$(SAMPLES),$(n)taxa/tree.$(i)))

$(NWK:.nwk=.kl_pdist_beast): %.kl_pdist_beast: %.nwk $$(BASE_DATA_DIR)/$$(basename $$*).beast.nwk
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

$(NWK:.nwk=.kl_pdist_mrbayes): %.kl_pdist_mrbayes: %.nwk $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes.nwk
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

BEAST_BEAGLE_KL   := $(foreach b,$(BASENAMES),$(b).beast-beagle.kl_pdist_beast)
MRBAYES_BEAGLE_KL := $(foreach b,$(BASENAMES),$(b).mrbayes-beagle.kl_pdist_mrbayes)

$(BEAST_BEAGLE_KL): %.beast-beagle.kl_pdist_beast: \
    $$(BASE_DATA_DIR)/%.beast-beagle.nwk $$(BASE_DATA_DIR)/%.beast.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

$(MRBAYES_BEAGLE_KL): %.mrbayes-beagle.kl_pdist_mrbayes: \
    $$(BASE_DATA_DIR)/%.mrbayes-beagle.nwk $$(BASE_DATA_DIR)/%.mrbayes.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

# Cross-method comparisons: MCMC methods measured against each other's reference too,
# so the "vs BEAST2" panel can show MrBayes/MrBayes+BEAGLE and vice versa.
MRBAYES_KL_BEAST        := $(foreach b,$(BASENAMES),$(b).mrbayes.kl_pdist_beast)
MRBAYES_BEAGLE_KL_BEAST := $(foreach b,$(BASENAMES),$(b).mrbayes-beagle.kl_pdist_beast)
BEAST_KL_MRBAYES        := $(foreach b,$(BASENAMES),$(b).beast.kl_pdist_mrbayes)
BEAST_BEAGLE_KL_MRBAYES := $(foreach b,$(BASENAMES),$(b).beast-beagle.kl_pdist_mrbayes)

$(MRBAYES_KL_BEAST): %.mrbayes.kl_pdist_beast: \
    $$(BASE_DATA_DIR)/%.mrbayes.nwk $$(BASE_DATA_DIR)/%.beast.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

$(MRBAYES_BEAGLE_KL_BEAST): %.mrbayes-beagle.kl_pdist_beast: \
    $$(BASE_DATA_DIR)/%.mrbayes-beagle.nwk $$(BASE_DATA_DIR)/%.beast.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

$(BEAST_KL_MRBAYES): %.beast.kl_pdist_mrbayes: \
    $$(BASE_DATA_DIR)/%.beast.nwk $$(BASE_DATA_DIR)/%.mrbayes.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

$(BEAST_BEAGLE_KL_MRBAYES): %.beast-beagle.kl_pdist_mrbayes: \
    $$(BASE_DATA_DIR)/%.beast-beagle.nwk $$(BASE_DATA_DIR)/%.mrbayes.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees --pairwise-dist-kl $< $(word 2,$^) 2>/dev/null | awk '/^Mean_pairwise_dist_KL:/ {print $$2}' > $@

$(KL_PDIST_BEAST_SUM): %.eval.all.kl_pdist_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).kl_pdist_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.kl_pdist_beast,$^) > $@

$(KL_PDIST_MRBAYES_SUM): %.eval.all.kl_pdist_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).kl_pdist_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.kl_pdist_mrbayes,$^) > $@

$(KL_PDIST_MRBAYES_VS_BEAST_SUM): %/mrbayes.eval.all.kl_pdist_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).mrbayes.kl_pdist_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.mrbayes.kl_pdist_beast,$^) > $@

$(KL_PDIST_MRBAYES_BEAGLE_VS_BEAST_SUM): %/mrbayes-beagle.eval.all.kl_pdist_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).mrbayes-beagle.kl_pdist_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.mrbayes-beagle.kl_pdist_beast,$^) > $@

$(KL_PDIST_BEAST_VS_MRBAYES_SUM): %/beast.eval.all.kl_pdist_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).beast.kl_pdist_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.beast.kl_pdist_mrbayes,$^) > $@

$(KL_PDIST_BEAST_BEAGLE_VS_MRBAYES_SUM): %/beast-beagle.eval.all.kl_pdist_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).beast-beagle.kl_pdist_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.beast-beagle.kl_pdist_mrbayes,$^) > $@

$(KL_PDIST_BEAST_BEAGLE_SUM): %/beast-beagle.eval.all.kl_pdist_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).beast-beagle.kl_pdist_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.kl_pdist_beast,$^) > $@

$(KL_PDIST_MRBAYES_BEAGLE_SUM): %/mrbayes-beagle.eval.all.kl_pdist_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).mrbayes-beagle.kl_pdist_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_pdist"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.mrbayes-beagle.kl_pdist_mrbayes,$^) > $@

# --- topology-based KL divergence (per-split Bernoulli KL, via compareTrees) ---

$(NWK:.nwk=.kl_topo_beast): %.kl_topo_beast: %.nwk $$(BASE_DATA_DIR)/$$(basename $$*).beast.nwk
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

$(NWK:.nwk=.kl_topo_mrbayes): %.kl_topo_mrbayes: %.nwk $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes.nwk
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

BEAST_BEAGLE_KL_TOPO   := $(foreach b,$(BASENAMES),$(b).beast-beagle.kl_topo_beast)
MRBAYES_BEAGLE_KL_TOPO := $(foreach b,$(BASENAMES),$(b).mrbayes-beagle.kl_topo_mrbayes)

$(BEAST_BEAGLE_KL_TOPO): %.beast-beagle.kl_topo_beast: \
    $$(BASE_DATA_DIR)/%.beast-beagle.nwk $$(BASE_DATA_DIR)/%.beast.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

$(MRBAYES_BEAGLE_KL_TOPO): %.mrbayes-beagle.kl_topo_mrbayes: \
    $$(BASE_DATA_DIR)/%.mrbayes-beagle.nwk $$(BASE_DATA_DIR)/%.mrbayes.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

# Cross-method comparisons: MCMC methods measured against each other's reference too,
# so the "vs BEAST2" panel can show MrBayes/MrBayes+BEAGLE and vice versa.
MRBAYES_KL_TOPO_BEAST        := $(foreach b,$(BASENAMES),$(b).mrbayes.kl_topo_beast)
MRBAYES_BEAGLE_KL_TOPO_BEAST := $(foreach b,$(BASENAMES),$(b).mrbayes-beagle.kl_topo_beast)
BEAST_KL_TOPO_MRBAYES        := $(foreach b,$(BASENAMES),$(b).beast.kl_topo_mrbayes)
BEAST_BEAGLE_KL_TOPO_MRBAYES := $(foreach b,$(BASENAMES),$(b).beast-beagle.kl_topo_mrbayes)

$(MRBAYES_KL_TOPO_BEAST): %.mrbayes.kl_topo_beast: \
    $$(BASE_DATA_DIR)/%.mrbayes.nwk $$(BASE_DATA_DIR)/%.beast.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

$(MRBAYES_BEAGLE_KL_TOPO_BEAST): %.mrbayes-beagle.kl_topo_beast: \
    $$(BASE_DATA_DIR)/%.mrbayes-beagle.nwk $$(BASE_DATA_DIR)/%.beast.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

$(BEAST_KL_TOPO_MRBAYES): %.beast.kl_topo_mrbayes: \
    $$(BASE_DATA_DIR)/%.beast.nwk $$(BASE_DATA_DIR)/%.mrbayes.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

$(BEAST_BEAGLE_KL_TOPO_MRBAYES): %.beast-beagle.kl_topo_mrbayes: \
    $$(BASE_DATA_DIR)/%.beast-beagle.nwk $$(BASE_DATA_DIR)/%.mrbayes.nwk
	mkdir -p $(dir $@)
	$(VINE_BIN1)/compareTrees $< $(word 2,$^) 2>/dev/null | awk '/^Mean_split_KL:/ {print $$2}' > $@

$(KL_TOPO_BEAST_SUM): %.eval.all.kl_topo_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).kl_topo_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.kl_topo_beast,$^) > $@

$(KL_TOPO_MRBAYES_SUM): %.eval.all.kl_topo_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).kl_topo_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.kl_topo_mrbayes,$^) > $@

$(KL_TOPO_MRBAYES_VS_BEAST_SUM): %/mrbayes.eval.all.kl_topo_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).mrbayes.kl_topo_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.mrbayes.kl_topo_beast,$^) > $@

$(KL_TOPO_MRBAYES_BEAGLE_VS_BEAST_SUM): %/mrbayes-beagle.eval.all.kl_topo_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).mrbayes-beagle.kl_topo_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.mrbayes-beagle.kl_topo_beast,$^) > $@

$(KL_TOPO_BEAST_VS_MRBAYES_SUM): %/beast.eval.all.kl_topo_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).beast.kl_topo_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.beast.kl_topo_mrbayes,$^) > $@

$(KL_TOPO_BEAST_BEAGLE_VS_MRBAYES_SUM): %/beast-beagle.eval.all.kl_topo_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).beast-beagle.kl_topo_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.beast-beagle.kl_topo_mrbayes,$^) > $@

$(KL_TOPO_BEAST_BEAGLE_SUM): %/beast-beagle.eval.all.kl_topo_beast.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).beast-beagle.kl_topo_beast)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.kl_topo_beast,$^) > $@

$(KL_TOPO_MRBAYES_BEAGLE_SUM): %/mrbayes-beagle.eval.all.kl_topo_mrbayes.txt: \
    $$(foreach i,$$(SAMPLES),$$*/tree.$$(i).mrbayes-beagle.kl_topo_mrbayes)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "kl_topo"} {print; sum+=$$1; n++} END {print "-----"; printf "%.6f\n", sum/n}' \
	    $(filter %.mrbayes-beagle.kl_topo_mrbayes,$^) > $@

$(NWK:.nwk=.dist): %.dist: %.var.dist.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast.dist.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast-beagle.dist.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes.dist.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes-beagle.dist.txt
	paste $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) > $@

$(NWK:.nwk=.var.rf.txt): %.var.rf.txt: %.nwk $$(BASE_DATA_DIR)/$$(basename $$*).true.nwk
	$(VINE_BIN)/evalTrees $< -t $(word 2,$^) > $@

$(NWK:.nwk=.rf): %.rf: %.var.rf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast.rf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast-beagle.rf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes.rf.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes-beagle.rf.txt
	awk '$$1 == "Mean:" {printf "%f\t", $$2}' $(filter %.rf.txt,$^) > $@
	echo >> $@

$(RF_SUM): %.eval.all.rf.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).rf)
	mkdir -p $(dir $@)
	nwk="$(BASE_DATA_DIR)/$(dir $*)tree.1.true.nwk"; \
	n=$$(grep -o ',' "$$nwk" | wc -l); n=$$((n + 1)); \
	div=$$((n - 3)); [ $$div -lt 1 ] && div=1; \
	cat $(filter %.rf,$^) | \
	awk -v d="$$div" 'BEGIN {print "vine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle"} \
	     {for(i=1;i<=NF;i++) {$$i=$$i/d; sum[i]+=$$i}; print} \
	     END {printf "-----\n"; for(i=1;i<=NF;i++) printf "%s%f", (i==1?"":"\t"), sum[i]/NR; printf "\n"}' > $@

# Expected normalized Kuhner-Felsenstein branch-score distance from each
# method's posterior sample to the true tree (matching hky_300sites/common.mk).
$(NWK:.nwk=.var.bsd.txt): %.var.bsd.txt: %.nwk $$(BASE_DATA_DIR)/$$(basename $$*).true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

$(foreach b,$(BASENAMES),$(foreach m,beast beast-beagle mrbayes mrbayes-beagle,$(b).$(m).bsd.txt)): \
%.bsd.txt: $$(BASE_DATA_DIR)/%.nwk $$(BASE_DATA_DIR)/$$(basename $$*).true.nwk
	$(VINE_BIN)/evalTrees $< -b $(word 2,$^) > $@

$(NWK:.nwk=.bsd): %.bsd: %.var.bsd.txt \
    $$(dir $$*)$$(basename $$(notdir $$*)).beast.bsd.txt \
    $$(dir $$*)$$(basename $$(notdir $$*)).beast-beagle.bsd.txt \
    $$(dir $$*)$$(basename $$(notdir $$*)).mrbayes.bsd.txt \
    $$(dir $$*)$$(basename $$(notdir $$*)).mrbayes-beagle.bsd.txt Makefile
	awk '/^Mean:/ {mean=$$2} /Reference tree length:/ {rl=$$NF} \
	     ENDFILE {printf "%s%.8f", (ARGIND == 1 ? "" : "\t"), (rl > 0 ? mean/rl : 0); mean=rl=0} \
	     END {print ""}' $(filter %.bsd.txt,$^) > $@

$(BSD_SUM): %.eval.all.bsd.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).bsd)
	mkdir -p $(dir $@)
	awk 'BEGIN {print "vine\tbeast\tbeast-beagle\tmrbayes\tmrbayes-beagle"} \
	     {for(i=1;i<=NF;i++) sum[i]+=$$i; print} \
	     END {print "-----"; for(i=1;i<=NF;i++) printf "%s%.8f", (i==1?"":"\t"), sum[i]/NR; print ""}' \
	    $(filter %.bsd,$^) > $@

$(NWK:.nwk=.var.ent.txt): %.var.ent.txt: %.nwk
	$(VINE_BIN)/evalTrees -e $< | awk '{printf "%f\t", $$NF} END {printf "\n"}' > $@

$(NWK:.nwk=.ent): %.ent: %.var.ent.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast.ent.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).beast-beagle.ent.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes.ent.txt \
    $$(BASE_DATA_DIR)/$$(basename $$*).mrbayes-beagle.ent.txt
	paste $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) > $@

$(DIST_SUM): %.eval.all.dist.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).dist)
	mkdir -p $(dir $@)
	cat $(filter %.dist,$^) | awk '{print $$1, $$2, $$3, $$4, $$6, $$7, $$8, $$9, $$11, $$12, $$13, $$14, $$16, $$17, $$18, $$19, $$21, $$22, $$23, $$24}' | \
	awk 'BEGIN {print "vine_r2 vine_RMSE vine_95CI vine_50CI beast_r2 beast_RMSE beast_95CI beast_50CI beast-beagle_r2 beast-beagle_RMSE beast-beagle_95CI beast-beagle_50CI mrbayes_r2 mrbayes_RMSE mrbayes_95CI mrbayes_50CI mrbayes-beagle_r2 mrbayes-beagle_RMSE mrbayes-beagle_95CI mrbayes-beagle_50CI"} \
	     {for(i=1;i<=NF;i++) sum[i]+=$$i; print} \
	     END {printf "-----\n"; for(i=1;i<=NF;i++) printf "%s%f", (i == 1 ? "" : "\t"), sum[i]/NR; printf "\n"}' > $@

$(ENT_SUM): %.eval.all.ent.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).ent)
	mkdir -p $(dir $@)
	cat $(filter %.ent,$^) | \
	awk 'BEGIN {print "vine_spl\tvine_top\tvine_br\tbeast_spl\tbeast_top\tbeast_br\tbeast-beagle_spl\tbeast-beagle_top\tbeast-beagle_br\tmrbayes_spl\tmrbayes_top\tmrbayes_br\tmrbayes-beagle_spl\tmrbayes-beagle_top\tmrbayes-beagle_br"} \
	     {for(i=1;i<=NF;i++) sum[i]+=$$i; print} \
	     END {printf "-----\n"; for(i=1;i<=NF;i++) printf "%s%f", (i == 1 ? "" : "\t"), sum[i]/NR; printf "\n"}' > $@

$(LNL_SUM): %.eval.all.lnl.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).lnl)
	mkdir -p $(dir $@)
	awk 'FNR==1 && NR==1 {print; next} \
	     FNR==2 {print; for(i=1;i<=NF;i++) sum[i]+=$$i; n++} \
	     END {if(n>0) {print "-----"; printf "all"; \
	          for(i=1;i<=NF;i++) printf "\t%.6f", sum[i]/n; printf "\n"}}' $(filter %.lnl,$^) > $@

$(MF_SUM): %.eval.all.mf.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).mf)
	mkdir -p $(dir $@)
	awk 'FNR==1 && NR==1 {print; next} \
	     FNR==2 {print; for(i=1;i<=NF;i++) sum[i]+=$$i; n++} \
	     END {if(n>0) {print "-----"; printf "all"; \
	          for(i=1;i<=NF;i++) printf "\t%.6f", sum[i]/n; printf "\n"}}' $(filter %.mf,$^) > $@

$(TIME_SUM): %.eval.all.time.txt: \
    $$(foreach i,$$(SAMPLES),$$(dir $$*)tree.$$(i).$$(notdir $$*).time_cmp)
	mkdir -p $(dir $@)
	awk 'FNR==1 && NR==1 {print; next} \
	     FNR==2 {print; for(i=2;i<=NF;i++) sum[i]+=$$i; n++} \
	     END {if(n>0) {print "-----"; printf "all"; \
	          for(i=2;i<=NF;i++) printf "\t%.6f", sum[i]/n; printf "\n"}}' $(filter %.time_cmp,$^) > $@

plots:
	mkdir -p $@

plots/hky300_summary.pdf: $(SUMMARIES) plot_posterior_quality.py | plots
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) plot_posterior_quality.py --taxa $(TAXA) --methods $(METHODS) --root .
	test -s $@

plots/hky300_posterior_quality_sup1.pdf: $(SUMMARIES) plot_posterior_quality_sup1.R | plots
	$(R_SCRIPT) plot_posterior_quality_sup1.R
	test -s $@

plots/hky300_posterior_quality_sup2.pdf: $(SUMMARIES) plot_posterior_quality_sup2.R | plots
	$(R_SCRIPT) plot_posterior_quality_sup2.R
	test -s $@

clean-vine:
	rm -f */tree.*.bootstrap_nj.nwk */tree.*.bootstrap_nj.nwk.log */tree.*.bootstrap_nj.time */tree.*.nj.time \
	      */tree.*.bootstrap_nj.varlnl */tree.*.bootstrap_nj.lnl */tree.*.bootstrap_nj.time_cmp \
	      */tree.*.bootstrap_nj.var.dist.txt */tree.*.bootstrap_nj.dist */tree.*.bootstrap_nj.var.ent.txt */tree.*.bootstrap_nj.ent \
	      */tree.*.bootstrap_nj.var.rf.txt */tree.*.bootstrap_nj.rf */tree.*.bootstrap_nj.var.bsd.txt */tree.*.bootstrap_nj.bsd \
	      */tree.*.bootstrap_nj.var.mf.txt */tree.*.bootstrap_nj.mf \
	      */tree.*.bootstrap_nj.kl_pdist_* */tree.*.bootstrap_nj.kl_topo_* \
	      */tree.*.vine.nwk */tree.*.vine.nwk.log */tree.*.vine.time */tree.*.vine.mean.nwk \
	      */tree.*.vine.varlnl */tree.*.vine.lnl */tree.*.vine.time_cmp \
	      */tree.*.vine.var.dist.txt */tree.*.vine.dist */tree.*.vine.var.ent.txt */tree.*.vine.ent \
	      */tree.*.vine.var.rf.txt */tree.*.vine.rf */tree.*.vine.var.bsd.txt */tree.*.vine.bsd \
	      */tree.*.vine.var.mf.txt */tree.*.vine.mf \
	      */tree.*.vine_flows.nwk */tree.*.vine_flows.nwk.log */tree.*.vine_flows.time */tree.*.vine_flows.mean.nwk \
	      */tree.*.vine_flows.varlnl */tree.*.vine_flows.lnl */tree.*.vine_flows.time_cmp \
	      */tree.*.vine_flows.var.dist.txt */tree.*.vine_flows.dist */tree.*.vine_flows.var.ent.txt */tree.*.vine_flows.ent \
	      */tree.*.vine_flows.var.rf.txt */tree.*.vine_flows.rf */tree.*.vine_flows.var.bsd.txt */tree.*.vine_flows.bsd \
	      */tree.*.vine_flows.var.mf.txt */tree.*.vine_flows.mf \
	      */tree.*.vine_mixture.nwk */tree.*.vine_mixture.nwk.log */tree.*.vine_mixture.time */tree.*.vine_mixture.mean.nwk \
	      */tree.*.vine_mixture.varlnl */tree.*.vine_mixture.lnl */tree.*.vine_mixture.time_cmp \
	      */tree.*.vine_mixture.var.dist.txt */tree.*.vine_mixture.dist */tree.*.vine_mixture.var.ent.txt */tree.*.vine_mixture.ent \
	      */tree.*.vine_mixture.var.rf.txt */tree.*.vine_mixture.rf */tree.*.vine_mixture.var.bsd.txt */tree.*.vine_mixture.bsd \
	      */tree.*.vine_mixture.var.mf.txt */tree.*.vine_mixture.mf \
	      */tree.*.vine_mixture_flows.nwk */tree.*.vine_mixture_flows.nwk.log */tree.*.vine_mixture_flows.time */tree.*.vine_mixture_flows.mean.nwk \
	      */tree.*.vine_mixture_flows.varlnl */tree.*.vine_mixture_flows.lnl */tree.*.vine_mixture_flows.time_cmp \
	      */tree.*.vine_mixture_flows.var.dist.txt */tree.*.vine_mixture_flows.dist */tree.*.vine_mixture_flows.var.ent.txt */tree.*.vine_mixture_flows.ent \
	      */tree.*.vine_mixture_flows.var.rf.txt */tree.*.vine_mixture_flows.rf */tree.*.vine_mixture_flows.var.bsd.txt */tree.*.vine_mixture_flows.bsd \
	      */tree.*.vine_mixture_flows.var.mf.txt */tree.*.vine_mixture_flows.mf \
	      */tree.*.vine*.pairwise_rf.* \
	      */tree.*.vine.pdist.txt */tree.*.vine.kl_pdist_beast */tree.*.vine.kl_pdist_mrbayes \
	      */tree.*.vine_flows.pdist.txt */tree.*.vine_flows.kl_pdist_beast */tree.*.vine_flows.kl_pdist_mrbayes \
	      */tree.*.vine_mixture.pdist.txt */tree.*.vine_mixture.kl_pdist_beast */tree.*.vine_mixture.kl_pdist_mrbayes \
	      */tree.*.vine_mixture_flows.pdist.txt */tree.*.vine_mixture_flows.kl_pdist_beast */tree.*.vine_mixture_flows.kl_pdist_mrbayes \
	      */tree.*.vine.kl_topo_beast */tree.*.vine.kl_topo_mrbayes \
	      */tree.*.vine_flows.kl_topo_beast */tree.*.vine_flows.kl_topo_mrbayes \
	      */tree.*.vine_mixture.kl_topo_beast */tree.*.vine_mixture.kl_topo_mrbayes \
	      */tree.*.vine_mixture_flows.kl_topo_beast */tree.*.vine_mixture_flows.kl_topo_mrbayes \
	      */tree.*.beast.bsd.txt */tree.*.beast-beagle.bsd.txt */tree.*.mrbayes.bsd.txt */tree.*.mrbayes-beagle.bsd.txt \
	      $(DIST_SUM) $(RF_SUM) $(ENT_SUM) $(LNL_SUM) $(TIME_SUM) $(MF_SUM) $(BSD_SUM) \
	      $(KL_PDIST_BEAST_SUM) $(KL_PDIST_MRBAYES_SUM) \
	      $(KL_TOPO_BEAST_SUM) $(KL_TOPO_MRBAYES_SUM)

clean-plots:
	rm -rf plots __pycache__

clean-rf:
	rm -f $(PAIRWISE_RF)

clean-mds:
	rm -f $(MDS_COORDS)

clean-varlnl:
	rm -rf */tree.*.varlnl
