# KL-only posterior comparison workflow.
#
# This file is intentionally isolated from `all`: only an explicit
# `make kl-metrics` (or `make kl-inputs`) activates it. Source posterior files
# are checked inside recipes rather than listed as prerequisites, preventing
# Make from invoking any inference/MCMC rule when an artifact is absent.

KL_TREE_COUNT ?= 1000
KL_METHODS := vine vine_flows dodonaphy geophy vaiphy beast mrbayes
KL_SAMPLES := $(shell seq 1 $(NSAMP))

KL_DOWNSAMPLE_NEXUS := $(PYTHON_SRC)/downsample_nexus_trees.py
KL_EXPORT_GEOPHY := $(PYTHON_SRC)/export_geophy_trees.py
KL_EXPORT_VAIPHY := $(PYTHON_SRC)/export_vaiphy_trees.py

KL_INPUTS := $(foreach i,$(KL_SAMPLES),$(foreach m,$(KL_METHODS),tree.$(i).$(m).kl.nwk))
KL_SUMMARIES := $(foreach m,$(KL_METHODS),\
	$(m).eval.all.kl_pdist_mrbayes.txt \
	$(m).eval.all.kl_pdist_beast.txt \
	$(m).eval.all.kl_topo_mrbayes.txt \
	$(m).eval.all.kl_topo_beast.txt)
KL_METHOD_TARGETS := $(addsuffix _kl,$(KL_METHODS))

.PHONY: kl-metrics kl-inputs kl-check-sources $(KL_METHOD_TARGETS)
kl-metrics: $(KL_SUMMARIES)
kl-inputs: $(KL_INPUTS)

# Per-method entry points. Each builds the four pairwise-distance/topology
# summaries against BEAST2 and MrBayes for one fitted method only.
define KL_METHOD_TARGET_RULE
$(1)_kl: $(filter $(1).%,$(KL_SUMMARIES))
endef
$(foreach m,$(KL_METHODS),$(eval $(call KL_METHOD_TARGET_RULE,$(m))))

# Read-only preflight. It reports missing fitted outputs and never asks Make to
# build them.
kl-check-sources:
	@set -eu; missing=0; \
	for f in $(GEOPHY_SIF) $(VAIPHY_SIF); do \
	  if [ ! -s "$$f" ]; then printf 'MISSING container %s\n' "$$f" >&2; missing=1; fi; \
	done; \
	for i in $(KL_SAMPLES); do \
	  for f in tree.$$i.var.nwk tree.$$i.dodonaphy.samples.t \
	      tree.$$i.geophy.latest.pt tree.$$i.geophy.latest.yaml \
	      tree.$$i.vaiphy.log \
	      tree.$$i/model_vaiphy_init_nj_phyml_samp_slantis_branch_ml_ng_0.1/S_128_seed_2/vaiphy_S_128_seed_2_object.pkl \
	      tree.$$i.beast-tree.trees tree.$$i.mrbayes.nex.t; do \
	    if [ ! -s "$$f" ]; then printf 'MISSING %s\n' "$$f" >&2; missing=1; fi; \
	  done; \
	done; \
	if [ "$$missing" -ne 0 ]; then exit 1; fi; \
	printf 'All fitted posterior sources are present; no inference was run.\n'

# Existing Vine posterior: retain at most the first 1000 nonempty Newick rows.
# Current runs contain exactly 1000, but the cap is enforced here.
tree.%.vine.kl.nwk:
	@set -eu; src=tree.$*.var.nwk; \
	test -s "$$src" || { printf 'Missing fitted Vine posterior: %s\n' "$$src" >&2; exit 1; }; \
	awk -v max=$(KL_TREE_COUNT) 'NF && n++ < max {print}' "$$src" > $@.tmp; \
	test -s $@.tmp; n=$$(awk 'NF{n++} END{print n+0}' $@.tmp); \
	test "$$n" -le $(KL_TREE_COUNT); mv $@.tmp $@

# This is the sole bridge from the isolated Vine + flows fit into KL metrics.
# Its raw prerequisite can only invoke the dedicated flow fit above.
tree.%.vine_flows.kl.nwk: tree.%.vine_flows.kl.raw.nwk
	@set -eu; \
	awk -v max=$(KL_TREE_COUNT) 'NF && n++ < max {print}' $< > $@.tmp; \
	n=$$(awk 'NF{n++} END{print n+0}' $@.tmp); \
	test "$$n" -eq $(KL_TREE_COUNT); mv $@.tmp $@

# Existing Dodonaphy posterior samples (Nexus); no burn-in applies to VI draws.
tree.%.dodonaphy.kl.nwk: $(KL_DOWNSAMPLE_NEXUS)
	@set -eu; src=tree.$*.dodonaphy.samples.t; \
	test -s "$$src" || { printf 'Missing fitted Dodonaphy posterior: %s\n' "$$src" >&2; exit 1; }; \
	trap 'rm -f $@.selected.nex $@.tmp' EXIT; \
	python3 $(KL_DOWNSAMPLE_NEXUS) "$$src" $@.selected.nex \
	  --max-trees $(KL_TREE_COUNT) --seed $*; \
	$(OTHER_BIN)/convertTrees -i nexus $@.selected.nex | \
	  sed 's/^\[&[^]]*\][[:space:]]*//' > $@.tmp; \
	n=$$(awk 'NF{n++} END{print n+0}' $@.tmp); test "$$n" -le $(KL_TREE_COUNT); \
	test "$$n" -gt 0; mv $@.tmp $@

# Draw from the saved GeoPhy variational state. This loads .latest.pt and does
# not execute run_gp.py or any optimizer/training code.
tree.%.geophy.kl.nwk: $(KL_EXPORT_GEOPHY)
	@set -eu; prefix=tree.$*.geophy; input=tree.$*.nex; state_input=$$prefix.nex; \
	test -s "$(GEOPHY_SIF)" || \
	  { printf 'Missing GeoPhy container: %s\n' '$(GEOPHY_SIF)' >&2; exit 1; }; \
	test -s "$$input" || \
	  { printf 'Missing existing GeoPhy alignment: %s\n' "$$input" >&2; exit 1; }; \
	test -s "$$prefix.latest.pt" && test -s "$$prefix.latest.yaml" || \
	  { printf 'Missing fitted GeoPhy checkpoint: %s.latest.{pt,yaml}\n' "$$prefix" >&2; exit 1; }; \
	trap 'rm -f "$$state_input" "$@.tmp"' EXIT; cp "$$input" "$$state_input"; \
	singularity exec --env PYTHONPATH=/opt/app \
	  --bind $(CURDIR):/mnt --bind $(MAIN_DIR):/repo --pwd /opt/app \
	  $(GEOPHY_SIF) python /repo/python/src/export_geophy_trees.py \
	  --state-prefix /mnt/$$prefix --output /mnt/$@.tmp \
	  --count $(KL_TREE_COUNT) --seed $*; \
	n=$$(awk 'NF{n++} END{print n+0}' $@.tmp); test "$$n" -eq $(KL_TREE_COUNT); mv $@.tmp $@

# Draw from the saved VaiPhy model object. This calls sample_trees on the fitted
# object and never calls train() or main.py.
tree.%.vaiphy.kl.nwk: $(KL_EXPORT_VAIPHY)
	@set -eu; \
	model=tree.$*/model_vaiphy_init_nj_phyml_samp_slantis_branch_ml_ng_0.1/S_128_seed_2/vaiphy_S_128_seed_2_object.pkl; \
	log=tree.$*.vaiphy.log; \
	alignment=tree.$*/tree.$*.nex; \
	test -s "$(VAIPHY_SIF)" || \
	  { printf 'Missing VaiPhy container: %s (build it from containers/vaiphy/vaiphy.def)\n' \
	      '$(VAIPHY_SIF)' >&2; exit 1; }; \
	test -s "$$model" && test -s "$$log" && test -s "$$alignment" || \
	  { printf 'Missing fitted VaiPhy model/log/alignment for tree.%s\n' '$*' >&2; exit 1; }; \
	singularity exec --env PYTHONPATH=/opt/vaiphy/src \
	  --bind $(CURDIR):/mnt --bind $(MAIN_DIR):/repo \
	  --bind $(MAIN_DIR)/bin/VaiPhy:/opt/vaiphy --pwd /opt/vaiphy/src \
	  $(VAIPHY_SIF) python /repo/python/src/export_vaiphy_trees.py \
	  --model /mnt/$$model --log /mnt/$$log --alignment /mnt/$$alignment \
	  --output /mnt/$@.tmp \
	  --count $(KL_TREE_COUNT) --seed $*; \
	n=$$(awk 'NF{n++} END{print n+0}' $@.tmp); test "$$n" -eq $(KL_TREE_COUNT); mv $@.tmp $@

# Extract existing MCMC samples with the configured 10% burn-in and a random
# cap of 1000. Raw logs are checked inside recipes so missing logs cannot
# trigger BEAST or MrBayes.
tree.%.beast.kl.nwk: $(KL_DOWNSAMPLE_NEXUS)
	@set -eu; src=tree.$*.beast-tree.trees; \
	test -s "$$src" || { printf 'Missing existing BEAST tree log: %s\n' "$$src" >&2; exit 1; }; \
	trap 'rm -f $@.selected.nex $@.tmp' EXIT; \
	python3 $(KL_DOWNSAMPLE_NEXUS) "$$src" $@.selected.nex \
	  --burnin-percent $(BURNIN_PCT) --max-trees $(KL_TREE_COUNT) --seed $*; \
	$(OTHER_BIN)/convertTrees -i nexus $@.selected.nex | \
	  sed 's/^\[&[^]]*\][[:space:]]*//' > $@.tmp; \
	n=$$(awk 'NF{n++} END{print n+0}' $@.tmp); test "$$n" -le $(KL_TREE_COUNT); \
	test "$$n" -gt 0; mv $@.tmp $@

tree.%.mrbayes.kl.nwk: $(KL_DOWNSAMPLE_NEXUS)
	@set -eu; src=tree.$*.mrbayes.nex.t; \
	test -s "$$src" || { printf 'Missing existing MrBayes tree log: %s\n' "$$src" >&2; exit 1; }; \
	trap 'rm -f $@.selected.nex $@.tmp' EXIT; \
	python3 $(KL_DOWNSAMPLE_NEXUS) "$$src" $@.selected.nex \
	  --burnin-percent $(BURNIN_PCT) --max-trees $(KL_TREE_COUNT) --seed $*; \
	$(OTHER_BIN)/convertTrees -i nexus $@.selected.nex | \
	  sed 's/^\[&[^]]*\][[:space:]]*//' > $@.tmp; \
	n=$$(awk 'NF{n++} END{print n+0}' $@.tmp); test "$$n" -le $(KL_TREE_COUNT); \
	test "$$n" -gt 0; mv $@.tmp $@

define KL_COMPARE_RULES
tree.$(1).$(2).kl_pdist_$(3): tree.$(1).$(2).kl.nwk tree.$(1).$(3).kl.nwk
	$$(VINE_BIN)/compareTrees --pairwise-dist-kl $$< $$(word 2,$$^) | \
	  awk '/^Mean_pairwise_dist_KL:/ {print $$$$2}' > $$@.tmp
	test -s $$@.tmp && mv $$@.tmp $$@

tree.$(1).$(2).kl_topo_$(3): tree.$(1).$(2).kl.nwk tree.$(1).$(3).kl.nwk
	$$(VINE_BIN)/compareTrees $$< $$(word 2,$$^) | \
	  awk '/^Mean_split_KL:/ {print $$$$2}' > $$@.tmp
	test -s $$@.tmp && mv $$@.tmp $$@
endef

$(foreach i,$(KL_SAMPLES),$(foreach r,beast mrbayes,\
  $(foreach m,$(filter-out $(r),$(KL_METHODS)),\
    $(eval $(call KL_COMPARE_RULES,$(i),$(m),$(r))))))

# A reference distribution compared with itself has KL divergence exactly 0.
KL_BEAST_SELF_PDIST := $(foreach i,$(KL_SAMPLES),tree.$(i).beast.kl_pdist_beast)
KL_BEAST_SELF_TOPO := $(foreach i,$(KL_SAMPLES),tree.$(i).beast.kl_topo_beast)
KL_MRBAYES_SELF_PDIST := $(foreach i,$(KL_SAMPLES),tree.$(i).mrbayes.kl_pdist_mrbayes)
KL_MRBAYES_SELF_TOPO := $(foreach i,$(KL_SAMPLES),tree.$(i).mrbayes.kl_topo_mrbayes)

$(KL_BEAST_SELF_PDIST): %.kl_pdist_beast: %.kl.nwk
	printf '0\n' > $@
$(KL_BEAST_SELF_TOPO): %.kl_topo_beast: %.kl.nwk
	printf '0\n' > $@
$(KL_MRBAYES_SELF_PDIST): %.kl_pdist_mrbayes: %.kl.nwk
	printf '0\n' > $@
$(KL_MRBAYES_SELF_TOPO): %.kl_topo_mrbayes: %.kl.nwk
	printf '0\n' > $@

define KL_SUMMARY_RULE
$(1).eval.all.kl_$(2)_$(3).txt: $(foreach i,$(KL_SAMPLES),tree.$(i).$(1).kl_$(2)_$(3))
	awk 'BEGIN {print "kl_$(2)"} {print; sum+=$$$$1; n++} \
	  END {print "-----"; if (n) printf "%.6f\n", sum/n}' $$^ > $$@.tmp
	test -s $$@.tmp && mv $$@.tmp $$@
endef

$(foreach m,$(KL_METHODS),$(foreach metric,pdist topo,$(foreach r,beast mrbayes,\
  $(eval $(call KL_SUMMARY_RULE,$(m),$(metric),$(r))))))
