export SHELL := /usr/bin/bash

# This benchmark reuses the simulations and comparison-method summaries from
# hky_300sites.  Only Vine with the alternate tree prior is run here.
.SECONDARY:
.PRECIOUS:

MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT := $(MAIN_DIR)/dna_sims/hky_300sites_treeprior
BASE_ROOT := $(MAIN_DIR)/dna_sims/hky_300sites
BASE_DIR := $(BASE_ROOT)/$(notdir $(CURDIR))
BIN := $(MAIN_DIR)/bin
VINE_BIN := $(BIN)/vine/bin
PYTHON_SRC := $(MAIN_DIR)/python/src

TREES := $(shell seq -f tree.%.0f.true.nwk 1 $(NSAMP))
VAR := $(patsubst %.true.nwk,%.var.nwk,$(TREES))
VARLOG := $(patsubst %.true.nwk,%.var.nwk.log,$(TREES))
VARTIME := $(patsubst %.true.nwk,%.var-time,$(TREES))
VARLNL := $(patsubst %.true.nwk,%.varlnl,$(TREES))
VARMF := $(patsubst %.true.nwk,%.var.mf.txt,$(TREES))
VARRF := $(patsubst %.true.nwk,%.var.rf.txt,$(TREES))
VARDIST := $(patsubst %.true.nwk,%.var.dist.txt,$(TREES))
VARENT := $(patsubst %.true.nwk,%.var.ent.txt,$(TREES))
VARBSD := $(patsubst %.true.nwk,%.var.bsd.txt,$(TREES))

SUMMARIES := eval.all.lnl.txt eval.all.time.txt eval.all.rf.txt \
             eval.all.mf.txt eval.all.dist.txt eval.all.ent.txt \
             eval.all.bsd.txt

.PHONY: all vine summaries clean_vine archive_vine

all: $(SUMMARIES)
vine: $(VAR) $(VARLOG)
summaries: $(SUMMARIES)

# Use the matching simulations directly from hky_300sites; no simulated data
# are copied into this alternate-method benchmark.
tree.%.var.nwk tree.%.var-time tree.%.var.nwk.log: $(BASE_DIR)/tree.%.fa
	/usr/bin/time -o tree.$*.var-time $(VINE_BIN)/vine $< \
	  -l tree.$*.var.nwk.log $(VAROPT) --mean tree.$*.mean.nwk \
	  > tree.$*.var.nwk

# Final unbiased MC estimate of E_q[lnL], falling back to LNL for older logs.
tree.%.varlnl: tree.%.var.nwk.log
	tail -1 $< | awk '{ll=""; mc=""; for(i=1;i<=NF;i++){if($$i=="LNL:")ll=$$(i+1); if($$i=="LNL_mc:")mc=$$(i+1)} v=(mc!=""?mc:ll); gsub(/,/,"",v); print v}' > $@

tree.%.var.mf.txt: tree.%.var.nwk $(BASE_DIR)/tree.%.heldout.fa tree.%.var.nwk.log
	kappa=$$(awk '{for(i=1;i<=NF;i++) if($$i=="kappa:") k=$$(i+1)} END{print k}' tree.$*.var.nwk.log); \
	$(VINE_BIN)/evalTrees tree.$*.var.nwk -f $(BASE_DIR)/tree.$*.heldout.fa -k "$$kappa" > $@.tmp && mv $@.tmp $@

tree.%.var.rf.txt: tree.%.var.nwk $(BASE_DIR)/tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -t $(BASE_DIR)/tree.$*.true.nwk > $@

tree.%.true.dist.txt: $(BASE_DIR)/tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< > $@

tree.%.var.dist.txt: tree.%.var.nwk tree.%.true.dist.txt
	$(VINE_BIN)/evalTrees $< > tmp.$*.var.dist
	python3 "$(PYTHON_SRC)/sumDists.py" tmp.$*.var.dist tree.$*.true.dist.txt | grep -v '^#' > $@
	rm -f tmp.$*.var.dist

tree.%.var.ent.txt: tree.%.var.nwk
	$(VINE_BIN)/evalTrees -e $< | awk '{printf "%f\t", $$NF} END{printf "\n"}' > $@

tree.%.var.bsd.txt: tree.%.var.nwk $(BASE_DIR)/tree.%.true.nwk
	$(VINE_BIN)/evalTrees $< -b $(BASE_DIR)/tree.$*.true.nwk > $@

# Each aggregate starts with the corresponding file from hky_300sites verbatim,
# then appends the new Vine-prior method block row by row.
eval.all.lnl.txt: $(BASE_DIR)/eval.all.lnl.txt $(VARLNL)
	{ echo vine-prior; for f in $(VARLNL); do cat "$$f"; done; awk '{s+=$$1} END{print "-----"; if(NR) printf "%f\n",s/NR}' $(VARLNL); } > $@.new
	paste $< $@.new > $@.tmp && mv $@.tmp $@
	rm -f $@.new

eval.all.time.txt: $(BASE_DIR)/eval.all.time.txt $(VARTIME)
	{ echo vine-prior; for f in $(VARTIME); do head -1 "$$f" | sed 's/user.*//'; done; awk '{s+=$$1} END{print "-----"; if(NR) printf "%.2f\n",s/NR}' $(VARTIME); } > $@.new
	paste $< $@.new > $@.tmp && mv $@.tmp $@
	rm -f $@.new

eval.all.rf.txt: $(BASE_DIR)/eval.all.rf.txt $(VARRF)
	n=$$(grep -o ',' $(BASE_DIR)/tree.1.true.nwk | wc -l); d=$$((n-2)); [ $$d -lt 1 ] && d=1; \
	{ echo 'vine-prior (sd)'; for f in $(VARRF); do awk -v d="$$d" '/Mean:/{m=$$2/d}/Std:/{printf "%f\t%f\n",m,$$2/d}' "$$f"; done; awk -v d="$$d" '/Mean:/{sm+=$$2/d}/Std:/{ss+=($$2/d)^2;n++} END{print "-----"; if(n) printf "%f\t%f\n",sm/n,sqrt(ss/n)}' $(VARRF); } > $@.new; \
	paste $< $@.new > $@.tmp && mv $@.tmp $@; rm -f $@.new

eval.all.mf.txt: $(BASE_DIR)/eval.all.mf.txt $(VARMF)
	{ echo 'vine-prior (sd)'; for f in $(VARMF); do awk '/Mean:/{m=$$2}/Std:/{printf "%f\t%f\n",m,$$2}' "$$f"; done; awk '/Mean:/{sm+=$$2}/Std:/{ss+=$$2^2;n++} END{print "-----"; if(n) printf "%f\t%f\n",sm/n,sqrt(ss/n)}' $(VARMF); } > $@.new
	paste $< $@.new > $@.tmp && mv $@.tmp $@
	rm -f $@.new

eval.all.dist.txt: $(BASE_DIR)/eval.all.dist.txt $(VARDIST)
	{ echo 'vine-prior_r2 vine-prior_RMSE vine-prior_95CI vine-prior_50CI'; cat $(VARDIST) | awk '{print $$1,$$2,$$3,$$4; for(i=1;i<=4;i++)s[i]+=$$i;n++} END{print "-----"; if(n) printf "%f\t%f\t%f\t%f\n",s[1]/n,s[2]/n,s[3]/n,s[4]/n}'; } > $@.new
	paste $< $@.new > $@.tmp && mv $@.tmp $@
	rm -f $@.new

eval.all.ent.txt: $(BASE_DIR)/eval.all.ent.txt $(VARENT)
	{ echo 'vine-prior_spl vine-prior_top vine-prior_br'; cat $(VARENT) | awk '{print $$1,$$2,$$3; for(i=1;i<=3;i++)s[i]+=$$i;n++} END{print "-----"; if(n) printf "%f\t%f\t%f\n",s[1]/n,s[2]/n,s[3]/n}'; } > $@.new
	paste $< $@.new > $@.tmp && mv $@.tmp $@
	rm -f $@.new

eval.all.bsd.txt: $(BASE_DIR)/eval.all.bsd.txt $(VARBSD)
	{ echo 'vine-prior (sd)'; for f in $(VARBSD); do awk '/Point.*BSD:/{p=$$NF}/Reference tree length:/{r=$$NF} END{print(r>0?p/r:0)}' "$$f"; done; awk '/Point.*BSD:/{p=$$NF}/Reference tree length:/{r=$$NF; if(r>0){x=p/r;s+=x;ss+=x*x;n++}} END{print "-----"; if(n){m=s/n;v=ss/n-m*m;if(v<0)v=0;printf "%f\t%f\n",m,sqrt(v)}}' $(VARBSD); } > $@.new
	paste $< $@.new > $@.tmp && mv $@.tmp $@
	rm -f $@.new

clean_vine:
	rm -f tree.*.var* tree.*.mean*.nwk $(SUMMARIES)
	rm -f tree.*.true.dist.txt

archive_vine:
	archive_dir="archive.vine_$$(date +%Y-%m-%d_%H:%M:%S)"; \
	mkdir "$$archive_dir"; \
	mv tree.*.var* tree.*.mean*.nwk $(SUMMARIES) "$$archive_dir"/
