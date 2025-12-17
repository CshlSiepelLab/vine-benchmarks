export SHELL=/usr/bin/bash

# Give beast2 and mrbayes access to beagle
export LD_LIBRARY_PATH := /local/storage/no-backup/vine-benchmarks/lib:$(LD_LIBRARY_PATH)

# Prevent Make from removing files it thinks it should clean up
#.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

# Edit for local structure; this is the only place absolute paths are used
MAIN_DIR := /local/storage/no-backup/vine-benchmarks
ROOT_SUFFIX := dna_sims/hky_10000sites_pilot_runs

ROOT := $(MAIN_DIR)/$(ROOT_SUFFIX)
BIN := $(MAIN_DIR)/bin
PHAST_BIN := $(BIN)/phast/bin
BEAST := $(BIN)/beast/bin/beast
MRBAYES := $(BIN)/mb
BEAST_TEMPLATE := $(ROOT)/beast_template.xml

TREES := $(shell seq -f tree.%.0f.true.nwk 1 $(NSAMP))
FA := $(patsubst %.true.nwk,%.fa,$(TREES))

BEASTLOG := $(patsubst %.true.nwk,%.beast.chain1.log,$(TREES)) $(patsubst %.true.nwk,%.beast.chain2.log,$(TREES))
MRBAYESLOG := $(patsubst %.true.nwk,%.mrbayes.nex.run1.p,$(TREES)) $(patsubst %.true.nwk,%.mrbayes.nex.run2.p,$(TREES))

BEAST_MCMCSTATS := $(patsubst %.true.nwk,%.beast.mcmc_convergence_stats_convergence_points.txt,$(TREES))
MRBAYES_MCMCSTATS := $(patsubst %.true.nwk,%.mrbayes.mcmc_convergence_stats_convergence_points.txt,$(TREES))

TREE_CONVERGENCE := $(patsubst %.true.nwk,%.mcmc_convergence.txt,$(TREES))

all: eval.all.convergence.txt

simulate: $(TREES) $(FA)
infer: $(BEASTLOG) $(MRBAYESLOG)
stats: $(BEAST_MCMCSTATS) $(MRBAYES_MCMCSTATS)

tree.%.true.nwk: 
	$(BIN)/bdTree3 -b 1 -d 0.5 --oversample-k 3 --height 5 --min-edge 0.02 --expected-height $(EXPHEIGHT) --no-stem --ucln-sd 0.6 --target-stat median -n $(NTAXA) | sed 's/\[\&[UR]\] //' > $@

tree.%.fa: tree.%.true.nwk
	cp ../base-hky.mod tmp.mod
	echo -n "TREE: " >> tmp.mod
	cat $< >> tmp.mod
	$(PHAST_BIN)/base_evolve --nsites $(NSITES) tmp.mod > $@
	rm tmp.mod

# Run two chains for beast2 to check convergence later
tree.%.beast.chain1.xml tree.%.beast.chain2.xml:
	cp "$(BEAST_TEMPLATE)" tree.$*.beast.chain1.xml
	cp "$(BEAST_TEMPLATE)" tree.$*.beast.chain2.xml

tree.%.beast.chain1.term tree.%.beast-tree.chain1.trees tree.%.beast.chain1.log: tree.%.beast.chain1.xml tree.%.fa
	rm -f tree.$*.beast-tree.chain1.trees tree.$*.beast.chain1.log
	"$(BEAST)" -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$(BEAST_SAMPLEFREQ) -D printfreq=$(PRINTFREQ) $< > tree.$*.beast.chain1.term

tree.%.beast.chain2.term tree.%.beast-tree.chain2.trees tree.%.beast.chain2.log: tree.%.beast.chain2.xml tree.%.fa
	rm -f tree.$*.beast-tree.chain2.trees tree.$*.beast.chain2.log
	"$(BEAST)" -working -D fastapath=tree.$*.fa -D mcmclength=$(BEAST_MCMCLEN) -D samplefreq=$(BEAST_SAMPLEFREQ) -D printfreq=$(PRINTFREQ) $< > tree.$*.beast.chain2.term


# Assess beast2 convergence
tree.%.beast.mcmc_convergence_stats_convergence_points.txt: tree.%.beast.chain1.log tree.%.beast.chain2.log
	$(BIN)/mcmc_convergence_stats \
		--logfiles tree.$*.beast.chain1.log,tree.$*.beast.chain2.log  \
		--parameters "posterior,likelihood,Tree.Length,Tree.height" \
		--outputprefix tree.$*.beast.mcmc_convergence_stats \
		--ess_calc_freq $(ESS_CALCFREQ) \
		--rhat_calc_freq $(RHAT_CALCFREQ)

tree.%.nex: tree.%.fa
	$(BIN)/fa2nex $< $@

# MrBayes input file prep (convert fasta to nexus and add MrBayes block to the end of nexus to specify the model)
tree.%.mrbayes.nex: tree.%.nex
	$(BIN)/addMrbayesModelToNex --in_nexus tree.$*.nex --out_nexus tree.$*.mrbayes.nex \
		--mcmc_length $(MRBAYES_MCMCLEN) --model HKY --nruns 2 --sample_freq $(MRBAYES_SAMPLEFREQ) \
		--print_freq $(PRINTFREQ) --diagn_freq $(PRINTFREQ) --use_beagle

# Run MrBayes with two chains to check convergence later
tree.%.mrbayes.term tree.%.mrbayes.nex.run1.p tree.%.mrbayes.nex.run1.t tree.%.mrbayes.nex.run2.p tree.%.mrbayes.nex.run2.t: tree.%.mrbayes.nex
	$(MRBAYES) tree.$*.mrbayes.nex > tree.$*.mrbayes.term

# Assess MrBayes convergence
tree.%.mrbayes.mcmc_convergence_stats_convergence_points.txt: tree.%.mrbayes.nex.run1.p tree.%.mrbayes.nex.run2.p
	$(BIN)/mcmc_convergence_stats \
		--logfiles tree.$*.mrbayes.nex.run1.p,tree.$*.mrbayes.nex.run2.p \
		--parameters "lnLike,TL" \
		--outputprefix tree.$*.mrbayes.mcmc_convergence_stats \
		--ess_calc_freq $(ESS_CALCFREQ) \
		--rhat_calc_freq $(RHAT_CALCFREQ)

# Summarize all convergence stats into single files
tree.%.mcmc_convergence.txt: tree.%.beast.mcmc_convergence_stats_convergence_points.txt tree.%.mrbayes.mcmc_convergence_stats_convergence_points.txt
	@beast_val="$$(tail -n 1 tree.$*.beast.mcmc_convergence_stats_convergence_points.txt | cut -f2)"; \
	mb_val="$$(tail -n 1 tree.$*.mrbayes.mcmc_convergence_stats_convergence_points.txt | cut -f2)"; \
	printf "beast\tmrbayes\n%s\t%s\n" "$$beast_val" "$$mb_val" > tree.$*.mcmc_convergence.txt

eval.all.convergence.txt: $(TREE_CONVERGENCE)
	echo -e "samp\tbeast\tmrbayes" > tmp
	for file in $^ ; do \
		samp="$${file#tree.}"; samp="$${samp%%.*}"; \
		tail -n 1 $${file} | awk -v s="$$samp" '{print s "\t" $$1 "\t" $$2}' >> tmp ;\
	done
	awk 'NR>1 {b+=$$2; m+=$$3; n++} END {printf "avg\t%.0f\t%.0f\n", b/n, m/n}' tmp >> tmp
	mv tmp $@

clean:
	rm -rf $(TREES) $(FA) tree.*.beast.* tree.*.mrbayes.* *.mcmc_convergence_stats* tree.*.nex

clean_beast:
	rm -rf tree.*.beast.*
	rm eval.all*

clean_mrbayes:
	rm -rf tree.*.mrbayes.*
	rm eval.all*

clean_stats:
	rm -rf *.mcmc_convergence_stats*

archive_beast:
	archive_dir="archive.beast_$$(date +%Y-%m-%d_%H:%M:%S)"; \
	mkdir $$archive_dir; \
	mv tree.*.beast* $$archive_dir/; \
	mv eval.all* $$archive_dir/

# Pilot runs without mcmc thinning use a lot of storage space, so this target removes log files from mcmc methods to free up space
release_storage:
	rm tree.*.beast.*.trees tree.*.beast.*.log tree.*.mrbayes.*.p tree.*.mrbayes.*.t


