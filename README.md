# Vine-Benchmarks

This *vine-benchmarks* repository contains the instructions and scripts for reproducing the data processing, benchmark simulations and figure generation for the manuscript available at https://pmc.ncbi.nlm.nih.gov/articles/PMC13042005/

## Dependencies

### Create Binary Directory

```sh
mkdir bin
cd bin
```

### Install Phast in Binary Directory

```sh
git clone git@github.com:CshlSiepelLab/phast.git
cd phast
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$(pwd)"
cmake --build build
cmake --install build
```

### Install Vine in Binary Directory

```sh
git clone git@github.com:CshlSiepelLab/vine.git
cd vine
cd bin/vine
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPHAST_ROOT=/path/to/bin/phast \
  -DCMAKE_INSTALL_PREFIX="$(pwd)"
cmake --build build
cmake --install build

```
Note: Can also install vine+phast with the below commands:
```sh
conda install -c conda-forge -c bioconda vine-phylo
```

### Install BEAST in Binary Directory

Follow instructions here: https://www.beast2.org/

```sh
wget https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.x86.tgz
tar -zxvf BEAST.v2.7.7.Linux.x86.tgz
```

### Install MrBayes in Binary Directory

Follow instructions here: https://github.com/NBISweden/MrBayes/blob/develop/INSTALL

```sh
git clone --depth=1 https://github.com/NBISweden/MrBayes.git
./configure
make
make install-exec bindir=/path/to/bin
```

### Install Helper Python Scripts

```sh
cd python/src
for f in addMrbayesModelToNex bdTree3 convertTrees fa2nex nex2nwk nwk2nex nex2fa; do
pyinstaller --onefile --distpath ../../bin --workpath ../../build --specpath ../../build "$f.py"
done
```

## Run Simulations under HKY Model for 300 site alignments

In dna_sims/hky_300sites/common.mk update `MAIN_DIR` to the absolute path of the cloned vine_benchmarks repository

Run make commands in common.mk

Examples:

To run vine, beast, mrbayes for the 10 taxa and 300 site alignments and to save the likelihoods and run times
```sh
cd 10taxa
make eval.all.lnl.txt
make eval.all.time.txt
```

## Plot Manuscript Figures

Intermediate data files for manuscript figures are already stored in the `graphs` directory, in subdirectories/files ending in `-data`, since regenerating all files from scratch using the necessary MCMC and ML methods is time-consuming. Note: If you have rerun simulations, extract results across taxa into the corresponding `graphs` `-data` location for plotting using the `graphs/extractGraphData.sh` script. Options for which simulation data to extract for plotting are described in `graphs/extractGraphData.sh`

For figure 2 run:
```sh
Rscript makeGraphs.R
```

For figure 3 run:
```sh
Rscript plot_posterior_quality_main.R
```

For figure 4 run:
```sh
Rscript makeGraphs.R
```
Note: Repository already contains the summary result files so plots can be generated without rerunning all of the crispr simulations. If all crispr simulations are completed, extract and summarize them across taxa into the graph folder with `graphs/figure4-crispr-testdata/extractGraphData.sh` after setting the `ROOT` variable at the top of the file

For figure 5:
Follow instructions in `dna_real_data/sars-cov2/README-download-instructions.txt` for downloading the nextstrain sars-cov2 dataset. Run experiments using the `dna_real_data/sars-cov2/Makefile` commands

In `graphs/figure5-sars-cov2/Makefile`, set `RSCRIPT` variable path and run make commands for building figure 5 panel plots