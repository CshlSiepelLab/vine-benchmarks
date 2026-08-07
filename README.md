# Vine-Benchmarks

This *vine-benchmarks* repository contains the instructions and scripts for reproducing the data processing, benchmark simulations and figure generation for the manuscript:

Siepel, A., Hassett, R., & Staklinski, S. J. (2026). VINE: Variational inference for scalable Bayesian reconstruction of species and cell-lineage phylogenies. *bioRxiv*. https://doi.org/10.64898/2025.12.24.696405 (under revision for Nature Methods)

## Dependencies

### Install Phast and Vine

Choose **one** of the two options below.

#### Option A: Install via Conda

```sh
conda install -c conda-forge -c bioconda vine-phylo
```

Note: The steps below still expect a `bin` directory at the repository root for installing BEAST, MrBayes, and the helper Python scripts, so create one if it doesn't already exist:
```sh
mkdir bin
```

The Makefiles throughout this repository expect Vine and Phast binaries at `bin/vine/bin` and `bin/phast/bin` (the layout produced by Option B). Rather than editing every Makefile's `VINE_BIN`/`PHAST_BIN` paths, symlink the conda-installed binaries into that layout instead:
```sh
mkdir -p bin/vine/bin bin/phast/bin
for f in vine evalTrees; do ln -s "$(command -v $f)" bin/vine/bin/; done
for f in phyloFit tree_doctor base_evolve; do ln -s "$(command -v $f)" bin/phast/bin/; done
```

#### Option B: Build from Source

Create a binary directory:
```sh
mkdir bin
cd bin
```

Install Phast in the binary directory:
```sh
git clone git@github.com:CshlSiepelLab/phast.git
cd phast
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$(pwd)"
cmake --build build
cmake --install build
```

Install Vine in the binary directory:
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

### Install BEAST in Binary Directory

Follow instructions here: https://www.beast2.org/

```sh
wget https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.x86.tgz
tar -zxvf BEAST.v2.7.7.Linux.x86.tgz
./beast/bin/packagemanager -list
./beast/bin/packagemanager -add feast
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
pip install -r requirements.txt
for f in addMrbayesModelToNex bdTree3 convertTrees fa2nex nex2nwk nwk2nex nex2fa; do
pyinstaller --onefile --distpath ../../bin --workpath ../../build --specpath ../../build "$f.py"
done
```

## Run Simulations under HKY Model for 300 site alignments

In `dna_sims/hky_300sites/common.mk` update `MAIN_DIR` to the absolute path of the cloned vine_benchmarks repository

Run make commands in common.mk

Examples:

To run vine, beast, mrbayes for the 10 taxa and 300 site alignments and to save the likelihoods and run times
```sh
cd 10taxa
make eval.all.lnl.txt
make eval.all.time.txt
```

### Run Simulations under JC69 Model for 300 site alignments

In `dna_sims/jc69_300sites/common.mk` update `MAIN_DIR` to the absolute path of the cloned vine_benchmarks repository

To run simulations for Dodonaphy, Geophy and Vaiphy models, you must build their singularity containers *.sif files in `containers/` from the existing *.def files

Run make commands in common.mk

## Plot Manuscript Figures

The R dependencies needed to run the plotting scripts below are listed in `graphs/environment.yaml`. Create and activate the conda environment with:
```sh
conda env create -f environment.yaml -n r_env
conda activate r_env
```

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