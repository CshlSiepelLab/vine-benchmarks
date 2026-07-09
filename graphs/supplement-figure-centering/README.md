# Supplement figure: centering of VINE's pairwise-distance estimates

Addresses the reviewer question: *is VINE's posterior mean recovered correctly (just with
narrow intervals) or is it mis-centered?* Answer: **well-centered, not mis-centered** — the
posterior-mean pairwise distances track the ground truth through the bulk of the range, with
only a mild, systematic underestimation of the longest distances (which trace back to
shrinkage of deep internal branches). See `VINE_DIST_BIAS_ANALYSIS.md` for the full write-up.

## Main figure
`dist_centering.pdf` / `.png` — three calibration panels (posterior-mean pairwise
distance vs. true distance) for **25 / 50 / 100 taxa** (panels A / B / C). Style, palette,
fonts, panel sizes, and A/B/C callouts match `../figure2-hky300/makeGraphs.R`
(Vine `#F28E2B`, BEAST2 `#59A14F`, MrBayes `#E15759`).

## How to build (two steps)
```bash
# 1. Extract data: reruns evalTrees on the sampled trees and joins to the ground truth.
./extractDistData.sh                 #  -> data/calibration_data.csv
# 2. Plot (needs R with ggplot2, patchwork, scales, dplyr — same as figure 2):
Rscript makeGraphs.R                 #  -> dist_centering.pdf/.png
```
Absolute lab paths (sims root, evalTrees/vine binaries) are set once at the top of
`extractDistData.sh`; edit there if running elsewhere.

## Files
| File | Purpose |
|---|---|
| `extractDistData.sh` | Reruns `evalTrees` on `tree.<rep>.{var,beast,mrbayes}.nwk`, writes `data/raw/*.dist`, then calls the join. |
| `buildCalibrationTable.py` | Joins per-pair estimates to `tree.<rep>.true.dist.txt` → `data/calibration_data.csv` (cols: ntaxa,rep,method,leaf1,leaf2,true,est). |
| `makeGraphs.R` | The figure (Fig-2 style). |
| `data/calibration_data.csv` | Tidy plotting table (committed). |
| `data/raw/` | Intermediate evalTrees output (regenerable). |
| `analysis/centering_stats.py` | Reproduces the centering numbers (slopes, per-bin bias, top-decile shrinkage) from the CSV. |
| `analysis/branch_level_bias.py` | Bipartition-matched branch bias, terminal vs internal (needs dendropy). |
| `analysis/ci_coverage_decomposition.py` | Splits 95% CI under-coverage into centering vs width (reruns evalTrees). Coverage is a *separate* paper topic. |
| `analysis/variance_reg_experiment.sh` | Reruns VINE across `-v` settings; shows centering is ~invariant to variance regularization. |
| `VINE_DIST_BIAS_ANALYSIS.md` | Full analysis log / findings. |

## Upstream
The tree samples (`tree.<rep>.{var,beast,mrbayes}.nwk`) and true trees are produced by the
simulation Makefile under `$SIMS_ROOT` (see that `Makefile` and `../../dna_sims/.../common.mk`);
this directory only re-derives distances and plots.
