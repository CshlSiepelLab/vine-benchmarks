# Centering of VINE's pairwise-distance estimates

Analysis behind the centering supplement figure. Everything needed to reproduce it lives
in this directory (see `README.md`); this file records the findings and the reasoning.

## Reviewer comment being addressed
> "The 95% CI inclusion analysis in Fig 3B suggests under-coverage, and Figs S2/S7 report
> worse RF distances of VINE. In these cases, is VINE able to recover the posterior mean
> correctly with narrower-than-true intervals, or also mis-centered? If mis-centered, the
> topologies recovered by VINE could be mostly different from the true distributions.
> Therefore, further reporting of biases of VINE's posterior mean relative to the
> simulation ground-truth could reflect this issue."

Scope here = **centering** (bias of the posterior mean vs. ground truth). Interval
width / coverage is a separate topic handled elsewhere in the paper (with `-v` settings
other than 0), so it is only summarized below.

## Answer (short)
VINE's posterior mean is **well-centered, not mis-centered.** Terminal branches and
short/medium distances are essentially unbiased (on par with BEAST2/MrBayes). The only
systematic effect is a **mild underestimation of the longest distances**, which traces to
shrinkage of **deep internal branches**. VINE still recovers ~85% of true clades
(BEAST2/MrBayes ~86–90%), so its topologies are mostly the correct ones — not "mostly
different from the true distribution." This bias is essentially invariant to the variance
regularizer, so it is a genuine point-estimate property, separate from the coverage issue.

## Key numbers (reproduce with `analysis/centering_stats.py`)

Calibration fit of estimated on true pairwise distance, `est = a + b·true`
(b < 1 ⇒ shrinkage of large distances):

| taxa | VINE slope | BEAST2 slope | MrBayes slope |
|---|---|---|---|
| 25  | 0.85 | 0.95 | 0.97 |
| 50  | 0.89 | 1.01 | 1.04 |
| 100 | 0.87 | 1.02 | 1.04 |

- VINE intercept is positive (+0.10 to +0.12): small distances slightly up, large
  distances down — classic shrinkage toward the middle.
- 30–50% of VINE's total MSE is removable by a single global linear recalibration
  (vs 1–26% for BEAST2/MrBayes) ⇒ most of the "error" is a scale effect, not scatter.
- R² is scale-invariant, so it stays comparable (~0.91 vs ~0.93–0.98); that is why RF and
  likelihood look fine while a raw RMSE looks worse than the centering warrants.
- Top-decile shrinkage `mean(est)/mean(true)`: VINE 0.82–0.92, BEAST2 0.92–1.02,
  MrBayes 0.95–1.05. BEAST2/MrBayes also underestimate the very longest distances, just
  less — an information limit at saturation (300 sites), same *kind* of effect.

## Branch-level mechanism (reproduce with `analysis/branch_level_bias.py`, needs dendropy)
Branches matched by bipartition between each method's tree sample and the true tree;
estimate = mean branch length over the sample (conditional on the clade being present).
Aggregate relative bias Σ(est−true)/Σtrue:

| taxa | VINE terminal | VINE internal | BEAST2/MrBayes (both) |
|---|---|---|---|
| 10 | −0.033 | **−0.108** | terminal ≈ 0, internal −0.02…+0.03 |
| 25 | −0.029 | **−0.088** | " |
| 50 | +0.004 | **−0.031** | " |

- Terminal branches: VINE ≈ unbiased. Internal (deep) branches: VINE shrinks them
  (more on smaller trees).
- Very short branches (<~0.05): all methods over-estimate (a floor/regularization
  effect); VINE is the mildest.
- VINE leaves ~8–10% of true internal splits fully absent from its sample vs ~0% for
  BEAST2/MrBayes (occasional dropped deep splits).

**Why it produces the long-distance bias:** a patristic distance between distant taxa
crosses many internal edges, so VINE's per-edge downward bias on deep branches
accumulates along long paths. Short (terminal-dominated) distances are unaffected. Likely
cause: the variational posterior shrinks weakly-constrained deep edges (cf. memory
`vine-prior-effect`).

## Coverage decomposition (context only; `analysis/ci_coverage_decomposition.py`)
The Fig 3B under-coverage is driven mostly by **too-narrow intervals**, not mis-centering:
VINE cov95 ≈ 0.05–0.11, median per-pair SD ≈ 0.012 vs BEAST2 ≈ 0.084 (~7× narrower). This
follows from the benchmark's `-v 0`, which disables VINE's variance-collapse regularizer
(`vine --var-reg` help: "A value of zero will disable regularization"). Sweeping
`-v 0/1/2/5` (`analysis/variance_reg_experiment.sh`) moves cov95 0.07→0.60 while the mean
error barely changes (−0.06→−0.10, slope 0.86–0.94) — confirming centering is independent
of the variance setting. Interval width is discussed separately in the main paper.

## Data provenance
- `eval.all.dist.txt` (in the sim dirs) holds only per-replicate summaries. Per-pair
  estimates are not retained by the sim Makefile, but the true pairwise distances
  (`tree.R.true.dist.txt`) and the sampled tree files (`tree.R.{var,beast,mrbayes}.nwk`)
  are. `extractDistData.sh` reruns `evalTrees` on those samples and joins to truth →
  `data/calibration_data.csv`. Method keys: `var`→vine, `beast`→beast, `mrbayes`→mrbayes.
