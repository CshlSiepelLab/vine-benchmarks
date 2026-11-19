# Supplement figure: Taylor-approximated ELBO vs. Monte Carlo

Addresses the reviewer request to compare the second-order Taylor approximation of
the ELBO against a fully sampled Monte Carlo estimate, at the fitted variational
optimum.

## Figure layout (4 panels, 2 per taxon count)

- **Top row (A: 25 taxa, B: 50 taxa)** — ELBO convergence trajectories for three
  replicate trees (1-3), each shown as a Taylor curve and a Monte Carlo curve
  (colored by method). The replicates converge to different optima and so separate
  into distinct bands; within every band the Taylor and 100-sample Monte Carlo
  curves converge at a similar rate to essentially the same optimum, and at the
  full scale of the ELBO the two are indistinguishable near convergence. Showing
  three trees (rather than one) guards against cherry-picking. (These are two
  separate optimization runs from the same initialization — a comparison of
  optimization dynamics, not a per-iteration same-(µ,Σ) estimator comparison.)

  The **gray "subsampled warm-up" band (iterations 1-40)** marks the phase in
  which the SGA optimizer uses site subsampling, which makes the ELBO estimate
  noticeably noisier. The scheduler grows the batch deterministically
  (256 sites -> 278 at iter 20 -> full 300 at iter 40, when the batch fraction
  crosses `tau_full = 0.95`), so subsampling switches off at iteration 40 in every
  run; the step-to-step ELBO jitter drops ~13x right there. After the band the two
  estimators are smooth and essentially superimposed.

  Note: in the original logs the `subsamp` column appeared stuck at 278 for the
  whole run — a logging bug (the field was only written on subsampling steps, so
  it went stale once full-batch mode kicked in). Fixed in the `vine` repo on
  branch `fix-subsamp-logging`; it is a logging-only change and does not affect any
  ELBO values here.
- **Bottom row (C: 25 taxa, D: 50 taxa)** — the ELBO at the maximum, Taylor − Monte
  Carlo at the *same* converged (µ,Σ), one point per replicate (n = 10), with the
  mean ± SD band. Same quantity as the top row, zoomed ~1000× so the residual is
  visible.

## What is compared

For each converged variational fit we compare, **at the identical (µ, Σ)**:

- `ELBO_taylor` — the Taylor/hybrid estimate VINE actually optimizes
  (`ll_µ + ½·tr(HΣ) − elbo_bias`, minus the analytic KLD). This is the `ELB`
  value on the `# Reverting …` line of the log.
- `ELBO_mc` — a **100-sample** Monte Carlo estimate of E_q[log p(x|z)] at the same
  parameters (the `LNL_mc` value on the same line), minus the same analytic KLD.

Because the KLD term is analytic and identical for both,
`ELBO_taylor − ELBO_mc = LNL − LNL_mc`, an apples-to-apples comparison of the two
estimators of the expected log likelihood at one point — not two separate
optimizations.

## Data provenance

The Taylor runs in `../../dna_sims/hky_300sites/{25taxa,50taxa}` used the default
`--batchsize 10`, so their internal `LNL_mc` was only a 10-sample estimate. To get
a fully sampled comparison the Taylor jobs were re-run with `--batchsize 100`
into sibling directories `taylor-25taxa` / `taylor-50taxa` (Makefiles identical to
the originals apart from `--batchsize 100`; same binary as the Monte Carlo runs,
`bin/vine/bin/vine`, v0.3.5). `extractElboData.sh` parses the `# Reverting …`
lines of those logs into `data/elbo_max.csv`.

Ten replicate datasets per size (25 and 50 taxa).

## Result

| size | Taylor − MC (mean ± SD) | mean \|Δ\| | max \|Δ\| | relative |
|------|-------------------------|-----------|----------|----------|
| 25 taxa | −1.22 ± 1.40 | 1.63 | 3.13 | 0.025% |
| 50 taxa | −2.17 ± 3.07 | 2.84 | 6.85 | 0.024% |

On ELBOs of 6,000–13,800, the Taylor estimate agrees with the 100-sample Monte
Carlo estimate to ~1–3 log-likelihood units (~0.02%). The small offset is
consistently **negative** — the Taylor estimate is mildly conservative (slightly
*under*-estimates E_q[lnL]), the opposite of what under-weighting low-likelihood
tail contributions would produce.

## Reproduce

```bash
./extractElboData.sh        # -> data/elbo_max.csv   (bottom-row points, all 10 reps)
python3 extractTrajData.py  # -> data/elbo_traj.csv  (top-row trajectories, rep 1)
/home/asiepel/miniforge3/envs/rcentering/bin/Rscript makeGraphs.R   # -> elbo_taylor_vs_mc.pdf/.png
```

To use a different representative replicate for the trajectory panels, edit `REP`
in `extractTrajData.py`.
