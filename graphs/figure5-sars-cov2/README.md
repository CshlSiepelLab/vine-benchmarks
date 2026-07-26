# Figure 5 — SARS-CoV-2 (VINE vs BEAST 2)

Generates the five panels of Figure 5 comparing VINE to BEAST 2 on a
SARS-CoV-2 data set. **This directory is figure generation only**; the data
analysis (VINE and BEAST 2 inference, raw posteriors) lives in the analysis
directory `../../dna_real_data/sars-cov2`. Nothing large is duplicated here:
raw posteriors are symlinked in from there at build time, and only the small
derived summaries the plots need are kept.

## Panels

| Output | Script | Shows |
|---|---|---|
| `vine-vs-beast.densitree.pdf` | `densitree-compare.R` | VINE vs BEAST 2 posterior tree clouds (364 taxa) |
| `tanglegram.compact.pdf` | `tanglegram.R` | MCC-tree tanglegram, VINE vs BEAST 2 |
| `scatter.pdf` | `scatter-plot.R` | pairwise patristic distances, VINE vs BEAST 2 |
| `vine-large.densitree.pdf` | `densitree-large.R` | VINE posterior tree cloud (1060 taxa) |
| `embedding-pca.pdf` | `embedding-pca.R` | PCA of the VINE latent embedding (1060 taxa) |

`tanglegram.R` also writes the full (non-compact) `vine-vs-beast.tanglegram.pdf`.

## Which VINE run

All panels use the **default** VINE run (`vine-small.nwk`, `vine-large.nwk` in
the analysis dir). Its point estimates match BEAST 2 tightly. The calibrated
`-v 1` run (`-S DIST -F -Z -v 1`) trades a little pairwise-distance accuracy for
wider, better-calibrated credible intervals and is used instead for the
**posterior-uncertainty supplement**, generated separately under
`../../dna_real_data/sars-cov2/graphs/uncertainty`. The `-v 1` renderings of
these same panels are kept for comparison in `archive.flows-v1-panels/`.

## Reproducing

```sh
make all
```

`make -n all` reports nothing to do when the panels are current. Individual
panels are their own targets (e.g. `make scatter.pdf`).

Build flow: raw posteriors are symlinked from the analysis dir; `make` converts
them to nexus, builds MCC trees (`treeannotator`) and pairwise-distance tables
(`evalTrees`), and renders each panel. The large full-posterior `.nex` files are
`.SECONDARY` intermediates — rebuilt on demand and not retained, so a re-render
with the panels already present does not regenerate them.

## R environment

Needs R with `ape`, `ggtree`, `treeio`, `phangorn`, `phytools`, `ggrastr`,
`ggpointdensity`, `uwot`, `tidyverse`, `cowplot`. `RSCRIPT` in the Makefile
points at a conda env (`rfig5`) providing these; override it for your setup:

```sh
make all RSCRIPT="/path/to/Rscript"
```

**Fonts:** panels are drawn with `cairo_pdf` and the default `sans` family. The
Makefile sets `FONTCONFIG_PATH=/etc/fonts` so `sans` resolves to Nimbus Sans
(the URW Helvetica clone, metrically identical to Helvetica). Without it, a
conda R env may embed DejaVu Sans instead. Every panel embeds
`NimbusSans-Regular`/`-Bold`.

## Layout

```
*.R                     plotting scripts (one per panel; embedding-umap.R is an
                        unused UMAP alternative to embedding-pca.R)
Makefile                build rules
python/time2subs.py     BEAST time-tree -> substitutions helper
*.pdf / *.png           rendered panels
*.mcc.nex, *.dist       small cached derived summaries
archive.flows-v1-panels/  the -v 1 renderings, for comparison
archive.preprint/       frozen preprint-version scripts (bulk data stripped)
```

The parameter exploration behind the default-vs-`-v 1` choice (regularization
sweep, seed replicates, convergence tests) is archived in
`../../dna_real_data/sars-cov2/archive.parameter-exploration.2026-07-24/` and
documented in `../../dna_real_data/sars-cov2/graphs/uncertainty/PROVENANCE.md`.
