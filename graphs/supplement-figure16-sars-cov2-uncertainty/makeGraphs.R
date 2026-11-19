#!/usr/bin/env Rscript
# Posterior uncertainty for the 364-taxon SARS-CoV-2 data set, BEAST 2 (x) vs Vine --posterior (y).
#   (A) Clade (split) posterior probabilities, for the splits supported by >= 1
#       expected substitution under BEAST 2 (posterior probability x mean
#       subtending branch length x alignment length), coloured by that number.
#       From splits.tsv (splits.py).
#   (B) Relative width of the 95% credible interval (CI width / posterior mean)
#       of every pairwise patristic distance, coloured by point density.
#       From {beast,vine_twostage}-small.x1000.dist (evalTrees on the Figure 5
#       trees with branch lengths x1000, for enough printed precision).
# Both use 1,000 Vine --posterior samples and the post-burn-in BEAST 2 samples.
# Usage: Rscript makeGraphs.R     (run via `make` in this folder)
# Output: plots/sars_cov2_posterior_uncertainty.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(MASS)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
fig5_dir <- file.path(dirname(script_dir), "figure5-sars-cov2")
out_dir <- file.path(script_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_support <- 1   # expected supporting substitutions under BEAST 2

# ---- (A) clade posterior probabilities ----
splits <- read.delim(file.path(script_dir, "splits.tsv"))
clades <- subset(splits, ref_support >= min_support)
clades <- clades[order(clades$ref_support), ]   # strongest drawn on top
message(nrow(clades), " splits with >= ", min_support,
        " expected supporting substitution(s) under BEAST 2")

paper_theme <- theme_classic(base_size = 7) +
  theme(
    plot.tag = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.height = unit(0.5, "lines"),
    legend.key.width = unit(0.5, "lines"),
    legend.margin = margin(1, 1, 1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA)
  )

panel_a <- ggplot(clades, aes(x = ref_pp, y = est_pp, color = ref_support)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray40") +
  geom_point(size = 1.1, alpha = 0.9) +
  scale_color_viridis_c(option = "viridis", trans = "log10", name = "Supporting\nsubstitutions") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
  coord_equal() +
  labs(x = "BEAST 2 posterior probability", y = "Vine --posterior posterior probability") +
  paper_theme +
  theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1))

# ---- (B) relative 95% CI width of pairwise distances ----
read_dist <- function(path) {
  d <- read.table(path, header = FALSE, comment.char = "#",
                  col.names = c("leaf1", "leaf2", "mean", "std", "med", "min", "max",
                                "low95", "high95", "low50", "high50"))
  d$key <- paste(d$leaf1, d$leaf2)
  d$rel_width <- (d$high95 - d$low95) / d$mean
  d[is.finite(d$rel_width), c("key", "rel_width")]
}
beast <- read_dist(file.path(script_dir, "beast-small.x1000.dist"))
vine <- read_dist(file.path(script_dir, "vine_twostage-small.x1000.dist"))
pairs <- merge(beast, vine, by = "key", suffixes = c("_beast", "_vine"))
message(nrow(pairs), " taxon pairs")

# point density from a 2-D kernel density estimate (as ggpointdensity's kde2d
# method), so the panel needs only MASS
dens <- kde2d(pairs$rel_width_beast, pairs$rel_width_vine, n = 300)
ix <- findInterval(pairs$rel_width_beast, dens$x, all.inside = TRUE)
iy <- findInterval(pairs$rel_width_vine, dens$y, all.inside = TRUE)
pairs$density <- dens$z[cbind(ix, iy)]
pairs <- pairs[order(pairs$density), ]
lim <- c(0, max(c(pairs$rel_width_beast, pairs$rel_width_vine)))

panel_b <- ggplot(pairs, aes(x = rel_width_beast, y = rel_width_vine, color = density)) +
  geom_point(size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray40") +
  scale_color_viridis_c(option = "magma", name = "Point density", labels = NULL) +
  scale_x_continuous(limits = lim, expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = lim, expand = expansion(mult = 0.02)) +
  coord_equal() +
  labs(x = "BEAST 2 relative 95% CI width", y = "Vine --posterior relative 95% CI width") +
  paper_theme +
  theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1))

figure <- (panel_a + panel_b) + plot_annotation(tag_levels = "A")
out_file <- file.path(out_dir, "sars_cov2_posterior_uncertainty.pdf")
ggsave(out_file, figure, width = 6, height = 3, units = "in", device = cairo_pdf)
message("Wrote ", out_file)
