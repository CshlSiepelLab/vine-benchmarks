#!/usr/bin/env Rscript
# Accuracy of VINE's posterior approximations on the
# simulated CRISPR data (n = 10 ... 1000), against the true trees and BEAM.
#   (A) 95% CI inclusion of the true pairwise distances (Vine, Vine --posterior, BEAM)
#   (B) pairwise-distance KL divergence from BEAM
#   (C) per-split topology KL divergence from BEAM
# For (A) and (B) VINE's branch lengths (expected edits per site) are
# multiplied by 100 to put them in days like the true and BEAM trees
# (crispr_sims/posterior_quality.mk, VINE_DIST_SCALE).
# Reads <method>.eval.all.{dist,kl_pdist_beam,kl_topo_refnorm_beam}.txt
# from each crispr_sims/<n>taxa folder (`make posterior-quality`).
#
# Usage: Rscript makeGraphs.R [crispr_root]   (default: ../../crispr_sims)
# Panel C uses compareTrees' Mean_split_KL_refnorm (normalized by the number
# of distinct splits in BEAM's sample).
# Output: plots/crispr_posterior_quality.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args)) args[1] else
  file.path(dirname(dirname(script_dir)), "crispr_sims")
out_dir <- file.path(script_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxa <- c(10, 25, 50, 100, 250, 500, 1000)
method_levels <- c("Vine", "Vine --posterior", "BEAM")
.
method_colors <- c(
  "Vine" = "#F28E2B",
  "Vine --posterior" = "#76B7B2",
  "BEAM" = "#59A14F"
)
method_prefix <- c("Vine" = "vine", "Vine --posterior" = "vine_twostage",
                   "BEAM" = "beam")

read_replicates <- function(path) {
  lines <- readLines(path, warn = FALSE)
  separator <- grep("^-+$", trimws(lines))
  if (length(separator)) lines <- lines[seq_len(separator[1] - 1)]
  read.table(text = lines, header = TRUE, check.names = FALSE)
}

mean_sd <- function(x) {
  x <- x[is.finite(x)]
  c(mean = mean(x), sd = if (length(x) > 1) sd(x) else 0)
}

collect <- function(methods, table, column) {
  d <- do.call(rbind, lapply(taxa, function(n) {
    do.call(rbind, lapply(methods, function(method) {
      path <- file.path(data_root, paste0(n, "taxa"),
                        paste0(method_prefix[[method]], ".eval.all.", table, ".txt"))
      z <- mean_sd(read_replicates(path)[[column]])
      data.frame(ntaxa = n, method = method, mean = z[["mean"]], sd = z[["sd"]])
    }))
  }))
  d$ntaxa <- factor(d$ntaxa, levels = taxa)
  d$method <- factor(d$method, levels = method_levels)
  d
}

paper_theme <- theme_minimal(base_size = 8) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray80", linewidth = 0.2),
    plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
    plot.tag.position = c(-0.02, 1.02),
    plot.margin = margin(8, 10, 5, 10)
  )

bar_plot <- function(d, y_label, upper = NA) {
  ggplot(d, aes(x = ntaxa, y = mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8, show.legend = TRUE) +
    geom_errorbar(
      aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
      position = position_dodge(width = 0.9), width = 0.18, linewidth = 0.35
    ) +
    scale_fill_manual(values = method_colors, breaks = method_levels,
                      limits = method_levels, drop = FALSE) +
    scale_y_continuous(limits = c(0, upper), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of Taxa (n)", y = y_label, fill = "Method") +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    paper_theme
}

vine_methods <- c("Vine", "Vine --posterior")
panel_a <- bar_plot(collect(method_levels, "dist", "CI95"), "95% CI Inclusion")
panel_b <- bar_plot(collect(vine_methods, "kl_pdist_beam", "kl_pdist"),
                    "Pairwise-distance KL divergence\nfrom BEAM")
panel_c <- bar_plot(collect(vine_methods, "kl_topo_refnorm_beam", "kl_topo_refnorm"),
                    "Per-split topology KL divergence\nfrom BEAM")

figure <- (panel_a + panel_b + panel_c) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical")

out_file <- file.path(out_dir, "crispr_posterior_quality.pdf")
ggsave(out_file, figure, width = 10, height = 3.2, units = "in", device = cairo_pdf)
message("Wrote ", out_file)
