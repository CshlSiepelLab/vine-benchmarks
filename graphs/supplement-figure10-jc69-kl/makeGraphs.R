#!/usr/bin/env Rscript
# Supplementary Figure S11: KL divergence of posterior approximations from the
# MCMC-based reference distributions on the JC69 300-site simulations
# (n = 10, 15, 20), including the VI methods Dodonaphy and GeoPhy. Same layout
# as Supplementary Fig. S10 (dna_sims/hky_300sites_posterior_quality/
# plot_posterior_quality_sup1.R):
#   (A) pairwise-distance KL from MrBayes   (B) pairwise-distance KL from BEAST2
#   (C) per-split topology KL from MrBayes  (D) per-split topology KL from BEAST2
# Reads <method>.eval.all.kl_{pdist,topo_refnorm}_{mrbayes,beast}.txt, built by
# `make kl-metrics` (dna_sims/jc69_300sites/kl_metrics.mk) in each size folder.
# Usage: Rscript makeGraphs.R [jc69_root]   (default: ../../dna_sims/jc69_300sites)
# Panels C-D use compareTrees' Mean_split_KL_refnorm (normalized by the number
# of distinct splits in the reference sample).
# Output: plots/jc69_kl_divergence.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args)) args[1] else
  file.path(dirname(dirname(script_dir)), "dna_sims", "jc69_300sites")
out_dir <- file.path(script_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxa <- c(10, 15, 20)
method_levels <- c("Vine", "Vine --posterior", "Dodonaphy", "GeoPhy",
                   "BEAST2", "MrBayes")
# Vine/BEAST2/MrBayes/two-stage colours match Fig. 3 and Supp. Figs. S8/S10;
# Dodonaphy/GeoPhy match the JC69 panels of Fig. 2.
method_colors <- c(
  "Vine" = "#F28E2B",
  "Vine --posterior" = "#76B7B2",
  "Dodonaphy" = "#B07AA1",
  "GeoPhy" = "#EDC948",
  "BEAST2" = "#59A14F",
  "MrBayes" = "#E15759"
)
method_prefix <- c(
  "Vine" = "vine", "Vine --posterior" = "vine_twostage",
  "Dodonaphy" = "dodonaphy", "GeoPhy" = "geophy",
  "BEAST2" = "beast", "MrBayes" = "mrbayes"
)

# Use only the replicate rows (drop the "-----" separator and mean row) so
# bars are replicate means and error bars sample SDs, as in Supp. Fig. S10.
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

collect_kl <- function(metric, reference) {
  column <- paste0("kl_", metric)
  reference_method <- if (reference == "mrbayes") "MrBayes" else "BEAST2"
  rows <- lapply(taxa, function(n) {
    do.call(rbind, lapply(method_levels, function(method) {
      if (method == reference_method) {
        # A reference compared with itself has KL exactly 0.
        return(data.frame(ntaxa = n, method = method, mean = 0, sd = 0))
      }
      path <- file.path(data_root, paste0(n, "taxa"),
                        paste0(method_prefix[[method]], ".eval.all.",
                               column, "_", reference, ".txt"))
      z <- mean_sd(read_replicates(path)[[column]])
      data.frame(ntaxa = n, method = method, mean = z[["mean"]], sd = z[["sd"]])
    }))
  })
  d <- do.call(rbind, rows)
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

kl_plot <- function(d, y_label) {
  ggplot(d, aes(x = ntaxa, y = mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    geom_errorbar(
      aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
      position = position_dodge(width = 0.9), width = 0.18, linewidth = 0.35
    ) +
    scale_fill_manual(values = method_colors, breaks = method_levels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of Taxa (n)", y = y_label, fill = "Method") +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    paper_theme
}

panel_a <- kl_plot(collect_kl("pdist", "mrbayes"),
                   "Pairwise-distance KL divergence\nfrom MrBayes")
panel_b <- kl_plot(collect_kl("pdist", "beast"),
                   "Pairwise-distance KL divergence\nfrom BEAST2")
panel_c <- kl_plot(collect_kl("topo_refnorm", "mrbayes"),
                   "Per-split topology KL divergence\nfrom MrBayes")
panel_d <- kl_plot(collect_kl("topo_refnorm", "beast"),
                   "Per-split topology KL divergence\nfrom BEAST2")

figure <- ((panel_a + panel_b) / (panel_c + panel_d)) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical")

out_file <- file.path(out_dir, "jc69_kl_divergence.pdf")
ggsave(out_file, figure, width = 8.5, height = 6.5, units = "in",
       device = cairo_pdf)
message("Wrote ", out_file)
