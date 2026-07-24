#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
root_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
out_dir <- file.path(root_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxa <- c(10, 25, 50, 100, 250, 500, 1000)
method_levels <- c("Vine", "Vine + flows", "BEAST2", "MrBayes")
method_colors <- c(
  "Vine" = "#F28E2B",
  "Vine + flows" = "#76B7B2",
  "BEAST2" = "#59A14F",
  "MrBayes" = "#E15759"
)

# Ignore each summary's dashed separator and precomputed final mean row; panel
# bars and error bars use the replicate mean and sample SD, as in the main plot.
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
  suffix <- paste0("eval.all.", column, "_", reference, ".txt")
  files <- if (reference == "mrbayes") {
    c("Vine" = paste0("vine.", suffix),
      "Vine + flows" = paste0("vine_flows.", suffix),
      "BEAST2" = paste0("beast.", suffix))
  } else {
    c("Vine" = paste0("vine.", suffix),
      "Vine + flows" = paste0("vine_flows.", suffix),
      "MrBayes" = paste0("mrbayes.", suffix))
  }
  reference_method <- if (reference == "mrbayes") "MrBayes" else "BEAST2"

  rows <- lapply(taxa, function(n) {
    observed <- lapply(names(files), function(method) {
      d <- read_replicates(file.path(root_dir, paste0(n, "taxa"), files[[method]]))
      z <- mean_sd(d[[column]])
      data.frame(ntaxa = n, method = method,
                 mean = z[["mean"]], sd = z[["sd"]])
    })
    observed[[length(observed) + 1]] <- data.frame(
      ntaxa = n, method = reference_method, mean = 0, sd = 0
    )
    do.call(rbind, observed)
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
    scale_fill_manual(values = method_colors, breaks = method_levels,
                      drop = FALSE) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of Taxa (n)", y = y_label, fill = "Method") +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    paper_theme
}

panel_a <- kl_plot(
  collect_kl("pdist", "mrbayes"),
  "Pairwise-distance KL divergence\nfrom MrBayes"
)
panel_b <- kl_plot(
  collect_kl("pdist", "beast"),
  "Pairwise-distance KL divergence\nfrom BEAST2"
)
panel_c <- kl_plot(
  collect_kl("topo", "mrbayes"),
  "Per-split topology KL divergence\nfrom MrBayes"
)
panel_d <- kl_plot(
  collect_kl("topo", "beast"),
  "Per-split topology KL divergence\nfrom BEAST2"
)

figure <- ((panel_a + panel_b) / (panel_c + panel_d)) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical")

out_file <- file.path(out_dir, "hky300_posterior_quality_sup1.pdf")
ggsave(out_file, figure, width = 8.5, height = 6.5, units = "in",
       device = cairo_pdf)
message("Wrote ", out_file)
