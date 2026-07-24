#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# First main posterior-quality figure. Panel A is intentionally reserved.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
root_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
out_dir <- file.path(root_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Staged copies of the raw per-replicate eval files, produced by
# `bash extractGraphData.sh hky300-postqual` (see graphs/extractGraphData.sh).
data_dir <- file.path(dirname(root_dir), "hky300-posterioranalysis-data")

taxa <- c(10, 25, 50, 100, 250, 500, 1000)
method_levels <- c("Vine", "Vine + flows", "BEAST2", "MrBayes")
method_colors <- c(
  "Vine" = "#F28E2B",
  "Vine + flows" = "#76B7B2",
  "BEAST2" = "#59A14F",
  "MrBayes" = "#E15759"
)

# Summary files end with a dashed separator and a precomputed mean row. Read
# only the replicate rows so that means and sample SDs are calculated here.
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

ci_rows <- lapply(taxa, function(n) {
  vine <- read_replicates(file.path(data_dir, paste0(n, "taxa"),
                                    "vine.eval.all.dist.txt"))
  flows <- read_replicates(file.path(data_dir, paste0(n, "taxa"),
                                     "vine_flows.eval.all.dist.txt"))
  values <- list(
    "Vine" = vine$vine_95CI,
    "Vine + flows" = flows$vine_95CI,
    "BEAST2" = vine$beast_95CI,
    "MrBayes" = vine$mrbayes_95CI
  )
  do.call(rbind, lapply(names(values), function(method) {
    z <- mean_sd(values[[method]])
    data.frame(ntaxa = n, method = method, mean = z[["mean"]], sd = z[["sd"]])
  }))
})
ci_data <- do.call(rbind, ci_rows)

kl_files <- c(
  "Vine" = "vine.eval.all.kl_pdist_mrbayes.txt",
  "Vine + flows" = "vine_flows.eval.all.kl_pdist_mrbayes.txt",
  "BEAST2" = "beast.eval.all.kl_pdist_mrbayes.txt"
)
kl_rows <- lapply(taxa, function(n) {
  rows <- lapply(names(kl_files), function(method) {
    d <- read_replicates(file.path(data_dir, paste0(n, "taxa"), kl_files[[method]]))
    z <- mean_sd(d$kl_pdist)
    data.frame(ntaxa = n, method = method, mean = z[["mean"]], sd = z[["sd"]])
  })
  rows[[length(rows) + 1]] <- data.frame(
    ntaxa = n, method = "MrBayes", mean = 0, sd = 0
  )
  do.call(rbind, rows)
})
kl_data <- do.call(rbind, kl_rows)

prepare_data <- function(d) {
  d$method <- factor(d$method, levels = method_levels)
  d$ntaxa <- factor(d$ntaxa, levels = taxa)
  d
}
ci_data <- prepare_data(ci_data)
kl_data <- prepare_data(kl_data)

# Panel A: 25-taxon, Euclidean-geometry embedding-dimension plot. This is a
# Euclidean-only, 25-taxa subset of the dual-axis panel built by
# make_dim_panel2() in graphs/supplement-figure6-embedding-dims/makeGraphs.R.
# Left axis (orange): Delta lnl = VINE - BEAST held-out log-likelihood, mean
# +/- SD over replicates. Right axis (green): VINE runtime, median with
# Q1-Q3 whiskers (robust to the heavy-tailed convergence time).
# dimSummary25.txt columns (by name): d, dlnlE, dlnlE_sd, timeE_med,
# timeE_q1, timeE_q3, dlnlH, dlnlH_sd, timeH_med, timeH_q1, timeH_q3.
dim_file <- file.path(dirname(root_dir), "hky300-data", "dimSummary25.txt")
dim_summary <- read.table(dim_file, header = TRUE, check.names = FALSE)
dim_value <- suppressWarnings(
  as.numeric(gsub(".*\\.D", "", as.character(dim_summary[[1]])))
)

d_euc <- dim_summary$dlnlE
d_euc_sd <- dim_summary$dlnlE_sd
te_m <- dim_summary$timeE_med
te_lo <- dim_summary$timeE_q1
te_hi <- dim_summary$timeE_q3

num_max <- max(abs(d_euc) + d_euc_sd, na.rm = TRUE)
den_max <- max(te_hi, na.rm = TRUE)
dim_scale_factor <- if (is.finite(num_max) && is.finite(den_max) &&
                        den_max > 0 && num_max > 0) num_max / den_max else 1

dim_data <- rbind(
  data.frame(
    dimension = dim_value, metric = "Delta lnl",
    y = d_euc, lo = d_euc - d_euc_sd, hi = d_euc + d_euc_sd
  ),
  data.frame(
    dimension = dim_value, metric = "time",
    y = te_m * dim_scale_factor,
    lo = te_lo * dim_scale_factor, hi = te_hi * dim_scale_factor
  )
)
dim_data$dimension <- factor(dim_data$dimension, levels = 2:8)
dim_data$metric <- factor(dim_data$metric, levels = c("Delta lnl", "time"))

dim_ll_lo <- min(c(d_euc - d_euc_sd, 0), na.rm = TRUE)
dim_ll_hi <- max(c(d_euc + d_euc_sd, 0), na.rm = TRUE)
dim_ll_breaks <- pretty(c(dim_ll_lo, dim_ll_hi), n = 6)
dim_t_breaks <- pretty(c(0, max(te_hi, 0, na.rm = TRUE)), n = 5)
dim_ymin <- min(dim_ll_breaks, min(dim_t_breaks) * dim_scale_factor, na.rm = TRUE)
dim_ymax <- max(dim_ll_breaks, max(dim_t_breaks) * dim_scale_factor, na.rm = TRUE)
# Merge in the (scaled) time breaks so major gridlines span the full panel,
# not just the Delta-lnl range -- make_dim_panel2() only grids the Delta-lnl
# breaks, leaving the time-bar region gridline-free.
dim_combined_breaks <- sort(unique(c(dim_ll_breaks, dim_t_breaks * dim_scale_factor)))

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

bar_plot <- function(d, y_label, limits = NULL) {
  ggplot(d, aes(x = ntaxa, y = mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    geom_errorbar(
      aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
      position = position_dodge(width = 0.9), width = 0.18, linewidth = 0.35
    ) +
    scale_fill_manual(values = method_colors, breaks = method_levels, drop = FALSE) +
    scale_y_continuous(limits = limits, expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of Taxa (n)", y = y_label, fill = "Method") +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    paper_theme
}

panel_a <- ggplot(dim_data, aes(x = dimension, y = y, fill = metric)) +
  geom_col(position = position_dodge(width = 0.55), width = 0.48) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    position = position_dodge(width = 0.55), width = 0.15, linewidth = 0.35
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  scale_fill_manual(values = c("Delta lnl" = "#F28E2B", "time" = "#59A14F")) +
  guides(fill = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
  scale_y_continuous(
    name = expression(Delta ~ lnl),
    limits = c(dim_ymin, dim_ymax), breaks = dim_combined_breaks,
    labels = function(x) ifelse(x <= 0, x, ""),
    sec.axis = sec_axis(~ . / dim_scale_factor, name = "time (s)",
                        breaks = dim_t_breaks)
  ) +
  labs(x = "Embedding Dimension (d)") +
  paper_theme +
  theme(
    legend.position = "none",
    axis.title.y = element_text(color = "#F28E2B"),
    axis.text.y = element_text(color = "#F28E2B"),
    axis.ticks.y = element_line(color = "#F28E2B"),
    axis.title.y.right = element_text(color = "#59A14F"),
    axis.text.y.right = element_text(color = "#59A14F"),
    axis.ticks.y.right = element_line(color = "#59A14F")
  )
panel_b <- bar_plot(ci_data, "95% CI Inclusion", limits = c(0, 1.05))
panel_c <- bar_plot(kl_data,
                    "Pairwise-distance KL divergence\nfrom MrBayes")

figure <- (panel_a + panel_b + panel_c) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", widths = c(0.85, 1.08, 1.08)) &
  theme(legend.position = "right", legend.direction = "vertical")

out_file <- file.path(out_dir, "hky300_posterior_quality_main.pdf")
ggsave(out_file, figure, width = 10.5, height = 3.35, units = "in",
       device = cairo_pdf)
message("Wrote ", out_file)
