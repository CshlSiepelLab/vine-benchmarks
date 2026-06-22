#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(grid))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 2) {
  stop(
    "Usage: Rscript plot_lnl.R <input.tsv> <output.pdf>",
    call. = FALSE
  )
}
tsv_path <- args_trailing[1]
outfile  <- args_trailing[2]

# Helper: save a plot as a compact PDF, using cairo if available
save_pdf <- function(plot, filename, width = 6, height = 4) {
  if (capabilities("cairo")) {
    ggsave(filename, plot = plot, width = width, height = height,
           units = "in", device = cairo_pdf)
  } else {
    pdf(filename, width = width, height = height, family = "Helvetica")
    print(plot)
    dev.off()
  }
}

# ------------ Theme ------------
theme_set(
  theme_minimal(base_size = 8) +
    theme(
      plot.title  = element_text(size = 9, face = "bold"),
      axis.title  = element_text(size = 9),
      axis.text   = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7),
      panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.2)
    )
)

# Colors
method_palette <- c(
  "LAML"        = "lightgrey",
  "BEAM-tree"        = "#029E73",
  "BEAM-migration" = "#5CBF9F",
  "VINE-tree"        = "#D55E00",
  "VINE-migration" = "#E69F73"
)

# ------------ Read data (skip last 2 lines) ------------
lines <- readLines(tsv_path, warn = FALSE)
if (length(lines) <= 2) stop("Input has <= 2 lines; cannot skip footer.", call. = FALSE)
lines2 <- lines[seq_len(length(lines) - 2)]

df <- suppressWarnings(read.table(
  text = lines2,
  header = TRUE,
  check.names = FALSE
))

df$beam_total <- df$beam_tree + df$beam_mig
df$vine_total <- df$vine_tree + df$vine_mig

means <- c(
  LAML = mean(df$laml,      na.rm = TRUE),
  "BEAM-tree" = mean(df$beam_tree, na.rm = TRUE),
  "BEAM-migration" = mean(df$beam_mig, na.rm = TRUE),
  "VINE-tree" = mean(df$vine_tree, na.rm = TRUE),
  "VINE-migration" = mean(df$vine_mig, na.rm = TRUE)
)

sds_total <- c(
  LAML = sd(df$laml,       na.rm = TRUE),
  "BEAM-tree" = sd(df$beam_tree, na.rm = TRUE),
  "BEAM-migration" = sd(df$beam_mig, na.rm = TRUE),
  "VINE-tree" = sd(df$vine_tree, na.rm = TRUE),
  "VINE-migration" = sd(df$vine_mig, na.rm = TRUE)
)

bar_labels <- c("LAML", "BEAM-tree", "BEAM-migration", "VINE-tree", "VINE-migration")

plot_df <- data.frame(
  method = factor(bar_labels, levels = bar_labels),
  total  = unname(means[bar_labels]),
  sd     = unname(sds_total[bar_labels])
)

p <- ggplot(plot_df, aes(x = method)) +
  geom_col(aes(y = total, fill = method), width = 0.9) +
  geom_errorbar(
    aes(ymin = total - sd, ymax = total + sd),
    width = 0.1,
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = method_palette,
    breaks = c("LAML", "BEAM-tree", "VINE-tree", "BEAM-migration", "VINE-migration")
  ) +
  labs(x = "Method", y = "Log-Likelihood", title = "") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)) +
  theme(legend.title = element_blank())

save_pdf(p, outfile, width = 4, height = 4)