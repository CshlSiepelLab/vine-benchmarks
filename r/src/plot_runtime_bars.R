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
    "Usage: Rscript plot_runtime_bars.R <input.tsv> <output.pdf>",
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
  LAML    = "lightgrey",
  Metient = "#0173B2",
  MACH2   = "#DE8F05",
  BEAM    = "#029E73",
  VINE    = "#D55E00"
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

# Drop improvement column if present
if ("beam_improvement" %in% names(df)) {
  df$beam_improvement <- NULL
}

# Totals for stacked methods (laml + method)
df$mach2_total   <- df$laml + df$mach2
df$metient_total <- df$laml + df$metient

mean_laml <- mean(df$laml, na.rm = TRUE)

means <- c(
  Metient = mean(df$metient_total, na.rm = TRUE),
  MACH2   = mean(df$mach2_total,   na.rm = TRUE),
  BEAM    = mean(df$beam,          na.rm = TRUE),
  VINE    = mean(df$vine,          na.rm = TRUE)
)

sds_total <- c(
  Metient = sd(df$metient_total, na.rm = TRUE),
  MACH2   = sd(df$mach2_total,   na.rm = TRUE),
  BEAM    = sd(df$beam,          na.rm = TRUE),
  VINE    = sd(df$vine,          na.rm = TRUE)
)

bar_labels <- c("Metient", "MACH2", "BEAM", "VINE")

plot_df <- data.frame(
  method = factor(bar_labels, levels = bar_labels),
  total  = unname(means[bar_labels]),
  sd     = unname(sds_total[bar_labels])
)

# ---- log-scale safe errorbars (cannot go <= 0) ----
eps <- min(plot_df$total[plot_df$total > 0], na.rm = TRUE) / 10
plot_df <- plot_df %>%
  mutate(
    ymin = pmax(total - sd, eps),
    ymax = total + sd
  )

# Compute the upper limit rounded up to nearest power of 10
y_max <- 10^ceiling(log10(max(plot_df$ymax, na.rm = TRUE)))

# LAML overlay df — only for Metient and MACH2
laml_df <- data.frame(
  method = factor(c("Metient", "MACH2"), levels = bar_labels),
  laml   = mean_laml
)

p <- ggplot(plot_df, aes(x = method, y = total, fill = method)) +
  geom_col(position = position_dodge(width = 0.9), width=0.9) +
  geom_col(data = laml_df, aes(x = method, y = laml, fill="LAML"), inherit.aes = FALSE, width=0.9) +
  geom_errorbar(aes(ymin = total - sd, ymax = total + sd),
                  width = 0.1,
                  linewidth = 0.4,
                  position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = method_palette,
    breaks = c("LAML", "Metient", "MACH2", "BEAM", "VINE")
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::label_comma(),
    limits = c(NA, y_max),
    expand = expansion(mult = c(0, 0))
  ) +
  annotation_logticks(sides = "l") +
  labs(x = "Method", y = "Time(s)", title = "") +
  theme(legend.title = element_blank())

save_pdf(p, outfile, width = 5, height = 4)