#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 2) {
  stop(
    "Usage: Rscript plot_f1_box.R <input.tsv> <output.pdf>",
    call. = FALSE
  )
}
tsv_path <- args_trailing[1]
outfile  <- args_trailing[2]

# Helper: save a plot as a compact PDF, using cairo if available
save_pdf <- function(plot, filename, width = 4, height = 4) {
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

pretty_method <- function(x) {
  ifelse(x == "METIENT", "Metient", x)
}

method_order <- c("Metient", "MACH2", "BEAM", "VINE")

# ------------ Read data ------------
# whitespace-delimited like sep=r"\s+"
df <- suppressWarnings(read.table(
  file = tsv_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
))

# Format and order methods
df$method <- pretty_method(df$method)
df$method <- factor(df$method, levels = method_order, ordered = TRUE)

# Keep only the ordered methods (avoids scale warnings if extra methods appear)
df <- df %>% filter(as.character(method) %in% method_order)

method_palette <- c(
  Metient = "#76B7B2",
  MACH2   = "#EDC948",
  BEAM    = "#E15759" ,
  VINE    = "#F28E2B"
)

# ------------ Plot ------------
p <- ggplot(df, aes(x = method, y = f1, fill = method)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.4, whisker.linewidth = 0.4, staplewidth = 0.1) +
  scale_fill_manual(values = method_palette, guide = "none") +
  labs(x = NULL, y = "F1 score") +
  coord_cartesian(ylim = c(0, 1.05)) +
  ylim(-0.01, 1.01) +
  labs(x = "Method")

save_pdf(p, outfile, width = 3, height = 4)