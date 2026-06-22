#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 2) {
  stop(
    "Usage: Rscript plot_precision_recall.R <input.tsv> <output.pdf>",
    call. = FALSE
  )
}
tsv_input  <- args_trailing[1]
pdf_output <- args_trailing[2]

# Helper: save a plot as a compact PDF, using cairo if available
save_pdf <- function(plot, filename, width = 5, height = 4) {
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

method_palette <- c(
  Metient = "#76B7B2",
  MACH2   = "#EDC948",
  BEAM    = "#E15759" ,
  VINE    = "#F28E2B"
)

method_labels <- c(
  Metient = "LAML - Metient",
  MACH2   = "LAML - MACH2",
  BEAM    = "BEAM",
  VINE    = "VINE"
)

label_map <- function(keys) {
  out <- method_labels[keys]
  out[is.na(out)] <- keys[is.na(out)]
  out
}

# ------------ Read & summarize ------------
df <- read_tsv(tsv_input, show_col_types = FALSE)

# Drop the f1 column if present (be forgiving)
if ("f1" %in% names(df)) {
  df <- df %>% select(-f1)
}

# Mean across sims for each method/threshold, then drop sim column
df_mean <- df %>%
  group_by(method, threshold) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

if ("sim" %in% names(df_mean)) {
  df_mean <- df_mean %>% select(-sim)
}

# Format method names and set ordering
df_mean <- df_mean %>%
  mutate(method = recode(method, "METIENT" = "Metient"))

method_order <- c("Metient", "MACH2", "BEAM", "VINE")
df_mean$method <- factor(df_mean$method, levels = method_order, ordered = TRUE)

# Keep only methods we know how to color (avoids scale warnings)
df_mean <- df_mean %>% filter(as.character(method) %in% names(method_palette))

# ------------ Plot ------------
p <- ggplot(df_mean, aes(x = recall, y = precision, color = method)) +
  geom_line(linewidth = 1.0) +
  labs(x = "Recall", y = "Precision", color = NULL) +
  scale_color_manual(
    values = method_palette[levels(df_mean$method)],
    breaks = levels(df_mean$method),
    labels = label_map(levels(df_mean$method))
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.0, 0.0))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.0, 0.0))
  ) +
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 6),
    legend.key.height = unit(0.8, "lines")
  )

# ------------ Save ------------
save_pdf(p, pdf_output, width = 4, height = 3)