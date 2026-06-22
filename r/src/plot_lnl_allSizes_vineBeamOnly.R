#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 2) {
  stop("Usage: Rscript plot_lnl_allSizes.R <input.tsv> <output.pdf>", call. = FALSE)
}
tsv_path <- args_trailing[1]
outfile  <- args_trailing[2]

save_pdf <- function(plot, filename, width = 12, height = 5) {
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
base_theme <- theme_minimal(base_size = 8) +
  theme(
    plot.title       = element_text(size = 9, face = "bold"),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 7),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray80", linewidth = 0.2)
  )
theme_set(base_theme)

method_palette <- c(
  "BEAM" = "#E15759",
  "VINE" = "#F28E2B"
)

method_levels <- c("BEAM", "VINE")

# ------------ Read & clean data ------------
df <- read.table(tsv_path, header = TRUE, check.names = FALSE)
df$size <- gsub("taxa", "", df$size)
size_levels <- df %>%
  pull(size) %>%
  as.integer() %>%
  sort() %>%
  unique() %>%
  as.character()

df <- df %>%
  select(-laml_tree) %>%
  mutate(
    BEAM = beam_tree + beam_mig,
    VINE = vine_tree + vine_mig
  ) %>%
  select(size, BEAM, VINE) %>%
  pivot_longer(-size, names_to = "method", values_to = "lnl") %>%
  mutate(
    method = factor(method, levels = method_levels),
    size   = factor(size,   levels = size_levels)
  )

# ------------ Summarise ------------
plot_df <- df %>%
  group_by(size, method) %>%
  summarise(mean = mean(lnl), sd = sd(lnl), .groups = "drop")

# ------------ Plot ------------
p <- ggplot(plot_df, aes(x = size, y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.9),
    width = 0.1, linewidth = 0.4
  ) +
  scale_fill_manual(values = method_palette, breaks = method_levels) +
  labs(x = "Number of Taxa (n)", y = "Log-Likelihood", fill = "Method") +
  theme(legend.title = element_blank())

save_pdf(p, outfile, width = 3, height = 3)