#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(grid))
suppressMessages(library(patchwork))

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
  "LAML-tree"           = "#59A14F",
  "BEAM-tree"      = "#E15759",
  "BEAM-migration" = "#E6706F",
  "VINE-tree"      = "#F28E2B",
  "VINE-migration" = "#F5A84A"
)

method_levels <- c("LAML-tree", "BEAM-tree", "BEAM-migration", "VINE-tree", "VINE-migration")

# ------------ Read & clean data ------------
df <- read.table(tsv_path, header = TRUE, check.names = FALSE)
df$size <- gsub("taxa", "", df$size)
size_levels <- df %>%
  pull(size) %>%
  as.integer() %>%
  sort() %>%
  unique() %>%
  as.character()

# ------------ Summarise per size group ------------
summary_df <- df %>%
  group_by(size) %>%
  summarise(
    mean_LAML_tree         = mean(laml_tree,      na.rm = TRUE),
    mean_BEAM_tree      = mean(beam_tree, na.rm = TRUE),
    mean_BEAM_migration = mean(beam_mig,  na.rm = TRUE),
    mean_VINE_tree      = mean(vine_tree, na.rm = TRUE),
    mean_VINE_migration = mean(vine_mig,  na.rm = TRUE),
    sd_LAML_tree             = sd(laml_tree,        na.rm = TRUE),
    sd_BEAM_tree        = sd(beam_tree,   na.rm = TRUE),
    sd_BEAM_migration   = sd(beam_mig,    na.rm = TRUE),
    sd_VINE_tree        = sd(vine_tree,   na.rm = TRUE),
    sd_VINE_migration   = sd(vine_mig,    na.rm = TRUE),
    .groups = "drop"
  )

mean_long <- summary_df %>%
  select(size, starts_with("mean_")) %>%
  pivot_longer(-size, names_to = "method", values_to = "mean") %>%
  mutate(method = sub("mean_", "", method),
         method = gsub("_", "-", method))

sd_long <- summary_df %>%
  select(size, starts_with("sd_")) %>%
  pivot_longer(-size, names_to = "method", values_to = "sd") %>%
  mutate(method = sub("sd_", "", method),
         method = gsub("_", "-", method))

plot_df <- left_join(mean_long, sd_long, by = c("size", "method")) %>%
  mutate(
    method    = factor(method, levels = method_levels),
    size      = factor(size,   levels = size_levels),
    component = case_when(
      grepl("migration", method) ~ "Migration log-likelihood",
      TRUE                       ~ "Tree log-likelihood"
    )
  )

# ------------ Left panel: original single-panel style ------------
p_left <- ggplot(plot_df, aes(x = size, y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.9),
    width = 0.1, linewidth = 0.4
  ) +
  scale_fill_manual(values = method_palette, breaks = method_levels) +
  labs(x = "Number of Taxa (n)", y = "Log-Likelihood", fill = "Method") +
  theme(legend.title = element_blank())

# ------------ Right panel: migration only ------------
p_right <- ggplot(plot_df %>% filter(grepl("migration", method)),
                  aes(x = size, y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.9),
    width = 0.1, linewidth = 0.4
  ) +
  scale_fill_manual(values = method_palette, breaks = method_levels, guide = "none") +
  labs(x = "Number of Taxa (n)", y = "Migration Log-Likelihood")

# ------------ Combine with patchwork ------------
# Collect legends so only one shared legend appears on the right
p_combined <- (p_left | p_right) +
  plot_layout(guides = "collect", widths = c(0.75, 0.25)) &
  theme(legend.position = "right")

save_pdf(p_combined, outfile, width = 9, height = 3)