#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(patchwork))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 2) {
  stop("Usage: Rscript plot_migration_stats_byCp.R <input.tsv> <output.pdf>", call. = FALSE)
}
tsv_path <- args_trailing[1]
outfile  <- args_trailing[2]

# ------------ Theme ------------
base_theme <- theme_minimal(base_size = 8) +
  theme(
    plot.title       = element_text(size = 9, face = "bold"),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 7),
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray80", linewidth = 0.2)
  )
theme_set(base_theme)

# Keep same general style/colors as before, but mapped to these methods
method_palette <- c(
  "Metient" = "#76B7B2",
  "BEAM"    = "#E15759",
  "VINE"    = "#F28E2B"
)
method_levels <- c("Metient", "BEAM", "VINE")

stat_labels <- c(
  mig       = "Number of Migrations",
  comig     = "Number of Co-migrations",
  seedsites = "Number of Seeding Sites"
)

# ------------ Read data ------------
raw_lines <- readLines(tsv_path)

# Drop blank lines and separator lines
clean_lines <- raw_lines[!grepl("^-+$", raw_lines) & nzchar(trimws(raw_lines))]

df_raw <- read.table(
  text = paste(clean_lines, collapse = "\n"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

# Keep only numeric cp rows (drops avg footer if present)
df_raw <- df_raw[grepl("^[0-9]+$", as.character(df_raw$cp)), ]
df_raw$cp <- as.integer(df_raw$cp)

cp_levels <- sort(unique(df_raw$cp))

# Compute avg row from numeric rows only
avg_row <- df_raw %>%
  summarise(across(-cp, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(cp = NA_integer_)

df_all <- bind_rows(df_raw, avg_row)

cp_levels_ext <- c(as.character(cp_levels), "ave")
df_all$cp_label <- ifelse(is.na(df_all$cp), "ave", as.character(df_all$cp))
df_all$cp_label <- factor(df_all$cp_label, levels = cp_levels_ext)

# ------------ Build one plot per statistic ------------
make_panel <- function(stat, ylabel, show_x_title = FALSE) {
  cols <- paste0(c("metient", "beam", "var"), "_", stat)

  plot_df <- df_all %>%
    select(cp_label, all_of(cols)) %>%
    rename(
      Metient = all_of(cols[1]),
      BEAM    = all_of(cols[2]),
      VINE    = all_of(cols[3])
    ) %>%
    pivot_longer(-cp_label, names_to = "method", values_to = "value") %>%
    mutate(
      method   = factor(method, levels = method_levels),
      cp_label = factor(cp_label, levels = cp_levels_ext)
    )

  ggplot(plot_df, aes(x = cp_label, y = value, fill = method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.85) +
    geom_vline(
      xintercept = length(cp_levels) + 0.5,
      color = "black",
      linewidth = 0.4
    ) +
    scale_x_discrete(expand = expansion(add = c(1.0, 0.5))) +
    scale_fill_manual(values = method_palette, breaks = method_levels) +
    labs(
      x = if (show_x_title) "Clonal Population (CP)" else NULL,
      y = ylabel,
      title = NULL
    ) +
    theme(
      legend.title = element_blank(),
      axis.title.x = if (show_x_title) element_text(size = 9) else element_blank()
    )
}

p_mig   <- make_panel("mig",       stat_labels["mig"],       show_x_title = FALSE)
p_comig <- make_panel("comig",     stat_labels["comig"],     show_x_title = FALSE)
p_seed  <- make_panel("seedsites", stat_labels["seedsites"], show_x_title = TRUE)

# ------------ Combine with patchwork ------------
combined <- p_mig / p_comig / p_seed +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Save: 9 wide x 9 tall (3 panels x 3 height each)
if (capabilities("cairo")) {
  ggsave(
    outfile,
    plot = combined,
    width = 9,
    height = 9,
    units = "in",
    device = cairo_pdf
  )
} else {
  pdf(outfile, width = 9, height = 9, family = "Helvetica")
  print(combined)
  dev.off()
}