#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(patchwork))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
hide_cass_cc <- "--hide-cass-cc" %in% args_trailing
args_trailing <- args_trailing[!args_trailing %in% "--hide-cass-cc"]
if (length(args_trailing) < 2) {
  stop("Usage: Rscript plot_stats_byCp.R <eval.all.stats.txt> <output.pdf> [--hide-cass-cc]", call. = FALSE)
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

method_palette <- c(
  "Cassiopeia-Greedy" = "#B07AA1",
  "LAML"       = "#59A14F",
  "BEAM"       = "#E15759",
  "VINE"       = "#F28E2B"
)
method_levels <- c("Cassiopeia-Greedy", "LAML", "BEAM", "VINE")

stat_labels <- c(
  cc  = "Cophenetic Correlation",
  pm  = "Parsimony Mutations",
  thi = "Tissue Homogeneity Index",
  pt  = "Parsimony Tissues"
)

# ------------ Read data ------------
# Expected columns (tab-separated, no header):
#   cp  cc_cass  pm_cass  thi_cass  pt_cass
#       cc_laml  pm_laml  thi_laml  pt_laml
#       cc_beam  pm_beam  thi_beam  pt_beam
#       cc_vine  pm_vine  thi_vine  pt_vine
# i.e. 1 + 4*4 = 17 columns, or with header we read it dynamically.
#
# The Makefile writes: cp \t cc_cass \t pm_cass \t thi_cass \t pt_cass
#                               \t cc_laml \t ...  (4 trees * 4 stats = 16 values)
# So columns: cp, then groups of 4 (cc,pm,thi,pt) for cass, laml, beam, vine

raw_lines <- readLines(tsv_path)
# Drop blank lines and separator lines
clean_lines <- raw_lines[!grepl("^-+$", raw_lines) & nzchar(trimws(raw_lines))]
df_raw <- read.table(text = paste(clean_lines, collapse = "\n"),
                     header = TRUE, sep = "\t", check.names = FALSE,
                     fill = TRUE)

# Keep only numeric cp rows (drops "avg" footer)
df_raw <- df_raw[grepl("^[0-9]+$", as.character(df_raw$cp)), ]
df_raw$cp <- as.integer(df_raw$cp)

cp_levels <- sort(unique(df_raw$cp))

# Compute avg row
avg_row <- df_raw %>%
  summarise(across(-cp, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(cp = NA_integer_)

df_all <- bind_rows(df_raw, avg_row)

# Optionally mask Cassiopeia cophenetic correlation
if (hide_cass_cc) df_all$cass_cc <- NA_real_

cp_levels_ext <- c(as.character(cp_levels), "ave")
df_all$cp_label <- ifelse(is.na(df_all$cp), "ave", as.character(df_all$cp))
df_all$cp_label <- factor(df_all$cp_label, levels = cp_levels_ext)

# ------------ Build one plot per statistic ------------
make_panel <- function(stat, ylabel) {
  # Gather the 4 method columns for this stat
  # File column order is method_stat (e.g. cass_cc, laml_cc, ...)
  cols <- paste0(c("cass", "laml", "beam", "vine"), "_", stat)

  plot_df <- df_all %>%
    select(cp_label, all_of(cols)) %>%
    rename(!!!setNames(cols, method_levels)) %>%
    pivot_longer(-cp_label, names_to = "method", values_to = "value") %>%
    mutate(
      method   = factor(method, levels = method_levels),
      cp_label = factor(cp_label, levels = cp_levels_ext)
    )

  ggplot(plot_df, aes(x = cp_label, y = value, fill = method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.85) +
    geom_vline(xintercept = length(cp_levels) + 0.5,
               color = "black", linewidth = 0.4) +
    scale_x_discrete(expand = expansion(add = c(1.0, 0.5))) +
    scale_fill_manual(values = method_palette, breaks = method_levels) +
    labs(x = "Clonal Population (CP)", y = ylabel, title = "") +
    theme(legend.title = element_blank())
}

p_cc  <- make_panel("cc",  stat_labels["cc"])
p_pm  <- make_panel("pm",  stat_labels["pm"])
p_thi <- make_panel("thi", stat_labels["thi"])
p_pt  <- make_panel("pt",  stat_labels["pt"])

# ------------ Combine with patchwork ------------
combined <- p_cc / p_pm / p_thi / p_pt +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Save: 9 wide x 12 tall (4 panels x 3 height each)
if (capabilities("cairo")) {
  ggsave(outfile, plot = combined, width = 9, height = 12,
         units = "in", device = cairo_pdf)
} else {
  pdf(outfile, width = 9, height = 12, family = "Helvetica")
  print(combined)
  dev.off()
}