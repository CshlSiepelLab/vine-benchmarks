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
  stop("Usage: Rscript plot_runtime_bars_allSizes.R <input.tsv> <output.pdf>", call. = FALSE)
}
tsv_path <- args_trailing[1]
outfile  <- args_trailing[2]

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
      plot.title       = element_text(size = 9, face = "bold"),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 7),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.2)
    )
)

method_palette <- c(
  LAML    = "lightgrey",
  Metient = "#0173B2",
  MACH2   = "#DE8F05",
  BEAM    = "#029E73",
  VINE    = "#D55E00"
)

# ------------ Read data ------------
df <- read.table(tsv_path, header = TRUE, check.names = FALSE)
df$size <- gsub("taxa", "", df$size)
if ("beam_improvement" %in% names(df)) df$beam_improvement <- NULL

# LEFT panel totals (include LAML baseline for Metient/MACH2)
df$mach2_total   <- df$laml + df$mach2
df$metient_total <- df$laml + df$metient

# ------------ Summarise per size group ------------
summary_df <- df %>%
  group_by(size) %>%
  summarise(
    # left totals
    mean_LAML          = mean(laml,          na.rm = TRUE),
    mean_Metient_total = mean(metient_total, na.rm = TRUE),
    mean_MACH2_total   = mean(mach2_total,   na.rm = TRUE),
    mean_BEAM          = mean(beam,          na.rm = TRUE),
    mean_VINE          = mean(vine,          na.rm = TRUE),
    sd_Metient_total   = sd(metient_total,   na.rm = TRUE),
    sd_MACH2_total     = sd(mach2_total,     na.rm = TRUE),
    sd_BEAM            = sd(beam,            na.rm = TRUE),
    sd_VINE            = sd(vine,            na.rm = TRUE),

    # right raw (no LAML added)
    mean_Metient_raw   = mean(metient,       na.rm = TRUE),
    mean_MACH2_raw     = mean(mach2,         na.rm = TRUE),
    mean_BEAM_raw      = mean(beam,          na.rm = TRUE),
    mean_VINE_raw      = mean(vine,          na.rm = TRUE),
    sd_Metient_raw     = sd(metient,         na.rm = TRUE),
    sd_MACH2_raw       = sd(mach2,           na.rm = TRUE),
    sd_BEAM_raw        = sd(beam,            na.rm = TRUE),
    sd_VINE_raw        = sd(vine,            na.rm = TRUE),

    .groups = "drop"
  )

# ------------ Pivot to long (LEFT) ------------
bar_levels <- c("Metient", "MACH2", "BEAM", "VINE")
size_levels <- df %>%
  pull(size) %>%
  as.integer() %>%
  sort() %>%
  unique() %>%
  as.character()

mean_left_long <- summary_df %>%
  select(size, mean_Metient_total, mean_MACH2_total, mean_BEAM, mean_VINE) %>%
  pivot_longer(-size, names_to = "method", values_to = "value") %>%
  mutate(method = sub("^mean_", "", method),
         method = sub("_total$", "", method))

sd_left_long <- summary_df %>%
  select(size, sd_Metient_total, sd_MACH2_total, sd_BEAM, sd_VINE) %>%
  pivot_longer(-size, names_to = "method", values_to = "sd") %>%
  mutate(method = sub("^sd_", "", method),
         method = sub("_total$", "", method))

laml_df <- summary_df %>%
  select(size, mean_LAML) %>%
  rename(laml = mean_LAML) %>%
  crossing(method = bar_levels) %>%
  mutate(laml = ifelse(method %in% c("BEAM", "VINE"), NA, laml),
         method = factor(method, levels = bar_levels),
         size   = factor(size, levels = size_levels))

plot_left_df <- left_join(mean_left_long, sd_left_long, by = c("size", "method")) %>%
  mutate(
    method = factor(method, levels = bar_levels),
    size   = factor(size,   levels = size_levels)
  )

# ------------ Pivot to long (RIGHT; raw values, no LAML added) ------------
mean_right_long <- summary_df %>%
  select(size, mean_Metient_raw, mean_MACH2_raw, mean_BEAM_raw, mean_VINE_raw) %>%
  pivot_longer(-size, names_to = "method", values_to = "value") %>%
  mutate(method = sub("^mean_", "", method),
         method = sub("_raw$", "", method))

sd_right_long <- summary_df %>%
  select(size, sd_Metient_raw, sd_MACH2_raw, sd_BEAM_raw, sd_VINE_raw) %>%
  pivot_longer(-size, names_to = "method", values_to = "sd") %>%
  mutate(method = sub("^sd_", "", method),
         method = sub("_raw$", "", method))

plot_right_df <- left_join(mean_right_long, sd_right_long, by = c("size", "method")) %>%
  mutate(
    method = factor(method, levels = bar_levels),
    size   = factor(size,   levels = size_levels)
  )

# ------------ Shared log-scale helpers (SAME SCALE BOTH PANELS) ------------
all_vals <- c(plot_left_df$value, plot_right_df$value, laml_df$laml)
eps <- min(all_vals[is.finite(all_vals) & all_vals > 0], na.rm = TRUE) / 10

all_tops <- c(plot_left_df$value + plot_left_df$sd,
              plot_right_df$value + plot_right_df$sd,
              laml_df$laml)
y_max <- 10^ceiling(log10(max(all_tops, na.rm = TRUE)))

dodge <- position_dodge(width = 0.9)

# ------------ LEFT PANEL ------------
p_left <- ggplot(plot_left_df, aes(x = size, y = value, fill = method)) +
  geom_col(position = dodge, width = 0.85) +
  geom_col(
    data = laml_df,
    aes(x = size, y = laml, fill = "LAML", group = method),
    inherit.aes = FALSE,
    position = dodge,
    width = 0.85
  ) +
  geom_errorbar(
    aes(ymin = pmax(value - sd, eps), ymax = value + sd),
    position = dodge,
    width = 0.1,
    linewidth = 0.4
  ) +
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
  labs(x = "Number of Taxa (n)", y = "Time (s)", title = "") +
  theme(legend.title = element_blank())

# ------------ RIGHT PANEL (raw; no totals / no LAML added) ------------
p_right <- ggplot(plot_right_df, aes(x = size, y = value, fill = method)) +
  geom_col(position = dodge, width = 0.85) +
  geom_errorbar(
    aes(ymin = pmax(value - sd, eps), ymax = value + sd),
    position = dodge,
    width = 0.1,
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = method_palette,
    breaks = c("Metient", "MACH2", "BEAM", "VINE"), guide = "none"
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::label_comma(),
    limits = c(NA, y_max),
    expand = expansion(mult = c(0, 0))
  ) +
  annotation_logticks(sides = "l") +
  labs(x = "Number of Taxa (n)", y = "Time (s)", title = "")

# ------------ Combine (2 columns, 1 row) ------------
p <- p_left + p_right + plot_layout(ncol = 2)

save_pdf(p, outfile, width = 9, height = 3)