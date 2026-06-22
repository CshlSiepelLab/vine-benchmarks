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
  LAML    = "#59A14F",
  Metient = "#76B7B2",
  MACH2   = "#EDC948",
  BEAM    = "#E15759" ,
  VINE    = "#F28E2B"
)

# ------------ Read data ------------
df <- read.table(tsv_path, header = TRUE, check.names = FALSE)
df$size <- gsub("taxa", "", df$size)
if ("beam_improvement" %in% names(df)) df$beam_improvement <- NULL

# Ppanel totals (include LAML baseline for Metient/MACH2)
df$mach2_total   <- df$laml + df$mach2
df$metient_total <- df$laml + df$metient

# ------------ Summarise per size group ------------
summary_df <- df %>%
  group_by(size) %>%
  summarise(
    mean_LAML          = mean(laml,          na.rm = TRUE),
    mean_Metient_total = mean(metient_total, na.rm = TRUE),
    mean_MACH2_total   = mean(mach2_total,   na.rm = TRUE),
    mean_BEAM          = mean(beam,          na.rm = TRUE),
    mean_VINE          = mean(vine,          na.rm = TRUE),
    sd_Metient_total   = sd(metient_total,   na.rm = TRUE),
    sd_MACH2_total     = sd(mach2_total,     na.rm = TRUE),
    sd_BEAM            = sd(beam,            na.rm = TRUE),
    sd_VINE            = sd(vine,            na.rm = TRUE),

    # Raw, no LAML added
    mean_Metient   = mean(metient,       na.rm = TRUE),
    mean_MACH2     = mean(mach2,         na.rm = TRUE),

    .groups = "drop"
  )

# Pivot to long
bar_levels <- c("Metient", "MACH2", "BEAM", "VINE")
size_levels <- df %>%
  pull(size) %>%
  as.integer() %>%
  sort() %>%
  unique() %>%
  as.character()

mean_long_raw <- summary_df %>%
  select(size, mean_Metient, mean_MACH2, mean_BEAM, mean_VINE) %>%
  pivot_longer(-size, names_to = "method", values_to = "value") %>%
  mutate(method = sub("^mean_", "", method))

plot_df_raw <- mean_long_raw %>%
  mutate(
    method = factor(method, levels = bar_levels),
    size   = factor(size,   levels = size_levels)
  )

mean_long <- summary_df %>%
  select(size, mean_Metient_total, mean_MACH2_total, mean_BEAM, mean_VINE) %>%
  pivot_longer(-size, names_to = "method", values_to = "value") %>%
  mutate(method = sub("^mean_", "", method),
         method = sub("_total$", "", method))

sd_long <- summary_df %>%
  select(size, sd_Metient_total, sd_MACH2_total, sd_BEAM, sd_VINE) %>%
  pivot_longer(-size, names_to = "method", values_to = "sd") %>%
  mutate(method = sub("^sd_", "", method),
         method = sub("_total$", "", method))

plot_df <- left_join(mean_long, sd_long, by = c("size", "method")) %>%
  mutate(
    method = factor(method, levels = bar_levels),
    size   = factor(size,   levels = size_levels)
  )

# ------------ Shared log-scale helpers (SAME SCALE BOTH PANELS) ------------
all_vals <- c(plot_df$value, plot_df_raw$value)
eps <- min(all_vals[is.finite(all_vals) & all_vals > 0], na.rm = TRUE) / 10

all_tops <- c(plot_df$value + plot_df$sd, plot_df_raw$value)
y_max <- 10^ceiling(log10(max(all_tops, na.rm = TRUE)))

dodge <- position_dodge(width = 0.9)


p <- ggplot(plot_df, aes(x = size, y = value, fill = "LAML", group=method)) +
  geom_col(position = dodge, width = 0.85) +
  geom_col(
    data = plot_df_raw,
    aes(x = size, y = value, fill = method, group = method),
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


save_pdf(p, outfile, width = 6, height = 3)