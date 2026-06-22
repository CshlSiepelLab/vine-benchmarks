#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(grid))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 2) {
  stop("Usage: Rscript plot_runtime_bars_perCP.R <input.tsv> <output.pdf>", call. = FALSE)
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
  BEAM    = "#E15759",
  VINE    = "#F28E2B"
)

# ------------ Read data ------------
# Read raw lines to find how many are valid data (skip separator/avg footer)
raw_lines <- readLines(tsv_path)
# Keep header + lines where the first field is an integer
is_data <- c(TRUE, sapply(raw_lines[-1], function(l) {
  first <- strsplit(trimws(l), "\\s+")[[1]][1]
  !is.na(suppressWarnings(as.integer(first)))
}))
clean_text <- paste(raw_lines[is_data], collapse = "\n")
df <- read.table(text = clean_text, header = TRUE, check.names = FALSE)
df$sim <- as.integer(df$sim)

# Compute Metient stacked total (Metient runtime + LAML baseline)
df$metient_total <- df$laml + df$metient

# ------------ Build long-format plot data ------------
# Two layers:
#   1. Background bar = LAML (shown via metient_total height, colored LAML)
#   2. Foreground bar = Metient raw (on top, colored Metient), BEAM, VINE

bar_levels <- c("Metient", "BEAM", "VINE")
sim_levels <- sort(unique(df$sim))

# Compute average row and append as a special "Avg" group
avg_row <- df %>%
  summarise(
    sim           = NA_integer_,
    vine          = mean(vine,          na.rm = TRUE),
    laml          = mean(laml,          na.rm = TRUE),
    metient       = mean(metient,       na.rm = TRUE),
    beam          = mean(beam,          na.rm = TRUE),
    metient_total = mean(metient_total, na.rm = TRUE)
  )
df_plot <- bind_rows(df, avg_row)
sim_levels_ext <- c(as.character(sim_levels), "ave")
df_plot$sim_label <- ifelse(is.na(df_plot$sim), "ave", as.character(df_plot$sim))
df_plot$sim_label <- factor(df_plot$sim_label, levels = sim_levels_ext)

# Raw values for BEAM, VINE, and Metient-only (foreground bars)
plot_df_raw <- df_plot %>%
  select(sim_label, metient, beam, vine) %>%
  rename(Metient = metient, BEAM = beam, VINE = vine) %>%
  pivot_longer(-sim_label, names_to = "method", values_to = "value") %>%
  mutate(method = factor(method, levels = bar_levels))

# Background LAML layer: full stacked height for Metient, 0 for BEAM/VINE
# (zero entries are needed so all methods occupy the same dodge slots)
plot_df_total <- bind_rows(
  df_plot %>% select(sim_label, metient_total) %>%
    rename(value = metient_total) %>%
    mutate(method = "Metient"),
  df_plot %>% select(sim_label) %>% mutate(method = "BEAM", value = NA_real_),
  df_plot %>% select(sim_label) %>% mutate(method = "VINE", value = NA_real_)
) %>%
  mutate(method = factor(method, levels = bar_levels))

# y-axis range
real_vals <- c(plot_df_raw$value, plot_df_total$value)
real_vals <- real_vals[is.finite(real_vals) & real_vals > 0]
y_min <- 10^floor(log10(min(real_vals, na.rm = TRUE)))
y_max <- 10^ceiling(log10(max(real_vals, na.rm = TRUE)))

dodge <- position_dodge(width = 0.9)

# ------------ Plot ------------
p <- ggplot() +
  # Background LAML layer (full stacked height for Metient bars only)
  geom_col(
    data = plot_df_total,
    aes(x = sim_label, y = value, fill = "LAML", group = method),
    position = dodge,
    width = 0.85
  ) +
  # Foreground: Metient (raw), BEAM, VINE
  geom_col(
    data = plot_df_raw,
    aes(x = sim_label, y = value, fill = method, group = method),
    position = dodge,
    width = 0.85
  ) +
  # Dashed separator before the Avg bar
  geom_vline(xintercept = length(sim_levels) + 0.5,
             color = "black", linewidth = 0.4) +
  scale_x_discrete(expand = expansion(add = c(1.5, 0.5))) +
  scale_fill_manual(
    values = method_palette,
    breaks = c("LAML", "Metient", "BEAM", "VINE")
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::label_comma(),
    limits = c(y_min, y_max),
    expand = expansion(mult = c(0, 0))
  ) +
  annotation_logticks(sides = "l") +
  labs(
    x     = "Clonal Population (CP)",
    y     = "Time (s)",
    title = ""
  ) +
  theme(
    legend.title  = element_blank(),
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
  )

save_pdf(p, outfile, width = 9, height = 3)