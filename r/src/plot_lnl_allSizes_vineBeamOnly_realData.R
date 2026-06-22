#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 2) {
  stop("Usage: Rscript plot_lnl_byCp.R <input.tsv> <output.pdf>", call. = FALSE)
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
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
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
raw <- readLines(tsv_path)
# Drop separator line (dashes) and blank lines
clean_lines <- raw[!grepl("^-+$", raw) & nzchar(trimws(raw))]
df <- read.table(text = paste(clean_lines, collapse = "\n"),
                 header = TRUE, check.names = FALSE)

# Keep only numeric cp rows (drop "all" summary row if present)
df <- df[grepl("^[0-9]+$", as.character(df$cp)), ]

cp_levels <- df %>%
  pull(cp) %>%
  as.integer() %>%
  sort() %>%
  unique() %>%
  as.character()

df <- df %>%
  mutate(
    BEAM = beam_tree + beam_mig,
    VINE = vine_tree + vine_mig
  )

# Compute avg row before pivoting
avg_row <- df %>%
  summarise(
    cp   = "ave",
    BEAM = mean(BEAM, na.rm = TRUE),
    VINE = mean(VINE, na.rm = TRUE)
  )

cp_levels_ext <- c(cp_levels, "ave")

df_long <- bind_rows(
    df %>% mutate(cp = as.character(cp)),
    avg_row
  ) %>%
  select(cp, BEAM, VINE) %>%
  pivot_longer(-cp, names_to = "method", values_to = "lnl") %>%
  mutate(
    method = factor(method, levels = method_levels),
    cp     = factor(cp, levels = cp_levels_ext)
  )

# ------------ Plot ------------
p <- ggplot(df_long, aes(x = cp, y = lnl, fill = method)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  # Dashed separator before the Avg bar
  geom_vline(xintercept = length(cp_levels) + 0.5,
             color = "black", linewidth = 0.4) +
  scale_x_discrete(expand = expansion(add = c(1.5, 0.5))) +
  scale_fill_manual(values = method_palette, breaks = method_levels) +
  labs(x = "Clonal Population (CP)", y = "Log-Likelihood", fill = "Method") +
  theme(legend.title = element_blank())

save_pdf(p, outfile, width = 9, height = 3)