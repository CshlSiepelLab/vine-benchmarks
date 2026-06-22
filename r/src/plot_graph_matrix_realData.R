#!/usr/bin/env Rscript

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 4) {
  stop(
    "Usage: Rscript plot_graph_matrix_realData.R <edge_file.csv> <threshold> <output.pdf> <primary_tissue>",
    call. = FALSE
  )
}
edge_file      <- args_trailing[1]
threshold      <- as.numeric(args_trailing[2])
outfile        <- args_trailing[3]
primary_tissue <- args_trailing[4]

# ------------ Theme ------------
base_theme <- theme_minimal(base_size = 8) +
  theme(
    plot.title       = element_text(size = 9, face = "bold"),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 7),
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA)
  )
theme_set(base_theme)

# ------------ Colors ------------
color_palette <- c("#59A14F", "#E15759", "#F28E2B", "#B07AA1", "#76B7B2", "#EDC948",
                   "#FF9DA7", "#9C755F", "#BAB0AC", "#4E79A7")

fill_palette <- c("#f7f7f7", "#f4a582", "#d6604d", "#b2182b")

# ------------ Read + parse ------------
edges_raw <- read_csv(
  edge_file,
  col_names = c("edge_id", "prob"),
  show_col_types = FALSE
) %>%
  mutate(
    source = str_extract(edge_id, "^[^_]+"),
    target = str_extract(edge_id, "(?<=_)[^_]+(?=_[^_]+$)")
  )

# Tissue ordering from the full file, before thresholding
all_tissues <- sort(unique(c(edges_raw$source, edges_raw$target)))
tissues <- c(primary_tissue, sort(setdiff(all_tissues, primary_tissue)))
tissue_pal <- setNames(rep(color_palette, length.out = length(tissues)), tissues)

# Now apply filtering for plotting
edges <- edges_raw %>%
  filter(prob > threshold) %>%
  filter(source != target)

if (nrow(edges) == 0) {
  p_empty <- ggplot() +
    annotate("text", x = 1, y = 1,
             label = paste0("No edges above threshold ", threshold),
             size = 3) +
    xlim(0.5, 1.5) +
    ylim(0.5, 1.5) +
    labs(x = NULL, y = NULL, title = "") +
    theme_void(base_size = 8)

  if (capabilities("cairo")) {
    ggsave(outfile, plot = p_empty, width = 3, height = 3,
           units = "in", device = cairo_pdf)
  } else {
    pdf(outfile, width = 3, height = 3, family = "Helvetica")
    print(p_empty)
    dev.off()
  }
  cat("Wrote:", outfile, "\n")
  quit(save = "no")
}

# ------------ Aggregate to full matrix ------------
agg <- edges %>%
  group_by(source, target) %>%
  summarise(
    n_edges   = n(),
    mean_prob = mean(prob),
    .groups   = "drop"
  ) %>%
  complete(
    source = tissues,
    target = tissues,
    fill = list(n_edges = 0, mean_prob = 0)
  ) %>%
  mutate(
    source = factor(source, levels = rev(tissues)),
    target = factor(target, levels = tissues),
    text_white = n_edges >= 0.5 * max(n_edges)
  )

# ------------ Plot ------------
p <- ggplot(agg, aes(x = target, y = source)) +
  geom_tile(
    aes(fill = n_edges),
    width = 0.92,
    height = 0.92,
    color = "grey60",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = ifelse(n_edges == 0, "", n_edges), color = text_white),
    size = 2.4,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "grey20"),
    guide = "none"
  ) +
  scale_fill_gradientn(
    colours = fill_palette,
    limits = c(0, 80),
    name = "Number of edges",
    guide = guide_colorbar(
      barwidth = 0.8,
      barheight = 6,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  scale_x_discrete(expand = c(0, 0), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Target tissue",
    y = "Source tissue",
    title = ""
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 1,
      face = "bold",
      color = tissue_pal[levels(agg$target)]
    ),
    axis.text.y = element_text(
      face = "bold",
      color = tissue_pal[levels(agg$source)]
    ),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(6, 6, 6, 6)
  )

# ------------ Save ------------
if (capabilities("cairo")) {
  ggsave(outfile, plot = p, width = 4, height = 3,
         units = "in", device = cairo_pdf)
} else {
  pdf(outfile, width = 4, height = 3, family = "Helvetica")
  print(p)
  dev.off()
}

cat("Wrote:", outfile, "\n")