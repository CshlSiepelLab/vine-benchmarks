#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(uwot)
  library(ggrepel)
  library(ggrastr)
  library(cowplot)
})

args      <- commandArgs(trailingOnly = TRUE)
emb_file  <- if (length(args) >= 1) args[1] else "~/sars-cov2/vine-large.emb.tsv"
meta_file <- if (length(args) >= 2) args[2] else "~/sars-cov2/small-subset.tsv"
out_file  <- if (length(args) >= 3) args[3] else "embedding-umap.pdf"

# colorbar settings matching densitree plots
YEAR_BREAKS <- as.Date(paste0(2019:2026, "-07-01"))
YEAR_LABELS <- 2019:2026

# ---- read embedding: col 1 = taxon name, remaining = embedding dimensions ----
emb    <- read_tsv(emb_file, col_names = FALSE, show_col_types = FALSE)
taxa   <- emb[[1]]
coords <- as.matrix(emb[, -1])

# ---- read metadata; join date by taxon name ----
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(!duplicated(strain)) %>%
  mutate(date = as.Date(date))

# ---- UMAP with 3 components ----
set.seed(1)
umap_xyz <- umap(coords,
                 n_components = 3L,
                 n_neighbors  = min(15L, nrow(coords) - 1L),
                 min_dist     = 0.1,
                 metric       = "euclidean",
                 verbose      = FALSE)

df <- tibble(taxon = taxa,
             UMAP1 = umap_xyz[, 1],
             UMAP2 = umap_xyz[, 2],
             UMAP3 = umap_xyz[, 3]) %>%
  left_join(meta %>% select(strain, date), by = c("taxon" = "strain"))

# ---- shared color scale ----
color_scale <- scale_color_viridis_c(
  name   = "Year",
  option = "plasma",
  limits = as.numeric(as.Date(c("2019-01-01", "2026-12-31"))),
  breaks = YEAR_BREAKS,
  labels = YEAR_LABELS,
  guide  = guide_colorbar(
    title.position = "top",
    title.hjust    = 0.5,
    barwidth       = unit(5, "cm"),
    barheight      = unit(0.4, "cm")
  )
)

panel_theme <- list(
  theme_classic(base_size = 7),
  theme(
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    axis.line         = element_line(linewidth = 0.3),
    legend.position   = "none"
  )
)

# ---- three panels: all pairs ----
p12 <- ggplot(df, aes(UMAP1, UMAP2, color = date)) +
  geom_point_rast(size = 1.2, alpha = 0.85, raster.dpi = 300) +
  color_scale + labs(x = "UMAP 1", y = "UMAP 2") + panel_theme

p13 <- ggplot(df, aes(UMAP1, UMAP3, color = date)) +
  geom_point_rast(size = 1.2, alpha = 0.85, raster.dpi = 300) +
  color_scale + labs(x = "UMAP 1", y = "UMAP 3") + panel_theme

p23 <- ggplot(df, aes(UMAP2, UMAP3, color = date)) +
  geom_point_rast(size = 1.2, alpha = 0.85, raster.dpi = 300) +
  color_scale + labs(x = "UMAP 2", y = "UMAP 3") + panel_theme

final <- plot_grid(p12, p13, p23, nrow = 1)

ggsave(out_file, final, width = 7.5, height = 2.5, device = cairo_pdf)
cat("Wrote:", out_file, "\n")
