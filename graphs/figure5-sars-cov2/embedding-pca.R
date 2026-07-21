#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrastr)
})

emb_file  <- "vine-large.emb.tsv"
meta_file <- "large-subset.tsv"
out_file  <- "embedding-pca.pdf"

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

# ---- PCA ----
pca <- prcomp(coords, center = TRUE, scale. = FALSE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

df <- tibble(taxon = taxa,
             PC1 = pca$x[, 1],
             PC2 = pca$x[, 2]) %>%
  left_join(meta %>% select(strain, date), by = c("taxon" = "strain"))

# ---- plot ----
p <- ggplot(df, aes(PC1, PC2, color = date)) +
  geom_point_rast(size = 1.2, alpha = 0.85, raster.dpi = 300) +
  scale_color_viridis_c(
    name   = "Year",
    option = "plasma",
    limits = as.numeric(as.Date(c("2019-01-01", "2026-12-31"))),
    breaks = YEAR_BREAKS,
    labels = YEAR_LABELS,
    guide  = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(4, "cm"),
      barheight      = unit(0.35, "cm")
    )
  ) +
  labs(x = sprintf("PC1 (%.1f%%)", pct[1]),
       y = sprintf("PC2 (%.1f%%)", pct[2])) +
  theme_classic(base_size = 6) +
  theme(
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    axis.line         = element_line(linewidth = 0.3),
    legend.position   = "none"
  )

ggsave(out_file, p, width = 2, height = 2, device = cairo_pdf)
cat("Wrote:", out_file, "\n")
