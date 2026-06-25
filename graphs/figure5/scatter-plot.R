#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpointdensity)
})

vine_file  <- "vine-small.dist"
beast_file <- "beast-small.dist"
out_file   <- "scatter.pdf"

vine.dist  <- read_table(vine_file,  show_col_types = FALSE)
beast.dist <- read_table(beast_file, show_col_types = FALSE)

df <- tibble(
  beast = beast.dist$mean,
  vine  = vine.dist$mean
)

r  <- cor(df$beast, df$vine, use = "complete.obs")
r2 <- r^2

p <- ggplot(df, aes(x = beast, y = vine)) +

        # density-colored points
        geom_pointdensity(size = 1.0) +

        # perceptually uniform color scale
        scale_color_viridis_c(
                option = "magma",
                name = "Point density"
        ) +

        # identity line (very important scientifically)
        geom_abline(
                slope = 1, intercept = 0,
                linetype = "dashed",
                linewidth = 0.4,
                color = "gray40"
        ) +
        labs(
                x = "BEAST 2 distance",
                y = "VINE distance"
        ) +
        coord_equal() +
        theme_classic(base_size = 7) +
        theme(
                # place legend inside panel
                legend.position = c(0.03, 0.97),
                legend.justification = c(0, 1),

                # smaller text
                legend.title = element_text(size = 6),
                legend.text = element_text(size = 5),

                # shrink color bar
                legend.key.height = unit(0.5, "lines"),
                legend.key.width = unit(0.5, "lines"),

                # reduce padding everywhere
                legend.spacing.y = unit(0.1, "lines"),
                legend.margin = margin(1, 1, 1, 1),
                legend.box.margin = margin(0, 0, 0, 0),

                # optional subtle background for readability
                legend.background = element_rect(
                        fill = scales::alpha("white", 0.8),
                        color = NA
                )
        ) +
        guides(
                color = guide_colorbar(
                        barheight = unit(10, "mm"),
                        barwidth = unit(2.5, "mm"),
                        ticks = FALSE
                )
        ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = sprintf("R\u00B2 = %.3f", r2),  # R²
    hjust = 1.05,
    vjust = -0.8,
    size = 2.5
  )

ggsave(out_file, p, width = 2, height = 2, device = cairo_pdf)
cat("Wrote:", out_file, "\n")
