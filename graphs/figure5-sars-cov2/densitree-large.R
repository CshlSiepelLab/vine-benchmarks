#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(grid)
  library(treeio)
  library(cowplot)
  library(phangorn)
})

# ---- hard compatibility shim: fixes "could not find function is.waive" ----
if (!exists("is.waive", mode = "function")) {
        is.waive <- function(x) inherits(x, "waiver")
}

process_metadata <- function(meta) {
  meta <- meta[!duplicated(meta$strain), ]
  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df$strain
  meta_df$date <- as.Date(meta_df$date)
  meta_df
}

read_mcc_phylo <- function(path) {
  tr <- read.nexus(path)
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  tr
}

root_on_clade <- function(tr, outgroup_tips) {
  og <- outgroup_tips[outgroup_tips %in% tr$tip.label]
  ape::root(tr, outgroup = og[1], resolve.root = TRUE)
}

prep_cloud_trees <- function(ts, tip_order, n_sub = 100, root_fn = phangorn::midpoint) {
        if (inherits(ts, "phylo")) ts <- structure(list(ts), class = "multiPhylo")
        n_sub <- min(n_sub, length(ts))
        idx <- sample.int(length(ts), n_sub)
        ts_sub <- ts[idx]

        ts_sub2 <- lapply(ts_sub, function(tr) {
                tr <- root_fn(tr)
                tr <- ape::rotateConstr(tr, tip_order)
                tr
        })
        class(ts_sub2) <- "multiPhylo"
        ts_sub2
}

make_cloud_mcc_panel <- function(ts_cloud, mcc_tree, meta_df,
                                 panel_title = NULL,
                                 cloud_alpha = 0.025,
                                 cloud_lwd = 0.2,
                                 mcc_lwd = 0.2,
                                 tip_size = 1.0,
                                 year_breaks = as.Date(paste0(2019:2026, "-07-01")),
                                 year_labels = 2019:2026) {
        p_cloud <- suppressWarnings(ggdensitree(ts_cloud, layout = "circular")) +
                theme_void() +
                theme(plot.margin = margin(0, 0, 0, 0), legend.position = "none")

        for (i in seq_along(p_cloud$layers)) {
                lyr <- p_cloud$layers[[i]]
                if (inherits(lyr$geom, "GeomSegment") || inherits(lyr$geom, "GeomPath") || inherits(lyr$geom, "GeomTree")) {
                        lyr$aes_params$alpha <- cloud_alpha
                        lyr$aes_params$colour <- "grey60"
                        if (!is.null(lyr$aes_params$linewidth)) {
                                lyr$aes_params$linewidth <- cloud_lwd
                        } else {
                                lyr$aes_params$size <- cloud_lwd
                        }
                        p_cloud$layers[[i]] <- lyr
                }
        }

        # MCC: use ggdensitree with a single tree so the coordinate space
        # matches the cloud exactly (ggtree and ggdensitree use opposite x conventions)
        p_mcc <- (suppressWarnings(ggdensitree(
                        structure(list(mcc_tree), class = "multiPhylo"),
                        layout = "circular"
                )) %<+% meta_df) +
                geom_tippoint(aes(color = date), size = tip_size) +
                scale_color_viridis_c(
                        name = "Year",
                        option = "plasma",
                        limits = as.numeric(as.Date(c("2019-01-01", "2026-12-31"))),
                        breaks = year_breaks,
                        labels = year_labels,
                        guide = guide_colorbar(
                                title.position = "top",
                                title.hjust = 0.5,
                                barwidth = unit(6, "cm"),
                                barheight = unit(0.4, "cm")
                        )
                ) +
          theme_void() +
          theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box.margin = margin(0, 0, 0, 0),
            legend.margin     = margin(-12, 0, 0, 0)
          )

        # Style MCC branches: black, fully opaque, thin
        for (i in seq_along(p_mcc$layers)) {
                lyr <- p_mcc$layers[[i]]
                if (inherits(lyr$geom, "GeomSegment") || inherits(lyr$geom, "GeomPath") || inherits(lyr$geom, "GeomTree")) {
                        lyr$aes_params$colour <- "black"
                        lyr$aes_params$alpha  <- 1.0
                        if (!is.null(lyr$aes_params$linewidth)) {
                                lyr$aes_params$linewidth <- mcc_lwd
                        } else {
                                lyr$aes_params$size <- mcc_lwd
                        }
                        p_mcc$layers[[i]] <- lyr
                }
        }

  # Both plots now share ggdensitree's coordinate convention — use identical scales
  tight_scales <- list(
    scale_x_continuous(expand = expansion(mult = 0, add = 0)),
    scale_y_continuous(expand = expansion(mult = 0, add = 0))
  )
  suppressMessages({
    p_cloud <- p_cloud + tight_scales
    p_mcc   <- p_mcc   + tight_scales
  })

  p_overlay <- ggdraw() +
    draw_plot(p_cloud, 0, 0, 1, 1) +
    draw_plot(p_mcc + theme(legend.position = "none"), 0, 0, 1, 1)

  if (!is.null(panel_title)) {
    p_overlay <- p_overlay +
      draw_label(panel_title,
                 x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
                 fontface = "bold", size = 12
                 )
  }

  list(panel = p_overlay, mcc_plot_for_legend = p_mcc)
}

# ---- inputs ----
vine_trees_file <- "vine-large.nex"
vine_mcc_file   <- "vine-large.mcc.nex"
meta_file       <- "large-subset.tsv"
out_pdf         <- "vine-large.densitree.pdf"

set.seed(1)
n_cloud <- 100

# ---- read metadata ----
meta <- read_tsv(meta_file, show_col_types = FALSE)
meta_df <- process_metadata(meta)
meta_2019 <- meta_df[format(meta_df$date, "%Y") == "2019", ]
tips_2019 <- rownames(meta_2019)[order(meta_2019$date)]  # earliest-first

# ---- read trees ----
vine_ts  <- read.nexus(vine_trees_file)
vine_mcc <- read_mcc_phylo(vine_mcc_file)

# ---- root MCC tree on 2019 Wuhan outgroup ----
vine_mcc_r <- root_on_clade(vine_mcc, tips_2019)

# ---- canonical tip order from MCC ----
canonical_order <- vine_mcc_r$tip.label

# ---- prep cloud trees ----
root_fn    <- function(tr) root_on_clade(tr, tips_2019)
vine_cloud <- prep_cloud_trees(vine_ts, canonical_order, n_sub = n_cloud, root_fn = root_fn)

# ---- build panel (no title) ----
panel_obj <- make_cloud_mcc_panel(vine_cloud, vine_mcc_r, meta_df)

# ---- scale bar ----
scale_w     <- 0.001
tree_height <- max(node.depth.edgelength(vine_mcc_r))

sb <- ggplot() +
  annotate("segment",
           x = 0, xend = scale_w,
           y = 0, yend = 0,
           linewidth = 0.6) +
  annotate("text",
           x = scale_w / 2,
           y = -0.25,
           label = paste0(scale_w, " substitutions/site"),
           size = 3) +
  coord_cartesian(
    xlim = c(0, tree_height),
    ylim = c(-1, 1),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void()

top <- ggdraw(panel_obj$panel) +
        draw_plot(sb,
                x = 0.12,
                y = 0.11,
                width = 0.5,   # single panel (vine_frac=1): sb_width = vine_frac/2 = 0.5
                height = 0.14
                )

# ---- legend ----
leg <- cowplot::get_legend(panel_obj$mcc_plot_for_legend)

final <- plot_grid(
  top,
  leg,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

ggsave(out_pdf, final, width = 8, height = 8, device = cairo_pdf)
cat("Wrote:", out_pdf, "\n")
