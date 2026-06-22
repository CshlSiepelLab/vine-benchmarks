#!/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(readr)
})


read_tree <- function(file) {
  tree <- tryCatch(
    read.newick(file = file),
    error = function(e) NULL
  )
  if (is.null(tree) || is.null(tree$tip.label)) {
    tree <- read.nexus(file = file)
    if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  }
  tree
}


# After plot.cophylo(), phytools stores tip coordinates in ape's .PlotPhyloEnv
# as last_plot.cophylo$left$xx/yy and $right$xx/yy.
# Rows 1:n_tips are tips in tip.label order (standard R phylo convention).
add_tip_points <- function(co, tc_left, tc_right, pt_cex) {
  env  <- ape:::.PlotPhyloEnv
  coph <- get("last_plot.cophylo", envir = env)
  n1   <- Ntip(co$trees[[1]])
  n2   <- Ntip(co$trees[[2]])
  L    <- coph$left
  R    <- coph$right
  points(L$xx[1:n1], L$yy[1:n1], pch = 20, col = tc_left,  cex = pt_cex)
  points(R$xx[1:n2], R$yy[1:n2], pch = 20, col = tc_right, cex = pt_cex)
}


# Place tissue legend just below the inner plot area (into the outer bottom margin).
# Must be called after plot.cophylo(); uses current par("usr") for coordinates.
add_tissue_legend <- function(leg_labels, leg_cols, ncols, cex_val, title_cex,
                              pt_cex_val = 1.1, ...) {
  usr <- par("usr")
  legend(x      = (usr[1] + usr[2]) / 2,
         y      = usr[3],
         legend = leg_labels,
         col    = leg_cols,
         pch    = 20, pt.cex = pt_cex_val,
         ncol   = ncols,
         bty    = "n", cex = cex_val,
         title  = "Tissue", title.cex = title_cex,
         xjust  = 0.5, yjust = 1,
         xpd    = NA, ...)
}


add_scale_bars <- function(co, bar_length_left = NULL, bar_length_right = NULL,
                           label_left = NULL, label_right = NULL,
                           cex_val = 0.45, lwd_val = 1.2) {
  env  <- ape:::.PlotPhyloEnv
  coph <- get("last_plot.cophylo", envir = env)
  usr  <- par("usr")

  L <- coph$left
  R <- coph$right

  lx_range <- range(L$xx)
  rx_range <- range(R$xx)

  # Real tree depths (in actual branch length units)
  depth_left  <- max(node.depth.edgelength(co$trees[[1]]))
  depth_right <- max(node.depth.edgelength(co$trees[[2]]))

  # Plotting span / real depth = pixels-per-unit scaling factor
  scale_left  <- diff(lx_range) / depth_left
  scale_right <- diff(rx_range) / depth_right

  # Auto-pick round bar lengths in REAL units if not specified
  if (is.null(bar_length_left))
    bar_length_left  <- signif(depth_left  * 0.25, 1)
  if (is.null(bar_length_right))
    bar_length_right <- signif(depth_right * 0.25, 1)

  # Convert real-unit bar lengths to plotting coordinates
  bar_plot_left  <- bar_length_left  * scale_left
  bar_plot_right <- bar_length_right * scale_right

  # Default labels
  if (is.null(label_left))  label_left  <- as.character(bar_length_left)
  if (is.null(label_right)) label_right <- as.character(bar_length_right)

  y_pos <- usr[3] - diff(c(usr[3], usr[4])) * 0.001

  # Left scale bar (anchored at left/root side)
  lx_start <- lx_range[1]
  lx_end   <- lx_start + bar_plot_left
  segments(lx_start, y_pos, lx_end, y_pos, xpd = NA, lwd = lwd_val)
  text((lx_start + lx_end) / 2, y_pos,
       labels = label_left, pos = 1, xpd = NA, cex = cex_val)

  # Right scale bar (anchored at right/root side)
  rx_start <- rx_range[2]
  rx_end   <- rx_start - bar_plot_right
  segments(rx_start, y_pos, rx_end, y_pos, xpd = NA, lwd = lwd_val)
  text((rx_start + rx_end) / 2, y_pos,
       labels = label_right, pos = 1, xpd = NA, cex = cex_val)
}


# ---- inputs ----
args <- commandArgs(trailingOnly = TRUE)
tree1_file <- args[1]
tree1_name <- args[2]
tree2_file <- args[3]
tree2_name  <- args[4]
meta_file  <- args[5]
out_full   <- args[6]
primary_tissue <- args[7]

# ---- metadata ----
# Expected columns: cell, tissue
meta <- read_csv(meta_file, col_names = c("cell", "tissue"), show_col_types = FALSE)
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$cell

# ---- read trees ----
tree1 <- read_tree(tree1_file)
tree2 <- read_tree(tree2_file)

# Remove node labels
tree1$node.label <- NULL
tree2$node.label <- NULL

# Resolve any polytomies
tree1 <- multi2di(tree1)
tree2 <- multi2di(tree2)

# ---- find common tips ----
common <- intersect(tree1$tip.label, tree2$tip.label)
if (length(common) == 0) stop("No shared tip labels between the two trees.")

# Prune both trees to common tips only
tree1 <- keep.tip(tree1, common)
tree2 <- keep.tip(tree2, common)

# ---- tip colors by tissue ----
tissues <- sort(unique(meta_df[common, "tissue"]))
n_tissues <- length(tissues)
color_palette <- c("#59A14F", "#E15759", "#F28E2B", "#B07AA1", "#76B7B2", "#EDC948",
                   "#FF9DA7", "#9C755F", "#BAB0AC", "#4E79A7")
tissue_pal <- setNames(color_palette[1:n_tissues], tissues)

tip_col_named <- tissue_pal[meta_df[common, "tissue"]]
names(tip_col_named) <- common

# Patch phytools' rotate.multi to skip NA or out-of-range nodes
safe_rotate_multi <- function(tree, nn) {
  for (i in seq_along(nn)) {
    node <- nn[i]
    if (is.na(node)) next
    if (node < 1 || node > (Ntip(tree) + Nnode(tree))) next
    tree <- tryCatch(
      rotate(tree, node),
      error = function(e) tree  # skip bad nodes silently
    )
    tree <- untangle(tree, "read.tree")
  }
  tree
}

# Overwrite the function in phytools' namespace
assignInNamespace("rotate.multi", safe_rotate_multi, ns = "phytools")

# Plot the tanglegram with auto-rotations to minimize crossings; tip colors by tissue.
co <- cophylo(tree1, tree2, rotate = TRUE, rotate.multi = TRUE)

# Colors in tip.label order (matches xx/yy[1:n_tips] in last_plot.cophylo)
tc_left  <- tip_col_named[co$trees[[1]]$tip.label]
tc_right <- tip_col_named[co$trees[[2]]$tip.label]

# ---- legend data ----
leg_labels <- c(primary_tissue, sort(setdiff(tissues, primary_tissue)))
leg_cols   <- tissue_pal[tissues]
ncols_leg  <- min(5, ceiling(n_tissues / 2))   # auto-wrap for many tissues


pdf(out_full, width = 3, height = 3.5)

# Outer margins keep labels/legend outside the tree area
par(oma = c(3.5, 0, 0.5, 0))

old_par <- par(col = adjustcolor("gray50", alpha.f = 0.5), lwd = 0.15)

plot(co,
     ftype     = "off",
     link.type = "curved",
     link.lwd  = 1.0,
     link.col  = adjustcolor("gray50", alpha.f = 0.5),
     link.lty  = "32",
     tip.lty   = "32",
     tip.col   = rep("gray50", Ntip(co$trees[[1]])),
     pts       = FALSE)
par(old_par)

add_tip_points(co, tc_left, tc_right, pt_cex = 0.8)
add_scale_bars(co)

# Panel labels in top outer margin
mtext(tree1_name, outer = TRUE, side = 3, at = 0.25, line = -0.5,
      font = 2, cex = 0.75, adj = 0.5)
mtext(tree2_name, outer = TRUE, side = 3, at = 0.75, line = -0.5,
      font = 2, cex = 0.75, adj = 0.5)

add_tissue_legend(leg_labels, leg_cols,
                  ncols     = ncols_leg,
                  cex_val   = 0.5,
                  title_cex = 0.5,
                  pt_cex_val = 0.8)

dev.off()
cat("Wrote:", out_full, "\n")
