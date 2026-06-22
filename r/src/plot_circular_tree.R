suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
tree_name <- args[2]
meta_file <- args[3]
out_pdf   <- args[4]
primary_tissue <- args[5]

tree <- read.newick(file = tree_file)
tree$node.label <- NULL
tree <- multi2di(tree)
if (tree_name == "CASSIOPEIA-GREEDY") {
  tree <- compute.brlen(tree, method = "Grafen", power = 0.5)  # Helps to even out the branch lengths since Cass trees have none and they are plotted with a bias in length for old splits to be very long
}

meta <- read_csv(meta_file, col_names = c("cell", "tissue"), show_col_types = FALSE)
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$cell

tissues_sorted <- sort(unique(meta_df$tissue))
tissues <- c(primary_tissue, sort(setdiff(tissues_sorted, primary_tissue)))
color_palette <- c("#59A14F", "#E15759", "#F28E2B", "#B07AA1", "#76B7B2", "#EDC948",
                   "#FF9DA7", "#9C755F", "#BAB0AC", "#4E79A7")
tissue_pal <- setNames(rep(color_palette, length.out = length(tissues)), tissues)
tip_cols   <- tissue_pal[meta_df[tree$tip.label, "tissue"]]
tip_cols[is.na(tip_cols)] <- "grey60"

pdf(out_pdf, width = 7, height = 7)

# Use a margin so the ring has room — mar controls whitespace around plot
par(mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))

plot.phylo(tree, type = "fan", show.tip.label = FALSE,
           edge.color = "grey50", edge.width = 1.0,
           no.margin = FALSE)

# grab tip coordinates from ape's internal state
pp      <- get("last_plot.phylo", envir = ape:::.PlotPhyloEnv)
tip_x   <- pp$xx[seq_len(Ntip(tree))]
tip_y   <- pp$yy[seq_len(Ntip(tree))]
tip_ang <- atan2(tip_y, tip_x)

# ring geometry
r_tip   <- sqrt(tip_x^2 + tip_y^2)
r_inner <- max(r_tip) * 1.24
r_outer <- max(r_tip) * 1.30

# expand the plot window to fit the ring
usr <- par("usr")
r_needed <- r_outer * 1.02
par(usr = c(-r_needed, r_needed, -r_needed, r_needed))

# angular slice per tip
n_tips     <- Ntip(tree)
slice_half <- (2 * pi / n_tips) / 2

for (i in seq_len(n_tips)) {
  a0  <- tip_ang[i] - slice_half
  a1  <- tip_ang[i] + slice_half
  arc <- seq(a0, a1, length.out = 8)
  px  <- c(r_inner * cos(arc), rev(r_outer * cos(arc)))
  py  <- c(r_inner * sin(arc), rev(r_outer * sin(arc)))
  polygon(px, py, col = tip_cols[i], border = NA)
}

legend("bottomleft", legend = tissues, col = tissue_pal[tissues],
       pch = 20, bty = "n", cex = 1.0, title = "Tissue")

mtext(tree_name, outer = TRUE, side = 3, at=0.5, line = -1.2, font = 2, cex = 1.5, adj = 0.5)

dev.off()
cat("Wrote:", out_pdf, "\n")