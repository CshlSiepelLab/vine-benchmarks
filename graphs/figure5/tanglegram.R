#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(readr)
})

# ---- helpers ----
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

# Place date legend just below the inner plot area (into the outer bottom margin).
# Must be called after plot.cophylo(); uses current par("usr") for coordinates.
add_date_legend <- function(leg_labels, leg_cols, ncols, cex_val, title_cex,
                             pt_cex_val = 1.1, ...) {
  usr <- par("usr")
  legend(x      = (usr[1] + usr[2]) / 2,
         y      = usr[3],
         legend = leg_labels,
         col    = leg_cols,
         pch    = 20, pt.cex = pt_cex_val,
         ncol   = ncols,
         bty    = "n", cex = cex_val,
         title  = "Collection date (2024)", title.cex = title_cex,
         xjust  = 0.5, yjust = 1,
         xpd    = NA, ...)
}

# ---- inputs ----
vine_mcc_file  <- "vine-small.mcc.nex"
beast_mcc_file <- "beast-small.mcc.nex"
meta_file      <- "small-subset.tsv"
out_full       <- "vine-vs-beast.tanglegram.pdf"
out_compact    <- "tanglegram.compact.pdf"

# ---- metadata ----
meta    <- read_tsv(meta_file, show_col_types = FALSE)
meta_df <- process_metadata(meta)

tips_2019 <- rownames(meta_df)[format(meta_df$date, "%Y") == "2019"]
tips_2024 <- rownames(meta_df)[format(meta_df$date, "%Y") == "2024"]
cat(sprintf("2019 outgroup tips: %d;  2024 tips in metadata: %d\n",
            length(tips_2019), length(tips_2024)))

# ---- read and root MCC trees ----
vine_mcc  <- read_mcc_phylo(vine_mcc_file)
beast_mcc <- read_mcc_phylo(beast_mcc_file)

vine_mcc_r  <- root_on_clade(vine_mcc,  tips_2019)
beast_mcc_r <- root_on_clade(beast_mcc, tips_2019)

# ---- prune to 2024 tips present in both trees ----
vine_2024  <- drop.tip(vine_mcc_r,  setdiff(vine_mcc_r$tip.label,  tips_2024))
beast_2024 <- drop.tip(beast_mcc_r, setdiff(beast_mcc_r$tip.label, tips_2024))

common <- intersect(vine_2024$tip.label, beast_2024$tip.label)
cat(sprintf("Common 2024 tips in both trees: %d\n", length(common)))

vine_2024  <- drop.tip(vine_2024,  setdiff(vine_2024$tip.label,  common))
beast_2024 <- drop.tip(beast_2024, setdiff(beast_2024$tip.label, common))

# ---- tip colors by collection date (plasma palette, Jan–Dec 2024) ----
tip_dates <- as.numeric(meta_df[common, "date"])
d_min <- as.numeric(as.Date("2024-01-01"))
d_max <- as.numeric(as.Date("2024-12-31"))
pal   <- hcl.colors(256, "plasma")
idx   <- pmin(256, pmax(1, round((tip_dates - d_min) / (d_max - d_min) * 255) + 1))
tip_col_named        <- pal[idx]
names(tip_col_named) <- common

# ---- build tanglegram once (rotate=TRUE minimises crossings) ----
co <- cophylo(vine_2024, beast_2024, rotate = TRUE)

# Colors in tip.label order (matches xx/yy[1:n_tips] in last_plot.cophylo)
tc_left  <- tip_col_named[co$trees[[1]]$tip.label]
tc_right <- tip_col_named[co$trees[[2]]$tip.label]

# ---- shared legend data ----
leg_dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "month")
leg_cols   <- pal[pmin(256, round((as.numeric(leg_dates) - d_min) /
                                  (d_max - d_min) * 255) + 1)]
leg_labels <- format(leg_dates, "%b")

# Every-other-month subset for compact legend (Jan Mar May Jul Sep Nov)
alt_idx    <- seq(1, 12, by = 2)
alt_labels <- leg_labels[alt_idx]
alt_cols   <- leg_cols[alt_idx]

# =========================================================
# FULL-PAGE VERSION (12 × 10")
# =========================================================
cairo_pdf(out_full, width = 12, height = 10)

# Outer margins keep labels/legend outside the tree area
par(oma = c(3.5, 0, 2, 0))

plot(co,
     ftype     = "i",
     fsize     = 0.4,
     link.type = "curved",
     link.lwd  = 0.4,
     link.col  = adjustcolor("gray50", alpha.f = 0.5),
     pts       = FALSE)

add_tip_points(co, tc_left, tc_right, pt_cex = 0.6)

# Panel labels in top outer margin (at = normalized device x, 0–1)
mtext("VINE",    outer = TRUE, side = 3, at = 0.25, line = 0.5,
      font = 2, cex = 1.2, adj = 0.5)
mtext("BEAST 2", outer = TRUE, side = 3, at = 0.75, line = 0.5,
      font = 2, cex = 1.2, adj = 0.5)

add_date_legend(leg_labels, leg_cols, ncols = 4, cex_val = 0.8,
                title_cex = 0.9, pt_cex_val = 1.4)

dev.off()
cat("Wrote:", out_full, "\n")

# =========================================================
# COMPACT VERSION (~3 × 3.5" panel — no sample names)
# =========================================================
cairo_pdf(out_compact, width = 3, height = 3.5)

# Top margin halved (0.9 lines) to reduce title-to-tree gap
par(oma = c(3.5, 0, 0.9, 0))

# Shared gray for all non-tip-point elements
gray_line <- adjustcolor("gray50", alpha.f = 0.5)
n_el <- nrow(co$trees[[1]]$edge)
n_er <- nrow(co$trees[[2]]$edge)

# par(col=gray_line): tip-extension lines inherit this color since phylogram
# draws them without an explicit col argument.  Explicitly colored elements
# (edge.col, link.col, tip points) are unaffected.
old_par <- par(col = gray_line, lwd = 0.15)  # tip-extension lines inherit both

plot(co,
     ftype     = "off",
     link.type = "curved",
     # lwd not passed → phylogram defaults to 1 (tree branches restored)
     link.lwd  = 1.0,              
     link.col  = gray_line,
     link.lty  = "32",
     tip.lty   = "32",
     edge.col  = list(left  = rep("black", n_el),
                      right = rep("black", n_er)),
     pts       = FALSE)

par(old_par)

add_tip_points(co, tc_left, tc_right, pt_cex = 1.4)

# Title-to-tree gap halved: line = 0.2 (was 0.4)
mtext("VINE",    outer = TRUE, side = 3, at = 0.25, line = 0.2,
      font = 2, cex = 0.75, adj = 0.5)
mtext("BEAST 2", outer = TRUE, side = 3, at = 0.75, line = 0.2,
      font = 2, cex = 0.75, adj = 0.5)

# All 12 months on one line, smaller font
add_date_legend(alt_labels, alt_cols, ncols = 6, cex_val = 0.42,
                title_cex = 0.5, pt_cex_val = 0.75)

dev.off()
cat("Wrote:", out_compact, "\n")
