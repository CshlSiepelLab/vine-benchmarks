#!/usr/bin/env Rscript
# Supplement figure: Convergence comparison: vine's log-likelihood and BEAST's # log-likelihood vs. wall-clock time (linear scale, full run length)
# Each panel contains: a zoomed-in view of the near-convergence region on the   # left, a view of the tail (including the true run endpoint) on the right, with
# "..." marks at the seam (top and bottom) flagging the break -- this lets
# the actual early convergence trajectory stay legible while still showing
# how long each run actually went. Two representative replicates per panel,
# shown as solid vs. dashed lines.
# vine's final log-likelihood is also extended as a light reference line
# across the whole panel (including the tail, where vine has no data of its
# own) so BEAST's eventual value can be compared against it directly.
#
# Usage: Rscript makeLoglikGraphs.R

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))

# Resolve script directory so we can build absolute paths
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()

data_dir <- file.path(script_dir, "data")
out_dir  <- script_dir

# Always compares against plain BEAST2 (no BEAGLE).
mcmc_label <- "BEAST2"

save_pdf <- function(plot, filename, width = 3, height = 3) {
  if (capabilities("cairo")) {
    ggsave(filename, plot = plot, width = width, height = height,
           units = "in", device = cairo_pdf)
  } else {
    pdf(filename, width = width, height = height, family = "Helvetica")
    print(plot); dev.off()
  }
}

theme_set(
  theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(
      plot.title  = element_text(size = 9, face = "bold"),
      axis.title  = element_text(size = 9),
      axis.text   = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7),
      panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.2)
    )
)

vine_col <- "#B35806"        # dark vine orange
vine_light_col <- "#FDBE85"  # light vine orange, for the reference line only
method_pal <- c("vine" = vine_col,
                "BEAST2" = "#1B7837")   # dark BEAST2 green
method_levels <- c("vine", mcmc_label)

TAXA_SUBSET <- c(50, 100, 250)

# ================= Convergence trajectories =================
traj <- read.csv(file.path(data_dir, "convergence_traj.csv"), stringsAsFactors = FALSE)
traj <- subset(traj, method %in% method_levels & ntaxa %in% TAXA_SUBSET)
traj$method <- factor(traj$method, levels = method_levels)

# Only two replicates, distinguished by linetype (solid vs. dashed) rather
# than by shading
reps_kept <- head(sort(unique(traj$rep)), 2)
traj <- subset(traj, rep %in% reps_kept)
rep_levels <- paste0("rep ", reps_kept)
traj$rep_label <- factor(paste0("rep ", traj$rep), levels = rep_levels)
rep_linetypes <- setNames(c("solid", "dashed")[seq_along(reps_kept)], rep_levels)

taxa_sizes <- sort(unique(traj$ntaxa))
tags <- LETTERS[seq_along(taxa_sizes)]

# Time (within a single replicate's own trajectory) at which it reaches ~97%
# of its total improvement
time_to_near_convergence <- function(d, tol = 0.03) {
  keys <- unique(d[, c("method", "rep")])
  t_conv <- numeric(nrow(keys))
  for (i in seq_len(nrow(keys))) {
    sub <- d[d$method == keys$method[i] & d$rep == keys$rep[i], ]
    sub <- sub[order(sub$time_sec), ]
    n <- nrow(sub)
    final <- mean(tail(sub$loglik, max(1, round(n * 0.05))))
    initial <- sub$loglik[1]
    thresh <- final - tol * (final - initial)
    idx <- which(sub$loglik >= thresh)[1]
    t_conv[i] <- if (is.na(idx)) max(sub$time_sec) else sub$time_sec[idx]
  }
  max(t_conv, na.rm = TRUE)
}

# Downsample each (method, rep) line to at most max_points evenly-spaced
# vertices, for DISPLAY ONLY (the underlying data/CSV is untouched -- this
# just affects what geom_line draws)
thin_for_display <- function(d, max_points = 150) {
  d %>%
    group_by(method, rep) %>%
    group_modify(function(sub, key) {
      sub <- sub[order(sub$time_sec), ]
      n <- nrow(sub)
      if (n <= max_points) return(sub)
      idx <- unique(round(seq(1, n, length.out = max_points)))
      sub[idx, ]
    }) %>%
    ungroup()
}

# "..." marks pinned to a panel's true edge (Inf/-Inf), not a data x-value,
# so they sit exactly at the seam and are never clipped. One pair near the
# x-axis (bottom) and one pair near the top of the panel, so the break reads
# as cutting through the whole plot, not just the axis.
mark_edge <- function(edge_x, hjust) {
  list(
    annotate("text", x = edge_x, y = -Inf, label = "⋯",
             hjust = hjust, vjust = -0.3, size = 3, color = "gray35"),
    annotate("text", x = edge_x, y = Inf, label = "⋯",
             hjust = hjust, vjust = 1.1, size = 3, color = "gray35")
  )
}

make_panel <- function(nt, tag, show_y = FALSE) {
  d <- subset(traj, ntaxa == nt)
  break_start <- time_to_near_convergence(d) * 1.15
  break_end <- max(d$time_sec) * 0.9
  y_range <- range(d$loglik)
  y_pad <- diff(y_range) * 0.04
  ylims <- c(y_range[1] - y_pad, y_range[2] + y_pad)

  # vine's final (converged) log-likelihood per replicate, extended as a
  # light reference line across the whole panel -- including the compressed
  # tail segment, where vine has no data of its own -- so BEAST's eventual
  # log-likelihood can be visually compared against where vine ended up.
  vine_final <- d %>%
    filter(method == "vine") %>%
    group_by(rep_label) %>%
    arrange(time_sec) %>%
    summarise(loglik = mean(tail(loglik, max(1, round(n() * 0.05)))), .groups = "drop")

  common <- list(
    geom_hline(data = vine_final, aes(yintercept = loglik, linetype = rep_label),
               color = vine_light_col, linewidth = 0.5, show.legend = FALSE,
               inherit.aes = FALSE),
    geom_line(linewidth = 0.4, alpha = 0.85),
    scale_color_manual(values = method_pal, breaks = method_levels, drop = FALSE),
    scale_linetype_manual(values = rep_linetypes, breaks = rep_levels, drop = FALSE),
    coord_cartesian(ylim = ylims),
    guides(color = guide_legend(override.aes = list(linewidth = 1)),
           linetype = guide_legend(override.aes = list(linewidth = 0.6)))
  )

  p_left <- ggplot(subset(d, time_sec <= break_start),
                    aes(time_sec, loglik, color = method, linetype = rep_label,
                        group = interaction(rep, method))) +
    common +
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    mark_edge(Inf, hjust = 1.1) +
    labs(tag = tag,
         x = NULL, y = if (show_y) "Log-likelihood" else NULL,
         color = NULL, linetype = NULL) +
    theme(plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
          plot.tag.position = c(0.01, 0.99),
          plot.margin = margin(t = 4, r = 8, b = 2, l = 2), legend.position = "none")

  p_right <- ggplot(thin_for_display(subset(d, time_sec >= break_end)),
                     aes(time_sec, loglik, color = method, linetype = rep_label,
                         group = interaction(rep, method))) +
    common +
    scale_x_continuous(breaks = pretty_breaks(n = 2)) +
    mark_edge(-Inf, hjust = -0.1) +
    labs(x = NULL, y = NULL, color = NULL, linetype = NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(t = 4, r = 6, b = 2, l = 8),
          legend.position = "none")

  (p_left + p_right) +
    plot_layout(widths = c(3, 1)) +
    plot_annotation(caption = "Wall-clock time (s)") &
    theme(plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 2)))
}

panels <- Map(function(nt, tag, i) {
  make_panel(nt, tag, show_y = (i %% 3 == 1))
}, taxa_sizes, tags, seq_along(taxa_sizes))

legend_ref <- ggplot(traj, aes(time_sec, loglik, color = method, linetype = rep_label)) +
  geom_line() +
  scale_color_manual(values = method_pal, breaks = method_levels) +
  scale_linetype_manual(values = rep_linetypes, breaks = rep_levels) +
  labs(color = "Method", linetype = "Replicate") +
  guides(color = guide_legend(order = 1, override.aes = list(linewidth = 1)),
         linetype = guide_legend(order = 2, override.aes = list(linewidth = 0.6))) +
  theme(legend.position = "right")
legend_grob <- cowplot::get_legend(legend_ref)

fig <- (wrap_plots(panels, ncol = 3) | wrap_elements(legend_grob)) +
  plot_layout(widths = c(1, 0.12))

out_stem <- "vine_vs_beast_loglik"
save_pdf(fig, file.path(out_dir, paste0(out_stem, ".pdf")), width = 10.5, height = 3.6)
ggsave(file.path(out_dir, paste0(out_stem, ".png")),
       plot = fig, width = 10.5, height = 3.6, units = "in", dpi = 220)
cat("wrote", out_stem, ".pdf / .png\n")
