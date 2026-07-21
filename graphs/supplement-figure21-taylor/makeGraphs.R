#!/usr/bin/env Rscript
# Supplement figure: accuracy of the Taylor-approximated ELBO vs. a 100-sample
# Monte Carlo estimate. Four panels, two per dataset size (25 / 50 taxa):
#   Top row  (A, B): ELBO convergence trajectory for one representative replicate,
#                    Taylor vs. Monte Carlo, showing similar convergence rate to
#                    essentially the same optimum (differences invisible at the
#                    full scale of the ELBO).
#   Bottom row (C, D): ELBO at the maximum, Taylor - Monte Carlo at the SAME
#                    converged (mu, Sigma), one point per replicate (n = 10), with
#                    the mean +/- SD band. Same quantity as the top row, zoomed
#                    ~1000x to make the residual visible.
# Style, palette, fonts, panel sizes, and A/B/C/D callouts match
# figure2-hky300/makeGraphs.R and supplement-figure-centering/makeGraphs.R.

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

save_pdf <- function(plot, filename, width = 3, height = 3) {
  if (capabilities("cairo")) {
    ggsave(filename, plot = plot, width = width, height = height,
           units = "in", device = cairo_pdf)
  } else {
    pdf(filename, width = width, height = height, family = "Helvetica")
    print(plot); dev.off()
  }
}

# ------------ Theme (identical to figure 2) ------------
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

taylor_col <- "#F28E2B"   # vine orange, matching figure 2
mc_col     <- "#4E79A7"   # blue, matching figure 2's NJ hue
method_pal <- c("Taylor" = taylor_col, "Monte Carlo" = mc_col)

# ================= Top row: convergence trajectories =================
traj <- read.csv(file.path(data_dir, "elbo_traj.csv"), stringsAsFactors = FALSE)
traj$method <- factor(traj$method, levels = c("Taylor", "Monte Carlo"))

# Site subsampling is active during the SGA warm-up: the scheduler's batch grows
# 256 -> 278 (iter 20) -> full data (iter 40), deterministically, so the ELBO is
# noticeably noisier over iterations 1-40 in every run. Shade that region gray.
WARMUP_END <- 40

make_traj <- function(nt, tag, show_y = FALSE) {
  d <- subset(traj, ntaxa == nt)
  ggplot(d, aes(iter, elbo, color = method,
                group = interaction(rep, method))) +
    annotate("rect", xmin = -Inf, xmax = WARMUP_END + 0.5,
             ymin = -Inf, ymax = Inf, fill = "gray55", alpha = 0.16) +
    geom_line(linewidth = 0.4, alpha = 0.85) +
    scale_color_manual(values = method_pal, breaks = c("Taylor", "Monte Carlo")) +
    scale_x_continuous(breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    labs(title = tag,
         x = "Iteration",
         y = if (show_y) "ELBO" else NULL,
         color = NULL) +
    guides(color = guide_legend(override.aes = list(linewidth = 1))) +
    theme(legend.position = c(0.98, 0.05), legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key.height = unit(9, "pt"))
}

# ================= Bottom row: ELBO-at-max differences =================
dat <- read.csv(file.path(data_dir, "elbo_max.csv"), stringsAsFactors = FALSE)
dat$diff <- dat$elbo_taylor - dat$elbo_mc   # Taylor - Monte Carlo, same (mu,Sigma)

summ <- dat %>%
  group_by(ntaxa) %>%
  summarise(mean = mean(diff), sd = sd(diff),
            absmean = mean(abs(diff)),
            rel = mean(abs(diff)) / mean(abs(elbo_mc)) * 100,
            .groups = "drop")

make_diff <- function(nt, tag, show_y = FALSE) {
  d <- subset(dat, ntaxa == nt)
  s <- subset(summ, ntaxa == nt)
  xr <- range(d$rep)
  lab <- sprintf("mean %+.2f ± %.2f  (%.3f%% of ELBO)", s$mean, s$sd, s$rel)
  ggplot(d, aes(rep, diff)) +
    annotate("rect", xmin = xr[1] - 0.5, xmax = xr[2] + 0.5,
             ymin = s$mean - s$sd, ymax = s$mean + s$sd,
             fill = taylor_col, alpha = 0.13) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray40", linewidth = 0.35) +
    geom_hline(yintercept = s$mean, color = taylor_col, linewidth = 0.5) +
    geom_point(color = taylor_col, size = 1.5, alpha = 0.9) +
    annotate("text", x = mean(xr), y = Inf, label = lab, hjust = 0.5, vjust = 1.4,
             size = 2.2, color = "gray25", family = "Helvetica") +
    scale_x_continuous(breaks = 1:10) +
    labs(title = tag,
         x = "Replicate",
         y = if (show_y) "ELBO at max:\nTaylor − Monte Carlo" else NULL)
}

pA <- make_traj(25, "A", show_y = TRUE)
pB <- make_traj(50, "B")
pC <- make_diff(25, "C", show_y = TRUE)
pD <- make_diff(50, "D")

fig <- (pA + pB + pC + pD) + plot_layout(ncol = 2)
fig <- fig & theme(
  plot.title = element_text(size = 12, face = "bold", family = "Helvetica",
                            hjust = 0, margin = margin(b = 4)),
  plot.margin = margin(t = 4, r = 6, b = 2, l = 2)
)

save_pdf(fig, file.path(out_dir, "elbo_taylor_vs_mc.pdf"), width = 6.6, height = 5.2)
ggsave(file.path(out_dir, "elbo_taylor_vs_mc.png"),
       plot = fig, width = 6.6, height = 5.2, units = "in", dpi = 220)
cat("wrote elbo_taylor_vs_mc.pdf / .png\n")
