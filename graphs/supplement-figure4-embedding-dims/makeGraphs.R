#!/usr/bin/env Rscript
# Supplement figure S4: effect of embedding dimension d on VINE accuracy
# (Delta lnl vs. BEAST) and runtime, for 25 taxa (A) and 50 taxa (B).
# Each panel overlays the Euclidean and hyperbolic geometries (previously
# shown as two separate figures / four panels).  Left axis (orange): Delta lnl
# = VINE - BEAST held-out log-likelihood; bars point down when VINE is worse.
# Right axis (green): VINE runtime.  Darker bars = Euclidean, lighter =
# hyperbolic.  Error bars are one SD over ten replicates.

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()

data_dir <- file.path(dirname(script_dir), "hky300-data")
out_dir <- script_dir

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

panel_tag_theme <- theme(
  plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
  plot.tag.position = c(0.01, 0.99)
)

# Palette: metric by hue (orange = Delta lnl, green = time),
# geometry by shade (dark = Euclidean, light = hyperbolic).
pal <- c(d_euc = "#F28E2B", d_hyp = "#F9C79A",
         t_euc = "#59A14F", t_hyp = "#A7D3A1")
lab <- c(d_euc = "Δ lnl, Euclidean",  d_hyp = "Δ lnl, hyperbolic",
         t_euc = "time, Euclidean",        t_hyp = "time, hyperbolic")
orange <- "#F28E2B"; green <- "#59A14F"

# ---- combined Euclidean + hyperbolic dual-axis panel ----
# Delta lnl (left, orange): mean +/- SD over replicates.
# time (right, green): median with IQR (Q1-Q3) whiskers -- robust to the
# heavy-tailed hyperbolic convergence time (a minority of replicates
# initialize with very large KLD and take many more iterations).
# dims_file columns:
#   d  dlnlE dlnlE_sd  timeE_med timeE_q1 timeE_q3  dlnlH dlnlH_sd  timeH_med timeH_q1 timeH_q3
make_dim_panel2 <- function(dims_file) {
  dims <- read.table(file.path(data_dir, dims_file),
                     header = TRUE, check.names = FALSE)
  dval <- suppressWarnings(as.numeric(gsub(".*\\.D", "", as.character(dims[[1]]))))

  d_euc <- dims$dlnlE; d_euc_sd <- dims$dlnlE_sd
  d_hyp <- dims$dlnlH; d_hyp_sd <- dims$dlnlH_sd
  te_m <- dims$timeE_med; te_lo <- dims$timeE_q1; te_hi <- dims$timeE_q3
  th_m <- dims$timeH_med; th_lo <- dims$timeH_q1; th_hi <- dims$timeH_q3

  have_hyp <- any(is.finite(d_hyp) & (d_hyp != 0 | th_m != 0))

  num_max <- max(abs(c(d_euc, if (have_hyp) d_hyp)) +
                 c(d_euc_sd, if (have_hyp) d_hyp_sd), na.rm = TRUE)
  den_max <- max(c(te_hi, if (have_hyp) th_hi), na.rm = TRUE)
  scale_factor <- if (is.finite(num_max) && is.finite(den_max) &&
                      den_max > 0 && num_max > 0) num_max / den_max else 1

  mk <- function(series, height, lo, hi)
    data.frame(d = dval, series = series, y = height, lo = lo, hi = hi)
  rows <- list(
    mk("d_euc", d_euc, d_euc - d_euc_sd, d_euc + d_euc_sd),
    mk("t_euc", te_m * scale_factor, te_lo * scale_factor, te_hi * scale_factor)
  )
  if (have_hyp) rows <- c(rows, list(
    mk("d_hyp", d_hyp, d_hyp - d_hyp_sd, d_hyp + d_hyp_sd),
    mk("t_hyp", th_m * scale_factor, th_lo * scale_factor, th_hi * scale_factor)
  ))
  plotdata <- bind_rows(rows)
  plotdata$series <- factor(plotdata$series, levels = c("d_euc","d_hyp","t_euc","t_hyp"))

  ll_lo <- min(c(d_euc - d_euc_sd, if (have_hyp) d_hyp - d_hyp_sd, 0), na.rm = TRUE)
  ll_hi <- max(c(d_euc + d_euc_sd, if (have_hyp) d_hyp + d_hyp_sd, 0), na.rm = TRUE)
  ll_breaks <- pretty(c(ll_lo, ll_hi), n = 6)
  t_hi_all <- max(c(te_hi, if (have_hyp) th_hi, 0), na.rm = TRUE)
  t_breaks <- pretty(c(0, t_hi_all), n = 5)
  axis_min <- min(ll_breaks, min(t_breaks) * scale_factor, na.rm = TRUE)
  axis_max <- max(ll_breaks, max(t_breaks) * scale_factor, na.rm = TRUE)

  dodge <- position_dodge(width = 0.8)
  ggplot(plotdata, aes(factor(d), y, fill = series)) +
    geom_col(width = 0.7, position = dodge) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  width = 0.25, linewidth = 0.3, position = dodge) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_fill_manual(values = pal, labels = lab, name = NULL,
                      breaks = names(pal), drop = FALSE) +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(
      name = expression(Delta ~ lnl),
      limits = c(axis_min, axis_max), breaks = ll_breaks,
      labels = function(y) ifelse(y <= 0, as.character(y), ""),
      sec.axis = sec_axis(~ . / scale_factor, name = "time (s)", breaks = t_breaks)
    ) +
    labs(x = "Embedding Dimension (d)") +
    theme(
      axis.title.y       = element_text(color = orange),
      axis.text.y        = element_text(color = orange),
      axis.ticks.y       = element_line(color = orange),
      axis.title.y.right = element_text(color = green),
      axis.text.y.right  = element_text(color = green),
      axis.ticks.y.right = element_line(color = green)
    )
}

p25 <- make_dim_panel2("dimSummary25.txt")
p50 <- make_dim_panel2("dimSummary50.txt")

fig <- (p25 + p50) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")
fig <- fig & panel_tag_theme & theme(legend.position = "right")

save_pdf(fig, file.path(out_dir, "dims_combined_panels.pdf"), width = 8, height = 3)
ggsave(file.path(out_dir, "dims_combined_panels.png"),
       plot = fig, width = 8, height = 3, units = "in", dpi = 220)
cat("wrote dims_combined_panels.pdf / .png\n")
