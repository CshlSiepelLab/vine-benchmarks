#!/usr/bin/env Rscript
# Supplement figure: centering of pairwise-distance estimates.
# Three calibration panels (posterior-mean distance vs. ground truth) for 25/50/100 taxa.
# Style, palette, fonts, panel sizes, and A/B/C callouts match figure2-hky300/makeGraphs.R.

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
  theme_minimal(base_size = 8) +
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

# ------------ Colors & Labels (identical to figure 2) ------------
method_palette <- c(
  NJ        = "#4E79A7",
  vine      = "#F28E2B",
  beast     = "#59A14F",
  "beast-beagle" = "#8BC184",
  mrbayes   = "#E15759",
  "mrbayes-beagle" = "#E98A8C"
)
method_labels <- c(
  NJ        = "NJ",
  vine      = "Vine",
  beast     = "BEAST2",
  "beast-beagle" = "BEAST2 + BEAGLE",
  mrbayes   = "MrBayes",
  "mrbayes-beagle" = "MrBayes + BEAGLE"
)
label_map <- function(keys) {
  out <- method_labels[keys]; out[is.na(out)] <- keys[is.na(out)]; out
}

# ------------ Data ------------
dat <- read.csv(file.path(data_dir, "calibration_data.csv"),
                stringsAsFactors = FALSE)
methods <- c("vine", "beast", "mrbayes")            # draw/legend order
dat$method <- factor(dat$method, levels = methods)

sizes <- c(25, 50, 100)
AXMAX <- 2.0                                        # data end at ~2.0 subs/site
NBIN  <- 12
POINTS_PER_GROUP <- 1500

set.seed(1)

# subsample points per (ntaxa, method) so the vector PDF stays light
pts <- dat %>%
  group_by(ntaxa, method) %>%
  slice_sample(n = POINTS_PER_GROUP) %>%
  ungroup()

# binned calibration curve: mean estimate per true-distance bin
brks <- seq(0, AXMAX, length.out = NBIN + 1)
cal <- dat %>%
  mutate(bin = cut(true, breaks = brks, include.lowest = TRUE)) %>%
  group_by(ntaxa, method, bin) %>%
  summarise(true = mean(true), est = mean(est), n = dplyr::n(), .groups = "drop") %>%
  filter(n > 3)

make_panel <- function(nt, tag, show_y = FALSE) {
  ggplot() +
    geom_point(data = subset(pts, ntaxa == nt),
               aes(true, est, color = method),
               size = 0.35, alpha = 0.12, stroke = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray40", linewidth = 0.3) +
    geom_line(data = subset(cal, ntaxa == nt),
              aes(true, est, color = method), linewidth = 0.7) +
    geom_point(data = subset(cal, ntaxa == nt),
               aes(true, est, color = method), size = 0.8) +
    scale_color_manual(values = method_palette[methods],
                       breaks = methods, labels = label_map(methods)) +
    scale_x_continuous(limits = c(0, AXMAX), breaks = pretty_breaks(n = 5)) +
    scale_y_continuous(limits = c(0, AXMAX), breaks = pretty_breaks(n = 5)) +
    coord_fixed() +
    labs(title = tag,
         x = "True distance",
         y = if (show_y) "Estimated distance" else NULL,
         color = "Method") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.4,
                                                    linewidth = 0.7)))
}

pA <- make_panel(sizes[1], "A", show_y = TRUE)
pB <- make_panel(sizes[2], "B")
pC <- make_panel(sizes[3], "C")

# A/B/C rendered as real panel titles (own strip above the panel) so they cannot
# overlap the axis; left-aligned and bold to match figure 2's callouts.
fig <- (pA + pB + pC) + plot_layout(ncol = 3, guides = "collect")
fig <- fig & theme(
  legend.position = "right",
  plot.title = element_text(size = 12, face = "bold", family = "Helvetica",
                            hjust = 0, margin = margin(b = 5)),
  plot.margin = margin(t = 4, r = 4, b = 2, l = 2)
)

save_pdf(fig, file.path(out_dir, "dist_centering.pdf"), width = 9, height = 3)
ggsave(file.path(out_dir, "dist_centering.png"),
       plot = fig, width = 9, height = 3, units = "in", dpi = 220)
cat("wrote dist_centering.pdf / .png\n")
