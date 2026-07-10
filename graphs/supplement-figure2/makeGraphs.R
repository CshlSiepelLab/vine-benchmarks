#!/usr/bin/env Rscript
# Supplement figure S2: accuracy vs. number of taxa for the HKY/300-site
# simulations.  Panel A: normalized Robinson-Foulds distance (topology).
# Panel B: normalized branch-score distance (branch lengths), using the point
# (posterior-mean-tree) BSD, which measures branch-length accuracy independent
# of posterior dispersion.  Both are means +/- SD over 10 replicate datasets.

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))

# Resolve script directory so we can build absolute paths
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}

graphs_dir <- dirname(script_dir)
data_dir <- file.path(graphs_dir, "hky300-data")

# Helper: save a plot as a compact PDF, using cairo if available
save_pdf <- function(plot, filename, width = 3, height = 3) {
  if (capabilities("cairo")) {
    ggsave(filename, plot = plot, width = width, height = height,
           units = "in", device = cairo_pdf)
  } else {
    pdf(filename, width = width, height = height, family = "Helvetica")
    print(plot)
    dev.off()
  }
}

# ------------ Theme ------------
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

# ------------ Colors & Labels ------------
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
  out <- method_labels[keys]
  out[is.na(out)] <- keys[is.na(out)]
  out
}

methods <- c("vine", "beast", "beast-beagle", "mrbayes", "mrbayes-beagle")

# ------------ Load a *Summary.txt table into long format ------------
# Columns: ntaxa NJ std vine std beast std beast-beagle std mrbayes std
#          mrbayes-beagle std   (NJ present but not plotted, matching prior S2)
load_summary <- function(path) {
  d <- read.table(path, header = TRUE)
  std_cols <- grep("^std", names(d))
  if (length(std_cols) < 6)
    stop("Expected >= 6 'std' columns in ", path, "; found ", length(std_cols))
  long <- rbind(
    data.frame(ntaxa = d$ntaxa, method = "vine",           mean = d$vine,           sd = d[[std_cols[2]]]),
    data.frame(ntaxa = d$ntaxa, method = "beast",          mean = d$beast,          sd = d[[std_cols[3]]]),
    data.frame(ntaxa = d$ntaxa, method = "beast-beagle",   mean = d$beast.beagle,   sd = d[[std_cols[4]]]),
    data.frame(ntaxa = d$ntaxa, method = "mrbayes",        mean = d$mrbayes,        sd = d[[std_cols[5]]]),
    data.frame(ntaxa = d$ntaxa, method = "mrbayes-beagle", mean = d$mrbayes.beagle, sd = d[[std_cols[6]]])
  )
  long$method <- factor(long$method, levels = methods)
  long$ymin <- pmax(long$mean - long$sd, 0)
  long$ymax <- long$mean + long$sd
  long
}

make_bar <- function(long, ylab, tag) {
  ggplot(long, aes(x = factor(ntaxa), y = mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2,
                  position = position_dodge(width = 0.9)) +
    labs(title = tag, x = "Number of Taxa", y = ylab, fill = "Method") +
    scale_fill_manual(
      values = method_palette[methods],
      breaks = methods,
      labels = label_map(methods)
    ) +
    scale_y_continuous(breaks = pretty_breaks(n = 8), limits = c(0, NA)) +
    guides(fill = guide_legend(override.aes = list(width = 0.6)))
}

rf_long  <- load_summary(file.path(data_dir, "rfSummary.txt"))
bsd_long <- load_summary(file.path(data_dir, "bsdSummary.txt"))

prf  <- make_bar(rf_long,  "Normalized RF distance", "A")
pbsd <- make_bar(bsd_long, "Normalized branch-score distance", "B")

fig <- (prf + pbsd) + plot_layout(ncol = 2, guides = "collect")
fig <- fig & theme(
  legend.position = "right",
  plot.title = element_text(size = 12, face = "bold", family = "Helvetica",
                            hjust = 0, margin = margin(b = 5))
)

save_pdf(fig, file.path(script_dir, "hky300_rf_bsd_bars.pdf"), width = 9, height = 3)
ggsave(file.path(script_dir, "hky300_rf_bsd_bars.png"),
       plot = fig, width = 9, height = 3, units = "in", dpi = 220)
cat("wrote hky300_rf_bsd_bars.pdf / .png\n")
