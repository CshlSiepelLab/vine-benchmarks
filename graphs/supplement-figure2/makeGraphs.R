#!/usr/bin/env Rscript

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

# ------------ Colors & Labels ------------
method_palette <- c(
  NJ        = "#4E79A7",
  vine      = "#F28E2B",
  beast     = "#59A14F",
  "beast-beagle" = "#8BC184",
  mrbayes   = "#E15759",
  "mrbayes-beagle" = "#E98A8C",
  dodonaphy = "#B07AA1",
  geophy    = "#EDC948",
  vaiphy    = "#76B7B2"
)

method_labels <- c(
  NJ        = "NJ",
  vine      = "Vine",
  beast     = "BEAST2",
  "beast-beagle" = "BEAST2 + BEAGLE",
  mrbayes   = "MrBayes",
  "mrbayes-beagle" = "MrBayes + BEAGLE",
  dodonaphy = "Dodonaphy",
  geophy    = "GeoPhy",
  vaiphy    = "VaiPhy"
)

label_map <- function(keys) {
  out <- method_labels[keys]
  out[is.na(out)] <- keys[is.na(out)]
  out
}

# ================================================================
# Robinson-Foulds distances
# ================================================================
rf <- read.table(file.path(data_dir, "rfSummary.txt"), header = TRUE)

## Identify the std columns in order (std, std.1, std.2, ...)
std_cols <- grep("^std", names(rf))
if (length(std_cols) < 4) {
  stop("Expected at least 4 'std' columns (for NJ, vine, beast, mrbayes). Found: ",
      length(std_cols))
}

# Long format (absolute RF values)
rf_long <- rbind(
  data.frame(ntaxa  = rf$ntaxa, method = "vine", mean   = rf$vine, sd     = rf[[std_cols[2]]]),
  data.frame(ntaxa  = rf$ntaxa, method = "beast", mean   = rf$beast, sd     = rf[[std_cols[3]]]),
  data.frame(ntaxa  = rf$ntaxa, method = "mrbayes", mean   = rf$mrbayes, sd     = rf[[std_cols[4]]])
)

rf_long$method <- factor(
  rf_long$method,
  levels = c("vine","beast","mrbayes")
)

# Error bars, truncated at zero on the lower end
rf_long$ymin <- pmax(rf_long$mean - rf_long$sd, 0)
rf_long$ymax <- rf_long$mean + rf_long$sd

# Absolute y-limits: min = 0, max = max(mean + sd)
y_max <- max(rf_long$ymax, na.rm = TRUE)
ylim  <- c(0, y_max)

# Plot (absolute RF distance; y starts at 0)
prf <- ggplot(rf_long, aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  labs(
    x = "Number of Taxa",
    y = "Normalized RF Distance",
    fill = "Method"
  ) +
  scale_fill_manual(
    values = method_palette[levels(rf_long$method)],
    breaks = levels(rf_long$method),
    labels = label_map(levels(rf_long$method))
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6))) +
  coord_cartesian(ylim = ylim)

# ================================================================
# Save
# ================================================================
save_pdf(prf,
         file.path(script_dir, "hky300_rf_bars.pdf"),
         width = 3, height = 3)
