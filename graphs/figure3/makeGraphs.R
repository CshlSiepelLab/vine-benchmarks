#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))

# Resolve script directory so we can build absolute paths
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}

data_dir <- file.path(dirname(script_dir), "hky300-data")
out_dir <- script_dir

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
# 95% CI Inclusion
# ================================================================
dist <- read.table(file.path(data_dir, "distSummary.txt"),
                   header = TRUE, fill = TRUE)

has_dist_sd <- all(
  c("vinedev", "beastdev", "vineflowsdev") %in% names(dist)
)
zero_sd <- rep(0, nrow(dist))

dist3_long <- rbind(
  data.frame(
    ntaxa = dist$ntaxa, method = "vine", value = dist$vine,
    sd = if (has_dist_sd) dist$vinedev else zero_sd
  ),
  data.frame(
    ntaxa = dist$ntaxa, method = "vine + flows",
    value = dist$vineflows,
    sd = if (has_dist_sd) dist$vineflowsdev else zero_sd
  ),
  data.frame(
    ntaxa = dist$ntaxa, method = "beast", value = dist$beast,
    sd = if (has_dist_sd) dist$beastdev else zero_sd
  )
)

dist3_long$method <- factor(
  dist3_long$method, levels = c("vine", "vine + flows", "beast")
)
dist3_long$ymin <- pmax(0, dist3_long$value - dist3_long$sd)
dist3_long$ymax <- pmin(1, dist3_long$value + dist3_long$sd)
y_top <- suppressWarnings(max(dist3_long$ymax, na.rm = TRUE))
if (!is.finite(y_top) || y_top <= 0) y_top <- 1
y_top <- y_top + 0.03

vine_col <- unname(method_palette["vine"])
beast_col <- unname(method_palette["beast"])
vine_flows_col <- "#E6550D"
beast_label <- unname(label_map("beast"))

pdist_flows <- ggplot(
  dist3_long, aes(x = factor(ntaxa), y = value, fill = method)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) +
  labs(
    x = "Number of Taxa (n)",
    y = "95% CI Inclusion",
    fill = "Method"
  ) +
  scale_y_continuous(
    limits = c(0, y_top),
    breaks = pretty_breaks(n = 6),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "vine" = vine_col,
      "vine + flows" = vine_flows_col,
      "beast" = beast_col
    ),
    breaks = c("vine", "vine + flows", "beast"),
    labels = c("Vine", "Vine + flows", beast_label)
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

save_pdf(pdist_flows,
         file.path(out_dir, "hky300_dist_flows_bars.pdf"),
         width = 3, height = 3)

# ================================================================
# Entropy
# ================================================================
entropy <- read.table(file.path(data_dir, "entropySummary.txt"),
                      header = TRUE, fill = TRUE)

has_entropy_sd <- all(
  c("vinedev", "beastdev", "vineflowsdev") %in% names(entropy)
)
zero_sd <- rep(0, nrow(entropy))

entropy_long <- rbind(
  data.frame(
    ntaxa = entropy$ntaxa, method = "vine",
    value = entropy$vine,
    sd = if (has_entropy_sd) entropy$vinedev else zero_sd
  ),
  data.frame(
    ntaxa = entropy$ntaxa, method = "vine + flows",
    value = entropy$vineflows,
    sd = if (has_entropy_sd) entropy$vineflowsdev else zero_sd
  ),
  data.frame(
    ntaxa = entropy$ntaxa, method = "beast",
    value = entropy$beast,
    sd = if (has_entropy_sd) entropy$beastdev else zero_sd
  )
)

entropy_long$method <- factor(
  entropy_long$method, levels = c("vine", "vine + flows", "beast")
)
entropy_long$ymin <- pmax(0, entropy_long$value - entropy_long$sd)
entropy_long$ymax <- entropy_long$value + entropy_long$sd
y_top <- suppressWarnings(max(entropy_long$ymax, na.rm = TRUE))
if (!is.finite(y_top) || y_top <= 0) y_top <- 1
y_top <- y_top + 0.03

pentropy <- ggplot(
  entropy_long, aes(x = factor(ntaxa), y = value, fill = method)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) +
  labs(
    x = "Number of Taxa (n)",
    y = "Entropy",
    fill = "Method"
  ) +
  scale_y_continuous(
    limits = c(0, y_top),
    breaks = pretty_breaks(n = 6),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "vine" = vine_col,
      "vine + flows" = vine_flows_col,
      "beast" = beast_col
    ),
    breaks = c("vine", "vine + flows", "beast"),
    labels = c("Vine", "Vine + flows", beast_label)
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

save_pdf(pentropy,
         file.path(out_dir, "hky300_entropy_bars.pdf"),
         width = 3, height = 3)
