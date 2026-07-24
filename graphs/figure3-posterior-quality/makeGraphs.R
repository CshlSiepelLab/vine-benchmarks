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

# ================================================================
# Embedding dimension (25 taxa) — Euclidean ("E") version
# ================================================================
# dimSummary25.txt columns (by name):
# d, dlnlE, dlnlE_sd, timeE_med, timeE_q1, timeE_q3,
# dlnlH, dlnlH_sd, timeH_med, timeH_q1, timeH_q3
dims25 <- read.table(file.path(data_dir, "dimSummary25.txt"),
                     header = TRUE, check.names = FALSE)

# Extract numeric embedding dimension from strings like "archive.D2"
dims25$d <- suppressWarnings(
  as.numeric(gsub(".*\\.D", "", as.character(dims25[[1]])))
)

delta_lnl <- dims25$dlnlE
sd_lnl    <- dims25$dlnlE_sd
time_med  <- dims25$timeE_med
time_q1   <- dims25$timeE_q1
time_q3   <- dims25$timeE_q3

# Shared scale factor for dual axis
num_max <- max(abs(delta_lnl), na.rm = TRUE)
den_max <- max(time_q3, na.rm = TRUE)
scale_factor <- if (is.finite(num_max) && is.finite(den_max) &&
                    den_max > 0) num_max / den_max else 1

# Long data; time uses q1/q3 as asymmetric error bounds (no time sd
# is reported), lnl uses mean +/- sd
plotdata <- rbind(
  data.frame(d = dims25$d, metric = "delta_lnl",
            mean_scaled = delta_lnl,
            ymin_scaled = delta_lnl - sd_lnl,
            ymax_scaled = delta_lnl + sd_lnl,
            stringsAsFactors = FALSE),
  data.frame(d = dims25$d, metric = "time",
            mean_scaled = time_med * scale_factor,
            ymin_scaled = time_q1 * scale_factor,
            ymax_scaled = time_q3 * scale_factor,
            stringsAsFactors = FALSE)
)
plotdata$metric <- factor(plotdata$metric, levels = c("delta_lnl", "time"))

# Axis breaks (only <=0 for Δlnl; only >=0 for time)
ll_min <- suppressWarnings(min(delta_lnl - sd_lnl, 0, na.rm = TRUE))
ll_breaks <- pretty(c(ll_min, 0), n = 6)
ll_breaks <- ll_breaks[ll_breaks <= 0]

t_max <- suppressWarnings(max(time_q3, 0, na.rm = TRUE))
t_breaks <- pretty(c(0, t_max), n = 6)
t_breaks <- t_breaks[t_breaks >= 0]

combined_breaks <- sort(unique(c(ll_breaks, t_breaks * scale_factor)))

dodge <- position_dodge(width = 0.6)
geom_palette <- c("delta_lnl" = "#F28E2B", "time" = "#59A14F")
orange <- unname(geom_palette["delta_lnl"])
green  <- unname(geom_palette["time"])

pdim25 <- ggplot(plotdata, aes(x = factor(d), y = mean_scaled, fill = metric)) +
  geom_col(width = 0.45, position = dodge) +
  geom_errorbar(aes(ymin = ymin_scaled, ymax = ymax_scaled),
                width = 0.15, position = dodge) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = geom_palette, guide = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
  scale_y_continuous(
    name = expression(Delta~lnl),
    breaks = combined_breaks,
    labels = function(y) ifelse(y <= 0, as.character(y), ""),
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "time (s)",
                        breaks = t_breaks)
  ) +
  labs(x = "Embedding Dimension (d)", title = NULL) +
  theme(
    axis.title.y      = element_text(color = orange),
    axis.text.y       = element_text(color = orange),
    axis.ticks.y      = element_line(color = orange),

    axis.title.y.right = element_text(color = green),
    axis.text.y.right  = element_text(color = green),
    axis.ticks.y.right = element_line(color = green)
  )

save_pdf(pdim25, file.path(out_dir, "hky300_dim25_bars.pdf"),
        width = 3, height = 3)

# ================================================================
# Embedding dimension (25 taxa) — Hyperbolic ("H") version
# ================================================================
delta_lnl <- dims25$dlnlH
sd_lnl    <- dims25$dlnlH_sd
time_med  <- dims25$timeH_med
time_q1   <- dims25$timeH_q1
time_q3   <- dims25$timeH_q3

num_max <- max(abs(delta_lnl), na.rm = TRUE)
den_max <- max(time_q3, na.rm = TRUE)
scale_factor <- if (is.finite(num_max) && is.finite(den_max) &&
                    den_max > 0) num_max / den_max else 1

plotdata <- rbind(
  data.frame(d = dims25$d, metric = "delta_lnl",
            mean_scaled = delta_lnl,
            ymin_scaled = delta_lnl - sd_lnl,
            ymax_scaled = delta_lnl + sd_lnl,
            stringsAsFactors = FALSE),
  data.frame(d = dims25$d, metric = "time",
            mean_scaled = time_med * scale_factor,
            ymin_scaled = time_q1 * scale_factor,
            ymax_scaled = time_q3 * scale_factor,
            stringsAsFactors = FALSE)
)
plotdata$metric <- factor(plotdata$metric, levels = c("delta_lnl", "time"))

# Axis limits/breaks include the error-bar extremes so bars are
# never clipped
ll_min_err <- suppressWarnings(min(delta_lnl - sd_lnl, 0, na.rm = TRUE))
ll_breaks <- pretty(c(ll_min_err, 0), n = 6)
if (length(ll_breaks) > 1 && min(ll_breaks) > ll_min_err) {
  step <- ll_breaks[2] - ll_breaks[1]
  ll_breaks <- c(min(ll_breaks) - step, ll_breaks)
}
ymin <- min(ll_breaks, na.rm = TRUE)

t_max_err <- suppressWarnings(max(time_q3, 0, na.rm = TRUE))
t_breaks <- pretty(c(0, t_max_err), n = 6)
if (length(t_breaks) > 1 && max(t_breaks) < t_max_err) {
  step <- t_breaks[2] - t_breaks[1]
  t_breaks <- c(t_breaks, max(t_breaks) + step)
}
ymax <- max(t_breaks, na.rm = TRUE) * scale_factor

combined_breaks <- sort(unique(c(ll_breaks, t_breaks * scale_factor)))

geom_palette <- c("delta_lnl" = "#F28E2B", "time" = "#59A14F")
orange <- unname(geom_palette["delta_lnl"])
green  <- unname(geom_palette["time"])

pdimH25 <- ggplot(plotdata, aes(x = factor(d), y = mean_scaled, fill = metric)) +
  geom_col(width = 0.45, position = dodge) +
  geom_errorbar(aes(ymin = ymin_scaled, ymax = ymax_scaled),
                width = 0.15, position = dodge) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = geom_palette, guide = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
  scale_y_continuous(
    name = expression(Delta~lnl),
    limits = c(ymin, ymax),
    breaks = combined_breaks,
    labels = function(y) ifelse(y <= 0, as.character(y), ""),
    expand = expansion(mult = c(0, 0)),
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "time (s)",
                        breaks = t_breaks)
  ) +
  labs(x = "Embedding Dimension (d)", title = NULL) +
  theme(
    axis.title.y      = element_text(color = orange),
    axis.text.y       = element_text(color = orange),
    axis.ticks.y      = element_line(color = orange),

    axis.title.y.right = element_text(color = green),
    axis.text.y.right  = element_text(color = green),
    axis.ticks.y.right = element_line(color = green)
  )

save_pdf(pdimH25, file.path(out_dir, "hky300_dimH25_bars.pdf"),
        width = 3, height = 3)
