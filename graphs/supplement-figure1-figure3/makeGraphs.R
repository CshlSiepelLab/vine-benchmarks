#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))

# Resolve script directory so we can build absolute paths
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}

graphs_dir <- dirname(script_dir)
data_dir <- file.path(graphs_dir, "hky-10k-graphs")

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

# ------------ Utilities ------------
melt_mean_sd <- function(df, methods, ntaxa_col = "ntaxa",
                         scale_by = NULL) {
  stopifnot(ntaxa_col %in% names(df))
  out <- list()
  for (m in methods) {
    if (!(m %in% names(df))) next
    m_idx <- which(names(df) == m)
    sd_idx <- if (m_idx < ncol(df) &&
                 identical(names(df)[m_idx + 1], "std")) {
      m_idx + 1
    } else {
      NA_integer_
    }
    mean_vals <- df[[m]]
    sd_vals <- if (!is.na(sd_idx)) df[[sd_idx]]
               else rep(NA_real_, length(mean_vals))
    if (!is.null(scale_by)) {
      denom <- abs(df[[scale_by]])
      mean_vals <- mean_vals / denom
      sd_vals <- sd_vals / denom
    }
    out[[m]] <- data.frame(
      ntaxa  = df[[ntaxa_col]],
      method = m,
      mean   = mean_vals,
      sd     = sd_vals,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# ================================================================
# 1) Log-likelihood deviations
# ================================================================
lnl <- read.table(file.path(data_dir, "lnlSummary.txt"),
                  header = TRUE, check.names = FALSE)

ll_methods <- intersect(c("NJ", "vine", "beast", "beast-beagle",
                           "mrbayes", "mrbayes-beagle",
                           "dodonaphy", "geophy", "vaiphy"),
                         names(lnl))

lnl_long <- melt_mean_sd(lnl, methods = ll_methods,
                         ntaxa_col = "ntaxa", scale_by = "ave")
lnl_long$method <- factor(lnl_long$method, levels = ll_methods)

lower <- with(lnl_long, mean - sd)
upper <- with(lnl_long, mean + sd)
L <- max(abs(c(lower, upper)), na.rm = TRUE)
if (!is.finite(L)) L <- 1
pad <- 0.02 * L
ylim_lnl <- c(-L - pad, L + pad)

plnl <- ggplot(lnl_long,
               aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.4,
                position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Number of Taxa (n)", y = "Δ Lnl", fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(lnl_long$method)],
    breaks = levels(lnl_long$method),
    labels = label_map(levels(lnl_long$method))
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 0.1),
    breaks = pretty_breaks(n = 8)
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6))) +
  coord_cartesian(ylim = ylim_lnl)

# ================================================================
# 2) Running times
# ================================================================
time <- read.table(file.path(data_dir, "timeSummary.txt"),
                   header = TRUE, check.names = FALSE)

time_methods <- intersect(c("vine", "beast", "beast-beagle",
                             "mrbayes", "mrbayes-beagle",
                             "dodonaphy", "geophy", "vaiphy"),
                           names(time))

time_long <- melt_mean_sd(time, methods = time_methods,
                          ntaxa_col = "ntaxa", scale_by = NULL)
time_long$method <- factor(time_long$method, levels = time_methods)

ptime <- ggplot(time_long,
                aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  labs(x = "Number of Taxa (n)", y = "Time (s)", fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(time_long$method)],
    breaks = levels(time_long$method),
    labels = label_map(levels(time_long$method))
  ) +
  scale_y_log10(
    limits = c(1, NA),
    breaks = scales::breaks_log(n = 6),
    labels = scales::label_number()
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

# ================================================================
# 3) Speedup
# ================================================================
spd <- read.table(file.path(data_dir, "speedSummary.txt"),
                  header = TRUE, check.names = FALSE)

speed_methods <- intersect(c("beast", "beast-beagle",
                             "mrbayes", "mrbayes-beagle"), names(spd))
speed_long <- melt_mean_sd(spd, methods = speed_methods,
                           ntaxa_col = "ntaxa", scale_by = NULL)
speed_long$method <- factor(speed_long$method, levels = speed_methods)

pspeed <- ggplot(speed_long,
                 aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  labs(x = "Number of Taxa (n)", y = "Speed Increase (x)",
       fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(speed_long$method)],
    breaks = levels(speed_long$method),
    labels = label_map(levels(speed_long$method))
  ) +
  guides(fill = "none") +
  theme(legend.position = "none")

# ================================================================
# 4) Model-fit deviations
# ================================================================
mf <- read.table(file.path(data_dir, "mfSummary.txt"),
                 header = TRUE, check.names = FALSE)
mf_methods <- intersect(
  c("NJ", "vine", "beast", "beast-beagle",
    "mrbayes", "mrbayes-beagle"),
  names(mf)
)
mf_long <- melt_mean_sd(mf, methods = mf_methods,
                        ntaxa_col = "ntaxa", scale_by = NULL)
true_by_ntaxa <- setNames(mf$true, mf$ntaxa)
mf_true <- true_by_ntaxa[as.character(mf_long$ntaxa)]
mf_long$mean <- (mf_long$mean - mf_true) / abs(mf_true)
mf_long$sd <- mf_long$sd / abs(mf_true)
mf_long$method <- factor(mf_long$method, levels = mf_methods)

lower <- with(mf_long, mean - sd)
upper <- with(mf_long, mean + sd)
L <- max(abs(c(lower, upper)))
pad <- 0.02 * L
ylim_mf <- c(-L - pad, L + pad)
fill_vals <- method_palette[levels(mf_long$method)]
fill_breaks <- levels(mf_long$method)
fill_labels <- label_map(levels(mf_long$method))

pmf <- ggplot(
  mf_long,
  aes(x = factor(ntaxa), y = mean, fill = method)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Number of Taxa (n)",
    y = "Held out Δ Lnl",
    fill = "Method"
  ) +
  scale_fill_manual(
    values = fill_vals,
    breaks = fill_breaks,
    labels = fill_labels
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 0.1),
    breaks = pretty_breaks(n = 8)
  ) +
  guides(fill = "none") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = ylim_mf)

pcombo <- (plnl + pmf) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "", tag_suffix = "")
pcombo <- pcombo & theme(
  legend.position = "right",
  plot.tag = element_text(
    size = 12, face = "bold", family = "Helvetica"
  ),
  plot.tag.position = c(0.01, 0.99)
)

# ================================================================
# Save
# ================================================================
save_pdf(
  pcombo,
  file.path(script_dir, "supplement-figure1.pdf"),
  width = 6, height = 3
)

ptime_speed <- (ptime + pspeed) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "", tag_suffix = "")
ptime_speed <- ptime_speed & theme(
  legend.position = "right",
  plot.tag = element_text(size = 12, face = "bold",
                          family = "Helvetica"),
  plot.tag.position = c(0.01, 0.99)
)
save_pdf(ptime_speed,
         file.path(script_dir, "supplement-figure3.pdf"),
         width = 6, height = 3)
