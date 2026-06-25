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

data_dir <- file.path(dirname(script_dir), "jc69-data")
out_dir <- script_dir

include_mf   <- FALSE
include_dist <- FALSE

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

ylim_lnl <- c(-0.023, 0.01)

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

exclude_rows_time <- time$ntaxa %in% c(25, 50, 100)
omit_methods_time <- intersect(
  c("beast-beagle", "mrbayes-beagle",
    "dodonaphy", "geophy", "vaiphy"),
  time_methods
)
for (m in omit_methods_time) {
  time[exclude_rows_time, m] <- ifelse(
    time[exclude_rows_time, m] == 0,
    NA, time[exclude_rows_time, m]
  )
  m_idx <- which(names(time) == m)
  if (length(m_idx) == 1 &&
      m_idx < ncol(time) &&
      identical(names(time)[m_idx + 1], "std")) {
    time[exclude_rows_time, m_idx + 1] <- NA
  }
}

time_long <- melt_mean_sd(time, methods = time_methods,
                          ntaxa_col = "ntaxa", scale_by = NULL)
time_long$method <- factor(time_long$method, levels = time_methods)

valid <- is.finite(time_long$mean) & is.finite(time_long$sd) &
  time_long$mean > 0 & time_long$sd >= 0
cv <- time_long$sd / time_long$mean
sigma_log <- sqrt(log1p(cv^2))
c_mult <- exp(sigma_log)
time_long$ymin_mult <- ifelse(valid, time_long$mean / c_mult, NA_real_)
time_long$ymax_mult <- ifelse(valid, time_long$mean * c_mult, NA_real_)

ptime <- ggplot(time_long,
                aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ymin_mult, ymax = ymax_mult),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  labs(x = "Number of Taxa (n)", y = "Time (s)", fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(time_long$method)],
    breaks = levels(time_long$method),
    labels = label_map(levels(time_long$method))
  ) +
  scale_y_log10() +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

ntaxa_sorted <- sort(unique(time_long$ntaxa))
k <- min(3, length(ntaxa_sorted))
first_ntaxa <- ntaxa_sorted[seq_len(k)]
last_ntaxa <- tail(ntaxa_sorted, k)

ptime_first <- ggplot(
  subset(time_long, ntaxa %in% first_ntaxa),
  aes(x = factor(ntaxa), y = mean, fill = method)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = ymin_mult, ymax = ymax_mult),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) +
  labs(x = "Number of Taxa (n)", y = "Time (s)", fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(time_long$method)],
    breaks = levels(time_long$method),
    labels = label_map(levels(time_long$method))
  ) +
  scale_y_log10() +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

ptime_last <- ggplot(
  subset(
    time_long,
    ntaxa %in% union(last_ntaxa, 10) &
      !(method %in% c("dodonaphy", "geophy", "vaiphy"))
  ),
  aes(x = factor(ntaxa), y = mean, fill = method)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) +
  labs(x = "Number of Taxa (n)", y = "Time (s)", fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(time_long$method)],
    breaks = levels(time_long$method),
    labels = label_map(levels(time_long$method))
  ) +
  guides(fill = "none")

ptime_split <- (ptime_first + ptime_last) +
  plot_layout(ncol = 2, guides = "collect")
ptime_split <- ptime_split & theme(legend.position = "right")

# -------------------- Lnl + Time panels (10/15/20) --------------------
sub_time <- subset(time_long, ntaxa %in% c(10, 15, 20))
data_max <- suppressWarnings(
  max(sub_time$ymax_mult, sub_time$mean, na.rm = TRUE)
)
if (!is.finite(data_max) || data_max <= 0) data_max <- 1
max_exp <- ceiling(log10(data_max))
min_exp <- -1
y_breaks_time <- 10^seq(min_exp, max_exp, by = 1)
half_breaks_time <- 10^seq(min_exp, max_exp - 1, by = 1) * sqrt(10)
labels_time <- function(x) {
  ifelse(x == 0.1, "0.1", as.character(as.integer(x)))
}
y_upper_time <- 10^max_exp

lnl_sub <- subset(lnl_long, ntaxa %in% c(10, 15, 20))
plnl_sub <- ggplot(
  lnl_sub,
  aes(x = factor(ntaxa), y = mean, fill = method)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.4,
    position = position_dodge(width = 0.9)
  ) +
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
  coord_cartesian(ylim = c(-0.023, 0.01))

ptime_sub <- ggplot(
  subset(time_long, ntaxa %in% c(10, 15, 20)),
  aes(x = factor(ntaxa), y = mean)
) +
  geom_hline(
    yintercept = half_breaks_time,
    color = "gray85",
    linewidth = 0.2
  ) +
  geom_linerange(
    aes(ymin = 0.1, ymax = mean, color = method),
    position = position_dodge(width = 0.9),
    linewidth = 2.0
  ) +
  geom_errorbar(
    aes(ymin = ymin_mult, ymax = ymax_mult, group = method),
    width = 0.2,
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  labs(x = "Number of Taxa (n)", y = "Time (s)", color = "Method") +
  scale_color_manual(
    values = method_palette[levels(time_long$method)],
    breaks = levels(time_long$method),
    labels = label_map(levels(time_long$method))
  ) +
  scale_y_log10(
    breaks = y_breaks_time,
    minor_breaks = half_breaks_time,
    expand = c(0, 0),
    labels = labels_time
  ) +
  coord_cartesian(ylim = c(0.1, y_upper_time)) +
  guides(color = "none")

lnl_time_panels <- (plnl_sub + ptime_sub) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    tag_levels = "A", tag_prefix = "", tag_suffix = ""
  )
lnl_time_panels <- lnl_time_panels & theme(
  legend.position = "right",
  plot.tag = element_text(
    size = 12, face = "bold", family = "Helvetica"
  ),
  plot.tag.position = c(0.01, 0.99)
)

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
  guides(fill = "none")

# -------------------- Model-fit deviations --------------------
if (include_mf) {
  mf <- read.table(file.path(data_dir, "mfSummary.txt"),
                   header = TRUE)
  cols <- c("NJ", "vine", "beast", "mrbayes")
  mf[cols] <- (mf[cols] - mf$true) / abs(mf$true)
  mf$`std.2` <- mf$`std.2` / abs(mf$true)
  mf$`std.3` <- mf$`std.3` / abs(mf$true)
  mf$`std.4` <- mf$`std.4` / abs(mf$true)

  mf_long <- rbind(
    data.frame(ntaxa = mf$ntaxa, method = "NJ",
               mean = mf$NJ, sd = mf$`std.1`),
    data.frame(ntaxa = mf$ntaxa, method = "vine",
               mean = mf$vine, sd = mf$`std.2`),
    data.frame(ntaxa = mf$ntaxa, method = "beast",
               mean = mf$beast, sd = mf$`std.3`),
    data.frame(ntaxa = mf$ntaxa, method = "mrbayes",
               mean = mf$mrbayes, sd = mf$`std.4`)
  )
  mf_long$method <- factor(
    mf_long$method,
    levels = c("NJ", "vine", "beast", "mrbayes")
  )

  desired_order <- c("NJ", "vine", "beast",
                     "mrbayes", "dodonaphy", "geophy")
  legend_methods <- desired_order[
    desired_order %in% union(
      levels(lnl_long$method), levels(mf_long$method)
    )
  ]
  fill_vals <- method_palette[legend_methods]
  fill_breaks <- legend_methods
  fill_labels <- label_map(legend_methods)

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
      labels = fill_labels,
      drop = TRUE
    ) +
    scale_y_continuous(
      labels = label_percent(accuracy = 0.1),
      breaks = pretty_breaks(n = 8)
    ) +
    guides(fill = "none") +
    theme(legend.position = "none")

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
}

# -------------------- 95% CIs --------------------
if (include_dist) {
  dist <- read.table(
    file.path(data_dir, "distSummary.txt"),
    header = TRUE
  )

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
      ntaxa = dist$ntaxa, method = "vine + flows", value = dist$vineflows,
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

  entropy <- read.table(
    file.path(data_dir, "entropySummary.txt"),
    header = TRUE
  )

  has_entropy_sd <- all(
    c("vinedev", "beastdev", "vineflowsdev") %in% names(entropy)
  )
  zero_sd <- rep(0, nrow(entropy))

  entropy_long <- rbind(
    data.frame(
      ntaxa = entropy$ntaxa, method = "vine", value = entropy$vine,
      sd = if (has_entropy_sd) entropy$vinedev else zero_sd
    ),
    data.frame(
      ntaxa = entropy$ntaxa, method = "vine + flows", value = entropy$vineflows,
      sd = if (has_entropy_sd) entropy$vineflowsdev else zero_sd
    ),
    data.frame(
      ntaxa = entropy$ntaxa, method = "beast", value = entropy$beast,
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

  vine_col <- unname(method_palette["vine"])
  beast_col <- unname(method_palette["beast"])
  vine_flows_col <- "#E6550D"
  beast_label <- unname(label_map("beast"))

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
}

# -------------------- Save --------------------
if (include_dist) {
  save_pdf(
    pdist_flows,
    file.path(out_dir, "jc69_dist_flows_bars.pdf"),
    width = 3, height = 3
  )
}

if (include_mf) {
  save_pdf(
    pcombo,
    file.path(out_dir, "jc69_lnl_mf_panels.pdf"),
    width = 6, height = 3
  )
}

save_pdf(
  lnl_time_panels,
  file.path(out_dir, "figure2_A_B.pdf"),
  width = 6, height = 3
)
