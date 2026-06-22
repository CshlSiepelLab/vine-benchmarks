#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))

# ------------ CLI ------------
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) < 1) {
  stop(
    "Usage: Rscript makeGraphs3.R [hky300|hky10k|jc69]",
    call. = FALSE
  )
}
raw_model <- args_trailing[1]
if (!(raw_model %in% c("hky300", "hky10k", "jc69"))) {
  stop(
    "Usage: Rscript makeGraphs3.R [hky300|hky10k|jc69]",
    call. = FALSE
  )
}
model <- raw_model
is_hky <- model %in% c("hky300", "hky10k")

# Resolve script directory so we can build absolute paths
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}

graphs_dir <- script_dir
data_dir <- switch(
  model,
  "hky300"   = file.path(graphs_dir, "hky-graphs"),
  "hky10k"   = file.path(graphs_dir, "hky-10k-graphs"),
  "jc69"     = file.path(graphs_dir, "jc69-graphs3")
)

# Feature switches for optional plots
# Set to FALSE when the corresponding data are unavailable
include_rf <- FALSE
include_mf   <- FALSE
include_dims <- FALSE
include_dist <- FALSE

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
  mrbayes   = "#E15759",
  dodonaphy = "#B07AA1",
  geophy    = "#EDC948",
  vaiphy    = "#76B7B2"
)

method_labels <- c(
  NJ        = "NJ",
  vine      = "Vine",
  beast     = "BEAST2",
  mrbayes   = "MrBayes",
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

ll_methods <- if (is_hky) {
  intersect(c("NJ", "vine", "beast", "mrbayes"), names(lnl))
} else {
  intersect(c("NJ", "vine", "beast", "mrbayes",
              "dodonaphy", "geophy", "vaiphy"),
            names(lnl))
}

lnl_long <- melt_mean_sd(lnl, methods = ll_methods,
                         ntaxa_col = "ntaxa", scale_by = "ave")
lnl_long$method <- factor(lnl_long$method, levels = ll_methods)

if (is_hky) {
  lower <- with(lnl_long, mean - sd)
  upper <- with(lnl_long, mean + sd)
  L <- max(abs(c(lower, upper)), na.rm = TRUE)
  if (!is.finite(L)) L <- 1
  pad <- 0.02 * L
  ylim_lnl <- c(-L - pad, L + pad)
} else {
  ylim_lnl <- c(-0.023, 0.01)
}

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

time_methods <- if (is_hky) {
  intersect(c("vine", "beast", "mrbayes"), names(time))
} else {
  intersect(c("vine", "beast", "mrbayes",
              "dodonaphy", "geophy", "vaiphy"),
            names(time))
}

if (model == "jc69") {
  exclude_rows_time <- time$ntaxa %in% c(25, 50, 100)
  omit_methods_time <- intersect(c("dodonaphy", "geophy", "vaiphy"),
                                 time_methods)
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
}

time_long <- melt_mean_sd(time, methods = time_methods,
                          ntaxa_col = "ntaxa", scale_by = NULL)
time_long$method <- factor(time_long$method, levels = time_methods)

if (model == "jc69" || model == "hky300") {
  valid <- is.finite(time_long$mean) & is.finite(time_long$sd) &
    time_long$mean > 0 & time_long$sd >= 0
  cv <- time_long$sd / time_long$mean
  sigma_log <- sqrt(log1p(cv^2))
  c_mult <- exp(sigma_log)
  time_long$ymin_mult <- ifelse(valid, time_long$mean / c_mult, NA_real_)
  time_long$ymax_mult <- ifelse(valid, time_long$mean * c_mult, NA_real_)
}

if (is_hky) {
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
} else {
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
}

if (model == "jc69") {
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
}

# -------------------- jc69: Lnl + Time panels (10/15/20) --------------------
if (model == "jc69") {
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
      linewidth = 3
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
} else if (model == "hky300") {
  sub_time <- time_long
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

  ptime_sub <- ggplot(
    time_long,
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
      linewidth = 3
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
}

# ================================================================
# 3) Speedup
# ================================================================
spd <- read.table(file.path(data_dir, "speedSummary.txt"),
                  header = TRUE, check.names = FALSE)

speed_methods <- intersect(c("beast", "mrbayes"), names(spd))
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
  guides(fill = guide_legend(override.aes = list(width = 0.6)))
pspeed <- pspeed + theme(legend.position = "none")
pspeed <- pspeed + guides(fill = "none")

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

  if (is_hky) {
    lower <- with(mf_long, mean - sd)
    upper <- with(mf_long, mean + sd)
    L <- max(abs(c(lower, upper)))
    pad <- 0.02 * L
    ylim_mf <- c(-L - pad, L + pad)
    fill_vals <- method_palette[levels(mf_long$method)]
    fill_breaks <- levels(mf_long$method)
    fill_labels <- label_map(levels(mf_long$method))
  } else {
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
  }

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
      drop = (model == "jc69")
    ) +
    scale_y_continuous(
      labels = label_percent(accuracy = 0.1),
      breaks = pretty_breaks(n = 8)
    ) +
    guides(fill = "none") +
    theme(legend.position = "none") +
    {
      if (is_hky) coord_cartesian(ylim = ylim_mf)
      else NULL
    }

  # Combined: lnl deviations (left) and model-fit deviations (right)
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

# -------------------- Dimensionality --------------------
if (include_dims && model == "hky300") {
  # Read wide-format data
  # Expected columns (by position):
  # 1:d, 2:vine, 3:std, 4:beast, 5:time, 6:std
  dims <- read.table(file.path(data_dir, "dimSummary25.txt"), header = TRUE, check.names = FALSE)

  # Extract numeric embedding dimension from strings like "archive.D2"
  dims$d <- suppressWarnings(as.numeric(gsub(".*\\.D", "", as.character(dims[[1]]))))

  # Build metrics & correctly paired SDs
  # Δlnl (primary axis)
  delta_lnl <- dims[[2]] - dims[[4]]    # vine - beast
  sd_lnl    <- dims[[3]]                # std immediately after vine

  # Time (secondary axis)
  time_mean <- dims[[5]]
  sd_time   <- dims[[6]]

  # Shared scale factor for dual axis
  num_max <- max(abs(delta_lnl), na.rm = TRUE)
  den_max <- max(time_mean, na.rm = TRUE)
  scale_factor <- if (is.finite(num_max) && is.finite(den_max) && den_max > 0) num_max / den_max else 1

  # Long data with error bars
  plotdata <- bind_rows(
    # Δlnl (primary axis; keep negative values)
    data.frame(d = dims$d, metric = "delta_lnl",
              mean_scaled = delta_lnl, sd_scaled = sd_lnl,
              stringsAsFactors = FALSE),
    # time (secondary axis; scale to primary for plotting)
    data.frame(d = dims$d, metric = "time",
              mean_scaled = time_mean * scale_factor,
              sd_scaled = sd_time * scale_factor,
              stringsAsFactors = FALSE)
  )
  plotdata$metric <- factor(plotdata$metric, levels = c("delta_lnl", "time"))

  # Axis breaks (only ≤0 for Δlnl; only ≥0 for time)
  # Left axis (primary, Δlnl):
  ll_min <- suppressWarnings(min(delta_lnl, 0, na.rm = TRUE))
  ll_breaks <- pretty(c(ll_min, 0), n = 6)
  ll_breaks <- ll_breaks[ll_breaks <= 0]

  # Right axis (secondary, time in seconds):
  t_max <- suppressWarnings(max(time_mean, 0, na.rm = TRUE))
  t_breaks <- pretty(c(0, t_max), n = 6)
  t_breaks <- t_breaks[t_breaks >= 0]

  combined_breaks <- sort(unique(c(ll_breaks, t_breaks * scale_factor)))

  # plotting
  dodge <- position_dodge(width = 0.6)
  geom_palette <- c("delta_lnl" = "#F28E2B", "time" = "#59A14F")
  orange <- unname(geom_palette["delta_lnl"])
  green  <- unname(geom_palette["time"])

  pdim <- ggplot(plotdata, aes(x = factor(d), y = mean_scaled, fill = metric)) +
    geom_col(width = 0.45, position = dodge) +
    geom_errorbar(aes(ymin = mean_scaled - sd_scaled, ymax = mean_scaled + sd_scaled),
                  width = 0.15, position = dodge) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = geom_palette, guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(
      name = expression(Delta~lnl),
      breaks = combined_breaks,
      labels = function(y) ifelse(y <= 0, as.character(y), ""),   # no labels for >0 on left axis
      sec.axis = sec_axis(~ . / scale_factor,
                          name = "time (s)",
                          breaks = t_breaks)
    ) +
    labs(x = "Embedding Dimension (d)", title = NULL) +
    theme(
      # Color the left (Δlnl) axis title, ticks, and labels orange
      axis.title.y      = element_text(color = orange),
      axis.text.y       = element_text(color = orange),
      axis.ticks.y      = element_line(color = orange),

      # Color the right (time) axis title, ticks, and labels green
      axis.title.y.right = element_text(color = green),
      axis.text.y.right  = element_text(color = green),
      axis.ticks.y.right = element_line(color = green)
    )

  # now do the same for 50 taxa
  # Read wide-format data
  # Expected columns (by position):
  # 1:d, 2:vine, 3:std, 4:beast, 5:time, 6:std
  dims <- read.table(file.path(data_dir, "dimSummary50.txt"), header = TRUE, check.names = FALSE)

  # Extract numeric embedding dimension from strings like "archive.D2"
  dims$d <- suppressWarnings(as.numeric(gsub(".*\\.D", "", as.character(dims[[1]]))))

  # Build metrics & correctly paired SDs
  # Δlnl (primary axis)
  delta_lnl <- dims[[2]] - dims[[4]]    # vine - beast
  sd_lnl    <- dims[[3]]                # std immediately after vine

  # Time (secondary axis)
  time_mean <- dims[[5]]
  sd_time   <- dims[[6]]

  # Shared scale factor for dual axis
  num_max <- max(abs(delta_lnl), na.rm = TRUE)
  den_max <- max(time_mean, na.rm = TRUE)
  scale_factor <- if (is.finite(num_max) && is.finite(den_max) && den_max > 0) num_max / den_max else 1

  # Long data with error bars
  plotdata <- bind_rows(
    # Δlnl (primary axis; keep negative values)
    data.frame(d = dims$d, metric = "delta_lnl",
              mean_scaled = delta_lnl, sd_scaled = sd_lnl,
              stringsAsFactors = FALSE),
    # time (secondary axis; scale to primary for plotting)
    data.frame(d = dims$d, metric = "time",
              mean_scaled = time_mean * scale_factor,
              sd_scaled = sd_time * scale_factor,
              stringsAsFactors = FALSE)
  )
  plotdata$metric <- factor(plotdata$metric, levels = c("delta_lnl", "time"))

  # Axis breaks (only ≤0 for Δlnl; only ≥0 for time)
  # Left axis (primary, Δlnl):
  ll_min <- suppressWarnings(min(delta_lnl, 0, na.rm = TRUE))
  ll_breaks <- pretty(c(ll_min, 0), n = 6)
  ll_breaks <- ll_breaks[ll_breaks <= 0]

  # Right axis (secondary, time in seconds):
  t_max <- suppressWarnings(max(time_mean, 0, na.rm = TRUE))
  t_breaks <- pretty(c(0, t_max), n = 6)
  t_breaks <- t_breaks[t_breaks >= 0]

  combined_breaks <- sort(unique(c(ll_breaks, t_breaks * scale_factor)))

  # plotting
  dodge <- position_dodge(width = 0.6)
  geom_palette <- c("delta_lnl" = "#F28E2B", "time" = "#59A14F")
  orange <- unname(geom_palette["delta_lnl"])
  green  <- unname(geom_palette["time"])

  pdim2 <- ggplot(plotdata, aes(x = factor(d), y = mean_scaled, fill = metric)) +
    geom_col(width = 0.45, position = dodge) +
    geom_errorbar(aes(ymin = mean_scaled - sd_scaled, ymax = mean_scaled + sd_scaled),
                  width = 0.15, position = dodge) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = geom_palette, guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(
      name = expression(Delta~lnl),
      breaks = combined_breaks,
      labels = function(y) ifelse(y <= 0, as.character(y), ""),   # no labels for >0 on left axis
      sec.axis = sec_axis(~ . / scale_factor,
                          name = "time (s)",
                          breaks = t_breaks)
    ) +
    labs(x = "Embedding Dimension (d)", title = NULL) +
    theme(
      # Color the left (Δlnl) axis title, ticks, and labels orange
      axis.title.y      = element_text(color = orange),
      axis.text.y       = element_text(color = orange),
      axis.ticks.y      = element_line(color = orange),

      # Color the right (time) axis title, ticks, and labels green
      axis.title.y.right = element_text(color = green),
      axis.text.y.right  = element_text(color = green),
      axis.ticks.y.right = element_line(color = green)
    )

  # Combined panels with A/B callouts
  pdim_panels <- (pdim + pdim2) +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      tag_prefix = "",
      tag_suffix = ""
    )
  pdim_panels <- pdim_panels & theme(
    plot.tag = element_text(
      size = 12, face = "bold", family = "Helvetica"
    ),
    plot.tag.position = c(0.01, 0.99)
  )
  save_pdf(pdim_panels, file.path(data_dir, "dims_panels.pdf"), width = 6, height = 3)
}

#--------------------- Dimensionality H Version ---------------------
if (include_dims && model == "hky300") {
    
  # Read wide-format data
  # Expected columns (by position):
  # 1:d, 2:vine, 3:std, 4:beast, 5:time, 6:std
  dims <- read.table(file.path(data_dir, "dimSummary25.txt"), header = TRUE, check.names = FALSE)

  # Extract numeric embedding dimension from strings like "archive.D2"
  dims$d <- suppressWarnings(as.numeric(gsub(".*\\.D", "", as.character(dims[[1]]))))

  # Build metrics & correctly paired SDs
  # Δlnl (primary axis)
  delta_lnl <- dims[[7]] - dims[[4]]    # vine - beast
  sd_lnl    <- dims[[8]]                # std immediately after vine

  # Time (secondary axis)
  time_mean <- dims[[9]]
  sd_time   <- dims[[10]]

  # Shared scale factor for dual axis
  num_max <- max(abs(delta_lnl), na.rm = TRUE)
  den_max <- max(time_mean, na.rm = TRUE)
  scale_factor <- if (is.finite(num_max) && is.finite(den_max) && den_max > 0) num_max / den_max else 1

  # Long data with error bars
  plotdata <- bind_rows(
    # Δlnl (primary axis; keep negative values)
    data.frame(d = dims$d, metric = "delta_lnl",
              mean_scaled = delta_lnl, sd_scaled = sd_lnl,
              stringsAsFactors = FALSE),
    # time (secondary axis; scale to primary for plotting)
    data.frame(d = dims$d, metric = "time",
              mean_scaled = time_mean * scale_factor,
              sd_scaled = sd_time * scale_factor,
              stringsAsFactors = FALSE)
  )
  plotdata$metric <- factor(plotdata$metric, levels = c("delta_lnl", "time"))

  # Axis limits/breaks include error bar extremes
  # Left axis (Δlnl) includes the lowest error bar tip
  ll_min_err <- suppressWarnings(
    min(delta_lnl - sd_lnl, 0, na.rm = TRUE)
  )
  ll_breaks <- pretty(c(ll_min_err, 0), n = 6)
  if (length(ll_breaks) > 1 && min(ll_breaks) > ll_min_err) {
    step <- ll_breaks[2] - ll_breaks[1]
    ll_breaks <- c(min(ll_breaks) - step, ll_breaks)
  }
  ymin <- min(ll_breaks, na.rm = TRUE)

  # Right axis (time) includes the highest error bar tip
  t_max_err <- suppressWarnings(
    max(time_mean + sd_time, 0, na.rm = TRUE)
  )
  t_breaks <- pretty(c(0, t_max_err), n = 6)
  if (length(t_breaks) > 1 && max(t_breaks) < t_max_err) {
    step <- t_breaks[2] - t_breaks[1]
    t_breaks <- c(t_breaks, max(t_breaks) + step)
  }
  ymax <- max(t_breaks, na.rm = TRUE) * scale_factor

  combined_breaks <- sort(
    unique(c(ll_breaks, t_breaks * scale_factor))
  )

  # plotting
  dodge <- position_dodge(width = 0.6)
  geom_palette <- c("delta_lnl" = "#F28E2B", "time" = "#59A14F")
  orange <- unname(geom_palette["delta_lnl"])
  green  <- unname(geom_palette["time"])

  pdim <- ggplot(plotdata, aes(x = factor(d), y = mean_scaled, fill = metric)) +
    geom_col(width = 0.45, position = dodge) +
    geom_errorbar(aes(ymin = mean_scaled - sd_scaled, ymax = mean_scaled + sd_scaled),
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
      # Color the left (Δlnl) axis title, ticks, and labels orange
      axis.title.y      = element_text(color = orange),
      axis.text.y       = element_text(color = orange),
      axis.ticks.y      = element_line(color = orange),

      # Color the right (time) axis title, ticks, and labels green
      axis.title.y.right = element_text(color = green),
      axis.text.y.right  = element_text(color = green),
      axis.ticks.y.right = element_line(color = green)
    )

  # now do the same for 50 taxa
  # Read wide-format data
  # Expected columns (by position):
  # 1:d, 2:vine, 3:std, 4:beast, 5:time, 6:std
  dims <- read.table(file.path(data_dir, "dimSummary50.txt"), header = TRUE, check.names = FALSE)

  # Extract numeric embedding dimension from strings like "archive.D2"
  dims$d <- suppressWarnings(as.numeric(gsub(".*\\.D", "", as.character(dims[[1]]))))

  # Build metrics & correctly paired SDs
  # Δlnl (primary axis)
  delta_lnl <- dims[[7]] - dims[[4]]    # vine - beast
  sd_lnl    <- dims[[8]]                # std immediately after vine

  # Time (secondary axis)
  time_mean <- dims[[9]]
  sd_time   <- dims[[10]]

  # Shared scale factor for dual axis
  num_max <- max(abs(delta_lnl), na.rm = TRUE)
  den_max <- max(time_mean, na.rm = TRUE)
  scale_factor <- if (is.finite(num_max) && is.finite(den_max) && den_max > 0) num_max / den_max else 1

  # Long data with error bars
  plotdata <- bind_rows(
    # Δlnl (primary axis; keep negative values)
    data.frame(d = dims$d, metric = "delta_lnl",
              mean_scaled = delta_lnl, sd_scaled = sd_lnl,
              stringsAsFactors = FALSE),
    # time (secondary axis; scale to primary for plotting)
    data.frame(d = dims$d, metric = "time",
              mean_scaled = time_mean * scale_factor,
              sd_scaled = sd_time * scale_factor,
              stringsAsFactors = FALSE)
  )
  plotdata$metric <- factor(plotdata$metric, levels = c("delta_lnl", "time"))

  # Axis limits/breaks include error bar extremes
  # Left axis (Δlnl) includes the lowest error bar tip
  ll_min_err <- suppressWarnings(
    min(delta_lnl - sd_lnl, 0, na.rm = TRUE)
  )
  ll_breaks <- pretty(c(ll_min_err, 0), n = 6)
  if (length(ll_breaks) > 1 && min(ll_breaks) > ll_min_err) {
    step <- ll_breaks[2] - ll_breaks[1]
    ll_breaks <- c(min(ll_breaks) - step, ll_breaks)
  }
  ymin <- min(ll_breaks, na.rm = TRUE)

  # Right axis (time) includes the highest error bar tip
  t_max_err <- suppressWarnings(
    max(time_mean + sd_time, 0, na.rm = TRUE)
  )
  t_breaks <- pretty(c(0, t_max_err), n = 6)
  if (length(t_breaks) > 1 && max(t_breaks) < t_max_err) {
    step <- t_breaks[2] - t_breaks[1]
    t_breaks <- c(t_breaks, max(t_breaks) + step)
  }
  ymax <- max(t_breaks, na.rm = TRUE) * scale_factor

  combined_breaks <- sort(
    unique(c(ll_breaks, t_breaks * scale_factor))
  )

  # plotting
  dodge <- position_dodge(width = 0.6)
  geom_palette <- c("delta_lnl" = "#F28E2B", "time" = "#59A14F")
  orange <- unname(geom_palette["delta_lnl"])
  green  <- unname(geom_palette["time"])

  pdim2 <- ggplot(plotdata, aes(x = factor(d), y = mean_scaled, fill = metric)) +
    geom_col(width = 0.45, position = dodge) +
    geom_errorbar(aes(ymin = mean_scaled - sd_scaled, ymax = mean_scaled + sd_scaled),
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
      # Color the left (Δlnl) axis title, ticks, and labels orange
      axis.title.y      = element_text(color = orange),
      axis.text.y       = element_text(color = orange),
      axis.ticks.y      = element_line(color = orange),

      # Color the right (time) axis title, ticks, and labels green
      axis.title.y.right = element_text(color = green),
      axis.text.y.right  = element_text(color = green),
      axis.ticks.y.right = element_line(color = green)
    )

  # Combined panels with A/B callouts
  pdim_panels <- (pdim + pdim2) +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      tag_prefix = "",
      tag_suffix = ""
    )
  pdim_panels <- pdim_panels & theme(
    plot.tag = element_text(
      size = 12, face = "bold", family = "Helvetica"
    ),
    plot.tag.position = c(0.01, 0.99)
  )
  save_pdf(pdim_panels, file.path(data_dir, "dimsH_panels.pdf"),
          width = 6, height = 3)
}
  # -------------------- Robinson-Foulds distances --------------------
if (include_rf) {
  if (model == "hky300") {
  rf <- read.table(file.path(data_dir, "rfSummary.txt"), header = TRUE)

  ## Identify the three std columns in order (std, std.1, std.2, ...)
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
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "Number of Taxa",
      y = "Robinson-Foulds Distance",
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

      save_pdf(prf,
            file.path(data_dir, "hky300_rf_bars.pdf"),
            width = 3, height = 3)
  }
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
if (model == "hky300") {
  if (include_dist) {
    save_pdf(
      pdist_flows,
      file.path(data_dir, "hky300_dist_flows_bars.pdf"),
      width = 3, height = 3
    )
    save_pdf(
      pentropy,
      file.path(data_dir, "hky300_entropy_bars.pdf"),
      width = 3, height = 3
    )
  }
  save_pdf(pspeed, file.path(data_dir, "hky300_speedup.pdf"),
          width = 3, height = 3)
  date_stamp <- format(Sys.Date(), "%m%d%y")
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
  if (include_mf) {
    save_pdf(
      pcombo,
      file.path(
        data_dir,
        paste0("hky300_lnl_mf_panels_", date_stamp, ".pdf")
      ),
      width = 6, height = 3
    )
  }
  save_pdf(ptime_speed,
           file.path(
             data_dir,
             paste0("hky300_time_speed_panels_",
                    date_stamp, ".pdf")
           ),
           width = 6, height = 3)
  if (include_mf) {
    plnl_mf_logtime <- (plnl + pmf + ptime_sub) +
      plot_layout(ncol = 3, guides = "collect") +
      plot_annotation(tag_levels = "A",
                      tag_prefix = "", tag_suffix = "")
    plnl_mf_logtime <- plnl_mf_logtime & theme(
      legend.position = "right",
      plot.tag = element_text(size = 12, face = "bold",
                              family = "Helvetica"),
      plot.tag.position = c(0.01, 0.99)
    )
    save_pdf(
      plnl_mf_logtime,
      file.path(
        data_dir,
        paste0("hky300_lnl_mf_logtime_panels_",
               date_stamp, ".pdf")
      ),
      width = 9, height = 3
    )
  }
           
} else if (model == "hky10k") {
  if (include_dist) {
    save_pdf(
      pdist, file.path(data_dir, "hky10k_dist_bars.pdf"),
      width = 3, height = 3
    )
    save_pdf(
      pdist_flows,
      file.path(data_dir, "hky10k_dist_flows_bars.pdf"),
      width = 3, height = 3
    )
  }
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
  if (include_mf) {
    save_pdf(
      pcombo,
      file.path(data_dir, "hky10k_lnl_mf_panels.pdf"),
      width = 6, height = 3
    )
  }
  save_pdf(ptime_speed,
           file.path(data_dir, "hky10k_time_speed_panels.pdf"),
           width = 6, height = 3)
} else {
  save_pdf(pspeed, file.path(data_dir, "jc69_speedup.pdf"),
           width = 3, height = 3)
  ptime_split_speed <- (ptime_first + ptime_last + pspeed) +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(tag_levels = "A",
                    tag_prefix = "", tag_suffix = "")
  ptime_split_speed <- ptime_split_speed & theme(
    legend.position = "right",
    plot.tag = element_text(size = 12, face = "bold",
                            family = "Helvetica"),
    plot.tag.position = c(0.01, 0.99)
  )
  save_pdf(ptime_split_speed,
           file.path(data_dir, "jc69_time_split_speed_panels.pdf"),
           width = 9, height = 3)
  if (include_dist) {
    save_pdf(
      pdist, file.path(data_dir, "jc69_dist_bars.pdf"),
      width = 3, height = 3
    )
    save_pdf(
      pdist_flows,
      file.path(data_dir, "jc69_dist_flows_bars.pdf"),
      width = 3, height = 3
    )
  }
  if (include_mf) {
    save_pdf(
      pcombo,
      file.path(data_dir, "jc69_lnl_mf_panels.pdf"),
      width = 6, height = 3
    )
  }
  save_pdf(
    lnl_time_panels,
    file.path(data_dir, "jc69_lnl_time_panels.pdf"),
    width = 6, height = 3
  )
}


