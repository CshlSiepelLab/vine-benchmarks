#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))

## -------------------- Styling suitable for 3"x3" panels --------------------
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

options(ggplot2.useDingbats = FALSE)
theme_set(theme_minimal(base_size = 8, base_family = "DejaVu Sans"))

## Helper: save a plot as a compact PDF, using cairo if available
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

## -------------------- Colors & Labels --------------------
# Palette compatible with your existing scheme (Tableau-ish)
method_palette <- c(
  vine     = "#F28E2B",  # orange
  laml    = "#59A14F",  # green
  beam    = "#E15759"  # red
)

# Legend display names (edit the values on the right as you like)
method_labels <- c(
  vine     = "Vine",
  laml    = "LAML",
  beam    = "BEAM"
)

label_map <- function(keys) {
  out <- method_labels[keys]
  out[is.na(out)] <- keys[is.na(out)]
  out
}

## -------------------- Utilities --------------------
# Given a data.frame that contains <method> columns and interleaved std columns,
# build a long data.frame with columns: ntaxa, method, mean, sd
# If scale_by is non-NULL (e.g., "ave"), divide mean/sd by abs(df[[scale_by]]).
melt_mean_sd <- function(df, methods, ntaxa_col = "ntaxa", scale_by = NULL) {
  stopifnot(ntaxa_col %in% names(df))
  out <- list()
  for (m in methods) {
    if (!(m %in% names(df))) next
    m_idx <- which(names(df) == m)
    # std column is assumed to be the immediate next column and to contain "std"
    sd_idx <- if (m_idx < ncol(df) && identical(names(df)[m_idx + 1], "std")) m_idx + 1 else NA_integer_
    mean_vals <- df[[m]]
    sd_vals   <- if (!is.na(sd_idx)) df[[sd_idx]] else rep(NA_real_, length(mean_vals))
    if (!is.null(scale_by)) {
      denom <- abs(df[[scale_by]])
      mean_vals <- mean_vals / denom
      sd_vals   <- sd_vals   / denom
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

## ================================================================
## 1) Log-likelihood deviations (scaled by |ave|); center around zero
## ================================================================
lnl <- read.table("lnlSummary.txt", header = TRUE, check.names = FALSE)

# Methods available (only include those actually present)
ll_methods <- intersect(c("vine","beam","laml"), names(lnl))

# Scale by |ave| (your inputs are true log-likelihoods; no sign flip)
lnl_long <- melt_mean_sd(lnl, methods = ll_methods, ntaxa_col = "ntaxa", scale_by = NULL)
lnl_long$method <- factor(lnl_long$method, levels = ll_methods)

ylim <- range(with(lnl_long, mean + c(-sd, sd)), na.rm = TRUE)

plnl <- ggplot(lnl_long, aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  labs(x = "Number of Taxa (n)", y = "Log-likelihood", fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(lnl_long$method)],
    breaks = levels(lnl_long$method),
    labels = label_map(levels(lnl_long$method))
  ) +
  scale_y_continuous(
    labels = function(x) formatC(x, format = "f", digits = 0),
    breaks = scales::pretty_breaks(n = 8)
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

## =====================================
## 2) Running times (seconds, raw scale)
## =====================================
time <- read.table("timeSummary.txt", header = TRUE, check.names = FALSE)

time_methods <- intersect(c("vine","beam","laml"), names(time))
time_long <- melt_mean_sd(time, methods = time_methods, ntaxa_col = "ntaxa", scale_by = NULL)
time_long$method <- factor(time_long$method, levels = time_methods)

## =====================================
## 2) Running times (seconds, baseline 0.01)
## =====================================
time <- read.table("timeSummary.txt", header = TRUE, check.names = FALSE)

time_methods <- intersect(c("vine","beam","laml"), names(time))
time_long <- melt_mean_sd(time, methods = time_methods, ntaxa_col = "ntaxa", scale_by = NULL)
time_long$method <- factor(time_long$method, levels = time_methods)

## ---- manual dodging geometry for geom_rect (bars start at 0.01) ----
baseline <- 0.01
bar_width <- 0.80                         # overall group width (like dodge ~0.9)
k <- length(levels(time_long$method))
step <- bar_width / k                      # width allotted per method within a group
rect_width <- step * 0.92                  # a little gap between bars

time_long$x <- as.numeric(factor(time_long$ntaxa))  # 1..G (groups)
time_long$m <- as.numeric(time_long$method)         # 1..k

# center methods within each group
time_long$x_center <- time_long$x + (time_long$m - (k + 1)/2) * step

time_long$xmin <- time_long$x_center - rect_width/2
time_long$xmax <- time_long$x_center + rect_width/2

# bars go from baseline up to mean (but never below baseline)
time_long$ymin_bar <- baseline
time_long$ymax_bar <- pmax(time_long$mean, baseline)

# error bars: clamp the bottom at baseline as well
time_long$ymin_err <- pmax(time_long$mean - time_long$sd, baseline)
time_long$ymax_err <- time_long$mean + time_long$sd

ptime <- ggplot(time_long, aes(fill = method)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin_bar, ymax = ymax_bar),
            color = NA) +
  geom_errorbar(aes(x = x_center, ymin = ymin_err, ymax = ymax_err),
                width = rect_width * 0.35) +
  geom_hline(yintercept = baseline, linewidth = 0.3) +
  scale_x_continuous(
    breaks = sort(unique(time_long$x)),
    labels = sort(unique(time_long$ntaxa))
  ) +
  labs(x = "Number of Taxa (n)", y = "Time (s)", fill = "Method") +
  scale_y_log10(
    labels = function(x) {
        ifelse(abs(x - 0.01) < 1e-8 | abs(x - 0.03) < 1e-8,
              sprintf("%.2f", x),
              sprintf("%.1f", x))
      },
    breaks = scales::breaks_log(10),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(baseline, NA)) +
  annotation_logticks(sides = "l") +
  scale_fill_manual(
    values = method_palette[levels(time_long$method)],
    breaks = levels(time_long$method),
    labels = label_map(levels(time_long$method))
  ) +
  guides(fill = "none") +
  theme(legend.position = "none")


## =========================================
## 3) Speedup (grouped bars from speed file)
## =========================================
# speedSummary.txt: like others but NO vine columns.
# Interpreted as: speedup of vine relative to each other method (means + sds).
spd <- read.table("speedSummary.txt", header = TRUE, check.names = FALSE)

# Use whichever of these methods are present in the speed file
speed_methods <- intersect(c("beam","laml"), names(spd))
speed_long <- melt_mean_sd(spd, methods = speed_methods, ntaxa_col = "ntaxa", scale_by = NULL)
speed_long$method <- factor(speed_long$method, levels = speed_methods)

pspeed <- ggplot(speed_long, aes(x = factor(ntaxa), y = mean, fill = method)) +
        geom_col(position = position_dodge(width = 0.9)) +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                position = position_dodge(width = 0.9)
        ) +
        labs(x = "Number of Taxa (n)", y = "Speed Increase (x)", fill = "Method") +
        scale_fill_manual(
                values = method_palette[levels(speed_long$method)],
                breaks = levels(speed_long$method),
                labels = label_map(levels(speed_long$method))
        ) +
        guides(fill = "none") +
        theme(legend.position = "none")
#           guide_legend(override.aes = list(width = 0.6)))

## =========================================
## Combined 3-panel figure with callouts A/B/C
## =========================================
pcombo <- (plnl + ptime + pspeed) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = ""
  )
pcombo <- pcombo & theme(
  plot.tag = element_text(
    size = 14, face = "bold", family = "Helvetica"
  ),
  plot.tag.position = c(0.01, 0.99)
)

## -------------------- Save --------------------
save_pdf(plnl,   "crispr_lnl_bars.pdf",  width = 3, height = 3)
save_pdf(ptime,  "crispr_time_bars.pdf", width = 3, height = 3)
save_pdf(pspeed, "crispr_speedup.pdf",   width = 3, height = 3)
save_pdf(pcombo, "crispr_sim_panels.pdf",    width = 9, height = 3)
