#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
use_relative_lnl <- "--relative" %in% args

suppressMessages(library(ggplot2))
suppressMessages(library(scales))

## -------------------- Styling (matches your prior scripts) --------------------
theme_set(
  theme_minimal(base_size = 8) +
    theme(
      plot.title   = element_text(size = 9, face = "bold"),
      axis.title   = element_text(size = 9),
      axis.text    = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7),
      panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.2)
    )
)

## Helper: save a plot as a compact PDF
save_pdf <- function(plot, filename, width = 6.8, height = 3.0) {
  if (capabilities("cairo")) {
    ggsave(filename, plot = plot, width = width, height = height,
           units = "in", device = cairo_pdf)
  } else {
    pdf(filename, width = width, height = height, family = "Helvetica")
    print(plot); dev.off()
  }
}

## -------------------- Colors & Labels (consistent with earlier) --------------------
method_palette <- c(
  NJ     = "#4E79A7",  # blue
  vine   = "#F28E2B",  # orange
  beast  = "#59A14F",  # green
  mrbayes  = "#E15759"  # red
)
method_labels <- c(
  NJ     = "NJ",
  vine   = "Vine",
  beast  = "BEAST2",
  mrbayes  = "MrBayes"
)
label_map <- function(keys) { out <- method_labels[keys]; out[is.na(out)] <- keys[is.na(out)]; out }

## -------------------- Utilities --------------------
# Read a table with possible dashed separator lines and trailing "all" row.
read_eval_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  stopifnot(length(lines) >= 2)
  header <- lines[1]
  rows   <- lines[-1]
  rows <- rows[!grepl("^[-]+$", rows)]                                     # drop "-----"
  rows <- rows[sapply(strsplit(rows, "\\s+"), function(x) length(x) >= 2)] # keep well-formed
  read.table(text = c(header, rows), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

# Build DS factor: ds1..dsn in input order (or numeric order if present), then "ave"
build_ds_levels <- function(ds_vec) {
  ds_vec2 <- ifelse(ds_vec == "all", "ave", ds_vec)
  # Keep original order for ds1..dsn if it looks like ds<number>, then append "ave"
  ds_main <- ds_vec2[ds_vec2 != "ave"]
  ds_levels <- unique(ds_main)
  if ("ave" %in% ds_vec2) ds_levels <- c(ds_levels, "ave")
  ds_levels
}

# Robust numeric conversion for fields like "-6866.92," (trailing comma)
numify <- function(x) as.numeric(gsub(",", "", x, fixed = TRUE))

## ====================================================
## 1) LOG-LIKELIHOODS (eval.all.lnl.txt) — NJ/vine/BEAST/RAxML
## ====================================================
lnl_raw <- read_eval_table("eval.all.lnl.txt")

# Normalize column names we care about; ignore 'ml'
# First col header is 's' per your file; rename it to 'ds' for consistency
names(lnl_raw)[1] <- "ds"

# Keep only ds + desired methods, in this order
keep_methods_lnl <- c("NJ","vine","beast","mrbayes")
# Column names in file are lowercase for methods except NJ spelled 'nj'? Let's normalize to match keep_methods_lnl
names(lnl_raw) <- sub("^nj$", "NJ", names(lnl_raw), ignore.case = TRUE)

present_lnl <- intersect(keep_methods_lnl, names(lnl_raw))
lnl_df <- lnl_raw[, c("ds", present_lnl), drop = FALSE]

# Rename "all" to "ave"
lnl_df$ds <- ifelse(lnl_df$ds == "all", "ave", lnl_df$ds)
ds_levels <- build_ds_levels(lnl_df$ds)

# Map ds labels to numeric indices (keep 'ave' literal)
ds_num_map <- {
  main <- ds_levels[ds_levels != "ave"]
  nums <- as.character(seq_along(main))
  out <- setNames(nums, main)
  if ("ave" %in% ds_levels) out["ave"] <- "ave"
  out
}

# Long format
lnl_long <- do.call(rbind, lapply(present_lnl, function(m) {
  data.frame(
    ds = lnl_df$ds,
    method = m,
    value = numify(lnl_df[[m]]),
    stringsAsFactors = FALSE
  )
}))
if (use_relative_lnl) {
  lnl_means <- tapply(lnl_long$value, lnl_long$ds, mean, na.rm = TRUE)
  lnl_long$ds_mean  <- lnl_means[lnl_long$ds]
  lnl_long$plot_y <- (lnl_long$value - lnl_long$ds_mean) / abs(lnl_long$ds_mean) * 100
} else {
  lnl_long$plot_y <- lnl_long$value
}

lnl_long$ds     <- factor(lnl_long$ds, levels = ds_levels)
lnl_long$method <- factor(lnl_long$method, levels = present_lnl)

# Highlight "ave"
ave_idx <- match("ave", levels(lnl_long$ds))
ave_bg  <- if (!is.na(ave_idx)) data.frame(xmin = ave_idx - 0.6, xmax = ave_idx + 0.6, ymin = -Inf, ymax = Inf) else NULL

if (use_relative_lnl) {
  y_label <- "ΔLnl (%)"
  y_scale <- scale_y_continuous(breaks = pretty_breaks(n = 8),
                                expand = expansion(mult = c(0.05, 0.05)))
} else {
  min_y <- min(lnl_long$plot_y, na.rm = TRUE)
  if (!is.finite(min_y) || min_y >= 0) min_y <- -1
  y_label <- "Log-likelihood"
  y_scale <- scale_y_continuous(limits = c(min_y, 0),
                                breaks = pretty_breaks(n = 8),
                                expand = c(0, 0))
}

p_lnl_ds <- ggplot(lnl_long, aes(x = ds, y = plot_y, fill = method)) +
  { if (!is.null(ave_bg))
      geom_rect(data = ave_bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "gray80", alpha = 0.3) } +
  { if (!is.na(ave_idx))
      geom_vline(xintercept = ave_idx - 0.7, color = unname(method_palette["NJ"]), linewidth = 0.5) } +
  geom_col(position = position_dodge(width = 0.7), width = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = "DS", y = y_label, fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(lnl_long$method)],
    breaks = levels(lnl_long$method),
    labels = label_map(levels(lnl_long$method))
  ) +
  y_scale +
  scale_x_discrete(
    labels = function(x) unname(ds_num_map[x]),
    expand = expansion(add = c(0.05, 0.25))
  ) +
  theme(
    axis.text.x = element_text(
      vjust = 1, hjust = 1, size = 7
    )
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

save_pdf(p_lnl_ds, "dna_lnl_by_ds.pdf", width = 6.8, height = 3.0)

## ===============================================
## 2) RUNNING TIMES (eval.all.time.txt) — vine/BEAST/RAxML
## ===============================================
time_raw <- read_eval_table("eval.all.time.txt")

# First col header is 'ds' per your file; ensure it is named 'ds'
names(time_raw)[1] <- "ds"

# Methods present (omit NJ if absent)
keep_methods_time <- c("vine","beast","mrbayes")
present_time <- intersect(keep_methods_time, names(time_raw))
time_df <- time_raw[, c("ds", present_time), drop = FALSE]

# Rename "all" to "ave"
time_df$ds <- ifelse(time_df$ds == "all", "ave", time_df$ds)
t_levels <- build_ds_levels(time_df$ds)

# Map ds labels to numeric indices (keep 'ave' literal)
t_num_map <- {
  main <- t_levels[t_levels != "ave"]
  nums <- as.character(seq_along(main))
  out <- setNames(nums, main)
  if ("ave" %in% t_levels) out["ave"] <- "ave"
  out
}

# Long format
time_long <- do.call(rbind, lapply(present_time, function(m) {
  data.frame(
    ds = time_df$ds,
    method = m,
    value = as.numeric(time_df[[m]]),
    stringsAsFactors = FALSE
  )
}))
time_long$ds     <- factor(time_long$ds, levels = t_levels)
time_long$method <- factor(time_long$method, levels = present_time)

# Y scale: 0 → max
max_y <- max(time_long$value, na.rm = TRUE); if (!is.finite(max_y) || max_y <= 0) max_y <- 1

# Highlight "ave"
t_ave_idx <- match("ave", levels(time_long$ds))
t_ave_bg  <- if (!is.na(t_ave_idx)) data.frame(xmin = t_ave_idx - 0.6, xmax = t_ave_idx + 0.6, ymin = -Inf, ymax = Inf) else NULL

p_time_ds <- ggplot(time_long, aes(x = ds, y = value, fill = method)) +
  { if (!is.null(t_ave_bg))
      geom_rect(data = t_ave_bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "gray80", alpha = 0.3) } +
  { if (!is.na(t_ave_idx))
      geom_vline(xintercept = t_ave_idx - 0.7, color = unname(method_palette["NJ"]), linewidth = 0.5) } +
  geom_col(position = position_dodge(width = 0.7), width = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = "DS", y = "Time (sec)", fill = "Method") +
  scale_fill_manual(
    values = method_palette[levels(time_long$method)],
    breaks = levels(time_long$method),
    labels = label_map(levels(time_long$method))
  ) +
  scale_y_continuous(
    limits = c(0, max_y),
    breaks = pretty_breaks(n = 8),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    labels = function(x) unname(t_num_map[x]),
    expand = expansion(add = c(0.05, 0.25))
  ) +
  theme(
    axis.text.x = element_text(
      vjust = 1, hjust = 1, size = 7
    )
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6)))

save_pdf(p_time_ds, "dna_time_by_ds.pdf", width = 6.8, height = 3.0)
