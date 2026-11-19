#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))

## -------------------- Colors --------------------
method_palette <- c(
  NJ   = "#4E79A7",
  vine = "#F28E2B",
  beast= "#59A14F"
)
fill_colors <- c(
  laml = unname(method_palette["beast"]),
  vine = unname(method_palette["vine"])
)
legend_labels <- c(laml = "LAML", vine = "Vine")

## -------------------- Read LNL data --------------------
all_lines <- readLines("summary.lnl.txt", warn = FALSE)
all_lines <- trimws(all_lines)
all_lines <- all_lines[nzchar(all_lines)]
header <- all_lines[1]
rows   <- all_lines[-1]
rows <- rows[!grepl("^[-]+$", rows)]
rows <- rows[sapply(strsplit(rows, "\\s+"), function(x) length(x) >= 4)]
raw <- read.table(text = c(header, rows), header = TRUE, stringsAsFactors = FALSE)
raw$cp <- ifelse(raw$cp == "all", "ave", raw$cp)

cp_numeric <- suppressWarnings(as.numeric(raw$cp))
is_num     <- !is.na(cp_numeric)
cp_levels  <- c(
  as.character(sort(unique(cp_numeric[is_num]))),
  if ("ave" %in% raw$cp) "ave" else character(0)
)

df_long <- rbind(
  data.frame(cp = raw$cp, method = "laml", value = as.numeric(raw$laml)),
  data.frame(cp = raw$cp, method = "vine", value = as.numeric(raw$vine))
)
df_long$cp     <- factor(df_long$cp, levels = cp_levels)
df_long$method <- factor(df_long$method, levels = c("laml", "vine"))

## -------------------- Read time data --------------------
t_lines <- readLines("summary.time.txt", warn = FALSE)
t_lines <- trimws(t_lines)
t_lines <- t_lines[nzchar(t_lines)]
t_header <- t_lines[1]
t_rows   <- t_lines[-1]
t_rows <- t_rows[!grepl("^[-]+$", t_rows)]
t_rows <- t_rows[sapply(strsplit(t_rows, "\\s+"), function(x) length(x) >= 4)]
time_raw <- read.table(text = c(t_header, t_rows), header = TRUE, stringsAsFactors = FALSE)
time_raw$cp <- ifelse(time_raw$cp == "all", "ave", time_raw$cp)

t_cp_num <- suppressWarnings(as.numeric(time_raw$cp))
t_is_num <- !is.na(t_cp_num)
t_levels <- c(
  as.character(sort(unique(t_cp_num[t_is_num]))),
  if ("ave" %in% time_raw$cp) "ave" else character(0)
)

time_long <- rbind(
  data.frame(cp = time_raw$cp, method = "laml", value = as.numeric(time_raw$laml)),
  data.frame(cp = time_raw$cp, method = "vine", value = as.numeric(time_raw$vine))
)
time_long$cp     <- factor(time_long$cp, levels = t_levels)
time_long$method <- factor(time_long$method, levels = c("laml", "vine"))
time_long <- subset(time_long, is.finite(value) & value > 0)

## -------------------- "ave" highlight geometry --------------------
# scale_y_discrete(limits=rev) reverses BOTH visual order AND the internal
# numeric coordinate mapping.  After reversal:
#   y=1  -> "ave"  (bottom)
#   y=2  -> "100"  (second from bottom)
#   ...
#   y=n  -> "4"    (top)
# So to shade "ave" and place the separator above it, use y=1.

# (highlight geometry defined inline via annotate() in each plot)

## -------------------- Common theme --------------------
base_sz <- 9
panel_theme <- theme_minimal(base_size = base_sz) +
  theme(
    axis.title.x     = element_text(size = 9),
    axis.title.y     = element_text(size = 9),
    axis.text.x      = element_text(size = 7),
    axis.text.y      = element_text(size = 5.5),
    panel.grid.major = element_line(color = "gray70", linewidth = 0.25),
    panel.grid.minor = element_line(color = "gray85", linewidth = 0.15),
    plot.margin      = margin(18, 6, 4, 4, "pt")
  )

## -------------------- Left panel: Log likelihood --------------------
min_lnl <- min(df_long$value, na.rm = TRUE)

p_lnl <- ggplot(df_long, aes(y = cp, x = value, fill = method)) +
  # Gray background for "ave" row
  annotate("rect", ymin = 0.5, ymax = 1.5, xmin = -Inf, xmax = Inf,
           fill = "gray80", alpha = 0.4) +
  # Separator line just above "ave"
  geom_hline(
    yintercept = 1.5,
    color = unname(method_palette["NJ"]), linewidth = 0.5
  ) +
  geom_col(position = position_dodge(width = 0.7), width = 0.55) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray40") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(
    limits = c(min_lnl * 1.02, 0),
    breaks = pretty_breaks(n = 5),
    expand = c(0, 0),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  scale_fill_manual(
    values = fill_colors, breaks = c("laml", "vine"), labels = legend_labels
  ) +
  labs(x = "Log likelihood", y = "CP", fill = NULL) +
  panel_theme

## -------------------- Right panel: Time --------------------
p_time <- ggplot(time_long, aes(y = cp, x = value, fill = method)) +
  # Gray background for "ave" row
  annotate("rect", ymin = 0.5, ymax = 1.5, xmin = -Inf, xmax = Inf,
           fill = "gray80", alpha = 0.4) +
  # Separator line just above "ave"
  geom_hline(
    yintercept = 1.5,
    color = unname(method_palette["NJ"]), linewidth = 0.5
  ) +
  geom_col(position = position_dodge(width = 0.7), width = 0.55) +
  scale_y_discrete(limits = rev) +
  scale_x_log10(
    breaks = scales::breaks_log(8),
    labels = scales::label_number(scale_cut = NULL, big.mark = "", decimal.mark = "."),
    expand = expansion(mult = c(0, 0.08))
  ) +
  annotation_logticks(sides = "b", size = 0.25) +
  scale_fill_manual(
    values = fill_colors, breaks = c("laml", "vine"), labels = legend_labels
  ) +
  labs(x = "Time (sec)", y = NULL, fill = NULL) +
  panel_theme +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.3, color = "gray40")
  )

## -------------------- Combine with patchwork --------------------
combined <- p_lnl + p_time +
  plot_layout(widths = c(1.1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.35, "cm"),
    legend.key.spacing.x = unit(4, "pt"),
    plot.tag          = element_text(size = 12, face = "bold", family = "Helvetica"),
    plot.tag.location = "panel",
    plot.tag.position = c(0, 1.025)
  )

## -------------------- Save --------------------
out_file <- "lnl_time_combined.pdf"
if (capabilities("cairo")) {
  ggsave(out_file, plot = combined,
         width = 7.5, height = 8.0, units = "in", device = cairo_pdf)
} else {
  pdf(out_file, width = 6, height = 8.0, family = "Helvetica")
  print(combined)
  dev.off()
}

cat("Saved:", out_file, "\n")
