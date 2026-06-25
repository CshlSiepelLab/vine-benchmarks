#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))

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

rf <- read.table("rfSummary.txt", header = TRUE)

## Identify the two std columns in order (std, std.1, ...)
std_cols <- grep("^std", names(rf))
if (length(std_cols) < 2) {
  stop("Expected at least 2 'std' columns (for vine, beam, laml). Found: ",
       length(std_cols))
}

# Long format (absolute RF values)
rf_long <- rbind(
  data.frame(
    ntaxa  = rf$ntaxa,
    method = "vine",
    mean   = rf$vine,
    sd     = rf[[std_cols[1]]]
  ),
  data.frame(
    ntaxa  = rf$ntaxa,
    method = "beam",
    mean   = rf$beam,
    sd     = rf[[std_cols[2]]]
  ),
  data.frame(
    ntaxa  = rf$ntaxa,
    method = "laml",
    mean   = rf$laml,
    sd     = rf[[std_cols[3]]]
  )
)

rf_long$method <- factor(rf_long$method, levels = c("vine", "beam", "laml"))

# Error bars, truncated at zero on the lower end
rf_long$ymin <- pmax(rf_long$mean - rf_long$sd, 0)
rf_long$ymax <- rf_long$mean + rf_long$sd

# Absolute y-limits: min = 0, max = max(mean + sd)
y_max <- max(rf_long$ymax, na.rm = TRUE)
ylim  <- c(0, y_max)

# Plot (absolute RF distance; y starts at 0)
pmf <- ggplot(rf_long, aes(x = factor(ntaxa), y = mean, fill = method)) +
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

## -------------------- Save --------------------
save_pdf(pmf, "rf_bars.pdf", width = 3, height = 3)
