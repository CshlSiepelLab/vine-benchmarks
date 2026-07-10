#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()

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

load_summary <- function(filename) {
  d <- read.table(file.path(script_dir, filename), header = TRUE)
  std_cols <- grep("^std", names(d))
  if (length(std_cols) < 3) stop("Expected three std columns in ", filename)
  out <- rbind(
    data.frame(ntaxa=d$ntaxa, method="vine", mean=d$vine, sd=d[[std_cols[1]]]),
    data.frame(ntaxa=d$ntaxa, method="laml", mean=d$laml, sd=d[[std_cols[2]]]),
    data.frame(ntaxa=d$ntaxa, method="beam", mean=d$beam, sd=d[[std_cols[3]]])
  )
  out$method <- factor(out$method, levels=c("vine", "beam", "laml"))
  out$ymin <- pmax(out$mean - out$sd, 0)
  out$ymax <- out$mean + out$sd
  out
}

make_bar <- function(long, ylab, tag) {
  methods <- levels(long$method)
  ggplot(long, aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  labs(
    title = tag,
    x = "Number of Taxa",
    y = ylab,
    fill = "Method"
  ) +
  scale_fill_manual(
    values = method_palette[methods],
    breaks = methods,
    labels = label_map(methods)
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6))) +
  scale_y_continuous(breaks = pretty_breaks(n = 8), limits = c(0, NA))
}

## -------------------- Save --------------------
prf <- make_bar(load_summary("rfSummary.txt"), "Normalized RF distance", "A")
pbsd <- make_bar(load_summary("bsdSummary.txt"),
                 "Normalized branch-score distance", "B")
fig <- (prf + pbsd) + plot_layout(ncol = 2, guides = "collect")
fig <- fig & theme(legend.position = "right",
                   plot.title = element_text(size = 12, face = "bold", hjust = 0))
output_file <- file.path(script_dir, sprintf("rf_bars_%s.pdf", format(Sys.Date(), "%m%d%y")))
save_pdf(fig, output_file, width = 9, height = 3)
cat("wrote", output_file, "\n")
