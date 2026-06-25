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

panel_tag_theme <- theme(
  plot.tag = element_text(
    size = 12, face = "bold", family = "Helvetica"
  ),
  plot.tag.position = c(0.01, 0.99)
)

# ================================================================
# Dimensionality helper
# ================================================================
make_dim_panel <- function(dims_file, lnl_cols, time_cols) {
  dims <- read.table(file.path(data_dir, dims_file),
                     header = TRUE, check.names = FALSE)
  dims$d <- suppressWarnings(
    as.numeric(gsub(".*\\.D", "", as.character(dims[[1]])))
  )

  delta_lnl <- dims[[lnl_cols[1]]] - dims[[lnl_cols[2]]]
  sd_lnl    <- dims[[lnl_cols[3]]]
  time_mean <- dims[[time_cols[1]]]
  sd_time   <- dims[[time_cols[2]]]

  num_max <- max(abs(delta_lnl), na.rm = TRUE)
  den_max <- max(time_mean, na.rm = TRUE)
  scale_factor <- if (is.finite(num_max) && is.finite(den_max) &&
                      den_max > 0 && num_max > 0) {
    num_max / den_max
  } else {
    1
  }

  plotdata <- bind_rows(
    data.frame(d = dims$d, metric = "delta_lnl",
               mean_scaled = delta_lnl, sd_scaled = sd_lnl,
               stringsAsFactors = FALSE),
    data.frame(d = dims$d, metric = "time",
               mean_scaled = time_mean * scale_factor,
               sd_scaled = sd_time * scale_factor,
               stringsAsFactors = FALSE)
  )
  plotdata$metric <- factor(plotdata$metric,
                            levels = c("delta_lnl", "time"))

  ll_min_err <- suppressWarnings(
    min(delta_lnl - sd_lnl, 0, na.rm = TRUE)
  )
  ll_breaks <- pretty(c(ll_min_err, 0), n = 6)
  if (length(ll_breaks) > 1 && min(ll_breaks) > ll_min_err) {
    step <- ll_breaks[2] - ll_breaks[1]
    ll_breaks <- c(min(ll_breaks) - step, ll_breaks)
  }
  ymin <- min(ll_breaks, na.rm = TRUE)

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

  dodge <- position_dodge(width = 0.6)
  geom_palette <- c("delta_lnl" = "#F28E2B", "time" = "#59A14F")
  orange <- unname(geom_palette["delta_lnl"])
  green  <- unname(geom_palette["time"])

  ggplot(plotdata,
         aes(x = factor(d), y = mean_scaled, fill = metric)) +
    geom_col(width = 0.45, position = dodge) +
    geom_errorbar(
      aes(ymin = mean_scaled - sd_scaled,
          ymax = mean_scaled + sd_scaled),
      width = 0.15, position = dodge
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = geom_palette, guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(
      name = expression(Delta ~ lnl),
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
      axis.title.y       = element_text(color = orange),
      axis.text.y        = element_text(color = orange),
      axis.ticks.y       = element_line(color = orange),
      axis.title.y.right = element_text(color = green),
      axis.text.y.right  = element_text(color = green),
      axis.ticks.y.right = element_line(color = green)
    )
}

# ================================================================
# Dimensionality (standard): 25 taxa (A) + 50 taxa (B)
# Columns: 2=vine, 3=std, 4=beast, 5=time, 6=std
# ================================================================
pdim_25 <- make_dim_panel("dimSummary25.txt",
                          lnl_cols = c(2, 4, 3),
                          time_cols = c(5, 6))
pdim_50 <- make_dim_panel("dimSummary50.txt",
                          lnl_cols = c(2, 4, 3),
                          time_cols = c(5, 6))

pdim_panels <- (pdim_25 + pdim_50) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "", tag_suffix = "")
pdim_panels <- pdim_panels & panel_tag_theme

save_pdf(pdim_panels, file.path(out_dir, "dims_panels.pdf"),
         width = 6, height = 3)

# ================================================================
# Dimensionality H version: 25 taxa (A) + 50 taxa (B)
# Columns: 7=vineH, 8=std, 4=beast, 9=timeH, 10=std
# ================================================================
pdimH_25 <- make_dim_panel("dimSummary25.txt",
                           lnl_cols = c(7, 4, 8),
                           time_cols = c(9, 10))
pdimH_50 <- make_dim_panel("dimSummary50.txt",
                           lnl_cols = c(7, 4, 8),
                           time_cols = c(9, 10))

pdimH_panels <- (pdimH_25 + pdimH_50) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "", tag_suffix = "")
pdimH_panels <- pdimH_panels & panel_tag_theme

save_pdf(pdimH_panels, file.path(out_dir, "dimsH_panels.pdf"),
         width = 6, height = 3)
