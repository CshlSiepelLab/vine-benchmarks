#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
root_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
out_dir <- file.path(root_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxa <- c(10, 25, 50, 100, 250, 500, 1000)
all_method_levels <- c(
  "Vine", "Vine + flows", "BEAST2", "BEAST2 + BEAGLE",
  "MrBayes", "MrBayes + BEAGLE"
)
core_method_levels <- c("Vine", "Vine + flows", "BEAST2", "MrBayes")
method_colors <- c(
  "Vine" = "#F28E2B",
  "Vine + flows" = "#76B7B2",
  "BEAST2" = "#59A14F",
  "BEAST2 + BEAGLE" = "#8BC184",
  "MrBayes" = "#E15759",
  "MrBayes + BEAGLE" = "#E98A8C"
)

read_replicates <- function(path) {
  lines <- readLines(path, warn = FALSE)
  separator <- grep("^-+$", trimws(lines))
  if (length(separator)) lines <- lines[seq_len(separator[1] - 1)]
  read.table(text = lines, header = TRUE, check.names = FALSE)
}

mean_sd <- function(x) {
  x <- x[is.finite(x)]
  c(mean = mean(x), sd = if (length(x) > 1) sd(x) else 0)
}

read_eval_sd <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hit <- grep("^Std:", lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(strsplit(trimws(hit[1]), "[[:space:]]+")[[1]][2])
}

metric_spec <- list(
  time = list(file = "time", key = "plain"),
  mf = list(file = "mf", key = "percent_true"),
  rf = list(file = "rf", key = "plain"),
  bsd = list(file = "bsd", key = "plain"),
  entropy = list(file = "ent", key = "entropy")
)

metric_values <- function(d, method_key, transform) {
  column <- if (transform == "entropy") paste0(method_key, "_top") else method_key
  values <- d[[column]]
  if (transform == "percent_true") {
    values <- 100 * (values - d$true) / abs(d$true)
  }
  values
}

collect_metric <- function(metric, include_beagle = FALSE) {
  spec <- metric_spec[[metric]]
  methods <- if (include_beagle) all_method_levels else core_method_levels
  rows <- lapply(taxa, function(n) {
    vine <- read_replicates(file.path(
      root_dir, paste0(n, "taxa"), paste0("vine.eval.all.", spec$file, ".txt")
    ))
    flows <- read_replicates(file.path(
      root_dir, paste0(n, "taxa"), paste0("vine_flows.eval.all.", spec$file, ".txt")
    ))
    sources <- list(
      "Vine" = list(data = vine, key = "vine"),
      "Vine + flows" = list(data = flows, key = "vine"),
      "BEAST2" = list(data = vine, key = "beast"),
      "BEAST2 + BEAGLE" = list(data = vine, key = "beast-beagle"),
      "MrBayes" = list(data = vine, key = "mrbayes"),
      "MrBayes + BEAGLE" = list(data = vine, key = "mrbayes-beagle")
    )
    do.call(rbind, lapply(methods, function(method) {
      source <- sources[[method]]
      z <- mean_sd(metric_values(source$data, source$key, spec$key))
      data.frame(ntaxa = n, method = method,
                 mean = z[["mean"]], sd = z[["sd"]])
    }))
  })
  d <- do.call(rbind, rows)
  d$ntaxa <- factor(d$ntaxa, levels = taxa)
  d$method <- factor(d$method, levels = all_method_levels)
  d
}

# Match the model-fit processing used by graphs/figure2-hky300: average the
# raw likelihoods over simulation replicates before normalizing by the mean
# true-tree likelihood. Error bars summarize within-posterior likelihood
# variation (RMS posterior SD), rather than between-replicate variation.
collect_model_fit <- function() {
  base_dir <- file.path(dirname(root_dir), "hky_300sites")
  rows <- lapply(taxa, function(n) {
    vine <- read_replicates(file.path(
      root_dir, paste0(n, "taxa"), "vine.eval.all.mf.txt"
    ))
    flows <- read_replicates(file.path(
      root_dir, paste0(n, "taxa"), "vine_flows.eval.all.mf.txt"
    ))
    true_mean <- mean(vine$true)
    sources <- list(
      "Vine" = list(data = vine, key = "vine", file_key = "vine",
                    local = TRUE),
      "Vine + flows" = list(data = flows, key = "vine",
                            file_key = "vine_flows", local = TRUE),
      "BEAST2" = list(data = vine, key = "beast", file_key = "beast",
                      local = FALSE),
      "MrBayes" = list(data = vine, key = "mrbayes", file_key = "mrbayes",
                       local = FALSE)
    )
    do.call(rbind, lapply(core_method_levels, function(method) {
      source <- sources[[method]]
      eval_dir <- if (source$local) root_dir else base_dir
      suffix <- if (source$local) ".var.mf.txt" else ".mf.txt"
      posterior_sd <- vapply(seq_len(nrow(source$data)), function(i) {
        read_eval_sd(file.path(
          eval_dir, paste0(n, "taxa"),
          paste0("tree.", i, ".", source$file_key, suffix)
        ))
      }, numeric(1))
      posterior_sd <- posterior_sd[is.finite(posterior_sd)]
      data.frame(
        ntaxa = n,
        method = method,
        mean = 100 * (mean(source$data[[source$key]]) - true_mean) /
          abs(true_mean),
        sd = 100 * sqrt(mean(posterior_sd^2)) / abs(true_mean)
      )
    }))
  })
  d <- do.call(rbind, rows)
  d$ntaxa <- factor(d$ntaxa, levels = taxa)
  d$method <- factor(d$method, levels = all_method_levels)
  d
}

paper_theme <- theme_minimal(base_size = 8) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.spacing.y = grid::unit(1.5, "pt"),
    panel.grid.major = element_line(color = "gray60", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray80", linewidth = 0.2),
    plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
    plot.tag.position = c(-0.02, 1.02),
    plot.margin = margin(8, 10, 5, 10)
  )

metric_plot <- function(d, y_label, log_y = FALSE, zero_floor = TRUE) {
  p <- ggplot(d, aes(x = ntaxa, y = mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8,
             show.legend = TRUE) +
    geom_errorbar(
      aes(ymin = if (log_y) pmax(mean - sd, .Machine$double.eps)
          else if (zero_floor) pmax(0, mean - sd) else mean - sd,
          ymax = mean + sd),
      position = position_dodge(width = 0.9), width = 0.18, linewidth = 0.35,
      show.legend = TRUE
    ) +
    scale_fill_manual(values = method_colors, breaks = all_method_levels,
                      drop = FALSE) +
    labs(x = "Number of Taxa (n)", y = y_label, fill = "Method") +
    guides(fill = guide_legend(
      ncol = 1, byrow = TRUE,
      override.aes = list(linewidth = 0, linetype = 0)
    )) +
    paper_theme
  if (log_y) {
    p + scale_y_log10(
      breaks = 10^(0:5),
      labels = function(x) format(x, scientific = FALSE, trim = TRUE)
    )
  } else {
    p + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  }
}

panel_a <- metric_plot(
  collect_metric("time", include_beagle = TRUE), "Runtime (seconds)", log_y = TRUE
) + labs(tag = "A") +
  theme(axis.title.y = element_text(margin = margin(r = -2)))
panel_b <- metric_plot(
  collect_model_fit(), "% change from true\nheld-out log likelihood",
  zero_floor = FALSE
) + labs(tag = "B")
panel_c <- metric_plot(
  collect_metric("rf"), "Normalized Robinson-Foulds\ndistance to true tree"
) + labs(tag = "C") + coord_cartesian(ylim = c(0, 0.8)) +
  theme(axis.title.y = element_text(margin = margin(r = -2)))
panel_d <- metric_plot(
  collect_metric("bsd"), "Normalized branch-score\ndistance to true tree"
) + labs(tag = "D")

# Use the same native collected-legend layout as Supplementary Figure 1 so
# patchwork allocates only the width actually needed by the legend.
figure <- ((panel_a + panel_b) / (panel_c + panel_d)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical")

out_file <- file.path(out_dir, "hky300_posterior_quality_sup2.pdf")
ggsave(out_file, figure, width = 8.5, height = 6.5, units = "in",
       device = cairo_pdf)
message("Wrote ", out_file)
