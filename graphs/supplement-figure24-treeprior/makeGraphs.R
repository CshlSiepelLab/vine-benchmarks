#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(patchwork))

# Resolve script directory so we can build absolute paths
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}

graphs_dir <- dirname(script_dir)
base_dir <- file.path(graphs_dir, "hky300-data")
tp_dir   <- file.path(graphs_dir, "hky300-treeprior-data")

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
  NJ           = "#4E79A7",
  vine_noprior = "#FED8B1",
  vine_prior   = "#F28E2B",
  beast        = "#59A14F",
  mrbayes      = "#E15759"
)

method_labels <- c(
  NJ           = "NJ",
  vine_noprior = "Vine (No Prior)",
  vine_prior   = "Vine (Tree Prior)",
  beast        = "BEAST2",
  mrbayes      = "MrBayes"
)

label_map <- function(keys) {
  out <- method_labels[keys]
  out[is.na(out)] <- keys[is.na(out)]
  out
}

methods_keep <- c("NJ", "vine_noprior", "vine_prior", "beast", "mrbayes")

# ------------ Utilities ------------
melt_mean_sd <- function(df, methods, ntaxa_col = "ntaxa",
                         scale_by = NULL) {
  out <- list()
  for (m in methods) {
    if (!(m %in% names(df))) next
    m_idx  <- which(names(df) == m)
    sd_idx <- if (m_idx < ncol(df) &&
                  identical(names(df)[m_idx + 1], "std")) {
      m_idx + 1
    } else {
      NA_integer_
    }
    mean_vals <- df[[m]]
    sd_vals   <- if (!is.na(sd_idx)) df[[sd_idx]]
                 else rep(NA_real_, length(mean_vals))
    if (!is.null(scale_by)) {
      denom     <- abs(df[[scale_by]])
      mean_vals <- mean_vals / denom
      sd_vals   <- sd_vals / denom
    }
    out[[m]] <- data.frame(
      ntaxa = df[[ntaxa_col]],
      method = m,
      mean   = mean_vals,
      sd     = sd_vals,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# ================================================================
# Read data
# ================================================================
max_taxa <- 100

lnl_base <- read.table(file.path(base_dir, "lnlSummary.txt"),
                       header = TRUE, check.names = FALSE)
lnl_base <- lnl_base[lnl_base$ntaxa <= max_taxa, ]
lnl_tp   <- read.table(file.path(tp_dir, "lnlSummary.txt"),
                       header = TRUE, check.names = FALSE)
lnl_tp   <- lnl_tp[lnl_tp$ntaxa <= max_taxa, ]

rf_base  <- read.table(file.path(base_dir, "rfSummary.txt"),
                       header = TRUE, check.names = FALSE)
rf_base  <- rf_base[rf_base$ntaxa <= max_taxa, ]
rf_tp    <- read.table(file.path(tp_dir, "rfSummary.txt"),
                       header = TRUE, check.names = FALSE)
rf_tp    <- rf_tp[rf_tp$ntaxa <= max_taxa, ]

no_prior_methods <- c("NJ", "beast", "mrbayes")

# ================================================================
# A) Log-likelihood deviations
# ================================================================
lnl_long_np <- melt_mean_sd(lnl_base, methods = no_prior_methods,
                            ntaxa_col = "ntaxa", scale_by = "ave")

lnl_long_vine_noprior <- melt_mean_sd(lnl_base, methods = "vine",
                                      ntaxa_col = "ntaxa", scale_by = "ave")
lnl_long_vine_noprior$method <- "vine_noprior"

lnl_long_vine_prior <- melt_mean_sd(lnl_tp, methods = "vine-prior",
                                    ntaxa_col = "ntaxa", scale_by = "ave")
lnl_long_vine_prior$method <- "vine_prior"

lnl_long <- rbind(lnl_long_np, lnl_long_vine_noprior, lnl_long_vine_prior)
lnl_long$method <- factor(lnl_long$method, levels = methods_keep)

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
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Number of Taxa", y = expression(Delta ~ "Lnl"),
       fill = "Method") +
  scale_fill_manual(
    values = method_palette[methods_keep],
    breaks = methods_keep,
    labels = label_map(methods_keep)
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 0.1),
    breaks = pretty_breaks(n = 8)
  ) +
  guides(fill = guide_legend(override.aes = list(width = 0.6))) +
  coord_cartesian(ylim = ylim_lnl)

# ================================================================
# B) Robinson-Foulds distances
# ================================================================
rf_long_np <- melt_mean_sd(rf_base, methods = no_prior_methods,
                           ntaxa_col = "ntaxa")

rf_long_vine_noprior <- melt_mean_sd(rf_base, methods = "vine",
                                     ntaxa_col = "ntaxa")
rf_long_vine_noprior$method <- "vine_noprior"

rf_long_vine_prior <- melt_mean_sd(rf_tp, methods = "vine-prior",
                                   ntaxa_col = "ntaxa")
rf_long_vine_prior$method <- "vine_prior"

rf_long <- rbind(rf_long_np, rf_long_vine_noprior, rf_long_vine_prior)
rf_long$method <- factor(rf_long$method, levels = methods_keep)

rf_long$ymin <- pmax(rf_long$mean - rf_long$sd, 0)
rf_long$ymax <- rf_long$mean + rf_long$sd
y_max_rf <- max(rf_long$ymax, na.rm = TRUE)

prf <- ggplot(rf_long,
              aes(x = factor(ntaxa), y = mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  labs(x = "Number of Taxa", y = "Normalized RF Distance",
       fill = "Method") +
  scale_fill_manual(
    values = method_palette[methods_keep],
    breaks = methods_keep,
    labels = label_map(methods_keep)
  ) +
  scale_y_continuous(breaks = pretty_breaks(n = 8)) +
  guides(fill = guide_legend(override.aes = list(width = 0.6))) +
  coord_cartesian(ylim = c(0, y_max_rf))

# ================================================================
# Combine and save
# ================================================================
pcombo <- (plnl | prf) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "", tag_suffix = "")
pcombo <- pcombo & theme(
  legend.position = "right",
  plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
  plot.tag.position = c(0.01, 0.99)
)

save_pdf(pcombo,
         file.path(script_dir, "treeprior_lnl_rf.pdf"),
         width = 9, height = 3.5)
