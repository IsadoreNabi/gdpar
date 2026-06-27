## report_hit_rate_multi.R
##
## Read inst/benchmarks/results/cp_ncp_hit_rate_multi.csv, print
## hit-rate summary tables (per component and per scenario), and write
## a faceted ggplot to inst/benchmarks/results/cp_ncp_hit_rate_multi.png
## with one panel per (scenario, component) showing hit/miss per coord
## and a side panel with n_div_pred vs n_div_alt per scenario.
##
## Usage (from the package root):
##   Rscript inst/benchmarks/scripts/report_hit_rate_multi.R
##
## Optional argument: the path to a custom CSV. If omitted, the
## canonical inst/benchmarks/results/cp_ncp_hit_rate_multi.csv is used.

suppressPackageStartupMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
})

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1L) args[1L]
            else "inst/benchmarks/results/cp_ncp_hit_rate_multi.csv"
if (!file.exists(csv_path)) {
  stop(sprintf("CSV not found: %s\nRun calibrate_cp_ncp_multi.R first.",
               csv_path))
}

results <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

cat("=== Hit rate per component (truths only) ===\n")
sub <- results[!is.na(results$regime_truth) & !is.na(results$hit), ]
if (nrow(sub) > 0L) {
  hits_comp <- aggregate(hit ~ component, data = sub,
                          FUN = function(x) mean(as.logical(x)))
  hits_comp$hit <- round(hits_comp$hit, 3)
  print(hits_comp, row.names = FALSE)
} else {
  cat("(no rows with truth defined)\n")
}

cat("\n=== Hit rate per scenario (truths only) ===\n")
if (nrow(sub) > 0L) {
  hits_sc <- aggregate(hit ~ scenario, data = sub,
                        FUN = function(x) mean(as.logical(x)))
  hits_sc$hit <- round(hits_sc$hit, 3)
  print(hits_sc, row.names = FALSE)
}

cat("\n=== Borderline (no truth) regime predictions ===\n")
bd <- results[is.na(results$regime_truth), ]
if (nrow(bd) > 0L) {
  print(bd[, c("scenario", "k", "component", "regime_pred",
                "decision_reason", "t_info_cp", "t_info_ncp")],
        row.names = FALSE)
} else {
  cat("(none)\n")
}

cat("\n=== n_div pred vs alt (per scenario, deduplicated) ===\n")
div_tbl <- unique(results[, c("scenario", "p", "n_div_pred",
                               "n_div_alt", "ebfmi_min")])
print(div_tbl, row.names = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\nNOTE: ggplot2 not installed; skipping plot generation.\n")
  cat("Install with: install.packages('ggplot2') (in an interactive R session).\n")
  quit(save = "no", status = 0L)
}

out_dir <- dirname(csv_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
plot_path <- file.path(out_dir, "cp_ncp_hit_rate_multi.png")

results$hit_label <- ifelse(
  is.na(results$hit), "borderline",
  ifelse(as.logical(results$hit), "hit", "miss")
)
results$hit_label <- factor(results$hit_label,
                             levels = c("hit", "miss", "borderline"))
results$component <- factor(results$component, levels = c("a", "W"))

p_hits <- ggplot2::ggplot(
  results,
  ggplot2::aes(x = factor(k), fill = hit_label)
) +
  ggplot2::geom_bar(width = 0.7) +
  ggplot2::facet_grid(scenario ~ component, scales = "free_x",
                       space = "free_x") +
  ggplot2::scale_fill_manual(values = c(
    "hit" = "#2c7bb6", "miss" = "#d7191c",
    "borderline" = "#fdae61"
  )) +
  ggplot2::labs(
    title = "preflight_parametrization_multi: per-(scenario, k, component) hit/miss",
    subtitle = sprintf("Read from %s", basename(csv_path)),
    x = "k (coordinate)", y = "count", fill = "verdict"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))

ggplot2::ggsave(plot_path, p_hits, width = 8, height = 10, dpi = 110)
cat(sprintf("\nWrote plot: %s\n", plot_path))

div_long <- rbind(
  data.frame(scenario = div_tbl$scenario, p = div_tbl$p,
             which = "n_div_pred", count = div_tbl$n_div_pred,
             stringsAsFactors = FALSE),
  data.frame(scenario = div_tbl$scenario, p = div_tbl$p,
             which = "n_div_alt", count = div_tbl$n_div_alt,
             stringsAsFactors = FALSE)
)
p_div <- ggplot2::ggplot(
  div_long,
  ggplot2::aes(x = scenario, y = count, fill = which)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                     width = 0.6) +
  ggplot2::scale_fill_manual(values = c(
    "n_div_pred" = "#2c7bb6", "n_div_alt" = "#d7191c"
  )) +
  ggplot2::labs(
    title = "Contrastive divergences: predicted vs alternative parametrization",
    subtitle = "Lower for predicted means the preflight chose the better regime",
    x = "scenario", y = "n_divergent (over 400 transitions)",
    fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

ggplot2::ggsave(file.path(out_dir, "cp_ncp_div_pred_vs_alt.png"),
                 p_div, width = 8, height = 4, dpi = 110)
cat(sprintf("Wrote plot: %s\n",
            file.path(out_dir, "cp_ncp_div_pred_vs_alt.png")))

cat("Done.\n")
