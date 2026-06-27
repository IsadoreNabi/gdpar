## Block 6 reporter: reads the long-tidy CSV emitted by
## run_synthetic_hard_multi.R and produces three artifacts in
## inst/benchmarks/results/:
##   - synthetic_hard_multi_plot.png:    facet_wrap(metric ~ scenario)
##                                         bar plot, one bar per method,
##                                         with error bars for elpd_loo.
##   - synthetic_hard_multi_wall.png:    wall time per method per scenario.
##   - synthetic_hard_multi_pareto.png:  max Pareto-k per method per
##                                         scenario (PSIS-LOO health).
##   - synthetic_hard_multi_summary.md:  markdown verdict table per
##                                         scenario summarising best
##                                         elpd_loo competitor + caveats.
##
## Usage:
##   Rscript inst/benchmarks/scripts/report_synthetic_hard_multi.R
##
## Reads from:
##   inst/benchmarks/results/synthetic_hard_multi_results.csv

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())) # nolint

repo_root <- normalizePath(".")
results_dir <- file.path(repo_root, "inst", "benchmarks", "results")
csv_path <- file.path(results_dir, "synthetic_hard_multi_results.csv")
if (!file.exists(csv_path)) {
  stop(sprintf(
    "Results CSV not found at %s. Run inst/benchmarks/run_synthetic_hard_multi.R first.",
    csv_path
  ))
}
results <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for plotting.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required for summaries.")
}

main_metrics <- c("rmse", "mae", "crps", "auc", "log_loss", "brier",
                    "elpd_loo")
plot_df <- subset(results, results$metric %in% main_metrics)

ordered_methods <- c("gdpar", "brms", "brms_rescor", "mgcv", "inla",
                       "rstanarm")
plot_df$method <- factor(plot_df$method,
                           levels = intersect(ordered_methods,
                                                unique(plot_df$method)))
plot_df$scenario_id <- factor(plot_df$scenario_id,
                                levels = sort(unique(plot_df$scenario_id)))
plot_df$metric <- factor(plot_df$metric, levels = main_metrics)

p_main <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = .data$method, y = .data$value, fill = .data$method)
) +
  ggplot2::geom_col(
    ggplot2::aes(group = .data$coord),
    position = ggplot2::position_dodge(width = 0.8), width = 0.7
  ) +
  ggplot2::geom_errorbar(
    data = subset(plot_df, !is.na(plot_df$se_value)),
    ggplot2::aes(ymin = .data$value - 2 * .data$se_value,
                  ymax = .data$value + 2 * .data$se_value,
                  group = .data$coord),
    position = ggplot2::position_dodge(width = 0.8), width = 0.2
  ) +
  ggplot2::facet_wrap(~ metric + scenario_id, scales = "free_y",
                       ncol = length(unique(plot_df$scenario_id))) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  ggplot2::labs(x = NULL, y = NULL,
                  title = "Block 6 synthetic adversarial multi: metrics",
                  subtitle = paste("Per-coord and aggregated values;",
                                    "error bars on elpd_loo are 2 * se."))
ggplot2::ggsave(file.path(results_dir, "synthetic_hard_multi_plot.png"),
                  p_main, width = 16, height = 12, dpi = 110)

wall_df <- subset(results, results$metric == "wall")
wall_df$method <- factor(wall_df$method,
                           levels = intersect(ordered_methods,
                                                unique(wall_df$method)))
p_wall <- ggplot2::ggplot(
  wall_df,
  ggplot2::aes(x = .data$method, y = .data$value, fill = .data$method)
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ scenario_id) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                          hjust = 1)) +
  ggplot2::labs(x = NULL, y = "seconds",
                  title = "Block 6 wall time per scenario per method")
ggplot2::ggsave(file.path(results_dir, "synthetic_hard_multi_wall.png"),
                  p_wall, width = 14, height = 10, dpi = 110)

pareto_df <- subset(results,
                      results$metric == "elpd_loo" &
                        !is.na(results$max_pareto_k))
p_pareto <- ggplot2::ggplot(
  pareto_df,
  ggplot2::aes(x = .data$method, y = .data$max_pareto_k,
                fill = .data$method)
) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed") +
  ggplot2::facet_wrap(~ scenario_id) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                          hjust = 1)) +
  ggplot2::labs(x = NULL, y = "max Pareto-k",
                  title = paste("Block 6 PSIS-LOO health:",
                                  "max Pareto-k per scenario per method"))
ggplot2::ggsave(file.path(results_dir, "synthetic_hard_multi_pareto.png"),
                  p_pareto, width = 14, height = 10, dpi = 110)

summary_df <- subset(results,
                       results$metric == "elpd_loo" &
                         !is.na(results$value))
summary_df <- summary_df[order(summary_df$scenario_id,
                                  -summary_df$value), ]

md_lines <- c(
  "# Block 6 — synthetic adversarial multi: verdict summary",
  "",
  paste("Generated at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "| scenario | best (elpd_loo) | runner-up | gdpar rank | caveats |",
  "|----------|------------------|-----------|------------|---------|"
)
for (sc in unique(summary_df$scenario_id)) {
  rows <- summary_df[summary_df$scenario_id == sc, ]
  ordered_meth <- rows$method
  best <- ordered_meth[1L]
  runner <- if (length(ordered_meth) > 1L) ordered_meth[2L] else "-"
  gdpar_rank <- match("gdpar", ordered_meth)
  high_pk <- rows[rows$n_pareto_k_above_07 > 0L &
                     !is.na(rows$n_pareto_k_above_07), ]
  caveats <- if (nrow(high_pk) > 0L) {
    paste0("Pareto-k > 0.7 in: ",
            paste(unique(high_pk$method), collapse = ", "))
  } else "-"
  md_lines <- c(md_lines,
                  sprintf("| %s | %s (%.1f) | %s | %s | %s |",
                          sc, best, rows$value[1L],
                          as.character(runner),
                          as.character(gdpar_rank),
                          caveats))
}
md_path <- file.path(results_dir, "synthetic_hard_multi_summary.md")
writeLines(md_lines, md_path)

message(sprintf("Wrote:\n  - %s\n  - %s\n  - %s\n  - %s",
                  file.path(results_dir, "synthetic_hard_multi_plot.png"),
                  file.path(results_dir, "synthetic_hard_multi_wall.png"),
                  file.path(results_dir, "synthetic_hard_multi_pareto.png"),
                  md_path))
