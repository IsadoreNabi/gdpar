## Block RG, RG.1.c analysis: the cheap offline pass over the signals collected
## by run_rg_calibration.R. Reconstructs the classifier's decision input from
## each stored cell, sweeps the thresholds of gdpar_geometry_thresholds(), and
## reports MEASURED error rates (the honest replacement for the "8/8 diagonal of
## one replica" smoke of B9.22).
##
## Threshold calibration = the union chosen in B9.23 (Q3): a data-driven
## proposal (Youden-optimal cut per gating signal) is REGULARISED by a small
## interpretable grid via coordinate ascent on a CALIBRATION fold, and the final
## thresholds are REPORTED on a held-out VALIDATION fold (out-of-sample, immune
## to the "tuned and tested on the same diagonal" critique). The classifier
## stays the transparent rule-based one in R/geometry_diagnostic.R; only its
## numeric thresholds move.
##
## All of this is a pure function of the stored signals: no resampling, runs in
## seconds. Source under devtools::load_all (uses gdpar:::.gdpar_geom_classify).

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

## --- collect: stored RDS cells -> one classification-input row per cell ----
## Reproduces the aggregation gdpar_geometry_diagnostic uses to classify: the
## signals of the largest-n pilot, with multimodality maximised across pilots.
rg_collect <- function(cells_dir) {
  rds <- list.files(cells_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds) == 0L) stop("No RDS cells in ", cells_dir)
  rows <- vector("list", length(rds))
  for (i in seq_along(rds)) {
    obj <- readRDS(rds[[i]])
    s <- obj$signals
    if (is.null(s) || all(s$failed)) {
      agg <- as.list(stats::setNames(rep(NA_real_, 10L), c(
        "divergent_rate", "ebfmi_min", "treedepth_sat_rate", "condition_number",
        "step_scale_ratio", "lambda_max_cov", "mean_leapfrog", "multimodality",
        "heavy_kurtosis", "boundary_proximity")))
    } else {
      last <- s[nrow(s), , drop = FALSE]
      agg <- as.list(last[, c(
        "divergent_rate", "ebfmi_min", "treedepth_sat_rate", "condition_number",
        "step_scale_ratio", "lambda_max_cov", "mean_leapfrog", "multimodality",
        "heavy_kurtosis", "boundary_proximity")])
      agg$multimodality <- suppressWarnings(max(s$multimodality, na.rm = TRUE))
      if (!is.finite(agg$multimodality)) agg$multimodality <- NA_real_
    }
    rows[[i]] <- data.frame(
      cell_id = obj$cell_id, cell_base = obj$cell_base,
      target_id = obj$target_id, diff_level = obj$diff_level,
      pilot_level = obj$pilot_level, rep = obj$rep,
      ground_truth = obj$ground_truth, slope = obj$slope %||% NA_real_,
      status = obj$status %||% "ok",
      as.data.frame(agg, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$diff_level <- factor(out$diff_level, levels = c("low", "mid", "high"))
  out
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

## --- reclassify a row under a threshold set (cheap, pure) ------------------
## culprit does not enter the classification rules (it is reported separately),
## so an empty culprit is passed.
rg_predict <- function(tab, th) {
  empty <- gdpar:::.gdpar_geom_empty_culprit()
  vapply(seq_len(nrow(tab)), function(i) {
    r <- tab[i, , drop = FALSE]
    sig <- list(
      divergent_rate = r$divergent_rate, ebfmi_min = r$ebfmi_min,
      treedepth_sat_rate = r$treedepth_sat_rate,
      condition_number = r$condition_number,
      step_scale_ratio = r$step_scale_ratio, lambda_max_cov = r$lambda_max_cov,
      mean_leapfrog = r$mean_leapfrog, multimodality = r$multimodality,
      heavy_kurtosis = r$heavy_kurtosis, boundary_proximity = r$boundary_proximity)
    ncurve <- list(
      slope = r$slope,
      grows_with_n = is.finite(r$slope) && r$slope > th$nslope_grows)
    gdpar:::.gdpar_geom_classify(sig, ncurve, empty, th)$pathology
  }, character(1))
}

## --- metrics ---------------------------------------------------------------
rg_classes <- function() c(
  "isotropic", "anisotropic", "funnel", "heavy_tails", "quasi_deterministic",
  "multimodal", "boundary", "flat_direction")

## Wilson score interval for a binomial proportion (better than Wald near 0/1).
rg_wilson <- function(x, n, conf = 0.95) {
  if (n == 0L) return(c(lo = NA_real_, hi = NA_real_, half = NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- x / n
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  halfw <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d
  c(lo = centre - halfw, hi = centre + halfw, half = halfw)
}

## Balanced accuracy = mean over ground-truth classes of the per-class recall.
rg_balanced_accuracy <- function(truth, pred) {
  cls <- intersect(rg_classes(), unique(truth))
  recalls <- vapply(cls, function(c) {
    idx <- truth == c
    if (!any(idx)) return(NA_real_)
    mean(pred[idx] == c)
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

rg_confusion <- function(truth, pred) {
  cls <- rg_classes()
  table(factor(truth, levels = cls), factor(pred, levels = cls))
}

## Design-uniform (macro) accuracy: because the adaptive Stage-2 top-up gives
## hard cell_bases MORE replicas, pooling over replicas biases the raw accuracy
## toward the hard subset. The honest population estimate weights each
## cell_base (target x difficulty x pilot) equally: compute the per-cell_base
## correct rate, then average -- overall and per ground-truth class.
rg_macro_accuracy <- function(tab, th) {
  pred <- rg_predict(tab, th)
  ok <- pred == tab$ground_truth
  by_base_rate <- tapply(ok, tab$cell_base, mean)
  base_class <- tapply(tab$ground_truth, tab$cell_base, function(z) z[1L])
  by_class <- tapply(by_base_rate, base_class[names(by_base_rate)], mean)
  list(macro_accuracy = mean(by_base_rate),
       macro_by_class = round(by_class, 3))
}

## --- data-driven proposal: Youden-optimal cuts -----------------------------
## For a "high" gate (signal >= cut => positive class), the optimal cut
## maximises Youden's J = TPR - FPR over candidate midpoints.
rg_youden_cut <- function(pos, neg, side = c("ge", "le")) {
  side <- match.arg(side)
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  if (length(pos) == 0L || length(neg) == 0L) return(NA_real_)
  cand <- sort(unique(c(pos, neg)))
  cuts <- (utils::head(cand, -1L) + utils::tail(cand, -1L)) / 2
  if (length(cuts) == 0L) cuts <- cand
  j <- vapply(cuts, function(c) {
    if (side == "ge") {
      tpr <- mean(pos >= c); fpr <- mean(neg >= c)
    } else {
      tpr <- mean(pos <= c); fpr <- mean(neg <= c)
    }
    tpr - fpr
  }, numeric(1))
  cuts[which.max(j)]
}

## Map each gating threshold to (signal, positive-class set, side). These are
## the thresholds that move a decision in .gdpar_geom_classify.
rg_threshold_specs <- function() {
  list(
    boundary_prox_high  = list(sig = "boundary_proximity",
                               pos = "boundary", side = "ge"),
    multimodal_high     = list(sig = "multimodality",
                               pos = "multimodal", side = "ge"),
    flat_var_high       = list(sig = "lambda_max_cov",
                               pos = "flat_direction", side = "ge"),
    nslope_grows        = list(sig = "slope",
                               pos = "quasi_deterministic", side = "ge"),
    condition_high      = list(sig = "condition_number",
                               pos = c("anisotropic", "quasi_deterministic"),
                               side = "ge"),
    heavy_kurtosis_high = list(sig = "heavy_kurtosis",
                               pos = c("funnel", "heavy_tails"), side = "ge"),
    funnel_ebfmi_low    = list(sig = "ebfmi_min", pos = "funnel", side = "le"),
    heavy_cond_max      = list(sig = "condition_number",
                               pos = "heavy_tails", side = "le")
  )
}

rg_propose_datadriven <- function(tab, th0) {
  specs <- rg_threshold_specs()
  th <- th0
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    pos <- tab[[sp$sig]][tab$ground_truth %in% sp$pos]
    neg <- tab[[sp$sig]][!(tab$ground_truth %in% sp$pos)]
    cut <- rg_youden_cut(pos, neg, sp$side)
    if (is.finite(cut)) th[[nm]] <- cut
  }
  th
}

## --- grid regularisation: coordinate ascent on the calibration fold --------
## Candidate grid per threshold = the default, the data-driven proposal, and
## interpretable multiplicative steps around the default. Coordinate ascent
## keeps the value maximising calibration balanced accuracy; a few passes.
rg_threshold_grid <- function(nm, th_def, th_dd) {
  base <- th_def[[nm]]
  prop <- th_dd[[nm]]
  steps <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3) * base
  cand <- sort(unique(c(steps, prop, base)))
  ## slope thresholds may legitimately be small/zero; others must be positive.
  if (identical(nm, "nslope_grows")) cand[is.finite(cand)]
  else cand[is.finite(cand) & cand > 0]
}

rg_grid_regularize <- function(cal, th_def, th_dd, passes = 3L) {
  th <- th_dd
  specs <- names(rg_threshold_specs())
  best <- rg_balanced_accuracy(cal$ground_truth, rg_predict(cal, th))
  for (p in seq_len(passes)) {
    improved <- FALSE
    for (nm in specs) {
      grid <- rg_threshold_grid(nm, th_def, th_dd)
      for (v in grid) {
        th_try <- th; th_try[[nm]] <- v
        ba <- rg_balanced_accuracy(cal$ground_truth, rg_predict(cal, th_try))
        if (isTRUE(ba > best + 1e-9)) { best <- ba; th <- th_try; improved <- TRUE }
      }
    }
    if (!improved) break
  }
  list(thresholds = th, calibration_balanced_accuracy = best)
}

## --- held-out fold split (deterministic, by replica parity) ----------------
rg_fold_split <- function(tab) {
  cal <- tab[tab$rep %% 2L == 1L, , drop = FALSE]   # odd reps
  val <- tab[tab$rep %% 2L == 0L, , drop = FALSE]   # even reps
  list(cal = cal, val = val)
}

## --- adaptive-allocation helper: which cell_bases sit on the decision
## boundary (widest Wilson CI for the default-threshold correct rate) and so
## warrant Stage-2 top-up. Minimax-robust: spend replicas where p ~ 0.5. -----
rg_boundary_cells <- function(tab, th, lo = 0.2, hi = 0.8) {
  pred <- rg_predict(tab, th)
  ok <- pred == tab$ground_truth
  agg <- tapply(ok, tab$cell_base, function(z) mean(z))
  names(agg)[is.finite(agg) & agg >= lo & agg <= hi]
}

## --- top-level analysis -----------------------------------------------------
rg_analyze <- function(cells_dir) {
  th_def <- gdpar_geometry_thresholds()
  tab <- rg_collect(cells_dir)
  split <- rg_fold_split(tab)
  th_dd <- rg_propose_datadriven(split$cal, th_def)
  reg <- rg_grid_regularize(split$cal, th_def, th_dd)
  th_final <- reg$thresholds

  ## evaluate default and tuned thresholds, in-sample and held-out.
  eval_one <- function(data, th) {
    pred <- rg_predict(data, th)
    list(
      balanced_accuracy = rg_balanced_accuracy(data$ground_truth, pred),
      accuracy = mean(pred == data$ground_truth),
      confusion = rg_confusion(data$ground_truth, pred),
      pred = pred)
  }
  list(
    tab = tab, folds = split,
    thresholds_default = th_def,
    thresholds_datadriven = th_dd,
    thresholds_final = th_final,
    calibration_balanced_accuracy = reg$calibration_balanced_accuracy,
    eval_default_all      = eval_one(tab, th_def),
    eval_default_heldout  = eval_one(split$val, th_def),
    eval_final_heldout    = eval_one(split$val, th_final),
    eval_final_all        = eval_one(tab, th_final),
    macro_default         = rg_macro_accuracy(tab, th_def),
    macro_final           = rg_macro_accuracy(tab, th_final))
}

## --- report helpers --------------------------------------------------------
## Per-class recall with a Wilson CI, on a given data set under given thresholds.
rg_perclass_table <- function(data, th) {
  pred <- rg_predict(data, th)
  cls <- intersect(rg_classes(), unique(data$ground_truth))
  rows <- lapply(cls, function(c) {
    idx <- data$ground_truth == c
    n <- sum(idx); x <- sum(pred[idx] == c)
    w <- rg_wilson(x, n)
    data.frame(class = c, n = n, correct = x, recall = round(x / n, 3),
               ci_lo = round(w["lo"], 3), ci_hi = round(w["hi"], 3),
               row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## Accuracy by ground-truth class x difficulty level (the error-vs-severity
## curve).
rg_severity_table <- function(tab, th) {
  pred <- rg_predict(tab, th)
  ok <- pred == tab$ground_truth
  acc <- tapply(ok, list(tab$ground_truth, tab$diff_level), mean)
  round(acc, 2)
}

## Accuracy by pilot budget (the accuracy-vs-pilot-cost curve).
rg_cost_table <- function(tab, th) {
  pred <- rg_predict(tab, th)
  ok <- pred == tab$ground_truth
  data.frame(
    pilot_level = names(tapply(ok, tab$pilot_level, mean)),
    accuracy = round(as.numeric(tapply(ok, tab$pilot_level, mean)), 3),
    balanced = round(as.numeric(tapply(seq_along(ok), tab$pilot_level,
      function(i) rg_balanced_accuracy(tab$ground_truth[i], pred[i]))), 3),
    n = as.integer(tapply(ok, tab$pilot_level, length)),
    row.names = NULL, stringsAsFactors = FALSE)
}

## Markdown for a confusion matrix (rows = truth, cols = predicted).
rg_confusion_md <- function(cm) {
  cn <- colnames(cm); ab <- substr(cn, 1, 5)
  hdr <- paste0("| truth \\\\ pred | ", paste(ab, collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", length(ab) + 1L), collapse = "|"), "|")
  body <- vapply(seq_len(nrow(cm)), function(i) {
    paste0("| ", rownames(cm)[i], " | ",
           paste(cm[i, ], collapse = " | "), " |")
  }, character(1))
  c(hdr, sep, body)
}

rg_num_md <- function(x) paste(format(unlist(x), digits = 4), collapse = ", ")

rg_write_report <- function(an, cells_dir, path) {
  tab <- an$tab
  n_stage1 <- sum(tab$rep <= 5L); n_stage2 <- sum(tab$rep > 5L)
  thd <- an$thresholds_default; thf <- an$thresholds_final
  changed <- names(thd)[vapply(names(thd), function(k)
    !isTRUE(all.equal(thd[[k]], thf[[k]])), logical(1))]
  L <- c(
    "# Block RG / RG.1.c -- posterior-geometry diagnostic calibration",
    "",
    sprintf("Generated %s by `inst/benchmarks/rg_calibration_analysis.R` over the cells in `%s`.",
            format(Sys.Date()), basename(cells_dir)),
    "",
    "This report replaces the **8/8 diagonal of one replica** smoke of B9.22 with **measured error rates**: the `gdpar_geometry_diagnostic` classifier run over a `gdpar_geometry_suite` grid of three axes -- difficulty (severity ladder), pilot budget (cheap 100/100 vs default 150/150), and replicas -- with an adaptive, minimax-robust replica allocation and an out-of-sample (held-out) threshold calibration.",
    "",
    "## Methodology",
    "",
    "- **Signals are collected once** by cheap NUTS pilots through cmdstan; the threshold sweep is a pure function of the stored signals (no resampling). Size-invariant signals only -- divergences, E-BFMI, tree depth, condition number, step/scale -- **never R-hat/ESS on short runs** (lesson B9.20/B9.21).",
    "- **Adaptive allocation (Q2 mix 3+4):** Stage 1 gives every cell_base R0=5 replicas; Stage 2 tops up to R=25 only the cell_bases whose Wilson confidence interval for the correct rate is *unresolved* (overlaps 0.5, or rate 0). Minimax-robust: replicas are spent where the verdict is statistically undecided; settled cells keep R0.",
    "- **Threshold calibration (Q3 union):** a data-driven Youden-optimal cut per gating signal *proposes* thresholds; a coordinate-ascent grid around the interpretable defaults *regularises* them on a CALIBRATION fold (odd replicas); the final thresholds are *reported* on a held-out VALIDATION fold (even replicas). The classifier stays the transparent rule-based one; only its numeric thresholds move.",
    "- **Honesty (no overreach):** these are measured rates, not a curated diagonal. gdpar's diagnostic *can* and *does* fail targets -- mild pathologies look benign, and some pathologies are mutually confusable. Reported below as-is.",
    "",
    sprintf("## Coverage: %d cells (%d Stage 1 + %d Stage 2), %d cell_bases.",
            nrow(tab), n_stage1, n_stage2, length(unique(tab$cell_base))),
    "",
    "## Headline accuracy",
    "",
    "**Design-uniform (macro) accuracy** is the honest population estimate: each cell_base (target x difficulty x pilot) is weighted equally, so the adaptive Stage-2 top-up (which gives hard cells more replicas) does not bias it. The replica-pooled figures are the tuning objective and are reported alongside but are pessimistic on the adaptive sample (they over-weight the hard, over-sampled cells).",
    "",
    "| thresholds | macro accuracy (design-uniform) | pooled accuracy | pooled balanced |",
    "|---|---|---|---|",
    sprintf("| default | %.3f | %.3f | %.3f |", an$macro_default$macro_accuracy,
            an$eval_default_all$accuracy, an$eval_default_all$balanced_accuracy),
    sprintf("| tuned | %.3f | %.3f | %.3f |", an$macro_final$macro_accuracy,
            an$eval_final_all$accuracy, an$eval_final_all$balanced_accuracy),
    "",
    "Out-of-sample (held-out validation fold) balanced accuracy -- the overfitting guard:",
    "",
    "| thresholds | held-out accuracy | held-out balanced |",
    "|---|---|---|",
    sprintf("| default | %.3f | %.3f |", an$eval_default_heldout$accuracy,
            an$eval_default_heldout$balanced_accuracy),
    sprintf("| tuned | %.3f | %.3f |", an$eval_final_heldout$accuracy,
            an$eval_final_heldout$balanced_accuracy),
    "",
    sprintf("Calibration-fold balanced accuracy at the tuned thresholds: %.3f.",
            an$calibration_balanced_accuracy),
    "",
    "Design-uniform per-class accuracy (tuned):",
    "",
    paste0("`", paste(sprintf("%s=%.2f", names(an$macro_final$macro_by_class),
                              an$macro_final$macro_by_class), collapse = "  "),
           "`"),
    "",
    "## Held-out confusion matrix (tuned thresholds)",
    "",
    rg_confusion_md(an$eval_final_heldout$confusion),
    "",
    "Class abbreviations are the first five letters (isotr, aniso, funne, heavy, quasi, multi, bound, flat_).",
    "",
    "## Held-out per-class recall (tuned thresholds), Wilson 95% CI",
    "",
    "| class | n | correct | recall | ci_lo | ci_hi |",
    "|---|---|---|---|---|---|")
  pct <- rg_perclass_table(an$folds$val, thf)
  for (i in seq_len(nrow(pct))) {
    L <- c(L, sprintf("| %s | %d | %d | %.3f | %.3f | %.3f |",
                      pct$class[i], pct$n[i], pct$correct[i], pct$recall[i],
                      pct$ci_lo[i], pct$ci_hi[i]))
  }
  sev <- rg_severity_table(tab, thf)
  L <- c(L, "",
         "## Error-vs-severity (accuracy by ground-truth class x difficulty)",
         "", "Difficulty rank low/mid/high = mild/default/severe pathology intensity (mid ~ the suite default).",
         "",
         paste0("| class | ", paste(colnames(sev), collapse = " | "), " |"),
         paste0("|", paste(rep("---", ncol(sev) + 1L), collapse = "|"), "|"))
  for (i in seq_len(nrow(sev))) {
    L <- c(L, paste0("| ", rownames(sev)[i], " | ",
                     paste(ifelse(is.na(sev[i, ]), "-", sev[i, ]),
                           collapse = " | "), " |"))
  }
  cst <- rg_cost_table(tab, thf)
  L <- c(L, "",
         "## Accuracy-vs-pilot-cost (tuned thresholds)",
         "",
         "| pilot budget | accuracy | balanced | n |",
         "|---|---|---|---|")
  for (i in seq_len(nrow(cst))) {
    L <- c(L, sprintf("| %s | %.3f | %.3f | %d |", cst$pilot_level[i],
                      cst$accuracy[i], cst$balanced[i], cst$n[i]))
  }
  L <- c(L, "",
         "## Thresholds: default vs tuned",
         "",
         "| threshold | default | tuned | changed |",
         "|---|---|---|---|")
  for (k in names(thd)) {
    L <- c(L, sprintf("| %s | %s | %s | %s |", k,
                      format(thd[[k]], digits = 4), format(thf[[k]], digits = 4),
                      if (k %in% changed) "yes" else ""))
  }
  L <- c(L, "",
         "## Anchors (no overreach)",
         "",
         "Girolami & Calderhead (2011, RMHMC); Betancourt (SoftAbs, E-BFMI, pathology taxonomy); Neal (the funnel); Lu et al. (relativistic Monte Carlo); Pennec / Bhattacharya-Patrangenaru (Riemannian statistics, the user's `ON_RIEMMANIAN_STATISTICS`). gdpar integrates and orchestrates these with its structural knowledge; it does not claim to invent the methods, nor that its classifier is unique or infallible.",
         "")
  writeLines(L, path)
  invisible(path)
}

