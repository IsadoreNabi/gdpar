# Tests for the S3 class gdpar_preflight_report (Phase H.0).
# Pure tests: no Stan, no cmdstanr, no fit. All defensive paths
# covered: shape, NA-resilience, type-stability, empty edges,
# aggregation invariants, accessor contracts, S3 method outputs.

# ---- helpers ----------------------------------------------------

make_per_dim_full <- function(p = 2L,
                              decisions_a = NULL,
                              decisions_W = NULL) {
  if (is.null(decisions_a)) decisions_a <- rep("NCP", p)
  if (is.null(decisions_W)) decisions_W <- rep("NCP", p)
  stopifnot(length(decisions_a) == p, length(decisions_W) == p)
  data.frame(
    component       = c(rep("a", p), rep("W", p)),
    dim             = c(seq_len(p), seq_len(p)),
    decision        = c(decisions_a, decisions_W),
    decision_reason = c(rep("filter_info_high", p),
                        rep("filter_info_low",  p)),
    n_divergent     = rep(2L, 2L * p),
    div_pct         = rep(0.005, 2L * p),
    ebfmi_min       = rep(0.85, 2L * p),
    t_attribution   = rep(NA_real_, 2L * p),
    t_info_cp       = c(rep(3.5, p), rep(-1.2, p)),
    t_info_ncp      = c(rep(7.1, p), rep(-3.4, p)),
    stringsAsFactors = FALSE
  )
}


# ---- aggregator unit tests --------------------------------------

test_that("aggregate_preflight_one_component: any_ncp returns CP only when all CP", {
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("CP", "CP", "CP"),
                                              "any_ncp"),
    "CP"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("CP", "NCP", "CP"),
                                              "any_ncp"),
    "NCP"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("NCP", "NCP", "NCP"),
                                              "any_ncp"),
    "NCP"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("CP", "absent"),
                                              "any_ncp"),
    "CP"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("absent", "absent"),
                                              "any_ncp"),
    "absent"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(character(0), "any_ncp"),
    "absent"
  )
})

test_that("aggregate_preflight_one_component: majority breaks ties toward NCP", {
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("CP", "CP", "NCP"),
                                              "majority"),
    "CP"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("CP", "NCP"),
                                              "majority"),
    "NCP"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("NCP", "NCP", "CP"),
                                              "majority"),
    "NCP"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("CP"), "majority"),
    "CP"
  )
})

test_that("aggregate_preflight_one_component: per_k returns sentinel or absent", {
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("CP", "NCP"), "per_k"),
    "per_k"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(c("absent", "absent"),
                                              "per_k"),
    "absent"
  )
  expect_identical(
    gdpar:::aggregate_preflight_one_component(character(0), "per_k"),
    "absent"
  )
})

test_that("aggregate_preflight_one_component rejects unknown methods", {
  expect_error(
    gdpar:::aggregate_preflight_one_component(c("CP", "NCP"),
                                              "majority_lite"),
    class = "gdpar_input_error"
  )
})


# ---- agreement unit tests ---------------------------------------

test_that("compute_preflight_agreement matches frequency under each strategy", {
  expect_equal(
    gdpar:::compute_preflight_agreement(c("CP", "CP", "NCP"), "CP"),
    2 / 3
  )
  expect_equal(
    gdpar:::compute_preflight_agreement(c("CP", "CP", "NCP"), "NCP"),
    1 / 3
  )
  expect_equal(
    gdpar:::compute_preflight_agreement(c("CP", "NCP", "NCP"),
                                        "per_k"),
    2 / 3
  )
  expect_true(is.na(
    gdpar:::compute_preflight_agreement(character(0), "CP")
  ))
  expect_true(is.na(
    gdpar:::compute_preflight_agreement(c("absent", "absent"), "CP")
  ))
  expect_true(is.na(
    gdpar:::compute_preflight_agreement(c("CP", "NCP"), "absent")
  ))
})


# ---- constructor + validators -----------------------------------

test_that("new_gdpar_preflight_report constructs with default any_ncp", {
  per_dim <- make_per_dim_full(p = 2L,
                               decisions_a = c("CP", "CP"),
                               decisions_W = c("NCP", "NCP"))
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  expect_s3_class(rep, "gdpar_preflight_report")
  expect_identical(rep$method, "any_ncp")
  expect_identical(rep$global$component, c("a", "W"))
  expect_identical(rep$global$global_decision, c("CP", "NCP"))
  expect_equal(rep$global$agreement, c(1, 1))
  expect_identical(rep$global$method, c("any_ncp", "any_ncp"))
})

test_that("new_gdpar_preflight_report aggregates correctly with majority and per_k", {
  per_dim <- make_per_dim_full(p = 3L,
                               decisions_a = c("CP", "CP", "NCP"),
                               decisions_W = c("CP", "NCP", "NCP"))
  rep_maj <- gdpar:::new_gdpar_preflight_report(per_dim,
                                                method = "majority")
  expect_identical(rep_maj$global$global_decision, c("CP", "NCP"))
  expect_equal(rep_maj$global$agreement, c(2 / 3, 2 / 3))

  rep_pk <- gdpar:::new_gdpar_preflight_report(per_dim,
                                               method = "per_k")
  expect_identical(rep_pk$global$global_decision,
                   c("per_k", "per_k"))
  expect_equal(rep_pk$global$agreement, c(2 / 3, 2 / 3))
})

test_that("new_gdpar_preflight_report handles edge case p=1 cleanly", {
  per_dim <- make_per_dim_full(p = 1L,
                               decisions_a = "CP",
                               decisions_W = "NCP")
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  expect_identical(rep$global$global_decision, c("CP", "NCP"))
  expect_equal(rep$global$agreement, c(1, 1))
})

test_that("new_gdpar_preflight_report handles all-absent component", {
  per_dim <- make_per_dim_full(p = 2L,
                               decisions_a = c("CP", "NCP"),
                               decisions_W = c("absent", "absent"))
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  expect_identical(rep$global$global_decision, c("NCP", "absent"))
  expect_equal(rep$global$agreement[1], 0.5)
  expect_true(is.na(rep$global$agreement[2]))
})

test_that("new_gdpar_preflight_report rejects unknown method", {
  per_dim <- make_per_dim_full(p = 1L)
  expect_error(
    gdpar:::new_gdpar_preflight_report(per_dim, method = "weird"),
    class = "gdpar_input_error"
  )
})

test_that("new_gdpar_preflight_report rejects non-list settings", {
  per_dim <- make_per_dim_full(p = 1L)
  expect_error(
    gdpar:::new_gdpar_preflight_report(per_dim, settings = "x"),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects missing columns", {
  per_dim <- make_per_dim_full(p = 2L)
  per_dim$decision <- NULL
  expect_error(
    gdpar:::validate_preflight_per_dim(per_dim),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects non-data-frame", {
  expect_error(
    gdpar:::validate_preflight_per_dim(list(a = 1)),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects invalid component values", {
  per_dim <- make_per_dim_full(p = 1L)
  per_dim$component[1] <- "b"
  expect_error(
    gdpar:::validate_preflight_per_dim(per_dim),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects invalid decision values", {
  per_dim <- make_per_dim_full(p = 1L)
  per_dim$decision[1] <- "MAYBE"
  expect_error(
    gdpar:::validate_preflight_per_dim(per_dim),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects non-integer dim", {
  per_dim <- make_per_dim_full(p = 1L)
  per_dim$dim <- as.numeric(per_dim$dim) + 0.5
  expect_error(
    gdpar:::validate_preflight_per_dim(per_dim),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects NA dim", {
  per_dim <- make_per_dim_full(p = 1L)
  per_dim$dim[1] <- NA_integer_
  expect_error(
    gdpar:::validate_preflight_per_dim(per_dim),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects duplicate (component, dim)", {
  per_dim <- make_per_dim_full(p = 2L)
  per_dim$dim <- c(1L, 1L, 1L, 1L)
  expect_error(
    gdpar:::validate_preflight_per_dim(per_dim),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_per_dim rejects non-numeric stats", {
  per_dim <- make_per_dim_full(p = 1L)
  per_dim$t_info_cp <- as.character(per_dim$t_info_cp)
  expect_error(
    gdpar:::validate_preflight_per_dim(per_dim),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_global rejects out-of-range agreement", {
  bad <- data.frame(
    component       = "a",
    global_decision = "CP",
    agreement       = 1.5,
    method          = "any_ncp",
    stringsAsFactors = FALSE
  )
  expect_error(
    gdpar:::validate_preflight_global(bad),
    class = "gdpar_input_error"
  )
})

test_that("validate_preflight_global tolerates NA agreement", {
  ok <- data.frame(
    component       = "a",
    global_decision = "absent",
    agreement       = NA_real_,
    method          = "any_ncp",
    stringsAsFactors = FALSE
  )
  expect_invisible(gdpar:::validate_preflight_global(ok))
})


# ---- accessors --------------------------------------------------

test_that("preflight_per_dim returns the per-dim data frame", {
  per_dim <- make_per_dim_full(p = 2L)
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  expect_identical(preflight_per_dim(rep), per_dim)
})

test_that("preflight_global_decision returns the global summary", {
  per_dim <- make_per_dim_full(p = 2L,
                               decisions_a = c("CP", "NCP"),
                               decisions_W = c("NCP", "NCP"))
  rep <- gdpar:::new_gdpar_preflight_report(per_dim,
                                            method = "majority")
  g <- preflight_global_decision(rep)
  expect_s3_class(g, "data.frame")
  expect_identical(g$global_decision, c("NCP", "NCP"))
  expect_identical(g$method, c("majority", "majority"))
})

test_that("preflight_per_dim and preflight_global_decision reject non-reports", {
  expect_error(preflight_per_dim(list()),
               class = "gdpar_input_error")
  expect_error(preflight_global_decision(list()),
               class = "gdpar_input_error")
})


# ---- S3 methods -------------------------------------------------

test_that("print.gdpar_preflight_report level='global' renders compact output", {
  per_dim <- make_per_dim_full(p = 2L,
                               decisions_a = c("CP", "CP"),
                               decisions_W = c("NCP", "NCP"))
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  out <- capture.output(print(rep))
  expect_true(any(grepl("<gdpar_preflight_report>", out)))
  expect_true(any(grepl("aggregation", out)))
  expect_true(any(grepl("global decisions:", out)))
  expect_true(any(grepl("\\bCP\\b", out)))
  expect_true(any(grepl("\\bNCP\\b", out)))
  expect_true(any(grepl("uniform", out)))
})

test_that("print.gdpar_preflight_report level='dim' renders per-dim table", {
  per_dim <- make_per_dim_full(p = 2L)
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  out <- capture.output(print(rep, level = "dim"))
  expect_true(any(grepl("per-coordinate decisions:", out)))
  expect_true(any(grepl("decision_reason", out)))
  expect_true(any(grepl("t_info_cp", out)))
})

test_that("print.gdpar_preflight_report level='both' shows both blocks", {
  per_dim <- make_per_dim_full(p = 2L)
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  out <- capture.output(print(rep, level = "both"))
  expect_true(any(grepl("global decisions:", out)))
  expect_true(any(grepl("per-coordinate decisions:", out)))
})

test_that("print.gdpar_preflight_report rejects unknown level", {
  rep <- gdpar:::new_gdpar_preflight_report(make_per_dim_full(p = 1L))
  expect_error(print(rep, level = "weird"))
})

test_that("print.gdpar_preflight_report returns invisibly", {
  rep <- gdpar:::new_gdpar_preflight_report(make_per_dim_full(p = 1L))
  out <- withVisible(print(rep))
  expect_false(out$visible)
  expect_identical(out$value, rep)
})

test_that("summary.gdpar_preflight_report contracts", {
  per_dim <- make_per_dim_full(p = 3L,
                               decisions_a = c("CP", "CP", "NCP"),
                               decisions_W = c("NCP", "NCP", "NCP"))
  rep <- gdpar:::new_gdpar_preflight_report(per_dim,
                                            method = "majority")
  s <- summary(rep)
  expect_named(s, c("n_components", "n_dims", "method",
                    "per_component", "overall_agreement",
                    "per_dim_counts", "settings"))
  expect_identical(s$n_components, 2L)
  expect_identical(s$n_dims, 3L)
  expect_identical(s$method, "majority")
  expect_s3_class(s$per_component, "data.frame")
  expect_equal(s$overall_agreement,
               mean(s$per_component$agreement))
  expect_identical(s$per_dim_counts,
                   c(CP = 2L, NCP = 4L, absent = 0L))
})

test_that("summary.gdpar_preflight_report handles all-NA agreement", {
  per_dim <- make_per_dim_full(p = 1L,
                               decisions_a = "absent",
                               decisions_W = "absent")
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  s <- summary(rep)
  expect_true(is.na(s$overall_agreement))
})

test_that("as.data.frame.gdpar_preflight_report returns per_dim", {
  per_dim <- make_per_dim_full(p = 2L)
  rep <- gdpar:::new_gdpar_preflight_report(per_dim)
  df <- as.data.frame(rep)
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), 4L)
  expect_identical(colnames(df), gdpar:::.preflight_per_dim_cols)
})

test_that("format.gdpar_preflight_report returns one-line string", {
  per_dim <- make_per_dim_full(p = 2L)
  rep <- gdpar:::new_gdpar_preflight_report(per_dim,
                                            method = "any_ncp")
  s <- format(rep)
  expect_type(s, "character")
  expect_length(s, 1L)
  expect_match(s, "p=2")
  expect_match(s, "aggregation=any_ncp")
  expect_match(s, "components=a,W")
})


# ---- type-stability + invariants --------------------------------

test_that("agreement is bounded in [0, 1] across all aggregation methods", {
  per_dim <- make_per_dim_full(p = 4L,
                               decisions_a = c("CP", "CP", "NCP", "NCP"),
                               decisions_W = c("NCP", "NCP", "NCP", "CP"))
  for (m in c("any_ncp", "majority", "per_k")) {
    rep <- gdpar:::new_gdpar_preflight_report(per_dim, method = m)
    expect_true(all(rep$global$agreement >= 0 - 1e-12 &
                    rep$global$agreement <= 1 + 1e-12))
  }
})

test_that("aggregator is deterministic (idempotent on identical input)", {
  per_dim <- make_per_dim_full(p = 3L,
                               decisions_a = c("CP", "NCP", "CP"),
                               decisions_W = c("NCP", "CP", "NCP"))
  rep1 <- gdpar:::new_gdpar_preflight_report(per_dim,
                                             method = "majority")
  rep2 <- gdpar:::new_gdpar_preflight_report(per_dim,
                                             method = "majority")
  expect_identical(rep1$global, rep2$global)
})

test_that("aggregator monotonicity: any_ncp <= majority for CP rate", {
  # If majority says NCP, any_ncp must also say NCP (any_ncp is
  # strictly more conservative). The contrapositive: any_ncp = CP
  # implies majority = CP.
  cases <- list(
    c("CP", "CP", "CP"),
    c("CP", "CP", "NCP"),
    c("CP", "NCP", "NCP"),
    c("NCP", "NCP", "NCP"),
    c("CP", "CP", "NCP", "NCP")
  )
  for (decs in cases) {
    a_any <- gdpar:::aggregate_preflight_one_component(decs, "any_ncp")
    a_maj <- gdpar:::aggregate_preflight_one_component(decs, "majority")
    if (identical(a_any, "CP")) expect_identical(a_maj, "CP")
  }
})

test_that("aggregator with a single CP coordinate yields any_ncp = CP", {
  expect_identical(
    gdpar:::aggregate_preflight_one_component("CP", "any_ncp"),
    "CP"
  )
})


# ---- settings list flows through the constructor ----------------

test_that("settings list is preserved by the constructor", {
  per_dim <- make_per_dim_full(p = 1L)
  s <- list(tau_cp = 5, tau_ncp = 2, n_chains = 2L,
            n_warmup = 200L, n_sampling = 200L,
            adapt_delta = 0.95, max_treedepth = 10L)
  rep <- gdpar:::new_gdpar_preflight_report(per_dim, settings = s)
  expect_identical(rep$settings, s)
  expect_identical(summary(rep)$settings, s)
})
