test_that("gdpar_family_multi() builds homogeneous gaussian from name", {
  fam <- gdpar_family_multi("gaussian", p = 2L)
  expect_s3_class(fam, "gdpar_family_multi")
  expect_equal(fam$p, 2L)
  expect_true(fam$homogeneous)
  expect_equal(fam$name, "gaussian")
  expect_equal(fam$link, "identity")
  expect_equal(fam$stan_id, 1L)
  expect_true(fam$has_dispersion)
  expect_equal(fam$did_status, "holds")
  expect_length(fam$families, 2L)
  for (k in seq_len(2L)) {
    expect_s3_class(fam$families[[k]], "gdpar_family")
    expect_equal(fam$families[[k]]$name, "gaussian")
  }
})

test_that("gdpar_family_multi() builds homogeneous poisson from name", {
  fam <- gdpar_family_multi("poisson", p = 3L)
  expect_equal(fam$p, 3L)
  expect_equal(fam$name, "poisson")
  expect_equal(fam$link, "log")
  expect_equal(fam$stan_id, 2L)
  expect_false(fam$has_dispersion)
  expect_length(fam$families, 3L)
})

test_that("gdpar_family_multi() builds homogeneous neg_binomial_2 from name", {
  fam <- gdpar_family_multi("neg_binomial_2", p = 2L)
  expect_equal(fam$stan_id, 3L)
  expect_true(fam$has_dispersion)
})

test_that("gdpar_family_multi() builds homogeneous bernoulli from name", {
  fam <- gdpar_family_multi("bernoulli", p = 4L)
  expect_equal(fam$stan_id, 4L)
  expect_equal(fam$link, "logit")
  expect_false(fam$has_dispersion)
})

test_that("gdpar_family_multi() accepts non-default link via name string", {
  fam <- gdpar_family_multi("gaussian", p = 2L, link = "log")
  expect_equal(fam$link, "log")
  expect_equal(fam$stan_id, 1L)
})

test_that("gdpar_family_multi() accepts a gdpar_family object", {
  base_fam <- gdpar_family("gaussian")
  fam <- gdpar_family_multi(base_fam, p = 2L)
  expect_equal(fam$name, "gaussian")
  expect_equal(fam$stan_id, 1L)
  expect_length(fam$families, 2L)
})

test_that("gdpar_family_multi() rejects link when family is a gdpar_family", {
  base_fam <- gdpar_family("gaussian")
  expect_error(
    gdpar_family_multi(base_fam, p = 2L, link = "log"),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_multi() accepts a homogeneous list of gdpar_family", {
  fams <- replicate(2L, gdpar_family("poisson"), simplify = FALSE)
  fam <- gdpar_family_multi(fams, p = 2L)
  expect_equal(fam$name, "poisson")
  expect_true(fam$homogeneous)
})

test_that("gdpar_family_multi() rejects heterogeneous list (deferred feature)", {
  fams <- list(gdpar_family("gaussian"), gdpar_family("poisson"))
  expect_error(
    gdpar_family_multi(fams, p = 2L),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_family_multi() rejects heterogeneous links (deferred feature)", {
  fams <- list(
    gdpar_family("gaussian", link = "identity"),
    gdpar_family("gaussian", link = "log")
  )
  expect_error(
    gdpar_family_multi(fams, p = 2L),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_family_multi() rejects list of wrong length", {
  fams <- replicate(3L, gdpar_family("gaussian"), simplify = FALSE)
  expect_error(
    gdpar_family_multi(fams, p = 2L),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_multi() rejects list with non-gdpar_family entries", {
  fams <- list(gdpar_family("gaussian"), "poisson")
  expect_error(
    gdpar_family_multi(fams, p = 2L),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_multi() rejects invalid p (non-positive)", {
  expect_error(
    gdpar_family_multi("gaussian", p = 0L),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_family_multi("gaussian", p = -1L),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_multi() rejects invalid p (non-integer)", {
  expect_error(
    gdpar_family_multi("gaussian", p = 2.5),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_multi() rejects invalid family argument", {
  expect_error(
    gdpar_family_multi(42, p = 2L),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_multi() works with p = 1L (degenerate multi case)", {
  fam <- gdpar_family_multi("gaussian", p = 1L)
  expect_equal(fam$p, 1L)
  expect_length(fam$families, 1L)
  expect_s3_class(fam, "gdpar_family_multi")
})

test_that("print.gdpar_family_multi() emits expected lines", {
  fam <- gdpar_family_multi("gaussian", p = 2L)
  out <- capture.output(print(fam))
  expect_match(out[1], "gdpar_family_multi")
  expect_true(any(grepl("homogeneous", out)))
  expect_true(any(grepl("p\\s+:\\s+2", out)))
  expect_true(any(grepl("name\\s+:\\s+gaussian", out)))
})

test_that("print.gdpar_family_multi() returns x invisibly", {
  fam <- gdpar_family_multi("gaussian", p = 2L)
  ret <- withVisible(print(fam))
  expect_false(ret$visible)
  expect_identical(ret$value, fam)
})

test_that("gdpar_family_multi exposes param_specs_per_coord (decision 1C)", {
  fam <- gdpar_family_multi("gaussian", p = 3L)
  expect_true(!is.null(fam$param_specs_per_coord))
  expect_length(fam$param_specs_per_coord, 3L)
  for (k in seq_len(3L)) {
    expect_length(fam$param_specs_per_coord[[k]], 2L)
    expect_s3_class(fam$param_specs_per_coord[[k]][[1L]], "gdpar_param_spec")
    expect_equal(fam$param_specs_per_coord[[k]][[1L]]$scope, "per_observation")
    expect_equal(fam$param_specs_per_coord[[k]][[2L]]$scope, "population")
  }
})

test_that("param_specs_per_coord is consistent across homogeneous coordinates", {
  fam <- gdpar_family_multi("neg_binomial_2", p = 2L)
  expect_identical(
    fam$param_specs_per_coord[[1L]],
    fam$param_specs_per_coord[[2L]]
  )
})
