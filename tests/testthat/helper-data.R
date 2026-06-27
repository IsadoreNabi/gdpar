make_gaussian_data <- function(n = 200, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n, sd = 0.5)
  data.frame(x1 = x1, x2 = x2, y = y)
}

make_poisson_data <- function(n = 200, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 1 + 0.3 * x1 - 0.2 * x2
  y <- rpois(n, exp(eta))
  data.frame(x1 = x1, x2 = x2, y = y)
}

make_bernoulli_data <- function(n = 200, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 0.3 * x1 - 0.2 * x2
  prob <- plogis(eta)
  y <- rbinom(n, size = 1, prob = prob)
  data.frame(x1 = x1, x2 = x2, y = y)
}

make_neg_binomial_data <- function(n = 200, seed = NULL,
                                   phi = 5) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 1 + 0.3 * x1 - 0.2 * x2
  mu <- exp(eta)
  y <- rnbinom(n, mu = mu, size = phi)
  data.frame(x1 = x1, x2 = x2, y = y)
}

skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  cs <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
  if (is.null(cs) || !nzchar(cs)) {
    testthat::skip("cmdstan not configured")
  }
}
