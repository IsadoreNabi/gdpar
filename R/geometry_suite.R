# Synthetic posterior-geometry suite of known difficulty (Block RG, RG.1).
#
# Each target is provided in dual form: a Stan program (so the RG.1 diagnostic
# can run cheap NUTS pilots through cmdstan) and an R log-density closure with
# an analytic gradient on the unconstrained scale (so the RG.2 R-native
# geometric engine can consume the same target, and so the analytic gradient
# can be cross-checked against Stan's grad_log_prob in a rigour-three-layers
# numerical test). Every target carries a ground-truth pathology label, the
# geometry level that remedies it (the remedy column of the Block RG taxonomy),
# the culprit parameter(s) where applicable, and whether the difficulty grows
# with the size knob n (the quasi-determinism vs structural-pathology
# discriminator and the 28-GB nuance of the charter).

#' Catalogue of synthetic posterior geometries of known difficulty
#'
#' Build the Block RG calibration suite: a registry of synthetic target
#' distributions, each engineered to exhibit one row of the posterior-geometry
#' pathology taxonomy with a \emph{known} ground truth. The suite is the
#' falsifiable backbone against which \code{\link{gdpar_geometry_diagnostic}}
#' is calibrated and its error rates are measured.
#'
#' Each target is returned in dual representation so that the same geometry can
#' be exercised by the cmdstan NUTS pilots of RG.1 and by the R-native
#' geometric engine of RG.2: a Stan program string plus an R log-density
#' closure carrying an analytic gradient on the unconstrained scale.
#'
#' @details
#' The eight targets cover the taxonomy of the Block RG charter:
#' \describe{
#'   \item{\code{G0_isotropic}}{Standard normal. Easy negative control;
#'     remedy is the Euclidean diagonal metric (the current default).}
#'   \item{\code{G1_anisotropic}}{Diagonal normal with a fixed, large condition
#'     number that does \emph{not} grow with \code{n}. A straight canyon;
#'     remedy is a constant (dense Euclidean) metric.}
#'   \item{\code{G2_funnel}}{Neal's funnel: a log-scale variable \code{v}
#'     governs the spread of the remaining coordinates, so curvature varies with
#'     position. Remedy is a position-dependent (Riemannian) metric.}
#'   \item{\code{G3_heavy_tails}}{Independent Student-t with a small number of
#'     degrees of freedom. Directional heaviness; remedy is a non
#'     inner-product (Finsler / relativistic) kinetic energy.}
#'   \item{\code{G4_quasi_deterministic}}{A near-degenerate canyon in which
#'     \code{d - 1} directions are pinned with variance \code{1 / n} while one
#'     direction stays free. The condition number \emph{grows} with \code{n}:
#'     the posterior collapses to a lower-dimensional manifold. This is the
#'     eBird count case; remedy is a sub-Riemannian (distribution-of-directions)
#'     treatment.}
#'   \item{\code{G5_multimodal}}{Equal-weight mixture of two separated normals.
#'     Remedy is tempering / a general metric space.}
#'   \item{\code{G6_boundary}}{One parameter bounded on (0, 1) with mass pinned
#'     against the lower bound, plus free nuisance coordinates. Singular
#'     curvature at the edge; the analogue of the Tweedie shape parameter
#'     \code{p} hugging its bound.}
#'   \item{\code{G7_flat_direction}}{A reparametrisation redundancy
#'     (\code{a + b} identified, the contrast \code{a - b} flat) with a wide
#'     prior, yielding a near-zero Hessian eigenvalue. The exact analogue of the
#'     declared-but-unused \code{sigma_a_k} scales of session B9.21; remedy is to
#'     reparametrise or eliminate (Option A).}
#' }
#'
#' Each registry entry is a list of class \code{gdpar_geometry_target} with the
#' static fields \code{id}, \code{label}, \code{pathology},
#' \code{geometry_remedy}, \code{culprit} (parameter names or \code{NA}),
#' \code{difficulty_scales_with_n} (the ground-truth answer the difficulty-vs-n
#' curve must recover), \code{bounds} (named list of constrained-scale bounds,
#' or \code{NULL}), \code{default_n}, \code{default_difficulty},
#' \code{n_grid} (a suggested size sweep), and a constructor \code{make(n,
#' difficulty)} returning an \emph{instance}: a list with \code{stan_code},
#' \code{stan_data}, \code{log_prob(theta)}, \code{grad_log_prob(theta)},
#' \code{dim}, and \code{param_names}.
#'
#' The R closures operate on the unconstrained scale (matching the scale on
#' which Hamiltonian samplers act). Absolute log-density values may differ from
#' Stan's by an additive constant; the gradients agree exactly and are what the
#' geometric machinery uses.
#'
#' @param which Optional character vector selecting a subset of target ids.
#'   Defaults to all eight targets.
#'
#' @return A named list of \code{gdpar_geometry_target} objects.
#'
#' @seealso \code{\link{gdpar_geometry_diagnostic}}.
#'
#' @examples
#' suite <- gdpar_geometry_suite()
#' names(suite)
#' funnel <- suite$G2_funnel$make(n = 1, difficulty = 3)
#' funnel$dim
#' funnel$log_prob(c(v = 0, x = rep(0, funnel$dim - 1)))
#'
#' @export
gdpar_geometry_suite <- function(which = NULL) {
  registry <- list(
    G0_isotropic           = .gdpar_geom_G0(),
    G1_anisotropic         = .gdpar_geom_G1(),
    G2_funnel              = .gdpar_geom_G2(),
    G3_heavy_tails         = .gdpar_geom_G3(),
    G4_quasi_deterministic = .gdpar_geom_G4(),
    G5_multimodal          = .gdpar_geom_G5(),
    G6_boundary            = .gdpar_geom_G6(),
    G7_flat_direction      = .gdpar_geom_G7()
  )
  if (is.null(which)) {
    return(registry)
  }
  if (!is.character(which) || !all(which %in% names(registry))) {
    gdpar_abort(
      sprintf(
        "Argument 'which' must be a subset of: %s.",
        paste(names(registry), collapse = ", ")
      ),
      class = "gdpar_input_error"
    )
  }
  registry[which]
}

# Internal constructor helper: stamps the common class and validates the make()
# contract lightly so a malformed target surfaces early.
.gdpar_geom_target <- function(id, label, pathology, geometry_remedy, culprit,
                               difficulty_scales_with_n, bounds, default_n,
                               default_difficulty, n_grid, make) {
  obj <- list(
    id = id, label = label, pathology = pathology,
    geometry_remedy = geometry_remedy, culprit = culprit,
    difficulty_scales_with_n = difficulty_scales_with_n,
    bounds = bounds, default_n = default_n,
    default_difficulty = default_difficulty, n_grid = n_grid,
    make = make
  )
  class(obj) <- c("gdpar_geometry_target", "list")
  obj
}

# ---------------------------------------------------------------------------
# G0 -- isotropic standard normal (easy negative control).
# ---------------------------------------------------------------------------
.gdpar_geom_G0 <- function() {
  stan_code <- paste(
    "data { int<lower=1> d; }",
    "parameters { vector[d] theta; }",
    "model { theta ~ std_normal(); }",
    sep = "\n"
  )
  make <- function(n = 1, difficulty = 5) {
    d <- as.integer(difficulty)
    list(
      stan_code = stan_code,
      stan_data = list(d = d),
      log_prob = function(theta) -0.5 * sum(theta^2),
      grad_log_prob = function(theta) -theta,
      dim = d,
      param_names = paste0("theta[", seq_len(d), "]")
    )
  }
  .gdpar_geom_target(
    id = "G0_isotropic", label = "Isotropic standard normal (easy control)",
    pathology = "isotropic", geometry_remedy = "euclidean_diagonal",
    culprit = NA_character_, difficulty_scales_with_n = FALSE, bounds = NULL,
    default_n = 1, default_difficulty = 5, n_grid = c(1, 4, 16), make = make
  )
}

# ---------------------------------------------------------------------------
# G1 -- anisotropic normal, fixed condition number (straight canyon).
# difficulty == condition number kappa (variance ratio). n-invariant.
# ---------------------------------------------------------------------------
.gdpar_geom_G1 <- function() {
  stan_code <- paste(
    "data { int<lower=1> d; vector<lower=0>[d] sigma; }",
    "parameters { vector[d] theta; }",
    "model { theta ~ normal(0, sigma); }",
    sep = "\n"
  )
  make <- function(n = 1, difficulty = 100) {
    d <- 5L
    kappa <- as.numeric(difficulty)
    vars <- exp(seq(0, log(kappa), length.out = d))
    sigma <- sqrt(vars)
    inv_var <- 1 / vars
    list(
      stan_code = stan_code,
      stan_data = list(d = d, sigma = sigma),
      log_prob = function(theta) -0.5 * sum((theta^2) * inv_var),
      grad_log_prob = function(theta) -theta * inv_var,
      dim = d,
      param_names = paste0("theta[", seq_len(d), "]")
    )
  }
  .gdpar_geom_target(
    id = "G1_anisotropic", label = "Anisotropic normal, fixed condition number",
    pathology = "anisotropic", geometry_remedy = "euclidean_dense",
    culprit = NA_character_, difficulty_scales_with_n = FALSE, bounds = NULL,
    default_n = 1, default_difficulty = 100, n_grid = c(1, 4, 16), make = make
  )
}

# ---------------------------------------------------------------------------
# G2 -- Neal's funnel (variable curvature). difficulty == prior sd of v.
# ---------------------------------------------------------------------------
.gdpar_geom_G2 <- function() {
  stan_code <- paste(
    "data { int<lower=1> d; real<lower=0> scale_v; }",
    "parameters { real v; vector[d - 1] x; }",
    "model {",
    "  v ~ normal(0, scale_v);",
    "  x ~ normal(0, exp(v / 2));",
    "}",
    sep = "\n"
  )
  make <- function(n = 1, difficulty = 3) {
    d <- 10L
    scale_v <- as.numeric(difficulty)
    list(
      stan_code = stan_code,
      stan_data = list(d = d, scale_v = scale_v),
      log_prob = function(theta) {
        v <- theta[1]
        x <- theta[-1]
        -0.5 * (v / scale_v)^2 - 0.5 * (d - 1) * v -
          0.5 * sum(x^2) * exp(-v)
      },
      grad_log_prob = function(theta) {
        v <- theta[1]
        x <- theta[-1]
        gv <- -v / scale_v^2 - 0.5 * (d - 1) + 0.5 * sum(x^2) * exp(-v)
        gx <- -x * exp(-v)
        c(gv, gx)
      },
      dim = d,
      param_names = c("v", paste0("x[", seq_len(d - 1), "]"))
    )
  }
  .gdpar_geom_target(
    id = "G2_funnel", label = "Neal's funnel (variable curvature)",
    pathology = "funnel", geometry_remedy = "riemannian",
    culprit = "v", difficulty_scales_with_n = FALSE, bounds = NULL,
    default_n = 1, default_difficulty = 3, n_grid = c(1, 4, 16), make = make
  )
}

# ---------------------------------------------------------------------------
# G3 -- heavy tails (independent Student-t). difficulty == nu (low = heavy).
# ---------------------------------------------------------------------------
.gdpar_geom_G3 <- function() {
  stan_code <- paste(
    "data { int<lower=1> d; real<lower=0> nu; }",
    "parameters { vector[d] theta; }",
    "model { theta ~ student_t(nu, 0, 1); }",
    sep = "\n"
  )
  make <- function(n = 1, difficulty = 2) {
    d <- 5L
    nu <- as.numeric(difficulty)
    list(
      stan_code = stan_code,
      stan_data = list(d = d, nu = nu),
      log_prob = function(theta) {
        -((nu + 1) / 2) * sum(log1p(theta^2 / nu))
      },
      grad_log_prob = function(theta) {
        -(nu + 1) * theta / (nu + theta^2)
      },
      dim = d,
      param_names = paste0("theta[", seq_len(d), "]")
    )
  }
  .gdpar_geom_target(
    id = "G3_heavy_tails", label = "Heavy tails (independent Student-t)",
    pathology = "heavy_tails", geometry_remedy = "finsler_relativistic",
    culprit = NA_character_, difficulty_scales_with_n = FALSE, bounds = NULL,
    default_n = 1, default_difficulty = 2, n_grid = c(1, 4, 16), make = make
  )
}

# ---------------------------------------------------------------------------
# G4 -- quasi-deterministic canyon (the eBird count case). Condition number
# grows with n: d - 1 directions pinned with variance 1 / n. n is the knob.
# ---------------------------------------------------------------------------
.gdpar_geom_G4 <- function() {
  stan_code <- paste(
    "data { int<lower=1> d; real<lower=0> n_eff; }",
    "parameters { vector[d] theta; }",
    "model {",
    "  theta[1] ~ normal(0, 1);",
    "  theta[2:d] ~ normal(0, inv_sqrt(n_eff));",
    "}",
    sep = "\n"
  )
  make <- function(n = 100, difficulty = 5) {
    d <- as.integer(difficulty)
    n_eff <- as.numeric(n)
    inv_var <- c(1, rep(n_eff, d - 1))
    list(
      stan_code = stan_code,
      stan_data = list(d = d, n_eff = n_eff),
      log_prob = function(theta) -0.5 * sum((theta^2) * inv_var),
      grad_log_prob = function(theta) -theta * inv_var,
      dim = d,
      param_names = paste0("theta[", seq_len(d), "]")
    )
  }
  .gdpar_geom_target(
    id = "G4_quasi_deterministic",
    label = "Quasi-deterministic canyon (condition number grows with n)",
    pathology = "quasi_deterministic", geometry_remedy = "sub_riemannian",
    culprit = NA_character_, difficulty_scales_with_n = TRUE, bounds = NULL,
    default_n = 100, default_difficulty = 5, n_grid = c(10, 100, 1000),
    make = make
  )
}

# ---------------------------------------------------------------------------
# G5 -- multimodal (equal mixture of two separated normals). difficulty == sep.
# ---------------------------------------------------------------------------
.gdpar_geom_G5 <- function() {
  stan_code <- paste(
    "data { int<lower=1> d; real sep; }",
    "parameters { vector[d] theta; }",
    "model {",
    "  target += log_mix(0.5,",
    "                    normal_lpdf(theta | -sep, 1),",
    "                    normal_lpdf(theta |  sep, 1));",
    "}",
    sep = "\n"
  )
  make <- function(n = 1, difficulty = 4) {
    d <- 3L
    sep <- as.numeric(difficulty)
    list(
      stan_code = stan_code,
      stan_data = list(d = d, sep = sep),
      log_prob = function(theta) {
        l1 <- -0.5 * sum((theta + sep)^2)
        l2 <- -0.5 * sum((theta - sep)^2)
        m <- max(l1, l2)
        log(0.5) + m + log(exp(l1 - m) + exp(l2 - m))
      },
      grad_log_prob = function(theta) {
        l1 <- -0.5 * sum((theta + sep)^2)
        l2 <- -0.5 * sum((theta - sep)^2)
        m <- max(l1, l2)
        w1 <- exp(l1 - m) / (exp(l1 - m) + exp(l2 - m))
        w2 <- 1 - w1
        w1 * (-(theta + sep)) + w2 * (-(theta - sep))
      },
      dim = d,
      param_names = paste0("theta[", seq_len(d), "]")
    )
  }
  .gdpar_geom_target(
    id = "G5_multimodal", label = "Two well-separated normal modes",
    pathology = "multimodal", geometry_remedy = "tempering",
    culprit = NA_character_, difficulty_scales_with_n = FALSE, bounds = NULL,
    default_n = 1, default_difficulty = 4, n_grid = c(1, 4, 16), make = make
  )
}

# ---------------------------------------------------------------------------
# G6 -- boundary-pegged (Beta(alpha, 1) mass at the lower bound, alpha < 1)
# plus free nuisance coordinates. The Tweedie-p-at-edge analogue.
# R closure works on the unconstrained scale phi_b = logit(b).
# ---------------------------------------------------------------------------
.gdpar_geom_G6 <- function() {
  stan_code <- paste(
    "data { int<lower=1> d; real<lower=0> alpha; }",
    "parameters { real<lower=0, upper=1> b; vector[d - 1] z; }",
    "model {",
    "  b ~ beta(alpha, 1);",
    "  z ~ std_normal();",
    "}",
    sep = "\n"
  )
  make <- function(n = 1, difficulty = 0.3) {
    d <- 4L
    alpha <- as.numeric(difficulty)
    inv_logit <- function(phi) 1 / (1 + exp(-phi))
    list(
      stan_code = stan_code,
      stan_data = list(d = d, alpha = alpha),
      log_prob = function(theta) {
        phi <- theta[1]
        z <- theta[-1]
        b <- inv_logit(phi)
        alpha * log(b) + log1p(-b) - 0.5 * sum(z^2)
      },
      grad_log_prob = function(theta) {
        phi <- theta[1]
        z <- theta[-1]
        b <- inv_logit(phi)
        gphi <- alpha - (alpha + 1) * b
        c(gphi, -z)
      },
      dim = d,
      param_names = c("phi_b", paste0("z[", seq_len(d - 1), "]"))
    )
  }
  .gdpar_geom_target(
    id = "G6_boundary", label = "Mass pinned against a (0, 1) bound",
    pathology = "boundary", geometry_remedy = "boundary_reparam",
    culprit = "b", difficulty_scales_with_n = FALSE,
    bounds = list(b = c(0, 1)),
    default_n = 1, default_difficulty = 0.3, n_grid = c(1, 4, 16), make = make
  )
}

# ---------------------------------------------------------------------------
# G7 -- flat direction / non-identification (the sigma_a_k case). a + b is
# weakly identified, the contrast a - b is flat under a wide prior (tau),
# giving a near-zero Hessian eigenvalue. Structural: n-invariant.
# ---------------------------------------------------------------------------
.gdpar_geom_G7 <- function() {
  stan_code <- paste(
    "data { real<lower=0> tau; }",
    "parameters { real a; real b; }",
    "model {",
    "  target += normal_lpdf(0 | a + b, 1);",
    "  a ~ normal(0, tau);",
    "  b ~ normal(0, tau);",
    "}",
    sep = "\n"
  )
  make <- function(n = 1, difficulty = 100) {
    tau <- as.numeric(difficulty)
    inv_tau2 <- 1 / tau^2
    list(
      stan_code = stan_code,
      stan_data = list(tau = tau),
      log_prob = function(theta) {
        a <- theta[1]
        b <- theta[2]
        -0.5 * (a + b)^2 - 0.5 * inv_tau2 * (a^2 + b^2)
      },
      grad_log_prob = function(theta) {
        a <- theta[1]
        b <- theta[2]
        s <- a + b
        c(-s - inv_tau2 * a, -s - inv_tau2 * b)
      },
      dim = 2L,
      param_names = c("a", "b")
    )
  }
  .gdpar_geom_target(
    id = "G7_flat_direction", label = "Reparametrisation redundancy (flat a - b)",
    pathology = "flat_direction", geometry_remedy = "reparam_eliminate",
    culprit = c("a", "b"), difficulty_scales_with_n = FALSE, bounds = NULL,
    default_n = 1, default_difficulty = 100, n_grid = c(1, 4, 16), make = make
  )
}

#' Print method for gdpar_geometry_target objects
#'
#' @param x An object of class \code{gdpar_geometry_target}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geometry_target <- function(x, ...) {
  cat("<gdpar_geometry_target> ", x$id, "\n", sep = "")
  cat("  ", x$label, "\n", sep = "")
  cat("  pathology: ", x$pathology,
      " | remedy: ", x$geometry_remedy, "\n", sep = "")
  cat("  culprit: ",
      if (all(is.na(x$culprit))) "none" else paste(x$culprit, collapse = ", "),
      " | difficulty scales with n: ", x$difficulty_scales_with_n, "\n",
      sep = "")
  invisible(x)
}
