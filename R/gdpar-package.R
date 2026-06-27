#' gdpar: General Dynamic Parameter Models via Reference Anchoring
#'
#' Unified predictive framework in which individual parameters are
#' decomposed as theta_i = theta_ref + Delta(x_i, theta_ref), with
#' theta_ref a population reference and Delta an explicit deviation
#' function. The decomposition follows the Additive-Multiplicative-
#' Modulated canonical form
#' Delta(x, theta) = a(x) + b(x) * theta + W(theta) x
#' (where * denotes the elementwise product) and is estimated through
#' three complementary paths: hierarchical Bayesian inference via Stan
#' (Path 1), varying-coefficient models via penalized splines using
#' \pkg{mgcv} (Path 2), and amortized inference via hypernetworks in
#' \pkg{torch} (Path 3).
#'
#' @section Theoretical foundation:
#' The package is accompanied by ten theoretical vignettes (Blocks 0-9)
#' that establish the canonical form and identifiability conditions
#' (Block 1), the gnoseological validity of the population reference
#' (Block 2), the special-case subsumptions (Block 3), the asymptotic
#' theory for each path (Blocks 4-6), the Empirical Bayes versus Fully
#' Bayes comparison (Block 7), the positioning relative to causal
#' inference (Block 8), and the cognitive-science motivation (Block 9).
#' The implementation follows these blocks closely; methodological
#' decisions are documented in each function's Roxygen entry with
#' references to the relevant block.
#'
#' @section Three paths:
#' \describe{
#'   \item{Path 1 (hierarchical Bayesian)}{\code{\link{gdpar}} with
#'     \code{path = "bayes"}. Uses \pkg{cmdstanr} when available.
#'     Suitable for finite-dimensional parametric AMM components in
#'     this version of the package.}
#'   \item{Path 2 (varying-coefficient)}{\code{\link{gdpar}} with
#'     \code{path = "vcm"}. Uses \pkg{mgcv} for penalized spline
#'     fitting. Recommended for non-parametric components.}
#'   \item{Path 3 (hypernetwork, planned)}{\code{\link{gdpar}} with
#'     \code{path = "hyper"}. Planned to use the \pkg{torch} R package
#'     (LibTorch backend, no Python) for amortized inference and high
#'     expressive capacity. Not implemented in this release; \pkg{torch}
#'     is not currently a declared dependency of the package. Invoking
#'     \code{path = "hyper"} raises a structured
#'     \code{gdpar_unsupported_feature_error}.}
#' }
#'
#' @section Default path:
#' Path 1 is the default for finite-dimensional parametric specifications
#' because it admits closed-form identifiability results (Theorems 1A,
#' 1E of Block 1) and full Bayesian uncertainty quantification with
#' calibrated coverage in the limit (Theorem 4C of Block 4).
#'
#' @keywords internal
#' @aliases gdpar-package
"_PACKAGE"

#' @importFrom stats model.matrix terms update sd quantile median lm coef qnorm setNames qt
#' @importFrom methods is
NULL
