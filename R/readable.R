#' @include parameter_class.R generics.R ar1.R compound_symmetry.R
#' @include autoregressive.R scaled_matrix.R simplex.R
NULL

# The free vector is what a fit estimates and not what a reader reads. Nobody
# reads the hyperbolic arctangent of a partial autocorrelation, and the
# quantity behind it is not recoverable from the printed summary of a
# multivariate fit either, whose standard deviations and correlations describe
# the matrix rather than the family that produced it. A family therefore
# declares the quantities it is about, together with the derivatives a
# consumer needs to carry a variance matrix onto them and the scale each
# interval should be built on.


#' Quantities a Family Is About
#'
#' @description
#' The interpretable quantities of a parameter, with the Jacobian of the map
#' from the free vector and the scale each quantity's interval belongs on.
#'
#' @details
#' A consumer that has estimated the free vector and its variance matrix
#' \eqn{V} obtains standard errors for these quantities by the delta method,
#' \eqn{J V J'}, and builds each interval on the declared scale before mapping
#' it back, so that a variance stays positive and a correlation stays inside
#' its interval. The Jacobians here are closed form: for the families whose
#' coordinates are separate scalar links it is the diagonal of the inverse
#' link's derivative, and for \code{\link{autoregressive}} the autoregressive
#' coefficients carry their own derivatives out of the Levinson-Durbin
#' recursion, which already carries their derivatives.
#'
#' The base class declares nothing, so a family that has no quantity beyond
#' the matrix it builds loses nothing by saying so: the matrix itself is what
#' a consumer already reports.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Passed to methods.
#'
#' @return \code{NULL} when the family declares nothing, otherwise a list with
#'   \describe{
#'     \item{\code{value}}{a named numeric vector of the quantities;}
#'     \item{\code{jacobian}}{a matrix with one row per quantity and one column
#'       per free value;}
#'     \item{\code{transform}}{a character vector naming the scale each
#'       interval is built on, one of \code{"identity"}, \code{"log"},
#'       \code{"atanh"} or \code{"logit"};}
#'     \item{\code{label}}{a single string naming the group, for a consumer
#'       laying out a printed summary. The family supplies it because the
#'       family is what holds the reading.}
#'   }
#'
#' @seealso \code{\link{param_value}}, \code{\link{autoregressive}}
#'
#' @examples
#' param_readable(ar1(5), c(log(2), atanh(0.6)))$value
#' param_readable(autoregressive(8, 2), c(0, 0.9, -0.4))$value
#'
#' @export
param_readable <- S7::new_generic("param_readable", "s",
  function(s, eta, ...) S7::S7_dispatch()
)


#' No Declared Quantities
#'
#' @description
#' The base-class method, which declares nothing.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Ignored.
#'
#' @return \code{NULL}.
#'
#' @name param_readable.parameter
#' @keywords internal
S7::method(param_readable, parameter) <- function(s, eta, ...) NULL


#' Quantities That Are Separate Links of Separate Free Values
#'
#' @description
#' Assembles the declaration of a family whose quantities are one scalar link
#' each of one free value each, so that the Jacobian is the diagonal matrix of
#' the inverse links' first derivatives.
#'
#' @param links A list of \pkg{linkfunctions7} links, one per quantity, in the
#'   order of the free values they read.
#' @param eta A numeric vector of free values.
#' @param nm A character vector naming the quantities.
#' @param transform A character vector naming the scale each interval is built
#'   on.
#' @param label A single string naming the group.
#'
#' @return A list as described in \code{\link{param_readable}}.
#'
#' @keywords internal
readable_diagonal <- function(links, eta, nm, transform, label) {
  value <- vapply(seq_along(links), function(k) {
    linkfunctions7::linkinv(links[[k]], eta[k])
  }, numeric(1))
  jac <- matrix(0, length(nm), length(eta), dimnames = list(nm, NULL))
  for (k in seq_along(links)) {
    jac[k, k] <- linkfunctions7::dlinkinv(links[[k]], eta[k])
  }
  list(value = stats::setNames(value, nm), jacobian = jac,
       transform = stats::setNames(transform, nm), label = label)
}


#' The Scale and the Correlation of an AR(1)
#'
#' @description
#' Declares the marginal variance and the correlation at lag one.
#'
#' @param s An \code{\link{Ar1Param}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Ignored.
#'
#' @return A list as described in \code{\link{param_readable}}.
#'
#' @name param_readable.Ar1Param
#' @keywords internal
S7::method(param_readable, Ar1Param) <- function(s, eta, ...) {
  readable_diagonal(
    list(s@param_params$link_scale, s@param_params$link_rho), eta,
    c("scale", "rho"), c("log", "atanh"), "Autoregressive structure"
  )
}


#' The Scale and the Common Correlation of a Compound Symmetry
#'
#' @description
#' Declares the marginal variance and the correlation shared by every pair.
#'
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Ignored.
#'
#' @return A list as described in \code{\link{param_readable}}.
#'
#' @name param_readable.CompoundSymmetryParam
#' @keywords internal
S7::method(param_readable, CompoundSymmetryParam) <- function(s, eta, ...) {
  readable_diagonal(
    list(s@param_params$link_scale, s@param_params$link_rho), eta,
    c("scale", "rho"), c("log", "atanh"), "Compound symmetry structure"
  )
}


#' The Scale of a Fixed Matrix
#'
#' @description
#' Declares the multiplier, when the family carries one.
#'
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of at most one free value.
#' @param ... Ignored.
#'
#' @return A list as described in \code{\link{param_readable}}, or \code{NULL}
#'   when the matrix is fixed and carries no free value.
#'
#' @name param_readable.ScaledMatrixParam
#' @keywords internal
S7::method(param_readable, ScaledMatrixParam) <- function(s, eta, ...) {
  if (is.null(s@param_params$link)) return(NULL)
  readable_diagonal(list(s@param_params$link), eta, "scale", "log",
                    "Scale of the fixed matrix")
}


#' The Probabilities Behind an Additive Log-Ratio
#'
#' @description
#' Declares the probability vector, whose Jacobian is the one of the softmax
#' map: \eqn{\partial p_i/\partial\eta_j = p_i(\delta_{ij} - p_j)}, the last
#' coordinate contributing \eqn{-p_k p_j}.
#'
#' @param s A \code{\link{SimplexParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Ignored.
#'
#' @return A list as described in \code{\link{param_readable}}.
#'
#' @name param_readable.SimplexParam
#' @keywords internal
S7::method(param_readable, SimplexParam) <- function(s, eta, ...) {
  p <- as.numeric(param_value(s, eta))
  k <- length(p)
  nm <- paste0("p", seq_len(k))
  jac <- matrix(0, k, k - 1L, dimnames = list(nm, NULL))
  for (i in seq_len(k)) {
    for (j in seq_len(k - 1L)) {
      jac[i, j] <- p[i] * ((i == j) - p[j])
    }
  }
  list(value = stats::setNames(p, nm), jacobian = jac,
       transform = stats::setNames(rep("logit", k), nm),
       label = "Probabilities")
}


#' The Scale, the Partial Autocorrelations and the Coefficients
#'
#' @description
#' Declares the marginal variance, the partial autocorrelations that
#' parametrize the family, and the autoregressive coefficients they produce.
#'
#' @details
#' The coefficients are the quantity an autoregression is usually reported in
#' and are not obtainable from the covariance the fit prints. They come out of
#' the Levinson-Durbin recursion the family already propagates derivatives
#' through, so their derivatives with respect to every free value are read
#' off the first-order block of those arrays rather than computed again. The
#' stationary region
#' in the coefficients is not a box, so their intervals are built on the
#' identity scale: no scalar transformation expresses the constraint they are
#' under, and an interval that respected it would not be an interval.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Ignored.
#'
#' @return A list as described in \code{\link{param_readable}}.
#'
#' @name param_readable.AutoregressiveParam
#' @keywords internal
S7::method(param_readable, AutoregressiveParam) <- function(s, eta, ...) {
  q <- s@param_params$order
  n <- q + 1L
  tay <- ar_taylor(s, eta)
  scl <- c(
    linkfunctions7::linkinv(s@param_params$link_scale, eta[1L]),
    linkfunctions7::dlinkinv(s@param_params$link_scale, eta[1L])
  )
  rv <- vapply(seq_len(q), function(k) {
    c(
      linkfunctions7::linkinv(s@param_params$link_pacf, eta[k + 1L]),
      linkfunctions7::dlinkinv(s@param_params$link_pacf, eta[k + 1L])
    )
  }, numeric(2))
  rv <- matrix(rv, nrow = 2L)
  nm <- c("scale", paste0("pacf", seq_len(q)), paste0("phi", seq_len(q)))
  jac <- matrix(0, length(nm), s@n_free, dimnames = list(nm, NULL))
  jac[1L, 1L] <- scl[2L]
  for (k in seq_len(q)) jac[1L + k, 1L + k] <- rv[2L, k]
  for (j in seq_len(q)) jac[1L + q + j, ] <- tay$phi[j, 1L + seq_len(n)]
  list(
    value = stats::setNames(c(scl[1L], rv[1L, ], tay$phi[, 1L]), nm),
    jacobian = jac,
    transform = stats::setNames(
      c("log", rep("atanh", q), rep("identity", q)), nm
    ),
    label = "Autoregressive structure"
  )
}
