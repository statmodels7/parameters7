#' @include numerical_fallbacks.R
NULL


#' Scaled Fixed Matrix Parameter
#'
#' @description
#' The S7 class of a fixed symmetric positive semidefinite matrix carried by a
#' single scale. Constructed by \code{\link{scaled_matrix}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{ScaledMatrixParam}.
#'
#' @seealso \code{\link{scaled_matrix}}
#'
#' @examples
#' S7::S7_inherits(scaled_matrix(diag(3)), ScaledMatrixParam)
#'
#' @export
ScaledMatrixParam <- S7::new_class("ScaledMatrixParam", parent = matrix_parameter)


#' Construct a Scaled Fixed Matrix
#'
#' @description
#' \eqn{M(\eta) = h(\eta)\,P} for a fixed symmetric positive semidefinite
#' \eqn{P} and a positive scale \eqn{h(\eta)}, or the fixed \eqn{P} itself when
#' no link is given.
#'
#' @details
#' This is the commonest penalty in semiparametric regression, and the reason
#' the package admits rank-deficient matrices at all. \eqn{P} may be a Gram
#' matrix of a basis derivative, a difference penalty
#' \eqn{\Delta^\top \Delta}, a neighbourhood matrix, or the identity, which
#' makes the parameter a ridge.
#'
#' Everything is a constant times a function of the scale, so nothing needs
#' deriving. With the default log link, where \eqn{h(\eta) = e^{\eta}} and
#' \eqn{\lambda = h(\eta)},
#' \deqn{\partial_\eta M = M, \qquad \partial^2_\eta M = M,}
#' \deqn{\log|M|_+ = r\,\eta + \log|P|_+, \qquad
#'       \partial_\eta \log|M|_+ = r, \qquad
#'       \partial^2_\eta \log|M|_+ = 0,}
#' with \eqn{r} the rank and \eqn{\log|P|_+} computed once at construction.
#'
#' The derivative of the log pseudo-determinant being the rank is what makes
#' the scale estimable. Writing a penalty as a negative log prior,
#' \eqn{\frac{\lambda}{2}\beta^\top P \beta - \frac{r}{2}\log\lambda}, the
#' stationary point is \eqn{\lambda = r / (\beta^\top P \beta)}, whereas
#' dropping the second term leaves a derivative of one sign and sends the scale
#' to zero. That second term is exactly the normalising constant of the prior,
#' which is the reason this package keeps it.
#'
#' A rank-deficient \eqn{P} makes the corresponding gaussian improper, which is
#' what makes it a legitimate penalty and not a legitimate density. Both
#' readings of a multivariate gaussian require full rank: a singular covariance
#' is degenerate on a subspace, and a singular precision does not normalise.
#'
#' @param p A symmetric positive semidefinite matrix.
#' @param link A \pkg{linkfunctions7} link mapping the free scale to the
#'   positive scale, or \code{NULL} for a fully known matrix with no free
#'   value. Defaults to \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}. Defaults to
#'   \code{"precision"}, since the consumer of a scaled fixed matrix is a
#'   penalty.
#' @param tol The relative tolerance below which a singular value of \code{p}
#'   counts as zero when its rank is determined.
#'
#' @return An object of class \code{\link{ScaledMatrixParam}}.
#'
#' @seealso \code{\link{log_cholesky}}, \code{\link{param_null_basis}}
#'
#' @examples
#' # a ridge: the identity, scaled
#' r <- scaled_matrix(diag(4))
#' c(rank = r@rank, logdet = param_logdet(r, 0))
#'
#' # a second-difference penalty: rank deficient, and its null space is the
#' # constants and the straight lines, which the penalty leaves alone
#' d <- diff(diag(6), differences = 2)
#' s <- scaled_matrix(crossprod(d))
#' c(dimension = s@dimension, rank = s@rank)
#'
#' # the derivative of the log pseudo-determinant is the rank, at any scale
#' c(param_dlogdet(s, -4), param_dlogdet(s, 4))
#'
#' # and a fully known matrix has no free value at all
#' scaled_matrix(diag(3), link = NULL)@n_free
#'
#' @export
scaled_matrix <- function(p, link = linkfunctions7::log_link(),
                          role = c("precision", "covariance", "either"),
                          tol = 1e-10) {
  role <- match.arg(role)

  if (!is.matrix(p) || !is.numeric(p)) {
    stop("'p' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(p) != ncol(p)) stop("'p' must be square.", call. = FALSE)
  if (anyNA(p)) stop("'p' must not contain missing values.", call. = FALSE)
  if (max(abs(p - t(p))) > 1e-8 * max(1, max(abs(p)))) {
    stop("'p' must be symmetric.", call. = FALSE)
  }
  p <- (p + t(p)) / 2
  dim_p <- as.integer(nrow(p))

  ev <- eigen(p, symmetric = TRUE, only.values = TRUE)$values
  if (max(ev) <= 0) {
    stop("'p' must be positive semidefinite and not identically zero.",
      call. = FALSE
    )
  }
  if (min(ev) < -tol * max(ev)) {
    stop(sprintf(paste0(
      "'p' is not positive semidefinite: its smallest eigenvalue is %s\n",
      "  against a largest of %s."
    ), format(min(ev)), format(max(ev))), call. = FALSE)
  }

  ns <- param_null_basis(p, tol = tol)
  # The pseudo-determinant of P: a constant of the parameter, computed once.
  logdet_p <- sum(log(ev[seq_len(ns$rank)]))

  if (!is.null(link)) check_positive_link(link)

  ScaledMatrixParam(
    param_name = if (is.null(link)) "fixed" else "scaled",
    dimension = dim_p,
    n_free = if (is.null(link)) 0L else 1L,
    free_names = if (is.null(link)) character(0) else "scale",
    rank = ns$rank,
    null_basis = ns$null_basis,
    role = role,
    param_params = list(p = p, link = link, logdet_p = logdet_p)
  )
}


#' The Scale Behind a Free Vector, and Its Derivatives
#'
#' @description
#' The value of the link and its first two derivatives at \eqn{\eta}, or the
#' constant 1 when the parameter has no free value.
#'
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with \code{h}, \code{d1} and \code{d2}.
#'
#' @keywords internal
scaled_scale <- function(s, eta) {
  lk <- s@param_params$link
  if (is.null(lk)) return(list(h = 1, d1 = 0, d2 = 0))
  list(
    h = linkfunctions7::linkinv(lk, eta),
    d1 = linkfunctions7::dlinkinv(lk, eta),
    d2 = linkfunctions7::d2linkinv(lk, eta)
  )
}


#' @title Matrix of a Scaled Parameter
#' @name param_value.ScaledMatrixParam
#' @description The fixed matrix times the scale.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A symmetric positive semidefinite matrix.
#' @keywords internal
S7::method(param_value, ScaledMatrixParam) <- function(s, eta, ...) {
  name_dims(scaled_scale(s, eta)$h * s@param_params$p, s)
}


#' @title First Derivative of a Scaled Parameter
#' @name param_d1.ScaledMatrixParam
#' @description Closed form: \eqn{h'(\eta) P}.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list with one symmetric matrix.
#' @keywords internal
S7::method(param_d1, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(list(), character(0)))
  stats::setNames(
    list(name_dims(scaled_scale(s, eta)$d1 * s@param_params$p, s)),
    s@free_names
  )
}


#' @title Second Derivative of a Scaled Parameter
#' @name param_d2.ScaledMatrixParam
#' @description Closed form: \eqn{h''(\eta) P}.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list with one symmetric matrix.
#' @keywords internal
S7::method(param_d2, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(list(), character(0)))
  stats::setNames(
    list(name_dims(scaled_scale(s, eta)$d2 * s@param_params$p, s)),
    param_tuple_names(s)
  )
}


#' @title Log-Determinant of a Scaled Parameter
#' @name param_logdet.ScaledMatrixParam
#' @description
#' Closed form: \eqn{r \log h(\eta) + \log|P|_+}, with \eqn{r} the rank and the
#' second term a constant of the parameter. This is the log-determinant when
#' the family is of full rank and the log pseudo-determinant otherwise, in one
#' expression.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  s@rank * log(scaled_scale(s, eta)$h) + s@param_params$logdet_p
}


#' @title Log-Determinant Gradient of a Scaled Parameter
#' @name param_dlogdet.ScaledMatrixParam
#' @description
#' Closed form: \eqn{r\, h'(\eta)/h(\eta)}, which under the default log link is
#' the rank itself, whatever the scale.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(numeric(0), character(0)))
  sc <- scaled_scale(s, eta)
  stats::setNames(s@rank * sc$d1 / sc$h, s@free_names)
}


#' @title Log-Determinant Hessian of a Scaled Parameter
#' @name param_d2logdet.ScaledMatrixParam
#' @description
#' Closed form: \eqn{r\{h''/h - (h'/h)^2\}}, which under the default log link
#' is zero, the log pseudo-determinant being linear in the free value.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d2logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(numeric(0), character(0)))
  sc <- scaled_scale(s, eta)
  stats::setNames(
    s@rank * (sc$d2 / sc$h - (sc$d1 / sc$h)^2),
    param_tuple_names(s)
  )
}


#' @title Free Vector of a Scaled Parameter
#' @name param_free.ScaledMatrixParam
#' @description
#' The scale read off the ratio to the fixed matrix, after checking that the
#' matrix really is a multiple of it. The ratio is taken at the entry of
#' \eqn{P} of largest magnitude, which is where it is best determined.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param m A symmetric matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(param_free, ScaledMatrixParam) <- function(s, m, ...) {
  p <- s@param_params$p
  at <- which.max(abs(p))
  h <- m[at] / p[at]
  if (h <= 0) {
    stop("'m' is not a positive multiple of the parameter's fixed matrix.",
      call. = FALSE
    )
  }
  if (max(abs(m - h * p)) > 1e-8 * max(abs(m))) {
    stop(paste0(
      "'m' is not a multiple of the parameter's fixed matrix, so it is not\n",
      "  in the set this parameter parametrises."
    ), call. = FALSE)
  }
  if (!s@n_free) {
    if (abs(h - 1) > 1e-8) {
      stop("this parameter is fixed, and 'm' is not its matrix.", call. = FALSE)
    }
    return(stats::setNames(numeric(0), character(0)))
  }
  stats::setNames(
    linkfunctions7::linkfun(s@param_params$link, h), s@free_names
  )
}


#' @title Third Derivatives of a Scaled Parameter
#' @name param_d3.ScaledMatrixParam
#' @description Closed form: \eqn{h'''(\eta)\, M_0}.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector with one free value.
#' @param ... Unused.
#' @return A named list with one matrix.
#' @keywords internal
S7::method(param_d3, ScaledMatrixParam) <- function(s, eta, ...) {
  v <- linkfunctions7::d3linkinv(s@param_params$link, eta)
  stats::setNames(
    list(name_dims(v * s@param_params$p, s)),
    param_tuple_names(s, 3L)
  )
}

#' @title Fourth Derivatives of a Scaled Parameter
#' @name param_d4.ScaledMatrixParam
#' @description Closed form: \eqn{h''''(\eta)\, M_0}.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector with one free value.
#' @param ... Unused.
#' @return A named list with one matrix.
#' @keywords internal
S7::method(param_d4, ScaledMatrixParam) <- function(s, eta, ...) {
  v <- linkfunctions7::d4linkinv(s@param_params$link, eta)
  stats::setNames(
    list(name_dims(v * s@param_params$p, s)),
    param_tuple_names(s, 4L)
  )
}

scaled_dlog <- function(s, eta, order) {
  link <- s@param_params$link
  h <- linkfunctions7::linkinv(link, eta)
  u1 <- linkfunctions7::dlinkinv(link, eta) / h
  u2 <- linkfunctions7::d2linkinv(link, eta) / h
  u3 <- linkfunctions7::d3linkinv(link, eta) / h
  u4 <- linkfunctions7::d4linkinv(link, eta) / h
  if (order == 3L) return(u3 - 3 * u1 * u2 + 2 * u1^3)
  u4 - 4 * u1 * u3 - 3 * u2^2 + 12 * u1^2 * u2 - 6 * u1^4
}

#' @title Higher Log-Determinant Derivatives of a Scaled Parameter
#' @name param_d3logdet.ScaledMatrixParam
#' @description
#' Closed form: the log-(pseudo-)determinant is
#' \eqn{r \log h(\eta) + \log\mathrm{pdet}(M_0)} with \eqn{r} the rank, so
#' every derivative is \eqn{r} times the matching derivative of
#' \eqn{\log h}.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector with one free value.
#' @param ... Unused.
#' @return A named numeric vector with one entry.
#' @keywords internal
S7::method(param_d3logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  stats::setNames(s@rank * scaled_dlog(s, eta, 3L), param_tuple_names(s, 3L))
}

#' @title Fourth Log-Determinant Derivatives of a Scaled Parameter
#' @name param_d4logdet.ScaledMatrixParam
#' @description Closed form; see \code{\link{param_d3logdet.ScaledMatrixParam}}.
#' @param s A \code{\link{ScaledMatrixParam}} object.
#' @param eta A numeric vector with one free value.
#' @param ... Unused.
#' @return A named numeric vector with one entry.
#' @keywords internal
S7::method(param_d4logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  stats::setNames(s@rank * scaled_dlog(s, eta, 4L), param_tuple_names(s, 4L))
}
