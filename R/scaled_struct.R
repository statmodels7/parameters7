#' @include numerical_fallbacks.R
NULL


#' Scaled Fixed Matrix Structure
#'
#' @description
#' The S7 class of a fixed symmetric positive semidefinite matrix carried by a
#' single scale. Constructed by \code{\link{scaled_struct}}.
#'
#' @inheritParams covstruct
#'
#' @return An object of class \code{ScaledStruct}.
#'
#' @seealso \code{\link{scaled_struct}}
#'
#' @examples
#' S7::S7_inherits(scaled_struct(diag(3)), ScaledStruct)
#'
#' @export
ScaledStruct <- S7::new_class("ScaledStruct", parent = covstruct)


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
#' makes the structure a ridge.
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
#' @return An object of class \code{\link{ScaledStruct}}.
#'
#' @seealso \code{\link{log_cholesky}}, \code{\link{struct_null_basis}}
#'
#' @examples
#' # a ridge: the identity, scaled
#' r <- scaled_struct(diag(4))
#' c(rank = r@rank, logdet = struct_logdet(r, 0))
#'
#' # a second-difference penalty: rank deficient, and its null space is the
#' # constants and the straight lines, which the penalty leaves alone
#' d <- diff(diag(6), differences = 2)
#' s <- scaled_struct(crossprod(d))
#' c(dimension = s@dimension, rank = s@rank)
#'
#' # the derivative of the log pseudo-determinant is the rank, at any scale
#' c(struct_dlogdet(s, -4), struct_dlogdet(s, 4))
#'
#' # and a fully known matrix has no free value at all
#' scaled_struct(diag(3), link = NULL)@n_free
#'
#' @export
scaled_struct <- function(p, link = linkfunctions7::log_link(),
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

  ns <- struct_null_basis(p, tol = tol)
  # The pseudo-determinant of P: a constant of the structure, computed once.
  logdet_p <- sum(log(ev[seq_len(ns$rank)]))

  if (!is.null(link)) check_positive_link(link)

  ScaledStruct(
    struct_name = if (is.null(link)) "fixed" else "scaled",
    dimension = dim_p,
    n_free = if (is.null(link)) 0L else 1L,
    free_names = if (is.null(link)) character(0) else "scale",
    rank = ns$rank,
    null_basis = ns$null_basis,
    role = role,
    struct_params = list(p = p, link = link, logdet_p = logdet_p)
  )
}


#' The Scale Behind a Free Vector, and Its Derivatives
#'
#' @description
#' The value of the link and its first two derivatives at \eqn{\eta}, or the
#' constant 1 when the structure has no free value.
#'
#' @param s A \code{\link{ScaledStruct}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with \code{h}, \code{d1} and \code{d2}.
#'
#' @keywords internal
scaled_scale <- function(s, eta) {
  lk <- s@struct_params$link
  if (is.null(lk)) return(list(h = 1, d1 = 0, d2 = 0))
  list(
    h = linkfunctions7::linkinv(lk, eta),
    d1 = linkfunctions7::dlinkinv(lk, eta),
    d2 = linkfunctions7::d2linkinv(lk, eta)
  )
}


#' @title Matrix of a Scaled Structure
#' @name struct_matrix.ScaledStruct
#' @description The fixed matrix times the scale.
#' @param s A \code{\link{ScaledStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A symmetric positive semidefinite matrix.
#' @keywords internal
S7::method(struct_matrix, ScaledStruct) <- function(s, eta, ...) {
  name_dims(scaled_scale(s, eta)$h * s@struct_params$p, s)
}


#' @title First Derivative of a Scaled Structure
#' @name struct_dmatrix.ScaledStruct
#' @description Closed form: \eqn{h'(\eta) P}.
#' @param s A \code{\link{ScaledStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list with one symmetric matrix.
#' @keywords internal
S7::method(struct_dmatrix, ScaledStruct) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(list(), character(0)))
  stats::setNames(
    list(name_dims(scaled_scale(s, eta)$d1 * s@struct_params$p, s)),
    s@free_names
  )
}


#' @title Second Derivative of a Scaled Structure
#' @name struct_d2matrix.ScaledStruct
#' @description Closed form: \eqn{h''(\eta) P}.
#' @param s A \code{\link{ScaledStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list with one symmetric matrix.
#' @keywords internal
S7::method(struct_d2matrix, ScaledStruct) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(list(), character(0)))
  stats::setNames(
    list(name_dims(scaled_scale(s, eta)$d2 * s@struct_params$p, s)),
    struct_pair_names(s)
  )
}


#' @title Log-Determinant of a Scaled Structure
#' @name struct_logdet.ScaledStruct
#' @description
#' Closed form: \eqn{r \log h(\eta) + \log|P|_+}, with \eqn{r} the rank and the
#' second term a constant of the structure. This is the log-determinant when
#' the family is of full rank and the log pseudo-determinant otherwise, in one
#' expression.
#' @param s A \code{\link{ScaledStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(struct_logdet, ScaledStruct) <- function(s, eta, ...) {
  s@rank * log(scaled_scale(s, eta)$h) + s@struct_params$logdet_p
}


#' @title Log-Determinant Gradient of a Scaled Structure
#' @name struct_dlogdet.ScaledStruct
#' @description
#' Closed form: \eqn{r\, h'(\eta)/h(\eta)}, which under the default log link is
#' the rank itself, whatever the scale.
#' @param s A \code{\link{ScaledStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(struct_dlogdet, ScaledStruct) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(numeric(0), character(0)))
  sc <- scaled_scale(s, eta)
  stats::setNames(s@rank * sc$d1 / sc$h, s@free_names)
}


#' @title Log-Determinant Hessian of a Scaled Structure
#' @name struct_d2logdet.ScaledStruct
#' @description
#' Closed form: \eqn{r\{h''/h - (h'/h)^2\}}, which under the default log link
#' is zero, the log pseudo-determinant being linear in the free value.
#' @param s A \code{\link{ScaledStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(struct_d2logdet, ScaledStruct) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(numeric(0), character(0)))
  sc <- scaled_scale(s, eta)
  stats::setNames(
    s@rank * (sc$d2 / sc$h - (sc$d1 / sc$h)^2),
    struct_pair_names(s)
  )
}


#' @title Free Vector of a Scaled Structure
#' @name struct_free.ScaledStruct
#' @description
#' The scale read off the ratio to the fixed matrix, after checking that the
#' matrix really is a multiple of it. The ratio is taken at the entry of
#' \eqn{P} of largest magnitude, which is where it is best determined.
#' @param s A \code{\link{ScaledStruct}} object.
#' @param m A symmetric matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(struct_free, ScaledStruct) <- function(s, m, ...) {
  p <- s@struct_params$p
  at <- which.max(abs(p))
  h <- m[at] / p[at]
  if (h <= 0) {
    stop("'m' is not a positive multiple of the structure's fixed matrix.",
      call. = FALSE
    )
  }
  if (max(abs(m - h * p)) > 1e-8 * max(abs(m))) {
    stop(paste0(
      "'m' is not a multiple of the structure's fixed matrix, so it is not\n",
      "  in the set this structure parametrises."
    ), call. = FALSE)
  }
  if (!s@n_free) {
    if (abs(h - 1) > 1e-8) {
      stop("this structure is fixed, and 'm' is not its matrix.", call. = FALSE)
    }
    return(stats::setNames(numeric(0), character(0)))
  }
  stats::setNames(
    linkfunctions7::linkfun(s@struct_params$link, h), s@free_names
  )
}
