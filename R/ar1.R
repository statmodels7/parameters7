#' @include compound_symmetry.R
NULL


#' AR(1) Parameter
#'
#' @description
#' The S7 class of first-order autoregressive covariance matrices: equal
#' variances and a correlation falling geometrically with the lag.
#' Constructed by \code{\link{ar1}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{Ar1Param}.
#'
#' @seealso \code{\link{ar1}}
#'
#' @examples
#' S7::S7_inherits(ar1(3), Ar1Param)
#'
#' @export
Ar1Param <- S7::new_class("Ar1Param", parent = matrix_parameter)


#' Construct an AR(1) Parameter
#'
#' @description
#' The first-order autoregressive covariance
#' \deqn{M(\eta)_{ij} = \sigma^2 \rho^{\lvert i - j \rvert},}
#' with \eqn{\sigma^2} positive and \eqn{\rho \in (-1, 1)}, so two free values
#' whatever the dimension.
#'
#' @details
#' Unlike compound symmetry, the correlation is bounded only by
#' \eqn{\lvert\rho\rvert < 1} at every dimension, so it is carried by
#' \code{\link[linkfunctions7]{rhobit_link}} -- the inverse hyperbolic tangent
#' -- and every free value gives a positive definite matrix.
#'
#' Two quantities are closed form. The determinant of the correlation pattern
#' is \eqn{(1-\rho^2)^{p-1}}, so
#' \deqn{\log\lvert M \rvert = p\log\sigma^2 + (p-1)\log(1-\rho^2),}
#' again a sum of a function of one free value and a function of the other,
#' with every mixed derivative exactly zero. And the inverse is
#' \strong{tridiagonal}, which is the property the family is used for: an
#' AR(1) process is Markov, so its precision has no entries beyond the first
#' off-diagonal, and the solve is returned from that form rather than from a
#' factorization.
#'
#' The pattern is not linear in the correlation as compound symmetry's is: an
#' entry is \eqn{\rho^{m}} for the lag \eqn{m}, so its derivatives in the free
#' value are the composition of a power with the link, taken to fourth order
#' by \code{\link{compose4}}. This is also why an AR(1) covariance and an
#' AR(1) precision are different models -- the inverse is tridiagonal but not
#' AR(1) -- which is what the \code{role} label exists to record.
#'
#' @param dimension The side \eqn{p} of the matrix, at least 2.
#' @param link_scale A \pkg{linkfunctions7} link onto the positive scale.
#'   Defaults to \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{Ar1Param}}.
#'
#' @seealso \code{\link{compound_symmetry}}, \code{\link{correlation_matrix}}
#'
#' @examples
#' s <- ar1(4)
#' eta <- c(log(2), atanh(0.6))
#' round(param_value(s, eta), 4)
#'
#' # the inverse is tridiagonal: an AR(1) process is Markov
#' round(param_solve(s, eta), 4)
#'
#' # two free values whatever the dimension
#' c(p4 = ar1(4)@n_free, p20 = ar1(20)@n_free)
#'
#' @export
ar1 <- function(dimension,
                link_scale = linkfunctions7::log_link(),
                role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)
  if (p < 2L) {
    stop(paste0(
      "'dimension' must be at least 2: a one by one matrix has no\n",
      "  correlation, and the family would carry a free value with no effect."
    ), call. = FALSE)
  }
  check_positive_link(link_scale)
  link_rho <- linkfunctions7::rhobit_link()

  Ar1Param(
    param_name = "ar1",
    dimension = p,
    n_free = 2L,
    free_names = c(tagged_name(link_scale, "scale"),
                   tagged_name(link_rho, "rho")),
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(
      link_scale = link_scale,
      link_rho = link_rho
    )
  )
}


#' The Pattern of an AR(1) Parameter
#'
#' @description
#' \eqn{P(\rho)_{ij} = \rho^{\lvert i-j \rvert}} and its four derivatives in
#' the free value, entry by entry: the composition of a power with the link.
#'
#' @details
#' Each distinct lag is composed once and written into every entry that
#' carries it, the matrix having only \eqn{p} distinct values.
#'
#' @param s An \code{\link{Ar1Param}} object.
#' @param sc The scalars of \code{\link{econ_scalars}}.
#'
#' @return A list of five matrices.
#'
#' @keywords internal
ar1_pattern <- function(s, sc) {
  p <- s@dimension
  lag <- abs(outer(seq_len(p), seq_len(p), "-"))
  r <- sc$rho
  by_lag <- lapply(0:(p - 1L), function(m) {
    c(list(r[[1L]]^m), compose4(power_derivs(r[[1L]], m), r[-1L]))
  })
  lapply(1:5, function(k) {
    matrix(vapply(as.vector(lag), function(m) by_lag[[m + 1L]][[k]], numeric(1)),
           p, p)
  })
}


#' @title Value of an AR(1) Parameter
#' @name param_value.Ar1Param
#' @description \eqn{\sigma^2 \rho^{\lvert i-j \rvert}}.
#' @param s An \code{\link{Ar1Param}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return An AR(1) positive definite matrix.
#' @keywords internal
S7::method(param_value, Ar1Param) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta)
  name_dims(sc$scale[[1L]] * ar1_pattern(s, sc)[[1L]], s)
}


#' @title Free Vector of an AR(1) Parameter
#' @name param_free.Ar1Param
#' @description
#' The scale is the common diagonal entry and the correlation is the first
#' off-diagonal one divided by it, both exact; the rest of the matrix is then
#' checked against the pattern those two imply, and a matrix that does not
#' match is refused rather than fitted.
#' @param s An \code{\link{Ar1Param}} object.
#' @param m An AR(1) matrix.
#' @param ... Unused.
#' @return A named numeric vector of two free values.
#' @keywords internal
S7::method(param_free, Ar1Param) <- function(s, m, ...) {
  p <- s@dimension
  d <- diag(m)
  scl <- max(1, max(abs(m)))
  if (diff(range(d)) > 1e-8 * scl) {
    stop("'m' does not have a constant diagonal, so it is not AR(1).",
      call. = FALSE)
  }
  if (d[1L] <= 0) stop("'m' has a non-positive diagonal.", call. = FALSE)
  r <- m[2L, 1L] / d[1L]
  if (abs(r) >= 1) {
    stop("'m' implies a correlation outside (-1, 1), so it is not AR(1).",
      call. = FALSE)
  }
  want <- d[1L] * r^abs(outer(seq_len(p), seq_len(p), "-"))
  if (max(abs(m - want)) > 1e-8 * scl) {
    stop(paste0(
      "'m' does not follow the geometric pattern its first two entries\n",
      "  imply, so it is not AR(1). It is refused rather than fitted."
    ), call. = FALSE)
  }
  stats::setNames(
    c(linkfunctions7::linkfun(s@param_params$link_scale, d[1L]),
      linkfunctions7::linkfun(s@param_params$link_rho, r)),
    s@free_names
  )
}


#' @title Solve of an AR(1) Parameter
#' @name param_solve.Ar1Param
#' @description
#' Exact and tridiagonal. The precision of an AR(1) pattern is
#' \eqn{(1-\rho^2)^{-1}} times the matrix with \eqn{1} at the two corners of
#' the diagonal, \eqn{1+\rho^2} elsewhere on it and \eqn{-\rho} on the first
#' off-diagonals; no factorization is performed.
#' @param s An \code{\link{Ar1Param}} object.
#' @param eta A numeric vector of two free values.
#' @param b A numeric matrix with \code{s@dimension} rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(param_solve, Ar1Param) <- function(s, eta, b = NULL, ...) {
  p <- s@dimension
  sc <- econ_scalars(s, eta)
  v <- sc$scale[[1L]]
  r <- sc$rho[[1L]]
  t_mat <- diag(c(1, rep(1 + r^2, p - 2L), 1), nrow = p)
  off <- cbind(seq_len(p - 1L), seq.int(2L, p))
  t_mat[off] <- -r
  # drop = FALSE, because at p = 2 the index matrix has one row and would
  # otherwise collapse to a length-two vector, which R reads as two LINEAR
  # positions rather than as one row-column pair -- writing into [1, 1].
  t_mat[off[, c(2L, 1L), drop = FALSE]] <- -r
  (t_mat / (v * (1 - r^2))) %*% b
}


#' @title Derivatives of an AR(1) Parameter
#' @name param_d1.Ar1Param
#' @description
#' Closed form at every order: a component with \eqn{a} scale indices and
#' \eqn{b} correlation indices is the \eqn{a}-th derivative of the scale times
#' the \eqn{b}-th derivative of the pattern, and the pattern's derivatives are
#' powers composed with the link.
#' @param s An \code{\link{Ar1Param}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 1L, ar1_pattern)
}

#' @rdname param_d1.Ar1Param
#' @name param_d2.Ar1Param
#' @keywords internal
S7::method(param_d2, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 2L, ar1_pattern)
}

#' @rdname param_d1.Ar1Param
#' @name param_d3.Ar1Param
#' @keywords internal
S7::method(param_d3, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 3L, ar1_pattern)
}

#' @rdname param_d1.Ar1Param
#' @name param_d4.Ar1Param
#' @keywords internal
S7::method(param_d4, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 4L, ar1_pattern)
}


#' The Log-Determinant Terms of an AR(1) Parameter
#'
#' @description
#' The affine-logarithm terms of \eqn{q(\rho) = (p-1)\log(1-\rho^2)}, split as
#' \eqn{(p-1)\{\log(1-\rho) + \log(1+\rho)\}} so that both pieces are
#' logarithms of affine functions.
#'
#' @param s An \code{\link{Ar1Param}} object.
#'
#' @return A list of numeric triples.
#'
#' @keywords internal
ar1_logdet_terms <- function(s) {
  p <- s@dimension
  list(c(p - 1, 1, -1), c(p - 1, 1, 1))
}


#' @title Log-Determinant of an AR(1) Parameter
#' @name param_logdet.Ar1Param
#' @description
#' Closed form: \eqn{p\log\sigma^2 + (p-1)\log(1-\rho^2)}, the determinant of
#' the correlation pattern being \eqn{(1-\rho^2)^{p-1}}.
#' @param s An \code{\link{Ar1Param}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, Ar1Param) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta)
  p <- s@dimension
  p * log(sc$scale[[1L]]) + (p - 1) * log(1 - sc$rho[[1L]]^2)
}

#' @title Log-Determinant Derivatives of an AR(1) Parameter
#' @name param_dlogdet.Ar1Param
#' @description
#' Closed form at every order, with every mixed component exactly zero, the
#' log-determinant being a sum of a function of the scale and a function of
#' the correlation.
#' @param s An \code{\link{Ar1Param}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 1L, ar1_logdet_terms(s))
}

#' @rdname param_dlogdet.Ar1Param
#' @name param_d2logdet.Ar1Param
#' @keywords internal
S7::method(param_d2logdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 2L, ar1_logdet_terms(s))
}

#' @rdname param_dlogdet.Ar1Param
#' @name param_d3logdet.Ar1Param
#' @keywords internal
S7::method(param_d3logdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 3L, ar1_logdet_terms(s))
}

#' @rdname param_dlogdet.Ar1Param
#' @name param_d4logdet.Ar1Param
#' @keywords internal
S7::method(param_d4logdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 4L, ar1_logdet_terms(s))
}
