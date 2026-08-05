#' @include chain.R numerical_fallbacks.R
NULL

# Compound symmetry and AR(1) are the economical families: whatever p is, two
# free values, a scale and a correlation. They share the whole skeleton -- a
# positive scale through one link, a bounded correlation through another, and
# a value that is the scale times a pattern depending on the correlation alone
# -- so the skeleton is written here once and ar1.R adds the second family to
# it, the way transition_matrix.R builds on simplex.R.
#
# The economy is the point rather than a convenience. A random effect over
# twenty occasions is 210 free values unstructured and two here, and the two
# are what a model can actually estimate from a handful of subjects.


#' Compound Symmetry Parameter
#'
#' @description
#' The S7 class of compound-symmetric covariance matrices: equal variances and
#' one common correlation. Constructed by \code{\link{compound_symmetry}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{CompoundSymmetryParam}.
#'
#' @seealso \code{\link{compound_symmetry}}
#'
#' @examples
#' S7::S7_inherits(compound_symmetry(3), CompoundSymmetryParam)
#'
#' @export
CompoundSymmetryParam <- S7::new_class("CompoundSymmetryParam",
  parent = matrix_parameter)


#' Construct a Compound Symmetry Parameter
#'
#' @description
#' The exchangeable covariance
#' \deqn{M(\eta) = \sigma^2\{(1-\rho)I + \rho J\},}
#' with \eqn{\sigma^2} positive and \eqn{\rho} the common correlation, so two
#' free values whatever the dimension.
#'
#' @details
#' Positive definiteness bounds the correlation below as well as above:
#' the eigenvalues are \eqn{\sigma^2\{1 + (p-1)\rho\}} once and
#' \eqn{\sigma^2(1-\rho)} with multiplicity \eqn{p-1}, so the matrix is
#' definite exactly when \eqn{-1/(p-1) < \rho < 1}. The correlation is
#' therefore carried by \code{\link[linkfunctions7]{bounded_link}} onto that
#' interval and not by \code{\link[linkfunctions7]{rhobit_link}}, which would
#' let a caller build an indefinite matrix at a perfectly ordinary free value.
#'
#' Two quantities are closed form and cost nothing. The eigenvalues above give
#' \deqn{\log\lvert M \rvert = p\log\sigma^2 + \log\{1 + (p-1)\rho\}
#'   + (p-1)\log(1-\rho),}
#' a sum of a function of one free value and a function of the other, so every
#' mixed derivative of the log-determinant is exactly zero. And the inverse is
#' compound symmetric again, by the Sherman-Morrison identity, so it is
#' returned in closed form rather than factorized: the inverse of an
#' exchangeable covariance is an exchangeable precision, which is what makes
#' this family closed under the choice of side.
#'
#' The value is the scale times a pattern that is \emph{linear} in the
#' correlation, \eqn{I + \rho(J - I)}, so all four derivative orders follow
#' from the two links' own derivatives with no further algebra.
#'
#' @param dimension The side \eqn{p} of the matrix, at least 2.
#' @param link_scale A \pkg{linkfunctions7} link onto the positive scale.
#'   Defaults to \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{CompoundSymmetryParam}}.
#'
#' @seealso \code{\link{ar1}}, \code{\link{correlation_matrix}}
#'
#' @examples
#' s <- compound_symmetry(4)
#' c(n_free = s@n_free)
#'
#' eta <- c(log(2), 0.8)
#' round(param_value(s, eta), 4)
#'
#' # the inverse is compound symmetric too, and is not factorized
#' round(param_solve(s, eta), 4)
#'
#' # the round trip closes exactly
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' @export
compound_symmetry <- function(dimension,
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
  link_rho <- linkfunctions7::bounded_link(lwr = -1 / (p - 1), upr = 1)

  CompoundSymmetryParam(
    param_name = "compound_symmetry",
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


#' The Scale and the Correlation of an Economical Parameter
#'
#' @description
#' The two scalars a compound-symmetric or AR(1) parameter is built from, with
#' the first four derivatives of each in its own free value.
#'
#' @param s A \code{\link{CompoundSymmetryParam}} or \code{\link{Ar1Param}}
#'   object.
#' @param eta A numeric vector of two free values.
#'
#' @return A list with \code{scale} and \code{rho}, each a list of five
#'   numbers: the value and four derivatives.
#'
#' @keywords internal
econ_scalars <- function(s, eta) {
  ls <- s@param_params$link_scale
  lr <- s@param_params$link_rho
  grab <- function(link, e) {
    list(
      linkfunctions7::linkinv(link, e),
      linkfunctions7::dlinkinv(link, e),
      linkfunctions7::d2linkinv(link, e),
      linkfunctions7::d3linkinv(link, e),
      linkfunctions7::d4linkinv(link, e)
    )
  }
  list(scale = grab(ls, eta[1L]), rho = grab(lr, eta[2L]))
}


#' Derivatives of a Sum of Logarithms of Affine Functions
#'
#' @description
#' The first four derivatives at \eqn{r} of
#' \eqn{\sum_t c_t \log(a_t + b_t r)}, from
#' \eqn{\mathrm{d}^k \log(a + br)/\mathrm{d}r^k =
#' (-1)^{k-1}(k-1)!\, b^k/(a + br)^k}.
#'
#' @details
#' Both economical families have a log-determinant of this shape -- compound
#' symmetry from its two distinct eigenvalues, AR(1) from
#' \eqn{\lvert R \rvert = (1-\rho^2)^{p-1}} -- so the derivatives are written
#' once rather than transcribed twice.
#'
#' @param r The point.
#' @param terms A list of numeric triples \code{c(coefficient, a, b)}.
#'
#' @return A list of four numbers.
#'
#' @keywords internal
log_affine_derivs <- function(r, terms) {
  lapply(1:4, function(k) {
    sum(vapply(terms, function(t) {
      t[1L] * (-1)^(k - 1L) * factorial(k - 1L) * t[3L]^k / (t[2L] + t[3L] * r)^k
    }, numeric(1)))
  })
}


#' Derivative Components of an Economical Parameter
#'
#' @description
#' Assembles one derivative order of \eqn{M = \sigma^2 P(\rho)} from the
#' scale's derivatives and the pattern's.
#'
#' @details
#' The value is a product of a function of the first free value and a function
#' of the second, so a component with \eqn{a} scale indices and \eqn{b}
#' correlation indices is the \eqn{a}-th derivative of the scale times the
#' \eqn{b}-th derivative of the pattern. Nothing is approximated and no order
#' is special.
#'
#' @param s The parameter.
#' @param eta A numeric vector of two free values.
#' @param order The derivative order, 1 to 4.
#' @param pattern A function of the parameter and the scalars, returning the
#'   pattern and its four derivatives in the second free value.
#'
#' @return A named list of symmetric matrices.
#'
#' @keywords internal
econ_derivative <- function(s, eta, order, pattern) {
  sc <- econ_scalars(s, eta)
  pt <- pattern(s, sc)
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(t) {
    a <- sum(t == 1L)
    b <- sum(t == 2L)
    name_dims(sc$scale[[a + 1L]] * pt[[b + 1L]], s)
  })
  stats::setNames(out, param_tuple_names(s, order))
}


#' Log-Determinant Components of an Economical Parameter
#'
#' @description
#' Assembles one derivative order of
#' \eqn{\log\lvert M \rvert = p\log\sigma^2 + q(\rho)}.
#'
#' @details
#' The log-determinant is a sum of a function of one free value and a function
#' of the other, so every mixed component is exactly zero and the pure ones are
#' the two chains taken separately.
#'
#' @param s The parameter.
#' @param eta A numeric vector of two free values.
#' @param order The derivative order, 1 to 4.
#' @param terms The affine-logarithm terms of \eqn{q}; see
#'   \code{\link{log_affine_derivs}}.
#'
#' @return A named numeric vector.
#'
#' @keywords internal
econ_logdet_derivative <- function(s, eta, order, terms) {
  sc <- econ_scalars(s, eta)
  p <- s@dimension
  sd <- sc$scale
  rd <- sc$rho
  # p * log(scale), composed through the link
  a_log <- lapply(1:4, function(k) {
    (-1)^(k - 1L) * factorial(k - 1L) / sd[[1L]]^k
  })
  d_scale <- compose4(a_log, sd[-1L])
  d_rho <- compose4(log_affine_derivs(rd[[1L]], terms), rd[-1L])

  idx <- param_tuple_indices(s, order)
  out <- vapply(idx, function(t) {
    if (all(t == 1L)) return(p * d_scale[[order]])
    if (all(t == 2L)) return(d_rho[[order]])
    0
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, order))
}


#' The Pattern of a Compound-Symmetric Parameter
#'
#' @description
#' \eqn{P(\rho) = I + \rho(J - I)} and its four derivatives in the free value.
#' The pattern is linear in the correlation, so every derivative is the
#' matching derivative of the correlation times the constant \eqn{J - I}.
#'
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param sc The scalars of \code{\link{econ_scalars}}.
#'
#' @return A list of five matrices.
#'
#' @keywords internal
cs_pattern <- function(s, sc) {
  p <- s@dimension
  off <- matrix(1, p, p) - diag(p)
  r <- sc$rho
  c(list(diag(p) + r[[1L]] * off), lapply(2:5, function(k) r[[k]] * off))
}


#' @title Value of a Compound Symmetry Parameter
#' @name param_value.CompoundSymmetryParam
#' @description \eqn{\sigma^2\{(1-\rho)I + \rho J\}}.
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return A compound-symmetric positive definite matrix.
#' @keywords internal
S7::method(param_value, CompoundSymmetryParam) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta)
  name_dims(sc$scale[[1L]] * cs_pattern(s, sc)[[1L]], s)
}


#' @title Free Vector of a Compound Symmetry Parameter
#' @name param_free.CompoundSymmetryParam
#' @description
#' The scale is the common diagonal entry and the correlation is the common
#' off-diagonal one divided by it, both exact. A matrix whose diagonal is not
#' constant, or whose off-diagonal entries are not all equal, is refused
#' rather than averaged.
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param m A compound-symmetric matrix.
#' @param ... Unused.
#' @return A named numeric vector of two free values.
#' @keywords internal
S7::method(param_free, CompoundSymmetryParam) <- function(s, m, ...) {
  p <- s@dimension
  d <- diag(m)
  scl <- max(1, max(abs(m)))
  if (diff(range(d)) > 1e-8 * scl) {
    stop("'m' does not have a constant diagonal, so it is not compound symmetric.",
      call. = FALSE)
  }
  if (d[1L] <= 0) stop("'m' has a non-positive diagonal.", call. = FALSE)
  off <- m[lower.tri(m)]
  if (diff(range(off)) > 1e-8 * scl) {
    stop(paste0(
      "the off-diagonal entries of 'm' are not all equal, so it is not\n",
      "  compound symmetric. It is refused rather than averaged."
    ), call. = FALSE)
  }
  stats::setNames(
    c(linkfunctions7::linkfun(s@param_params$link_scale, d[1L]),
      linkfunctions7::linkfun(s@param_params$link_rho, off[1L] / d[1L])),
    s@free_names
  )
}


#' @title Solve of a Compound Symmetry Parameter
#' @name param_solve.CompoundSymmetryParam
#' @description
#' Exact, by Sherman-Morrison: the inverse of
#' \eqn{\sigma^2\{(1-\rho)I + \rho J\}} is
#' \eqn{\{\sigma^2(1-\rho)\}^{-1}[I - \rho J/\{1 + (p-1)\rho\}]}, compound
#' symmetric again, so no factorization is performed.
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param eta A numeric vector of two free values.
#' @param b A numeric matrix with \code{s@dimension} rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(param_solve, CompoundSymmetryParam) <- function(s, eta, b = NULL, ...) {
  p <- s@dimension
  sc <- econ_scalars(s, eta)
  v <- sc$scale[[1L]]
  r <- sc$rho[[1L]]
  inv <- (diag(p) - (r / (1 + (p - 1) * r)) * matrix(1, p, p)) / (v * (1 - r))
  inv %*% b
}


#' @title Derivatives of a Compound Symmetry Parameter
#' @name param_d1.CompoundSymmetryParam
#' @description
#' Closed form at every order: the value is the scale times a pattern linear
#' in the correlation, so a component with \eqn{a} scale indices and \eqn{b}
#' correlation indices is one derivative of each.
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 1L, cs_pattern)
}

#' @rdname param_d1.CompoundSymmetryParam
#' @name param_d2.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d2, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 2L, cs_pattern)
}

#' @rdname param_d1.CompoundSymmetryParam
#' @name param_d3.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d3, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 3L, cs_pattern)
}

#' @rdname param_d1.CompoundSymmetryParam
#' @name param_d4.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d4, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 4L, cs_pattern)
}


#' The Log-Determinant Terms of a Compound-Symmetric Parameter
#'
#' @description
#' The affine-logarithm terms of \eqn{q(\rho) = \log\{1 + (p-1)\rho\} +
#' (p-1)\log(1-\rho)}, the correlation's half of the log-determinant.
#'
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#'
#' @return A list of numeric triples.
#'
#' @keywords internal
cs_logdet_terms <- function(s) {
  p <- s@dimension
  list(c(1, 1, p - 1), c(p - 1, 1, -1))
}


#' @title Log-Determinant of a Compound Symmetry Parameter
#' @name param_logdet.CompoundSymmetryParam
#' @description
#' Closed form from the two distinct eigenvalues:
#' \eqn{p\log\sigma^2 + \log\{1+(p-1)\rho\} + (p-1)\log(1-\rho)}. No
#' factorization and no determinant is computed.
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta)
  p <- s@dimension
  r <- sc$rho[[1L]]
  p * log(sc$scale[[1L]]) + log(1 + (p - 1) * r) + (p - 1) * log(1 - r)
}

#' @title Log-Determinant Derivatives of a Compound Symmetry Parameter
#' @name param_dlogdet.CompoundSymmetryParam
#' @description
#' Closed form at every order. The log-determinant is a sum of a function of
#' the scale and a function of the correlation, so every mixed component is
#' exactly zero.
#' @param s A \code{\link{CompoundSymmetryParam}} object.
#' @param eta A numeric vector of two free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 1L, cs_logdet_terms(s))
}

#' @rdname param_dlogdet.CompoundSymmetryParam
#' @name param_d2logdet.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d2logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 2L, cs_logdet_terms(s))
}

#' @rdname param_dlogdet.CompoundSymmetryParam
#' @name param_d3logdet.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d3logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 3L, cs_logdet_terms(s))
}

#' @rdname param_dlogdet.CompoundSymmetryParam
#' @name param_d4logdet.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d4logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 4L, cs_logdet_terms(s))
}
