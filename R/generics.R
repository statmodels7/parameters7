#' @include parameter_class.R
NULL


#' The Value a Parameter Produces
#'
#' @description
#' Maps the free vector \eqn{\eta} to the constrained value it parametrizes:
#' a symmetric matrix for the \code{\link{matrix_parameter}} branch, a
#' probability vector for \code{\link{simplex}}, a row-stochastic matrix for
#' \code{\link{transition_matrix}}.
#'
#' @details
#' This is the only generic a parameter must implement. Everything else in the
#' package has a numerical method registered on the \code{\link{parameter}}
#' class and is therefore available from this one alone.
#'
#' The generic validates the free vector before dispatching, so every method,
#' including one written outside the package, rejects a vector of the wrong
#' length or one containing a non-finite value. The free scale has no boundary
#' to reach, so a non-finite entry is a defect in the caller rather than a
#' point of the domain.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return The constrained value, shaped as the family declares: a symmetric
#'   \code{s@dimension} by \code{s@dimension} matrix on the matrix branch, a
#'   vector or row-stochastic matrix otherwise.
#'
#' @seealso \code{\link{param_free}}, \code{\link{param_d1}},
#'   \code{\link{param_logdet}}
#'
#' @examples
#' s <- log_cholesky(2)
#' round(param_value(s, c(0, 0, 0.5)), 4)
#'
#' @export
param_value <- S7::new_generic("param_value", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' The Free Vector Behind a Value
#'
#' @description
#' The inverse map: given a value in the family's set, returns the free vector
#' that produces it.
#'
#' @details
#' Implemented exactly or not at all, never by optimization. An
#' optimization-based inverse would return a plausible \eqn{\eta} for a matrix
#' outside the set, indistinguishable from a correct answer.
#' A family whose map has no closed-form inverse signals an error, and the base
#' class does so on behalf of any parameter that does not implement this. A
#' value outside the set -- a matrix that is not positive definite, a vector
#' off the simplex -- is rejected rather than repaired.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param m A value of the family's shape: a symmetric matrix on the matrix
#'   branch, a probability vector or row-stochastic matrix otherwise.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of length \code{s@n_free}, named by
#'   \code{s@free_names}.
#'
#' @seealso \code{\link{param_value}}
#'
#' @examples
#' s <- log_cholesky(2)
#' eta <- c(0.3, -0.2, 0.5)
#' param_free(s, param_value(s, eta))
#'
#' @export
param_free <- S7::new_generic("param_free", "s", function(s, m, ...) {
  if (S7::S7_inherits(s, matrix_parameter)) m <- check_matrix(s, m)
  S7::S7_dispatch()
})


#' First Derivatives of a Parameter's Value
#'
#' @description
#' Returns \eqn{\partial V / \partial \eta_k} for every free value, as a
#' list of objects each shaped like the value.
#'
#' @details
#' A parameter that registers no method for this generic gets the numerical one
#' of the \code{\link{parameter}} class, which applies a single central
#' difference to \code{\link{param_value}} in each component.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A list of \code{s@n_free} symmetric matrices, named by
#'   \code{s@free_names}.
#'
#' @seealso \code{\link{param_d2}}, \code{\link{param_is_numerical}}
#'
#' @examples
#' s <- scaled_matrix(diag(2))
#' param_d1(s, 0)
#'
#' @export
param_d1 <- S7::new_generic("param_d1", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Second Derivatives of a Parameter's Value
#'
#' @description
#' Returns the \eqn{d(d+1)/2} distinct second derivatives
#' \eqn{\partial^2 V / \partial \eta_k \partial \eta_l}, each shaped like
#' the value.
#'
#' @details
#' The components are keyed by \code{\link{param_tuple_names}}, generated from
#' the same enumeration that produces the names rather than by taking a name
#' apart. Recovering an index by splitting a name on its separator is the
#' obvious route and it is wrong, because a free value whose own label contains
#' the separator splits into the wrong number of pieces.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A named list keyed as
#'   \code{param_tuple_names(s)}.
#'
#' @seealso \code{\link{param_d1}}, \code{\link{param_tuple_names}}
#'
#' @examples
#' s <- scaled_matrix(diag(2))
#' param_d2(s, 0)
#'
#' @export
param_d2 <- S7::new_generic("param_d2", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Log-Determinant of a Parameter's Matrix
#'
#' @description
#' Returns \eqn{\log|M|} when the family is of full rank, and the log
#' pseudo-determinant -- the sum of the logs of the non-zero eigenvalues --
#' when it is not.
#'
#' @details
#' One generic covers both because a consumer asks the same question of either:
#' what normalizing constant does this matrix contribute. The object records
#' which answer is the right one, through its declared rank.
#'
#' The sign in front of the result is the consumer's arithmetic. A gaussian
#' log-density written in the covariance carries \eqn{-\frac{1}{2}\log|\Sigma|}
#' and one written in the precision carries \eqn{+\frac{1}{2}\log|\Omega|};
#' the parameter answers what the log-determinant is, and the likelihood
#' decides where it goes.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A single number.
#'
#' @seealso \code{\link{param_dlogdet}}, \code{\link{param_solve}}
#'
#' @examples
#' param_logdet(log_cholesky(2), c(0, 0, 0))
#'
#' # a rank-deficient precision: the pseudo-determinant, over the 4 non-zero
#' # eigenvalues of a second-difference penalty on 6 coefficients
#' p <- crossprod(diff(diag(6), differences = 2))
#' s <- scaled_matrix(p)
#' c(rank = s@rank, logdet = param_logdet(s, 0))
#'
#' @export
param_logdet <- S7::new_generic("param_logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Gradient of the Log-Determinant
#'
#' @description
#' Returns \eqn{\partial \log|M| / \partial \eta_k} for every free value, or
#' the same derivative of the log pseudo-determinant when the family is rank
#' deficient.
#'
#' @details
#' The identity behind it is
#' \eqn{\partial_k \log|M| = \mathrm{tr}(M^{-1} \partial_k M)}, with the
#' Moore-Penrose inverse in place of \eqn{M^{-1}} in the rank-deficient case.
#' A closed form is therefore never an independent claim: it must agree with
#' \code{\link{param_d1}} through that identity, and
#' \code{\link{check_parameter}} compares the two routes.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of length \code{s@n_free}, named by
#'   \code{s@free_names}.
#'
#' @seealso \code{\link{param_logdet}}, \code{\link{param_d2logdet}}
#'
#' @examples
#' # for a scaled precision the derivative is the rank, whatever the scale
#' s <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
#' c(param_dlogdet(s, -3), param_dlogdet(s, 5), rank = s@rank)
#'
#' @export
param_dlogdet <- S7::new_generic("param_dlogdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Hessian of the Log-Determinant
#'
#' @description
#' Returns the distinct second derivatives of the log-determinant, or of the
#' log pseudo-determinant, keyed as \code{\link{param_tuple_names}}.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A named numeric vector, keyed as \code{param_tuple_names(s)}.
#'
#' @seealso \code{\link{param_dlogdet}}
#'
#' @examples
#' param_d2logdet(log_cholesky(2), c(0.2, -0.1, 0.4))
#'
#' @export
param_d2logdet <- S7::new_generic("param_d2logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Solve Through a Parameter's Matrix
#'
#' @description
#' Returns \eqn{M^{-1} B}, computed through a factor rather than through an
#' explicit inverse.
#'
#' @details
#' A rank-deficient parameter rejects rather than returning a pseudo-inverse.
#' What a consumer of an improper prior needs is the quadratic form and the
#' log pseudo-determinant -- the penalized normal equations invert
#' \eqn{X^\top X + \lambda P}, which is non-singular even when \eqn{P} is not,
#' and is assembled by the consumer -- so a pseudo-inverse would be a plausible
#' matrix answering a question nobody asked.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param b A numeric matrix or vector with \code{s@dimension} rows. Defaults to the
#'   identity, which returns the inverse.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with \code{s@dimension} rows.
#'
#' @seealso \code{\link{param_factor}}, \code{\link{param_logdet}}
#'
#' @examples
#' s <- log_cholesky(2)
#' round(param_solve(s, c(0, 0, 0.5)), 4)
#'
#' @export
param_solve <- S7::new_generic("param_solve", "s", function(s, eta, b = NULL, ...) {
  eta <- check_eta(s, eta)
  if (!S7::S7_inherits(s, matrix_parameter)) {
    stop(sprintf(
      "param_solve() is not defined for '%s': the value is not a symmetric matrix.",
      s@param_name
    ), call. = FALSE)
  }
  if (s@rank < s@dimension) {
    stop(sprintf(paste0(
      "'%s' is rank deficient (%d of %d), so it has no inverse. A consumer of\n",
      "  an improper prior needs the quadratic form and the log\n",
      "  pseudo-determinant, not a pseudo-inverse: assemble the matrix the\n",
      "  model actually inverts and solve that."
    ), s@param_name, s@rank, s@dimension), call. = FALSE)
  }
  if (is.null(b)) b <- diag(s@dimension)
  if (!is.matrix(b)) b <- as.matrix(b)
  if (nrow(b) != s@dimension) {
    stop(sprintf("'b' must have %d rows.", s@dimension), call. = FALSE)
  }
  S7::S7_dispatch()
})


#' A Factor of a Parameter's Matrix
#'
#' @description
#' Returns a lower triangular \eqn{L} with \eqn{M = L L^\top}.
#'
#' @details
#' Rejected for a rank-deficient family, for the reason given in
#' \code{\link{param_solve}}.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A lower triangular numeric matrix with \code{s@dimension} rows and
#'   columns.
#'
#' @seealso \code{\link{param_solve}}
#'
#' @examples
#' round(param_factor(log_cholesky(2), c(0.1, 0.2, -0.3)), 4)
#'
#' @export
param_factor <- S7::new_generic("param_factor", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  if (!S7::S7_inherits(s, matrix_parameter)) {
    stop(sprintf(
      "param_factor() is not defined for '%s': the value is not a symmetric matrix.",
      s@param_name
    ), call. = FALSE)
  }
  if (s@rank < s@dimension) {
    stop(sprintf(
      "'%s' is rank deficient (%d of %d), so it has no Cholesky factor.",
      s@param_name, s@rank, s@dimension
    ), call. = FALSE)
  }
  S7::S7_dispatch()
})


#' Third Derivatives of a Parameter's Value
#'
#' @description
#' Returns the distinct third derivatives
#' \eqn{\partial^3 V / \partial \eta_k \partial \eta_l \partial \eta_m},
#' keyed as \code{\link{param_tuple_names}(s, 3)}, each shaped like the
#' value.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A named list, keyed as \code{param_tuple_names(s, 3)}.
#'
#' @seealso \code{\link{param_d2}}, \code{\link{param_d4}}
#'
#' @examples
#' param_d3(scalar_matrix(2), 0.3)
#'
#' @export
param_d3 <- S7::new_generic("param_d3", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Fourth Derivatives of a Parameter's Value
#'
#' @description
#' Returns the distinct fourth derivatives, keyed as
#' \code{\link{param_tuple_names}(s, 4)}, each shaped like the value.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A named list, keyed as \code{param_tuple_names(s, 4)}.
#'
#' @seealso \code{\link{param_d3}}
#'
#' @examples
#' param_d4(scalar_matrix(2), 0.3)
#'
#' @export
param_d4 <- S7::new_generic("param_d4", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Third and Fourth Derivatives of the Log-Determinant
#'
#' @description
#' The higher derivatives of the log-(pseudo-)determinant, keyed as
#' \code{\link{param_tuple_names}} of the matching order.
#'
#' @param s An object inheriting from class \code{\link{matrix_parameter}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A named numeric vector.
#'
#' @seealso \code{\link{param_d2logdet}}
#'
#' @examples
#' param_d3logdet(log_cholesky(2), c(0.1, -0.2, 0.4))
#'
#' @export
param_d3logdet <- S7::new_generic("param_d3logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})

#' @rdname param_d3logdet
#' @export
param_d4logdet <- S7::new_generic("param_d4logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Names of the Distinct Derivative Components
#'
#' @description
#' The keys of the derivative lists: one per unordered tuple of free values.
#' At order 2, the keys of \code{\link{param_d2}} and
#' \code{\link{param_d2logdet}}, diagonal pairs first; at orders 3 and 4,
#' the keys of \code{\link{param_d3}} and \code{\link{param_d4}},
#' lexicographic combinations with repetition.
#'
#' @details
#' This exists so that nothing has to recover a tuple by splitting a key
#' apart. Generating the keys and the index tuples from one enumeration
#' cannot be fooled by a free name that contains the separator, which
#' splitting can. Use \code{\link{param_tuple_indices}} for the tuples
#' themselves.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param order The derivative order: 2 (default), 3 or 4.
#'
#' @return A character vector, one entry per distinct component.
#'
#' @seealso \code{\link{param_tuple_indices}}
#'
#' @examples
#' param_tuple_names(log_cholesky(2))
#' param_tuple_names(scalar_matrix(2), 3)
#'
#' @export
param_tuple_names <- function(s, order = 2L) {
  idx <- param_tuple_indices(s, order)
  nm <- s@free_names
  vapply(idx, function(t) paste(nm[t], collapse = ":"), character(1))
}


#' Index Tuples Behind the Derivative Component Names
#'
#' @description
#' The unordered index tuples \code{\link{param_tuple_names}} names, in
#' exactly the same order. At order 2 the diagonal pairs come first and then
#' the off-diagonal ones, because consumers index a Hessian that way; at
#' orders 3 and 4 the tuples are the lexicographic combinations with
#' repetition, matching the enumeration \pkg{distributions7} uses for its
#' higher derivatives.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param order The derivative order: 1 to 4.
#'
#' @return A list of integer vectors of length \code{order}.
#'
#' @seealso \code{\link{param_tuple_names}}
#'
#' @examples
#' param_tuple_indices(log_cholesky(2), 2)
#' length(param_tuple_indices(log_cholesky(2), 3))
#'
#' @export
param_tuple_indices <- function(s, order = 2L) {
  tuple_indices(s@n_free, order)
}


#' The Index Tuples of a Given Width
#'
#' @description
#' The enumeration behind \code{\link{param_tuple_indices}}, taken over a
#' number of variables rather than over a parameter, so that anything holding
#' derivatives over \eqn{d} variables -- a jet, for instance -- can share it.
#'
#' @param d The number of variables.
#' @param order The derivative order, 1 to 4.
#'
#' @return A list of integer vectors of length \code{order}.
#'
#' @seealso \code{\link{param_tuple_indices}}
#'
#' @keywords internal
tuple_indices <- function(d, order = 2L) {
  # One enumeration for the whole toolkit: delegating is what makes it
  # impossible for this copy to disagree with the one a jet is keyed by.
  numericals7::tuple_indices(d, order)
}
