#' @include covstruct_class.R
NULL


#' The Matrix a Structure Produces
#'
#' @description
#' Maps the free vector \eqn{\eta} to the matrix it parametrises.
#'
#' @details
#' This is the only generic a structure must implement. Everything else in the
#' package has a numerical method registered on the \code{\link{covstruct}}
#' class and is therefore available from this one alone.
#'
#' The generic validates the free vector before dispatching, so every method,
#' including one written outside the package, refuses a vector of the wrong
#' length or one containing a non-finite value. The free scale has no boundary
#' to reach, so a non-finite entry is a defect in the caller rather than a
#' point of the domain.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A symmetric numeric matrix with \code{s@dimension} rows and columns.
#'
#' @seealso \code{\link{struct_free}}, \code{\link{struct_dmatrix}},
#'   \code{\link{struct_logdet}}
#'
#' @examples
#' s <- log_cholesky(2)
#' round(struct_matrix(s, c(0, 0, 0.5)), 4)
#'
#' @export
struct_matrix <- S7::new_generic("struct_matrix", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' The Free Vector Behind a Matrix
#'
#' @description
#' The inverse map: given a matrix in the family's set, returns the free vector
#' that produces it.
#'
#' @details
#' Exact or refused, never obtained by optimisation. A generic
#' optimisation-based inverse would return a plausible \eqn{\eta} for a matrix
#' outside the set, which is a wrong answer wearing the shape of a right one.
#' A family whose map has no closed-form inverse refuses instead, and the base
#' class refuses on behalf of any structure that does not implement this.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param m A symmetric numeric matrix of the structure's dimension.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of length \code{s@n_free}, named by
#'   \code{s@free_names}.
#'
#' @seealso \code{\link{struct_matrix}}
#'
#' @examples
#' s <- log_cholesky(2)
#' eta <- c(0.3, -0.2, 0.5)
#' struct_free(s, struct_matrix(s, eta))
#'
#' @export
struct_free <- S7::new_generic("struct_free", "s", function(s, m, ...) {
  m <- check_matrix(s, m)
  S7::S7_dispatch()
})


#' First Derivatives of a Structure's Matrix
#'
#' @description
#' Returns \eqn{\partial M / \partial \eta_k} for every free value, as a list
#' of matrices.
#'
#' @details
#' A structure that registers no method for this generic gets the numerical one
#' of the \code{\link{covstruct}} class, which applies a single central
#' difference to \code{\link{struct_matrix}} in each component.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A list of \code{s@n_free} symmetric matrices, named by
#'   \code{s@free_names}.
#'
#' @seealso \code{\link{struct_d2matrix}}, \code{\link{struct_is_numerical}}
#'
#' @examples
#' s <- scaled_struct(diag(2))
#' struct_dmatrix(s, 0)
#'
#' @export
struct_dmatrix <- S7::new_generic("struct_dmatrix", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Second Derivatives of a Structure's Matrix
#'
#' @description
#' Returns the \eqn{d(d+1)/2} distinct second derivatives
#' \eqn{\partial^2 M / \partial \eta_k \partial \eta_l}.
#'
#' @details
#' The components are keyed by \code{\link{struct_pair_names}}, generated from
#' the same enumeration that produces the names rather than by taking a name
#' apart. Recovering an index by splitting a name on its separator is the
#' obvious route and it is wrong, because a free value whose own label contains
#' the separator splits into the wrong number of pieces.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A named list of symmetric matrices, keyed as
#'   \code{struct_pair_names(s)}.
#'
#' @seealso \code{\link{struct_dmatrix}}, \code{\link{struct_pair_names}}
#'
#' @examples
#' s <- scaled_struct(diag(2))
#' struct_d2matrix(s, 0)
#'
#' @export
struct_d2matrix <- S7::new_generic("struct_d2matrix", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Log-Determinant of a Structure's Matrix
#'
#' @description
#' Returns \eqn{\log|M|} when the family is of full rank, and the log
#' pseudo-determinant -- the sum of the logs of the non-zero eigenvalues --
#' when it is not.
#'
#' @details
#' One generic covers both because a consumer asks the same question of either:
#' what normalising constant does this matrix contribute. The object knows
#' which answer is the right one, through its declared rank.
#'
#' The sign in front of the result is the consumer's arithmetic. A gaussian
#' log-density written in the covariance carries \eqn{-\frac{1}{2}\log|\Sigma|}
#' and one written in the precision carries \eqn{+\frac{1}{2}\log|\Omega|};
#' the structure answers what the log-determinant is, and the likelihood
#' decides where it goes.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A single number.
#'
#' @seealso \code{\link{struct_dlogdet}}, \code{\link{struct_solve}}
#'
#' @examples
#' struct_logdet(log_cholesky(2), c(0, 0, 0))
#'
#' # a rank-deficient precision: the pseudo-determinant, over the 4 non-zero
#' # eigenvalues of a second-difference penalty on 6 coefficients
#' p <- crossprod(diff(diag(6), differences = 2))
#' s <- scaled_struct(p)
#' c(rank = s@rank, logdet = struct_logdet(s, 0))
#'
#' @export
struct_logdet <- S7::new_generic("struct_logdet", "s", function(s, eta, ...) {
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
#' \code{\link{struct_dmatrix}} through that identity, and
#' \code{\link{check_covstruct}} compares the two routes.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of length \code{s@n_free}, named by
#'   \code{s@free_names}.
#'
#' @seealso \code{\link{struct_logdet}}, \code{\link{struct_d2logdet}}
#'
#' @examples
#' # for a scaled precision the derivative is the rank, whatever the scale
#' s <- scaled_struct(crossprod(diff(diag(6), differences = 2)))
#' c(struct_dlogdet(s, -3), struct_dlogdet(s, 5), rank = s@rank)
#'
#' @export
struct_dlogdet <- S7::new_generic("struct_dlogdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Hessian of the Log-Determinant
#'
#' @description
#' Returns the distinct second derivatives of the log-determinant, or of the
#' log pseudo-determinant, keyed as \code{\link{struct_pair_names}}.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A named numeric vector, keyed as \code{struct_pair_names(s)}.
#'
#' @seealso \code{\link{struct_dlogdet}}
#'
#' @examples
#' struct_d2logdet(log_cholesky(2), c(0.2, -0.1, 0.4))
#'
#' @export
struct_d2logdet <- S7::new_generic("struct_d2logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Solve Through a Structure's Matrix
#'
#' @description
#' Returns \eqn{M^{-1} B}, computed through a factor rather than through an
#' explicit inverse.
#'
#' @details
#' A rank-deficient structure refuses rather than returning a pseudo-inverse.
#' What a consumer of an improper prior needs is the quadratic form and the
#' log pseudo-determinant -- the penalised normal equations invert
#' \eqn{X^\top X + \lambda P}, which is non-singular even when \eqn{P} is not,
#' and is assembled by the consumer -- so a pseudo-inverse would be a plausible
#' matrix answering a question nobody asked.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param b A numeric matrix or vector with \code{s@dimension} rows. Defaults to the
#'   identity, which returns the inverse.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with \code{s@dimension} rows.
#'
#' @seealso \code{\link{struct_factor}}, \code{\link{struct_logdet}}
#'
#' @examples
#' s <- log_cholesky(2)
#' round(struct_solve(s, c(0, 0, 0.5)), 4)
#'
#' @export
struct_solve <- S7::new_generic("struct_solve", "s", function(s, eta, b = NULL, ...) {
  eta <- check_eta(s, eta)
  if (s@rank < s@dimension) {
    stop(sprintf(paste0(
      "'%s' is rank deficient (%d of %d), so it has no inverse. A consumer of\n",
      "  an improper prior needs the quadratic form and the log\n",
      "  pseudo-determinant, not a pseudo-inverse: assemble the matrix the\n",
      "  model actually inverts and solve that."
    ), s@struct_name, s@rank, s@dimension), call. = FALSE)
  }
  if (is.null(b)) b <- diag(s@dimension)
  if (!is.matrix(b)) b <- as.matrix(b)
  if (nrow(b) != s@dimension) {
    stop(sprintf("'b' must have %d rows.", s@dimension), call. = FALSE)
  }
  S7::S7_dispatch()
})


#' A Factor of a Structure's Matrix
#'
#' @description
#' Returns a lower triangular \eqn{L} with \eqn{M = L L^\top}.
#'
#' @details
#' Refused for a rank-deficient family, for the reason given in
#' \code{\link{struct_solve}}.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Passed to methods.
#'
#' @return A lower triangular numeric matrix with \code{s@dimension} rows and
#'   columns.
#'
#' @seealso \code{\link{struct_solve}}
#'
#' @examples
#' round(struct_factor(log_cholesky(2), c(0.1, 0.2, -0.3)), 4)
#'
#' @export
struct_factor <- S7::new_generic("struct_factor", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  if (s@rank < s@dimension) {
    stop(sprintf(
      "'%s' is rank deficient (%d of %d), so it has no Cholesky factor.",
      s@struct_name, s@rank, s@dimension
    ), call. = FALSE)
  }
  S7::S7_dispatch()
})


#' Names of the Distinct Second-Derivative Components
#'
#' @description
#' The keys of \code{\link{struct_d2matrix}} and
#' \code{\link{struct_d2logdet}}: one per unordered pair of free values,
#' diagonal first and then the off-diagonal pairs in lexicographic order.
#'
#' @details
#' This exists so that nothing has to recover a pair by splitting a key apart.
#' Generating the keys and the index pairs from one enumeration cannot be
#' fooled by a free name that contains the separator, which splitting can.
#' Use \code{\link{struct_pair_indices}} for the pairs themselves.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#'
#' @return A character vector of length \eqn{d(d+1)/2}.
#'
#' @seealso \code{\link{struct_pair_indices}}
#'
#' @examples
#' struct_pair_names(log_cholesky(2))
#'
#' @export
struct_pair_names <- function(s) {
  idx <- struct_pair_indices(s)
  nm <- s@free_names
  vapply(idx, function(p) paste0(nm[p[1L]], ":", nm[p[2L]]), character(1))
}


#' Index Pairs Behind the Second-Derivative Names
#'
#' @description
#' The unordered index pairs \code{\link{struct_pair_names}} names, in exactly
#' the same order: the diagonal \eqn{(k, k)} first, then the off-diagonal
#' pairs \eqn{(k, l)} with \eqn{k < l} in lexicographic order.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#'
#' @return A list of integer vectors of length 2.
#'
#' @seealso \code{\link{struct_pair_names}}
#'
#' @examples
#' struct_pair_indices(log_cholesky(2))
#'
#' @export
struct_pair_indices <- function(s) {
  d <- s@n_free
  if (d == 0L) return(list())
  out <- lapply(seq_len(d), function(k) c(k, k))
  if (d > 1L) {
    for (k in seq_len(d - 1L)) {
      for (l in seq.int(k + 1L, d)) out[[length(out) + 1L]] <- c(k, l)
    }
  }
  lapply(out, as.integer)
}
