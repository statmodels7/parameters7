#' @include numerical_fallbacks.R
NULL

#' Identical Blocks of a Matrix Parameter
#'
#' @description
#' The S7 class of a block-diagonal matrix built from \eqn{m} identical
#' copies of an inner matrix parameter. Constructed by
#' \code{\link{kron_identity}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{KronIdentityParam}.
#'
#' @seealso \code{\link{kron_identity}}
#'
#' @examples
#' S7::S7_inherits(kron_identity(log_cholesky(2), 3), KronIdentityParam)
#'
#' @export
KronIdentityParam <- S7::new_class("KronIdentityParam", parent = matrix_parameter)


#' Construct a Block Replication of a Matrix Parameter
#'
#' @description
#' \eqn{M(\eta) = I_m \otimes S(\eta)} for an inner matrix parameter
#' \eqn{S} of dimension \eqn{d}: a block-diagonal matrix of dimension
#' \eqn{md} whose \eqn{m} diagonal blocks are the same matrix, sharing one
#' free vector.
#'
#' @details
#' This is the covariance (or precision) of \eqn{m} independent groups
#' whose within-group structure is common -- the shape a grouped
#' random-effect term needs, with \eqn{S} the per-group matrix over the
#' coefficients of one group. Every quantity of the contract is a linear
#' lift of the inner parameter's:
#' \deqn{\partial_k M = I_m \otimes \partial_k S, \qquad
#'       \log\lvert M \rvert_{+} = m \log\lvert S \rvert_{+}, \qquad
#'       M^{-1} = I_m \otimes S^{-1},}
#' so nothing is rederived, and the rank is \eqn{m} times the inner rank
#' with the null basis replicated blockwise.
#'
#' The free names are the inner parameter's own: the blocks share their
#' values by construction, which is what distinguishes this composition
#' from a general block-diagonal of distinct structures.
#'
#' @param structure An object inheriting from
#'   \code{\link{matrix_parameter}}: the per-block parameter.
#' @param m The number of blocks, a single integer of at least 1.
#'
#' @return An object of class \code{\link{KronIdentityParam}}.
#'
#' @seealso \code{\link{log_cholesky}}, \code{\link{diagonal_matrix}}
#'
#' @examples
#' s <- kron_identity(log_cholesky(2), 3)
#' c(dimension = s@dimension, n_free = s@n_free)
#' param_logdet(s, c(0.1, 0, -0.2))
#'
#' @export
kron_identity <- function(structure, m) {
  if (!S7::S7_inherits(structure, matrix_parameter)) {
    stop("'structure' must inherit from 'matrix_parameter'.", call. = FALSE)
  }
  if (length(m) != 1L || is.na(m) || m < 1 || m != round(m)) {
    stop("'m' must be a single integer of at least 1.", call. = FALSE)
  }
  m <- as.integer(m)
  KronIdentityParam(
    param_name = sprintf("kron(I%d, %s)", m, structure@param_name),
    dimension = m * structure@dimension,
    n_free = structure@n_free,
    free_names = structure@free_names,
    rank = m * structure@rank,
    null_basis = kronecker(diag(m), structure@null_basis),
    role = structure@role,
    param_params = list(inner = structure, m = m)
  )
}

.kron_inner <- function(s) s@param_params$inner
.kron_m <- function(s) s@param_params$m
.kron_lift <- function(m, x) kronecker(diag(m), x)

S7::method(param_value, KronIdentityParam) <- function(s, eta, ...) {
  .kron_lift(.kron_m(s), param_value(.kron_inner(s), eta))
}

S7::method(param_free, KronIdentityParam) <- function(s, m, ...) {
  inner <- .kron_inner(s)
  d <- inner@dimension
  blocks <- .kron_m(s)
  b1 <- m[seq_len(d), seq_len(d), drop = FALSE]
  eta <- param_free(inner, b1)
  rebuilt <- .kron_lift(blocks, param_value(inner, eta))
  if (max(abs(rebuilt - m)) > 1e-8 * max(1, max(abs(m)))) {
    stop(paste("'m' is not m identical diagonal copies of a matrix the inner",
               "parameter can represent."), call. = FALSE)
  }
  eta
}

S7::method(param_d1, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d1(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

S7::method(param_d2, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d2(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

S7::method(param_d3, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d3(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

S7::method(param_d4, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d4(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

S7::method(param_logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_logdet(.kron_inner(s), eta)
}

S7::method(param_dlogdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_dlogdet(.kron_inner(s), eta)
}

S7::method(param_d2logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_d2logdet(.kron_inner(s), eta)
}

S7::method(param_d3logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_d3logdet(.kron_inner(s), eta)
}

S7::method(param_d4logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_d4logdet(.kron_inner(s), eta)
}

S7::method(param_solve, KronIdentityParam) <- function(s, eta, b = NULL, ...) {
  inner <- .kron_inner(s)
  d <- inner@dimension
  out <- b
  for (j in seq_len(.kron_m(s))) {
    rows <- (j - 1L) * d + seq_len(d)
    out[rows, ] <- param_solve(inner, eta, b[rows, , drop = FALSE])
  }
  out
}

S7::method(param_factor, KronIdentityParam) <- function(s, eta, ...) {
  .kron_lift(.kron_m(s), param_factor(.kron_inner(s), eta))
}
