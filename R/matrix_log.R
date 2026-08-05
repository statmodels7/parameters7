#' @include numerical_fallbacks.R
NULL


#' Matrix Logarithm Parameter
#'
#' @description
#' The S7 class of unstructured symmetric positive definite matrices in the
#' matrix logarithm parametrization. Constructed by \code{\link{matrix_log}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{MatrixLogParam}.
#'
#' @seealso \code{\link{matrix_log}}
#'
#' @examples
#' S7::S7_inherits(matrix_log(2), MatrixLogParam)
#'
#' @export
MatrixLogParam <- S7::new_class("MatrixLogParam", parent = matrix_parameter)


#' Construct a Matrix Logarithm Parameter
#'
#' @description
#' The matrix logarithm parametrization of a symmetric positive definite
#' matrix: \eqn{M = \exp(S)} with \eqn{S} symmetric and genuinely free -- no
#' entry is transformed, the free values fill the lower triangle of \eqn{S}
#' directly, diagonal first and then below the diagonal column by column.
#'
#' @details
#' Two quantities are exact by construction and cost nothing. The
#' log-determinant is \eqn{\log|M| = \mathrm{tr}(S)}, the sum of the diagonal
#' free values, so it is linear and its higher derivatives vanish; and the
#' inverse is \eqn{M^{-1} = \exp(-S)}, through the same eigendecomposition
#' that evaluates the map, with no factorization.
#'
#' The derivatives of the map are the Frechet derivatives of the matrix
#' exponential, by the Daleckii-Krein representation: with
#' \eqn{S = Q \Lambda Q^\top} and directions rotated by \eqn{Q}, the
#' \eqn{k}-th derivative contracts the directions against divided differences
#' of \eqn{e^x} of order \eqn{k+1} at the eigenvalues, summed over the
#' orderings of the directions. The divided differences are computed by the
#' Opitz theorem -- the exponential of a small upper bidiagonal matrix, read
#' off its corner -- which stays exact under repeated and nearly repeated
#' eigenvalues, where the quotient formula cancels catastrophically.
#'
#' Next to \code{\link{log_cholesky}}: the log-Cholesky derivatives are
#' sparse products and cheaper, while here the log-determinant and the
#' inverse are the free quantities. The matrix logarithm is the chart to
#' reach for when both the covariance and the precision are wanted at once.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{MatrixLogParam}}.
#'
#' @references
#' Daleckii, J. L. and Krein, S. G. (1965). Integration and differentiation
#' of functions of Hermitian operators. \emph{American Mathematical Society
#' Translations} 47, 1-30.
#'
#' Opitz, G. (1964). Steigungsmatrizen. \emph{Zeitschrift fur Angewandte
#' Mathematik und Mechanik} 44, T52-T54.
#'
#' @seealso \code{\link{log_cholesky}}, \code{\link{param_value}}
#'
#' @examples
#' s <- matrix_log(2)
#' eta <- c(0.2, -0.3, 0.4)
#' round(param_value(s, eta), 4)
#'
#' # the log-determinant is the trace of S, linear in the free vector
#' c(param_logdet(s, eta), 0.2 - 0.3)
#'
#' # the round trip closes exactly
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' @export
matrix_log <- function(dimension, role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)
  pos <- chol_positions(p)
  nm <- paste0("S", pos$row, ".", pos$col)
  nm[pos$on_diagonal] <- paste0("S", pos$row[pos$on_diagonal])

  MatrixLogParam(
    param_name = "matrix_log",
    dimension = p,
    n_free = length(nm),
    free_names = nm,
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(positions = pos)
  )
}


#' The Symmetric Matrix Behind a Free Vector
#'
#' @description
#' Fills the lower triangle of \eqn{S} with the free values and mirrors it.
#'
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A symmetric numeric matrix.
#'
#' @keywords internal
mlog_s <- function(s, eta) {
  p <- s@dimension
  pos <- s@param_params$positions
  m <- matrix(0, p, p)
  m[cbind(pos$row, pos$col)] <- eta
  m[cbind(pos$col, pos$row)] <- eta
  m
}


#' The Basis Direction of One Free Value
#'
#' @description
#' The symmetric matrix \eqn{\partial S / \partial \eta_k}: a single diagonal
#' entry, or a symmetric pair below and above the diagonal.
#'
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param k The free-value index.
#'
#' @return A symmetric numeric matrix.
#'
#' @keywords internal
mlog_basis <- function(s, k) {
  p <- s@dimension
  pos <- s@param_params$positions
  b <- matrix(0, p, p)
  b[pos$row[k], pos$col[k]] <- 1
  b[pos$col[k], pos$row[k]] <- 1
  b
}


#' Exponential of a Small Upper Triangular Matrix
#'
#' @description
#' Scaling and squaring with a Taylor series, for the tiny bidiagonal
#' matrices the Opitz theorem builds. Deterministic and dependency free; the
#' matrices are at most 5 by 5.
#'
#' @param a A square numeric matrix.
#'
#' @return Its exponential.
#'
#' @keywords internal
mlog_expm_small <- function(a) {
  nrm <- max(abs(a))
  j <- max(0L, ceiling(log2(max(nrm, .Machine$double.eps) / 0.25)))
  b <- a / 2^j
  term <- diag(nrow(a))
  acc <- term
  for (i in 1:24) {
    term <- term %*% b / i
    acc <- acc + term
    if (max(abs(term)) < 1e-18 * max(1, max(abs(acc)))) break
  }
  for (i in seq_len(j)) acc <- acc %*% acc
  acc
}


#' Divided Differences of the Exponential
#'
#' @description
#' \eqn{e[\lambda_{i_1}, \ldots, \lambda_{i_m}]} by the Opitz theorem: the
#' exponential of the upper bidiagonal matrix with the arguments on the
#' diagonal and ones above it, read off the corner entry. Exact under
#' repeated arguments, where the recursive quotient cancels.
#'
#' @param lams The arguments, in any order.
#'
#' @return A single number.
#'
#' @keywords internal
dd_exp <- function(lams) {
  m <- length(lams)
  if (m == 1L) return(exp(lams))
  j <- diag(lams, nrow = m)
  j[cbind(seq_len(m - 1L), seq.int(2L, m))] <- 1
  mlog_expm_small(j)[1L, m]
}


#' The Eigendecomposition and Divided-Difference Tables at a Point
#'
#' @description
#' Everything the derivative contractions need: the eigendecomposition of
#' \eqn{S}, the rotated basis directions, and the divided-difference tables
#' of the orders asked for, computed once per free vector.
#'
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The highest derivative order wanted.
#'
#' @return A list with \code{q}, \code{lam}, \code{e} (rotated directions)
#'   and the tables \code{dd2} to \code{dd5} up to \code{order + 1} points.
#'
#' @keywords internal
mlog_tables <- function(s, eta, order) {
  p <- s@dimension
  es <- eigen(mlog_s(s, eta), symmetric = TRUE)
  q <- es$vectors
  lam <- es$values
  e <- lapply(seq_len(s@n_free), function(k) {
    crossprod(q, mlog_basis(s, k)) %*% q
  })

  out <- list(q = q, lam = lam, e = e)
  dd2 <- outer(seq_len(p), seq_len(p),
    Vectorize(function(i, j) dd_exp(lam[c(i, j)]))
  )
  out$dd2 <- dd2
  if (order >= 2L) {
    dd3 <- array(0, c(p, p, p))
    for (i in seq_len(p)) for (a in seq_len(p)) for (j in seq_len(p)) {
      dd3[i, a, j] <- dd_exp(lam[c(i, a, j)])
    }
    out$dd3 <- dd3
  }
  if (order >= 3L) {
    dd4 <- array(0, c(p, p, p, p))
    for (i in seq_len(p)) for (a in seq_len(p)) for (b in seq_len(p)) {
      for (j in seq_len(p)) dd4[i, a, b, j] <- dd_exp(lam[c(i, a, b, j)])
    }
    out$dd4 <- dd4
  }
  if (order >= 4L) {
    dd5 <- array(0, c(p, p, p, p, p))
    for (i in seq_len(p)) for (a in seq_len(p)) for (b in seq_len(p)) {
      for (cc in seq_len(p)) for (j in seq_len(p)) {
        dd5[i, a, b, cc, j] <- dd_exp(lam[c(i, a, b, cc, j)])
      }
    }
    out$dd5 <- dd5
  }
  out
}


#' One Derivative Component by the Daleckii-Krein Contraction
#'
#' @description
#' Contracts the rotated directions of one index tuple against the
#' divided-difference table of the matching order, summing over the
#' orderings of the directions, and rotates back.
#'
#' @param tb The tables of \code{\link{mlog_tables}}.
#' @param dirs The free-value indices of the tuple, possibly repeated.
#'
#' @return A symmetric numeric matrix.
#'
#' @keywords internal
mlog_contract <- function(tb, dirs) {
  p <- length(tb$lam)
  ord <- length(dirs)
  perms <- unique(combinat_perms(dirs))
  t_out <- matrix(0, p, p)

  for (pm in perms) {
    es <- tb$e[pm]
    if (ord == 1L) {
      t_out <- t_out + tb$dd2 * es[[1L]]
    } else if (ord == 2L) {
      for (i in seq_len(p)) for (j in seq_len(p)) {
        t_out[i, j] <- t_out[i, j] +
          sum(tb$dd3[i, , j] * es[[1L]][i, ] * es[[2L]][, j])
      }
    } else if (ord == 3L) {
      for (i in seq_len(p)) for (j in seq_len(p)) {
        acc <- 0
        for (a in seq_len(p)) {
          acc <- acc + es[[1L]][i, a] *
            sum(tb$dd4[i, a, , j] * es[[2L]][a, ] * es[[3L]][, j])
        }
        t_out[i, j] <- t_out[i, j] + acc
      }
    } else {
      for (i in seq_len(p)) for (j in seq_len(p)) {
        acc <- 0
        for (a in seq_len(p)) for (b in seq_len(p)) {
          acc <- acc + es[[1L]][i, a] * es[[2L]][a, b] *
            sum(tb$dd5[i, a, b, , j] * es[[3L]][b, ] * es[[4L]][, j])
        }
        t_out[i, j] <- t_out[i, j] + acc
      }
    }
  }
  # The multilinear form sums over ALL orderings of the directions; the loop
  # above runs the DISTINCT ones, and each of those occurs prod(mult!) times
  # among the k! permutations when directions repeat.
  t_out * prod(factorial(table(dirs)))
}


#' All Orderings of a Tuple
#'
#' @description
#' The distinct permutations of an index tuple, as a list.
#'
#' @param x An integer vector.
#'
#' @return A list of integer vectors.
#'
#' @keywords internal
combinat_perms <- function(x) {
  n <- length(x)
  if (n == 1L) return(list(x))
  out <- list()
  for (i in seq_len(n)) {
    rest <- combinat_perms(x[-i])
    for (r in rest) out[[length(out) + 1L]] <- c(x[i], r)
  }
  unique(out)
}


#' Third and Fourth Derivatives of a Matrix Exponential
#'
#' @description
#' The derivative components of orders three and four of \eqn{M = e^{S}} with
#' respect to the free entries of the symmetric matrix \eqn{S}.
#'
#' @details
#' The Frechet derivatives of the exponential contract chains of directions
#' against divided differences of \eqn{\exp} in the eigenvalues, and the sum
#' runs over every ordering of the directions rather than over the distinct
#' ones, so a tuple with repeated indices is counted with its multiplicity
#' rather than corrected for it afterwards. The divided differences come from
#' the Opitz representation, an exponential of a small bidiagonal matrix read
#' off its corner, which stays exact where the quotient recursion cancels
#' catastrophically under near-repeated eigenvalues.
#'
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order, 3 or 4.
#'
#' @return A named list of matrices, keyed by
#'   \code{\link{param_tuple_names}(s, order)}.
#'
#' @seealso \code{\link{matrix_log}}
#'
#' @keywords internal
mlog_higher <- function(s, eta, order) {
  tb <- mlog_tables(s, eta, order)
  idx <- param_tuple_indices(s, order)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, order)
  for (i in seq_along(idx)) {
    raw <- mlog_contract(tb, idx[[i]])
    m <- tb$q %*% raw %*% t(tb$q)
    out[[i]] <- name_dims((m + t(m)) / 2, s)
  }
  out
}


#' @title Value of a Matrix Logarithm Parameter
#' @name param_value.MatrixLogParam
#' @description \eqn{M = \exp(S)}, through the eigendecomposition of \eqn{S}.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A symmetric positive definite matrix.
#' @keywords internal
S7::method(param_value, MatrixLogParam) <- function(s, eta, ...) {
  es <- eigen(mlog_s(s, eta), symmetric = TRUE)
  name_dims(es$vectors %*% (exp(es$values) * t(es$vectors)), s)
}


#' @title Free Vector of a Matrix Logarithm Parameter
#' @name param_free.MatrixLogParam
#' @description
#' The matrix logarithm by eigendecomposition, exact for a symmetric positive
#' definite input and refused otherwise.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param m A symmetric positive definite matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(param_free, MatrixLogParam) <- function(s, m, ...) {
  es <- eigen(m, symmetric = TRUE)
  if (any(es$values <= 0)) {
    stop(paste0(
      "'m' is not positive definite, so it is not in the set matrix_log()\n",
      "  parametrizes."
    ), call. = FALSE)
  }
  sm <- es$vectors %*% (log(es$values) * t(es$vectors))
  pos <- s@param_params$positions
  stats::setNames(sm[cbind(pos$row, pos$col)], s@free_names)
}


#' @title First Derivatives of a Matrix Logarithm Parameter
#' @name param_d1.MatrixLogParam
#' @description
#' The Frechet derivative of the exponential by Daleckii-Krein: the rotated
#' direction weighted entrywise by first divided differences of \eqn{e^x} at
#' the eigenvalues.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, MatrixLogParam) <- function(s, eta, ...) {
  tb <- mlog_tables(s, eta, 1L)
  out <- lapply(seq_len(s@n_free), function(k) {
    m <- tb$q %*% (tb$dd2 * tb$e[[k]]) %*% t(tb$q)
    name_dims((m + t(m)) / 2, s)
  })
  stats::setNames(out, s@free_names)
}


#' @title Second Derivatives of a Matrix Logarithm Parameter
#' @name param_d2.MatrixLogParam
#' @description
#' The second Frechet derivative: chains of two rotated directions against
#' three-point divided differences, summed over the two orderings.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d2, MatrixLogParam) <- function(s, eta, ...) {
  mlog_higher(s, eta, 2L)
}


#' @title Third Derivatives of a Matrix Logarithm Parameter
#' @name param_d3.MatrixLogParam
#' @description Chains of three directions against four-point divided
#'   differences, over the six orderings.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d3, MatrixLogParam) <- function(s, eta, ...) {
  mlog_higher(s, eta, 3L)
}


#' @title Fourth Derivatives of a Matrix Logarithm Parameter
#' @name param_d4.MatrixLogParam
#' @description Chains of four directions against five-point divided
#'   differences, over the twenty-four orderings.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d4, MatrixLogParam) <- function(s, eta, ...) {
  mlog_higher(s, eta, 4L)
}


#' @title Log-Determinant of a Matrix Logarithm Parameter
#' @name param_logdet.MatrixLogParam
#' @description Closed form and linear: \eqn{\log|M| = \mathrm{tr}(S)}, the
#'   sum of the diagonal free values.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, MatrixLogParam) <- function(s, eta, ...) {
  sum(eta[s@param_params$positions$on_diagonal])
}

#' @title Log-Determinant Gradient of a Matrix Logarithm Parameter
#' @name param_dlogdet.MatrixLogParam
#' @description Closed form: 1 in each diagonal direction and 0 elsewhere.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, MatrixLogParam) <- function(s, eta, ...) {
  stats::setNames(
    ifelse(s@param_params$positions$on_diagonal, 1, 0), s@free_names
  )
}

#' @title Log-Determinant Hessian of a Matrix Logarithm Parameter
#' @name param_d2logdet.MatrixLogParam
#' @description Closed form: zero, the trace being linear.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector of zeros.
#' @keywords internal
S7::method(param_d2logdet, MatrixLogParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 2L)
  stats::setNames(rep(0, length(nm)), nm)
}

#' @title Higher Log-Determinant Derivatives of a Matrix Logarithm Parameter
#' @name param_d3logdet.MatrixLogParam
#' @description Closed form: zero at both orders, the trace being linear.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector of zeros.
#' @keywords internal
S7::method(param_d3logdet, MatrixLogParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 3L)
  stats::setNames(rep(0, length(nm)), nm)
}

#' @title Fourth Log-Determinant Derivatives of a Matrix Logarithm Parameter
#' @name param_d4logdet.MatrixLogParam
#' @description Closed form: zero.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector of zeros.
#' @keywords internal
S7::method(param_d4logdet, MatrixLogParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 4L)
  stats::setNames(rep(0, length(nm)), nm)
}

#' @title Solve of a Matrix Logarithm Parameter
#' @name param_solve.MatrixLogParam
#' @description Exact: \eqn{M^{-1} = \exp(-S)}, through the same
#'   eigendecomposition as the value, applied to \code{b}.
#' @param s A \code{\link{MatrixLogParam}} object.
#' @param eta A numeric vector of free values.
#' @param b A numeric matrix with \code{s@dimension} rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(param_solve, MatrixLogParam) <- function(s, eta, b = NULL, ...) {
  es <- eigen(mlog_s(s, eta), symmetric = TRUE)
  inv <- es$vectors %*% (exp(-es$values) * t(es$vectors))
  inv %*% b
}
