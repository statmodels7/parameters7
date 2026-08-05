#' @include numerical_fallbacks.R chain.R
NULL


#' Unstructured Positive Definite Parameter
#'
#' @description
#' The S7 class of unstructured symmetric positive definite matrices in the
#' log-Cholesky parametrization. Constructed by \code{\link{log_cholesky}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{LogCholeskyParam}. Use
#'   \code{\link{log_cholesky}} rather than calling the class directly.
#'
#' @seealso \code{\link{log_cholesky}}
#'
#' @examples
#' S7::S7_inherits(log_cholesky(2), LogCholeskyParam)
#'
#' @export
LogCholeskyParam <- S7::new_class("LogCholeskyParam", parent = matrix_parameter)


#' Positions of the Free Values in the Lower Triangle
#'
#' @description
#' The row and column of each free value of a log-Cholesky parameter, in the
#' order the free vector uses: the diagonal first, then the below-diagonal
#' entries column by column.
#'
#' @details
#' The ordering is fixed and is part of the contract, because
#' \code{free_names} depends on it and every consumer builds parameter tables
#' from those names.
#'
#' @param p The side of the matrix.
#'
#' @return A list with the integer vectors \code{row}, \code{col} and the
#'   logical \code{on_diagonal}.
#'
#' @keywords internal
chol_positions <- function(p) {
  row <- seq_len(p)
  col <- seq_len(p)
  on_diag <- rep(TRUE, p)
  if (p > 1L) {
    for (j in seq_len(p - 1L)) {
      for (i in seq.int(j + 1L, p)) {
        row <- c(row, i)
        col <- c(col, j)
        on_diag <- c(on_diag, FALSE)
      }
    }
  }
  list(row = as.integer(row), col = as.integer(col), on_diagonal = on_diag)
}


#' Construct an Unstructured Positive Definite Parameter
#'
#' @description
#' The log-Cholesky parametrization of a symmetric positive definite matrix:
#' \eqn{M = L L^\top} with \eqn{L} lower triangular and positive on the
#' diagonal, the free values being the logarithms of the diagonal entries of
#' \eqn{L} and the entries below it.
#'
#' @details
#' The parametrization is that of Pinheiro and Bates (1996), and it is the one
#' to reach for when nothing is known about the matrix. It is a bijection onto
#' the positive definite cone, smooth in both directions, with no boundary to
#' reach on the free scale, and unique because the Cholesky factor with a
#' positive diagonal is.
#'
#' The free vector runs the logarithms of the diagonal first and then the
#' below-diagonal entries column by column, so for \eqn{p = 3} it is
#' \eqn{(\log L_{11}, \log L_{22}, \log L_{33}, L_{21}, L_{31}, L_{32})}. The
#' ordering is fixed rather than incidental: \code{free_names} follows it and
#' every consumer builds its parameter tables from those names.
#'
#' The log-determinant is linear in the free vector,
#' \eqn{\log|M| = 2\sum_i \log L_{ii}}, so its gradient is the constant 2 in
#' the diagonal directions and 0 elsewhere, and its Hessian vanishes. Both are
#' supplied in closed form.
#'
#' The log on the diagonal is intrinsic to the parametrization and not a
#' swappable link, which is why it appears in the free names.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param role A label, one of \code{"either"} (the default),
#'   \code{"covariance"} or \code{"precision"}. Nothing computed depends on it;
#'   it records which side of a model the matrix parametrizes, since the family
#'   name does not say.
#'
#' @return An object of class \code{\link{LogCholeskyParam}}.
#'
#' @references
#' Pinheiro, J. C. and Bates, D. M. (1996). Unconstrained parametrizations for
#' variance-covariance matrices. \emph{Statistics and Computing} 6, 289-296.
#'
#' @seealso \code{\link{diagonal_matrix}}, \code{\link{scaled_matrix}},
#'   \code{\link{check_parameter}}
#'
#' @examples
#' s <- log_cholesky(3)
#' s
#'
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#' round(param_value(s, eta), 4)
#'
#' # the round trip closes exactly
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' # and the log-determinant is linear in the free vector
#' param_dlogdet(s, eta)
#'
#' @export
log_cholesky <- function(dimension, role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)

  pos <- chol_positions(p)
  nm <- paste0("L", pos$row, ".", pos$col)
  nm[pos$on_diagonal] <- paste0("log_L", pos$row[pos$on_diagonal])

  LogCholeskyParam(
    param_name = "log_cholesky",
    dimension = p,
    n_free = length(nm),
    free_names = nm,
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(positions = pos)
  )
}


#' The Cholesky Factor Behind a Free Vector
#'
#' @description
#' Assembles \eqn{L} from the free vector: the diagonal is the exponential of
#' the first \code{dimension} values, the rest are placed below it.
#'
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A lower triangular numeric matrix.
#'
#' @keywords internal
chol_assemble <- function(s, eta) {
  p <- s@dimension
  pos <- s@param_params$positions
  l <- matrix(0, p, p)
  # Not ifelse(): it evaluates both branches over the whole vector, so the
  # exponential would run on the below-diagonal values as well. Harmless here
  # and not harmless in param_free(), where the discarded branch is a
  # logarithm of numbers that are free to be negative.
  v <- eta
  v[pos$on_diagonal] <- exp(eta[pos$on_diagonal])
  l[cbind(pos$row, pos$col)] <- v
  l
}


#' @title Matrix of a Log-Cholesky Parameter
#' @name param_value.LogCholeskyParam
#' @description \eqn{M = L L^\top}, with \eqn{L} assembled from the free vector.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A symmetric positive definite matrix.
#' @keywords internal
S7::method(param_value, LogCholeskyParam) <- function(s, eta, ...) {
  l <- chol_assemble(s, eta)
  name_dims(tcrossprod(l), s)
}


#' @title Factor of a Log-Cholesky Parameter
#' @name param_factor.LogCholeskyParam
#' @description The factor is the parametrization, so it needs no computing.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A lower triangular numeric matrix.
#' @keywords internal
S7::method(param_factor, LogCholeskyParam) <- function(s, eta, ...) {
  chol_assemble(s, eta)
}


#' @title Free Vector of a Log-Cholesky Parameter
#' @name param_free.LogCholeskyParam
#' @description
#' The Cholesky factor of the matrix, with its diagonal logged. Exact, and the
#' inverse of \code{\link{param_value}} because the factor with a positive
#' diagonal is unique.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param m A symmetric positive definite matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(param_free, LogCholeskyParam) <- function(s, m, ...) {
  l <- chol_pd(m)
  if (is.null(l)) {
    stop(paste0(
      "'m' is not positive definite, so it is not in the set log_cholesky()\n",
      "  parametrizes. The verdict is spectral, not a failed factorization."
    ), call. = FALSE)
  }
  pos <- s@param_params$positions
  v <- l[cbind(pos$row, pos$col)]
  # The below-diagonal entries are free to be negative, so the logarithm is
  # taken on the diagonal alone rather than through ifelse(), which evaluates
  # both branches over the whole vector and warns about the NaNs it discards.
  v[pos$on_diagonal] <- log(v[pos$on_diagonal])
  stats::setNames(v, s@free_names)
}


#' @title First Derivatives of a Log-Cholesky Parameter
#' @name param_d1.LogCholeskyParam
#' @description
#' Closed form. Writing \eqn{L_k} for the derivative of the factor in the
#' \eqn{k}-th free value,
#' \eqn{\partial_k M = L_k L^\top + L L_k^\top}. The factor's derivative is
#' \eqn{L_{ii} E_{ii}} for a diagonal value, because the parametrization is
#' its logarithm, and \eqn{E_{ij}} for a value below the diagonal.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, LogCholeskyParam) <- function(s, eta, ...) {
  p <- s@dimension
  pos <- s@param_params$positions
  l <- chol_assemble(s, eta)

  out <- vector("list", s@n_free)
  names(out) <- s@free_names
  for (k in seq_len(s@n_free)) {
    lk <- matrix(0, p, p)
    lk[pos$row[k], pos$col[k]] <- if (pos$on_diagonal[k]) {
      l[pos$row[k], pos$col[k]]
    } else {
      1
    }
    a <- tcrossprod(lk, l)
    out[[k]] <- name_dims(a + t(a), s)
  }
  out
}


#' @title Second Derivatives of a Log-Cholesky Parameter
#' @name param_d2.LogCholeskyParam
#' @description
#' Closed form. Differentiating \eqn{\partial_k M = L_k L^\top + L L_k^\top}
#' again gives
#' \eqn{\partial_{kl} M = L_{kl} L^\top + L_k L_l^\top + L_l L_k^\top +
#' L L_{kl}^\top}, and the factor's second derivative \eqn{L_{kl}} is
#' non-zero only when \eqn{k = l} is a diagonal value, where it is
#' \eqn{L_{ii} E_{ii}} again.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d2, LogCholeskyParam) <- function(s, eta, ...) {
  p <- s@dimension
  pos <- s@param_params$positions
  l <- chol_assemble(s, eta)

  dfac <- lapply(seq_len(s@n_free), function(k) {
    lk <- matrix(0, p, p)
    lk[pos$row[k], pos$col[k]] <- if (pos$on_diagonal[k]) {
      l[pos$row[k], pos$col[k]]
    } else {
      1
    }
    lk
  })

  idx <- param_tuple_indices(s)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s)
  for (i in seq_along(idx)) {
    k <- idx[[i]][1L]
    l2 <- idx[[i]][2L]
    a <- tcrossprod(dfac[[k]], dfac[[l2]])
    m <- a + t(a)
    if (k == l2 && pos$on_diagonal[k]) {
      # The factor bends only in a diagonal direction, and only in its own.
      b <- tcrossprod(dfac[[k]], l)
      m <- m + b + t(b)
    }
    out[[i]] <- name_dims(m, s)
  }
  out
}


#' @title Log-Determinant of a Log-Cholesky Parameter
#' @name param_logdet.LogCholeskyParam
#' @description
#' Closed form and linear in the free vector:
#' \eqn{\log|M| = 2 \sum_i \log L_{ii}}, which is twice the sum of the free
#' values on the diagonal. No factorization and no determinant is computed.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, LogCholeskyParam) <- function(s, eta, ...) {
  2 * sum(eta[s@param_params$positions$on_diagonal])
}


#' @title Log-Determinant Gradient of a Log-Cholesky Parameter
#' @name param_dlogdet.LogCholeskyParam
#' @description
#' Closed form: 2 in each diagonal direction and 0 elsewhere, since the
#' log-determinant is linear in the free vector.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, LogCholeskyParam) <- function(s, eta, ...) {
  stats::setNames(
    ifelse(s@param_params$positions$on_diagonal, 2, 0),
    s@free_names
  )
}


#' @title Log-Determinant Hessian of a Log-Cholesky Parameter
#' @name param_d2logdet.LogCholeskyParam
#' @description Closed form: zero, the log-determinant being linear.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector of zeros.
#' @keywords internal
S7::method(param_d2logdet, LogCholeskyParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s)
  stats::setNames(rep(0, length(nm)), nm)
}


#' The Factor of a Log-Cholesky Parameter, and Its Derivatives
#'
#' @description
#' Returns \eqn{\partial^S L} for a multiset \eqn{S} of free-value indices, or
#' \code{NULL} when that derivative is identically zero.
#'
#' @details
#' A free value below the diagonal enters \eqn{L} linearly, so its second
#' derivative vanishes; a diagonal one enters through its exponential, so
#' every repeated derivative in that same value regenerates
#' \eqn{L_{ii} E_{ii}}. Everything else is zero, which is what makes the
#' Leibniz sum of \code{\link{chol_leibniz}} short.
#'
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param l The factor at the point, from \code{\link{chol_assemble}}.
#' @param ks A multiset of free-value indices, possibly empty; the empty one
#'   gives \eqn{L} itself.
#'
#' @return A numeric matrix, or \code{NULL}.
#'
#' @seealso \code{\link{chol_leibniz}}, \code{\link{leibniz_gram}}
#'
#' @keywords internal
chol_dfactor <- function(s, l, ks) {
  pos <- s@param_params$positions
  k1 <- ks[1L]
  if (length(ks) == 0L) return(l)
  if (length(ks) > 1L && (!pos$on_diagonal[k1] || any(ks != k1))) return(NULL)
  lk <- matrix(0, s@dimension, s@dimension)
  lk[pos$row[k1], pos$col[k1]] <- if (pos$on_diagonal[k1]) {
    l[pos$row[k1], pos$col[k1]]
  } else {
    1
  }
  lk
}

#' Derivative Components of a Log-Cholesky Parameter
#'
#' @description
#' Assembles one derivative order by the Leibniz rule on \eqn{M = LL^	op},
#' the factor's derivatives coming from \code{\link{chol_dfactor}}.
#'
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order.
#'
#' @return A named list of symmetric matrices.
#'
#' @seealso \code{\link{leibniz_gram}}
#'
#' @keywords internal
chol_leibniz <- function(s, eta, order) {
  l <- chol_assemble(s, eta)
  dfac <- function(ks) chol_dfactor(s, l, ks)
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(t) {
    name_dims(leibniz_gram(dfac, t, s@dimension), s)
  })
  stats::setNames(out, param_tuple_names(s, order))
}

#' @title Third Derivatives of a Log-Cholesky Parameter
#' @name param_d3.LogCholeskyParam
#' @description
#' Closed form by the Leibniz rule on \eqn{M = L L^\top}: each component
#' distributes its three differentiations over the two factors, and a factor
#' differentiated more than once survives only in a repeated diagonal
#' direction, where every derivative of \eqn{e^{\eta_k}} is itself.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d3, LogCholeskyParam) <- function(s, eta, ...) {
  chol_leibniz(s, eta, 3L)
}

#' @title Fourth Derivatives of a Log-Cholesky Parameter
#' @name param_d4.LogCholeskyParam
#' @description Closed form; see \code{\link{param_d3.LogCholeskyParam}}.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d4, LogCholeskyParam) <- function(s, eta, ...) {
  chol_leibniz(s, eta, 4L)
}

#' @title Higher Log-Determinant Derivatives of a Log-Cholesky Parameter
#' @name param_d3logdet.LogCholeskyParam
#' @description Closed form: zero, the log-determinant being linear in the
#'   free vector.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector of zeros.
#' @keywords internal
S7::method(param_d3logdet, LogCholeskyParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 3L)
  stats::setNames(rep(0, length(nm)), nm)
}

#' @title Fourth Log-Determinant Derivatives of a Log-Cholesky Parameter
#' @name param_d4logdet.LogCholeskyParam
#' @description Closed form: zero, as at third order.
#' @param s A \code{\link{LogCholeskyParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector of zeros.
#' @keywords internal
S7::method(param_d4logdet, LogCholeskyParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 4L)
  stats::setNames(rep(0, length(nm)), nm)
}
