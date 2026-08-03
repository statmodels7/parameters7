#' @include numerical_fallbacks.R
NULL


#' Unstructured Positive Definite Structure
#'
#' @description
#' The S7 class of unstructured symmetric positive definite matrices in the
#' log-Cholesky parametrisation. Constructed by \code{\link{log_cholesky}}.
#'
#' @inheritParams covstruct
#'
#' @return An object of class \code{LogCholeskyStruct}. Use
#'   \code{\link{log_cholesky}} rather than calling the class directly.
#'
#' @seealso \code{\link{log_cholesky}}
#'
#' @examples
#' S7::S7_inherits(log_cholesky(2), LogCholeskyStruct)
#'
#' @export
LogCholeskyStruct <- S7::new_class("LogCholeskyStruct", parent = covstruct)


#' Positions of the Free Values in the Lower Triangle
#'
#' @description
#' The row and column of each free value of a log-Cholesky structure, in the
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


#' Construct an Unstructured Positive Definite Structure
#'
#' @description
#' The log-Cholesky parametrisation of a symmetric positive definite matrix:
#' \eqn{M = L L^\top} with \eqn{L} lower triangular and positive on the
#' diagonal, the free values being the logarithms of the diagonal entries of
#' \eqn{L} and the entries below it.
#'
#' @details
#' The parametrisation is that of Pinheiro and Bates (1996), and it is the one
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
#' The log on the diagonal is intrinsic to the parametrisation and not a
#' swappable link, which is why it appears in the free names.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param role A label, one of \code{"either"} (the default),
#'   \code{"covariance"} or \code{"precision"}. Nothing computed depends on it;
#'   it records which side of a model the matrix parametrises, since the family
#'   name does not say.
#'
#' @return An object of class \code{\link{LogCholeskyStruct}}.
#'
#' @references
#' Pinheiro, J. C. and Bates, D. M. (1996). Unconstrained parametrizations for
#' variance-covariance matrices. \emph{Statistics and Computing} 6, 289-296.
#'
#' @seealso \code{\link{diag_struct}}, \code{\link{scaled_struct}},
#'   \code{\link{check_covstruct}}
#'
#' @examples
#' s <- log_cholesky(3)
#' s
#'
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#' round(struct_matrix(s, eta), 4)
#'
#' # the round trip closes exactly
#' max(abs(struct_free(s, struct_matrix(s, eta)) - eta))
#'
#' # and the log-determinant is linear in the free vector
#' struct_dlogdet(s, eta)
#'
#' @export
log_cholesky <- function(dimension, role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_struct_args(dimension, role)

  pos <- chol_positions(p)
  nm <- paste0("L", pos$row, ".", pos$col)
  nm[pos$on_diagonal] <- paste0("log_L", pos$row[pos$on_diagonal])

  LogCholeskyStruct(
    struct_name = "log_cholesky",
    dimension = p,
    n_free = length(nm),
    free_names = nm,
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    struct_params = list(positions = pos)
  )
}


#' The Cholesky Factor Behind a Free Vector
#'
#' @description
#' Assembles \eqn{L} from the free vector: the diagonal is the exponential of
#' the first \code{dimension} values, the rest are placed below it.
#'
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A lower triangular numeric matrix.
#'
#' @keywords internal
chol_assemble <- function(s, eta) {
  p <- s@dimension
  pos <- s@struct_params$positions
  l <- matrix(0, p, p)
  # Not ifelse(): it evaluates both branches over the whole vector, so the
  # exponential would run on the below-diagonal values as well. Harmless here
  # and not harmless in struct_free(), where the discarded branch is a
  # logarithm of numbers that are free to be negative.
  v <- eta
  v[pos$on_diagonal] <- exp(eta[pos$on_diagonal])
  l[cbind(pos$row, pos$col)] <- v
  l
}


#' @title Matrix of a Log-Cholesky Structure
#' @name struct_matrix.LogCholeskyStruct
#' @description \eqn{M = L L^\top}, with \eqn{L} assembled from the free vector.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A symmetric positive definite matrix.
#' @keywords internal
S7::method(struct_matrix, LogCholeskyStruct) <- function(s, eta, ...) {
  l <- chol_assemble(s, eta)
  name_dims(tcrossprod(l), s)
}


#' @title Factor of a Log-Cholesky Structure
#' @name struct_factor.LogCholeskyStruct
#' @description The factor is the parametrisation, so it needs no computing.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A lower triangular numeric matrix.
#' @keywords internal
S7::method(struct_factor, LogCholeskyStruct) <- function(s, eta, ...) {
  chol_assemble(s, eta)
}


#' @title Free Vector of a Log-Cholesky Structure
#' @name struct_free.LogCholeskyStruct
#' @description
#' The Cholesky factor of the matrix, with its diagonal logged. Exact, and the
#' inverse of \code{\link{struct_matrix}} because the factor with a positive
#' diagonal is unique.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param m A symmetric positive definite matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(struct_free, LogCholeskyStruct) <- function(s, m, ...) {
  l <- chol_pd(m)
  if (is.null(l)) {
    stop(paste0(
      "'m' is not positive definite, so it is not in the set log_cholesky()\n",
      "  parametrises. The verdict is spectral, not a failed factorisation."
    ), call. = FALSE)
  }
  pos <- s@struct_params$positions
  v <- l[cbind(pos$row, pos$col)]
  # The below-diagonal entries are free to be negative, so the logarithm is
  # taken on the diagonal alone rather than through ifelse(), which evaluates
  # both branches over the whole vector and warns about the NaNs it discards.
  v[pos$on_diagonal] <- log(v[pos$on_diagonal])
  stats::setNames(v, s@free_names)
}


#' @title First Derivatives of a Log-Cholesky Structure
#' @name struct_dmatrix.LogCholeskyStruct
#' @description
#' Closed form. Writing \eqn{L_k} for the derivative of the factor in the
#' \eqn{k}-th free value,
#' \eqn{\partial_k M = L_k L^\top + L L_k^\top}. The factor's derivative is
#' \eqn{L_{ii} E_{ii}} for a diagonal value, because the parametrisation is
#' its logarithm, and \eqn{E_{ij}} for a value below the diagonal.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(struct_dmatrix, LogCholeskyStruct) <- function(s, eta, ...) {
  p <- s@dimension
  pos <- s@struct_params$positions
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


#' @title Second Derivatives of a Log-Cholesky Structure
#' @name struct_d2matrix.LogCholeskyStruct
#' @description
#' Closed form. Differentiating \eqn{\partial_k M = L_k L^\top + L L_k^\top}
#' again gives
#' \eqn{\partial_{kl} M = L_{kl} L^\top + L_k L_l^\top + L_l L_k^\top +
#' L L_{kl}^\top}, and the factor's second derivative \eqn{L_{kl}} is
#' non-zero only when \eqn{k = l} is a diagonal value, where it is
#' \eqn{L_{ii} E_{ii}} again.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(struct_d2matrix, LogCholeskyStruct) <- function(s, eta, ...) {
  p <- s@dimension
  pos <- s@struct_params$positions
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

  idx <- struct_pair_indices(s)
  out <- vector("list", length(idx))
  names(out) <- struct_pair_names(s)
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


#' @title Log-Determinant of a Log-Cholesky Structure
#' @name struct_logdet.LogCholeskyStruct
#' @description
#' Closed form and linear in the free vector:
#' \eqn{\log|M| = 2 \sum_i \log L_{ii}}, which is twice the sum of the free
#' values on the diagonal. No factorisation and no determinant is computed.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(struct_logdet, LogCholeskyStruct) <- function(s, eta, ...) {
  2 * sum(eta[s@struct_params$positions$on_diagonal])
}


#' @title Log-Determinant Gradient of a Log-Cholesky Structure
#' @name struct_dlogdet.LogCholeskyStruct
#' @description
#' Closed form: 2 in each diagonal direction and 0 elsewhere, since the
#' log-determinant is linear in the free vector.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(struct_dlogdet, LogCholeskyStruct) <- function(s, eta, ...) {
  stats::setNames(
    ifelse(s@struct_params$positions$on_diagonal, 2, 0),
    s@free_names
  )
}


#' @title Log-Determinant Hessian of a Log-Cholesky Structure
#' @name struct_d2logdet.LogCholeskyStruct
#' @description Closed form: zero, the log-determinant being linear.
#' @param s A \code{\link{LogCholeskyStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector of zeros.
#' @keywords internal
S7::method(struct_d2logdet, LogCholeskyStruct) <- function(s, eta, ...) {
  nm <- struct_pair_names(s)
  stats::setNames(rep(0, length(nm)), nm)
}
