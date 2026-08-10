#' @include numerical_fallbacks.R
NULL

#' Distinct Blocks of a Matrix Parameter
#'
#' @description
#' The S7 class of a block-diagonal matrix built from several matrix
#' parameters, each carrying its own free values. Constructed by
#' \code{\link{block_diag}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{BlockDiagParam}.
#'
#' @seealso \code{\link{block_diag}}, \code{\link{kron_identity}}
#'
#' @examples
#' S7::S7_inherits(block_diag(log_cholesky(2), ar1(3)), BlockDiagParam)
#'
#' @export
BlockDiagParam <- S7::new_class("BlockDiagParam", parent = matrix_parameter)


#' Construct a Block Diagonal of Matrix Parameters
#'
#' @description
#' \eqn{M(\eta) = \mathrm{diag}(S_1(\eta_1), \ldots, S_B(\eta_B))} for matrix
#' parameters \eqn{S_1, \ldots, S_B}: a block-diagonal matrix whose blocks are
#' distinct families, each governed by its own stretch of the free vector.
#'
#' @details
#' This is the covariance (or precision) of several independent groups of
#' coefficients whose structures differ, which is what a model carrying more
#' than one random-effect term needs. It differs from
#' \code{\link{kron_identity}}, where the blocks are identical and share one
#' free vector; here the free vectors are concatenated, so the composite has
#' \eqn{\sum_b d_b} free values.
#'
#' Every quantity of the contract follows from the blocks without rederivation,
#' and the reason is that the free values of one block do not enter another:
#' \deqn{\partial_k M = \mathrm{diag}(0, \ldots, \partial_k S_b, \ldots, 0),
#'       \qquad \log\lvert M \rvert_{+} = \sum_b \log\lvert S_b \rvert_{+},
#'       \qquad M^{-1} = \mathrm{diag}(S_1^{-1}, \ldots, S_B^{-1}).}
#' A derivative whose indices do not all belong to one block is therefore
#' identically zero, at every order and for the log-determinant as well as for
#' the value. The rank is the sum of the blocks' ranks and the null basis is
#' their block diagonal, both read from the components rather than from an
#' assembled matrix.
#'
#' The free names are prefixed by the block's label, since two blocks of the
#' same family would otherwise report the same names and the class requires
#' them to be unique. Labels come from the names of the arguments where they
#' are given, and are \code{b1}, \code{b2}, ... otherwise.
#'
#' @param ... Two or more objects inheriting from
#'   \code{\link{matrix_parameter}}, or a single list of them. Named arguments
#'   supply the block labels.
#' @param role One of \code{"covariance"}, \code{"precision"} or
#'   \code{"either"}. Defaults to the blocks' common role, and to
#'   \code{"either"} when they disagree.
#'
#' @return An object of class \code{\link{BlockDiagParam}}.
#'
#' @seealso \code{\link{kron_identity}}, \code{\link{log_cholesky}}
#'
#' @examples
#' s <- block_diag(subject = log_cholesky(2), time = ar1(3))
#' c(dimension = s@dimension, n_free = s@n_free)
#' s@free_names
#' param_logdet(s, c(0, 0, 0, 0, 0))
#'
#' @export
block_diag <- function(..., role = NULL) {
  blocks <- list(...)
  if (length(blocks) == 1L && is.list(blocks[[1L]]) &&
    !S7::S7_inherits(blocks[[1L]], parameter)) {
    blocks <- blocks[[1L]]
  }
  if (length(blocks) < 1L) {
    stop("'block_diag' needs at least one block.", call. = FALSE)
  }
  ok <- vapply(blocks, S7::S7_inherits, logical(1), class = matrix_parameter)
  if (!all(ok)) {
    stop(sprintf("Every block must inherit from 'matrix_parameter'; block %s does not.",
                 paste(which(!ok), collapse = ", ")), call. = FALSE)
  }

  labels <- names(blocks)
  if (is.null(labels)) labels <- rep("", length(blocks))
  blank <- !nzchar(labels)
  labels[blank] <- paste0("b", seq_along(blocks))[blank]
  if (anyDuplicated(labels)) {
    stop("Block labels must be unique.", call. = FALSE)
  }

  dims <- vapply(blocks, function(b) b@dimension, integer(1))
  nfr <- vapply(blocks, function(b) b@n_free, integer(1))
  ranks <- vapply(blocks, function(b) b@rank, integer(1))
  p <- sum(dims)

  if (is.null(role)) {
    roles <- unique(vapply(blocks, function(b) b@role, character(1)))
    role <- if (length(roles) == 1L) roles else "either"
  }

  free_names <- unlist(Map(function(b, lab) paste0(lab, "_", b@free_names),
                           blocks, labels), use.names = FALSE)
  if (is.null(free_names)) free_names <- character(0)

  # the ranges each block occupies, in the matrix and in the free vector
  rows <- split_ranges(dims)
  free <- split_ranges(nfr)
  owner <- rep(seq_along(blocks), nfr)

  null_basis <- matrix(0, p, p - sum(ranks))
  col <- 0L
  for (j in seq_along(blocks)) {
    nb <- blocks[[j]]@null_basis
    if (ncol(nb)) {
      null_basis[rows[[j]], col + seq_len(ncol(nb))] <- nb
      col <- col + ncol(nb)
    }
  }

  BlockDiagParam(
    param_name = sprintf("blockdiag(%s)",
                         paste(vapply(blocks, function(b) b@param_name,
                                      character(1)), collapse = ", ")),
    dimension = as.integer(p),
    n_free = as.integer(sum(nfr)),
    free_names = free_names,
    rank = as.integer(sum(ranks)),
    null_basis = null_basis,
    role = role,
    param_params = list(blocks = blocks, labels = labels, rows = rows,
                        free = free, owner = owner)
  )
}


#' Consecutive Index Ranges of Given Widths
#'
#' @description
#' Turns the widths \eqn{n_1, \ldots, n_B} into the ranges they occupy when
#' laid end to end, as a list of integer vectors. A width of zero gives an
#' empty range rather than being dropped, so that the list stays aligned with
#' the blocks.
#'
#' @param widths An integer vector of widths.
#'
#' @return A list of integer vectors, one per width.
#'
#' @keywords internal
split_ranges <- function(widths) {
  ends <- cumsum(widths)
  starts <- ends - widths + 1L
  lapply(seq_along(widths), function(j) {
    if (widths[j] == 0L) integer(0) else seq.int(starts[j], ends[j])
  })
}


.bd <- function(s) s@param_params

#' The Blocks' Derivative Components, Keyed by Local Index Tuple
#'
#' @description
#' Fetches one block's derivatives of a given order and re-keys them by the
#' sorted local index tuple, so that a lookup from the composite's enumeration
#' needs no assumption about how the two orderings correspond.
#'
#' @param block A \code{\link{matrix_parameter}}.
#' @param eta The block's own stretch of the free vector.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of matrices, keyed as \code{"1"}, \code{"1,2"} and so
#'   on over the block's own free indices.
#'
#' @keywords internal
block_derivs_by_tuple <- function(block, eta, order) {
  d <- switch(order,
    param_d1(block, eta),
    param_d2(block, eta),
    param_d3(block, eta),
    param_d4(block, eta)
  )
  stats::setNames(d, vapply(param_tuple_indices(block, order),
                            function(t) paste(sort(t), collapse = ","),
                            character(1)))
}


#' Assemble a Block Diagonal's Derivatives of a Given Order
#'
#' @description
#' Places each block's own component in the rows and columns that block
#' occupies, and returns a zero matrix for any tuple whose indices are not all
#' owned by one block.
#'
#' @param s A \code{\link{BlockDiagParam}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of matrices keyed as \code{param_tuple_names(s, order)}.
#'
#' @keywords internal
block_diag_derivs <- function(s, eta, order) {
  b <- .bd(s)
  p <- s@dimension
  idx <- param_tuple_indices(s, order)
  zero <- matrix(0, p, p)
  cache <- vector("list", length(b$blocks))
  out <- lapply(idx, function(t) {
    own <- b$owner[t]
    j <- own[1L]
    if (any(own != j)) return(zero)
    if (is.null(cache[[j]])) {
      cache[[j]] <<- block_derivs_by_tuple(b$blocks[[j]],
                                           eta[b$free[[j]]], order)
    }
    local <- t - (b$free[[j]][1L] - 1L)
    v <- zero
    v[b$rows[[j]], b$rows[[j]]] <- cache[[j]][[paste(sort(local), collapse = ",")]]
    v
  })
  stats::setNames(out, param_tuple_names(s, order))
}


#' Assemble a Block Diagonal's Log-Determinant Derivatives
#'
#' @description
#' The log-determinant is the sum of the blocks', hence separable across
#' blocks, so a component is the owning block's own and zero for any tuple
#' spanning two blocks.
#'
#' @param s A \code{\link{BlockDiagParam}}.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named numeric vector keyed as \code{param_tuple_names(s, order)}.
#'
#' @keywords internal
block_diag_logdet_derivs <- function(s, eta, order) {
  b <- .bd(s)
  idx <- param_tuple_indices(s, order)
  cache <- vector("list", length(b$blocks))
  out <- vapply(idx, function(t) {
    own <- b$owner[t]
    j <- own[1L]
    if (any(own != j)) return(0)
    if (is.null(cache[[j]])) {
      e <- eta[b$free[[j]]]
      v <- switch(order,
        param_dlogdet(b$blocks[[j]], e),
        param_d2logdet(b$blocks[[j]], e),
        param_d3logdet(b$blocks[[j]], e),
        param_d4logdet(b$blocks[[j]], e)
      )
      cache[[j]] <<- stats::setNames(
        v, vapply(param_tuple_indices(b$blocks[[j]], order),
                  function(u) paste(sort(u), collapse = ","), character(1))
      )
    }
    local <- t - (b$free[[j]][1L] - 1L)
    cache[[j]][[paste(sort(local), collapse = ",")]]
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, order))
}


#' @title Value of a Block-Diagonal Parameter
#' @name param_value.BlockDiagParam
#' @description Assembles the blocks along the diagonal.
#' @param s A \code{\link{BlockDiagParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A symmetric matrix of side \code{s@dimension}.
#' @keywords internal
S7::method(param_value, BlockDiagParam) <- function(s, eta, ...) {
  b <- .bd(s)
  out <- matrix(0, s@dimension, s@dimension)
  for (j in seq_along(b$blocks)) {
    out[b$rows[[j]], b$rows[[j]]] <- param_value(b$blocks[[j]], eta[b$free[[j]]])
  }
  out
}

#' @title Free Vector of a Block-Diagonal Parameter
#' @name param_free.BlockDiagParam
#' @description
#' Inverts each diagonal block through its own parameter and concatenates the
#' results. A matrix whose off-diagonal blocks are not zero is rejected, that
#' being a matrix the family cannot represent.
#' @param s A \code{\link{BlockDiagParam}} object.
#' @param m A symmetric matrix of side \code{s@dimension}.
#' @param ... Unused.
#' @return A numeric vector of length \code{s@n_free}.
#' @keywords internal
S7::method(param_free, BlockDiagParam) <- function(s, m, ...) {
  b <- .bd(s)
  scale <- max(1, max(abs(m)))
  eta <- numeric(s@n_free)
  for (j in seq_along(b$blocks)) {
    r <- b$rows[[j]]
    eta[b$free[[j]]] <- param_free(b$blocks[[j]], m[r, r, drop = FALSE])
  }
  off <- m
  for (j in seq_along(b$blocks)) {
    r <- b$rows[[j]]
    off[r, r] <- 0
  }
  if (max(abs(off)) > 1e-8 * scale) {
    stop("'m' is not block diagonal in the blocks of this parameter.",
         call. = FALSE)
  }
  eta
}

#' @title Derivatives of a Block-Diagonal Parameter
#' @name param_d1.BlockDiagParam
#' @description
#' Each block's own derivatives, placed in the rows and columns that block
#' occupies. A component whose indices span two blocks is exactly zero, the
#' free values of one block not entering another.
#' @param s A \code{\link{BlockDiagParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, BlockDiagParam) <- function(s, eta, ...) {
  stats::setNames(block_diag_derivs(s, eta, 1L), s@free_names)
}

#' @rdname param_d1.BlockDiagParam
#' @name param_d2.BlockDiagParam
#' @keywords internal
S7::method(param_d2, BlockDiagParam) <- function(s, eta, ...) {
  block_diag_derivs(s, eta, 2L)
}

#' @rdname param_d1.BlockDiagParam
#' @name param_d3.BlockDiagParam
#' @keywords internal
S7::method(param_d3, BlockDiagParam) <- function(s, eta, ...) {
  block_diag_derivs(s, eta, 3L)
}

#' @rdname param_d1.BlockDiagParam
#' @name param_d4.BlockDiagParam
#' @keywords internal
S7::method(param_d4, BlockDiagParam) <- function(s, eta, ...) {
  block_diag_derivs(s, eta, 4L)
}

#' @title Log-Determinant of a Block-Diagonal Parameter
#' @name param_logdet.BlockDiagParam
#' @description
#' The sum of the blocks' log-(pseudo-)determinants, which is what the
#' determinant of a block-diagonal matrix is.
#' @param s A \code{\link{BlockDiagParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, BlockDiagParam) <- function(s, eta, ...) {
  b <- .bd(s)
  sum(vapply(seq_along(b$blocks),
             function(j) param_logdet(b$blocks[[j]], eta[b$free[[j]]]),
             numeric(1)))
}

#' @title Log-Determinant Derivatives of a Block-Diagonal Parameter
#' @name param_dlogdet.BlockDiagParam
#' @description
#' The blocks' own, concatenated at the first order and placed by owner above
#' it, every cross-block component being zero because the log-determinant is a
#' sum over the blocks.
#' @param s A \code{\link{BlockDiagParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, BlockDiagParam) <- function(s, eta, ...) {
  stats::setNames(block_diag_logdet_derivs(s, eta, 1L), s@free_names)
}

#' @rdname param_dlogdet.BlockDiagParam
#' @name param_d2logdet.BlockDiagParam
#' @keywords internal
S7::method(param_d2logdet, BlockDiagParam) <- function(s, eta, ...) {
  block_diag_logdet_derivs(s, eta, 2L)
}

#' @rdname param_dlogdet.BlockDiagParam
#' @name param_d3logdet.BlockDiagParam
#' @keywords internal
S7::method(param_d3logdet, BlockDiagParam) <- function(s, eta, ...) {
  block_diag_logdet_derivs(s, eta, 3L)
}

#' @rdname param_dlogdet.BlockDiagParam
#' @name param_d4logdet.BlockDiagParam
#' @keywords internal
S7::method(param_d4logdet, BlockDiagParam) <- function(s, eta, ...) {
  block_diag_logdet_derivs(s, eta, 4L)
}

#' @title Solve and Factor of a Block-Diagonal Parameter
#' @name param_solve.BlockDiagParam
#' @description
#' Both are blockwise: the inverse of a block-diagonal matrix is the block
#' diagonal of the inverses, and the same holds of a triangular factor.
#' @param s A \code{\link{BlockDiagParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param b A matrix with \code{s@dimension} rows, or \code{NULL}.
#' @param ... Unused.
#' @return A matrix.
#' @keywords internal
S7::method(param_solve, BlockDiagParam) <- function(s, eta, b = NULL, ...) {
  bd <- .bd(s)
  out <- b
  for (j in seq_along(bd$blocks)) {
    r <- bd$rows[[j]]
    out[r, ] <- param_solve(bd$blocks[[j]], eta[bd$free[[j]]],
                            b[r, , drop = FALSE])
  }
  out
}

#' @rdname param_solve.BlockDiagParam
#' @name param_factor.BlockDiagParam
#' @keywords internal
S7::method(param_factor, BlockDiagParam) <- function(s, eta, ...) {
  bd <- .bd(s)
  out <- matrix(0, s@dimension, s@dimension)
  for (j in seq_along(bd$blocks)) {
    r <- bd$rows[[j]]
    out[r, r] <- param_factor(bd$blocks[[j]], eta[bd$free[[j]]])
  }
  out
}
