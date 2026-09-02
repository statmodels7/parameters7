#' @include block_diag.R
NULL

#' The Inverse of a Matrix Parameter
#'
#' @description
#' The S7 class of the family whose value is the inverse of another family's.
#' [inverse_of()] builds one. It carries the inner family's free vector
#' unchanged, so its coordinates are read as those of the matrix it inverts.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `InverseParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `inner`, the family being inverted.
#'
#' @seealso [inverse_of()], the constructor, and [ar1_inv()] and
#'   [autoregressive_inv()], the two families whose inverse has a structure of
#'   its own and which are written out rather than composed.
#'
#' @examples
#' # The value is the inner family's inverse, and the free vector is the
#' # inner family's own.
#' s <- inverse_of(correlation_matrix(3))
#' s@free_names
#' eta <- c(0.4, -0.3, 0.8)
#' max(abs(param_value(s, eta) %*% param_value(correlation_matrix(3), eta) -
#'         diag(3)))
#'
#' @export
InverseParam <- S7::new_class("InverseParam", parent = matrix_parameter)


#' Construct the Inverse of a Matrix Parameter
#'
#' @description
#' Returns the family \eqn{N(\eta) = S(\eta)^{-1}} for a matrix family \eqn{S}:
#' the same free vector, the inverted value. It is the composition that says
#' "the other side" to a consumer that fixes one, as
#' [penalties7::structured_penalty()] fixes the precision.
#'
#' @details
#' # What it is for
#'
#' A family and its inverse are the same set of matrices only when the family
#' is closed under inversion, which [log_cholesky()], [matrix_log()],
#' [diagonal_matrix()], [scalar_matrix()], [compound_symmetry()] and
#' [dr_prod()] are and [correlation_matrix()], [ar1()] and [autoregressive()]
#' are not. For a family that is not closed the two sides are different models,
#' and this wrapper is how the one the family does not name is written:
#' `inverse_of(correlation_matrix(3))` is the matrix whose inverse is a
#' correlation matrix, which no chart of the package produces directly.
#'
#' Where the inverse has a structure worth exploiting the package writes it out
#' instead: [ar1_inv()] and [autoregressive_inv()] give the same values as
#' `inverse_of(ar1())` and `inverse_of(autoregressive())` and reach them in
#' \eqn{O(p)} entries rather than through a factorization.
#'
#' # The derivatives
#'
#' Differentiating \eqn{N = S^{-1}} repeatedly gives a sum over the **ordered**
#' set partitions of the differentiated positions,
#'
#' \deqn{\partial_I N = \sum_{(B_1, \ldots, B_q)} (-1)^q\,
#'   N\,(\partial_{B_1} S)\,N \cdots N\,(\partial_{B_q} S)\,N,}
#'
#' the blocks being ordered because matrices do not commute. At first order it
#' is \eqn{-N (\partial_k S) N} and at second
#' \eqn{N(\partial_k S\,N\,\partial_l S + \partial_l S\,N\,\partial_k S)N -
#' N(\partial_{kl}S)N}, which are the two expressions written by hand
#' elsewhere in the toolkit. The number of terms is the Fubini number of the
#' order: 1, 3, 13 and 75.
#'
#' Every factor is the inner family's own closed form, so nothing is
#' differenced: the arithmetic is exact wherever the inner family's is.
#'
#' # The log-determinant
#'
#' \eqn{\log\lvert N\rvert = -\log\lvert S\rvert}, so the value and all four
#' derivative orders are the inner family's negated, and no determinant of the
#' inverted matrix is taken.
#'
#' # Rank
#'
#' The inner family must be full rank. A singular matrix has no inverse, so
#' there is no family to build, and the constructor rejects it naming the rank
#' and the dimension.
#'
#' @section Notation:
#' \eqn{S} is the inner family's map, \eqn{N = S^{-1}} the one this builds,
#' \eqn{I} a multiset of differentiated positions and \eqn{(B_1,\ldots,B_q)} an
#' ordered partition of it into non-empty blocks.
#'
#' @param structure An object inheriting from [matrix_parameter()], of full
#'   rank.
#'
#' @return An object of class [InverseParam()], with `dimension`, `rank`,
#'   `n_free` and `free_names` the inner family's, `null_basis` a `dimension`
#'   by 0 matrix, and `param_name` `"inverse_of(inner)"`.
#'
#' @seealso [ar1_inv()] and [autoregressive_inv()] for the two written-out
#'   inverses, [block_diag()], [kron_identity()], [dr_prod()] and
#'   [sum_struct()] for the other compositions, and [param_solve()], which a
#'   family with a closed inverse implements and which this reads.
#'
#' @examples
#' # The inverse of an AR(1) correlation is tridiagonal, so the family whose
#' # value is that inverse is a different model from ar1() itself.
#' s <- inverse_of(ar1(5))
#' eta <- c(log(2), atanh(0.6))
#' round(param_value(s, eta), 4)
#'
#' # It is the inner family's inverse, exactly.
#' max(abs(param_value(s, eta) - solve(param_value(ar1(5), eta))))
#'
#' # The log-determinant is the inner one negated.
#' c(inverse = param_logdet(s, eta), inner = param_logdet(ar1(5), eta))
#'
#' # The first derivative is -N (d S) N.
#' N <- param_value(s, eta)
#' A <- param_d1(ar1(5), eta)
#' max(abs(param_d1(s, eta)[[1]] - (-N %*% A[[1]] %*% N)))
#'
#' # The free vector is the inner family's, so the round trip closes on it.
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' @export
inverse_of <- function(structure) {
  if (!S7::S7_inherits(structure, matrix_parameter)) {
    stop("'structure' must inherit from 'matrix_parameter'.", call. = FALSE)
  }
  if (structure@rank < structure@dimension) {
    stop(sprintf(paste0(
      "'%s' has rank %d of %d, so it has no inverse and there is no family to\n",
      "  build. A rank-deficient structure describes a matrix that is singular\n",
      "  at every free value, and inverting it is not a reparametrization but\n",
      "  an operation that does not exist."
    ), structure@param_name, structure@rank, structure@dimension),
    call. = FALSE)
  }
  p <- structure@dimension
  InverseParam(
    param_name = sprintf("inverse_of(%s)", structure@param_name),
    n_free = structure@n_free,
    free_names = structure@free_names,
    dimension = p,
    rank = p,
    null_basis = matrix(numeric(0), p, 0L),
    param_params = list(inner = structure)
  )
}


#' The Inner Family of an Inverse Parameter
#'
#' @description
#' Reads `param_params$inner`, written once so the methods do not repeat the
#' accessor.
#'
#' @param s An [InverseParam()] object.
#'
#' @return The inner [matrix_parameter()].
#'
#' @keywords internal
.inv_inner <- function(s) s@param_params$inner


#' Ordered Set Partitions of the First k Positions
#'
#' @description
#' Every way of splitting `1:k` into non-empty blocks **and ordering the
#' blocks**, which is what the derivative of an inverse sums over: the factors
#' are matrices and do not commute, so two orderings of the same split are two
#' terms.
#'
#' @details
#' Built from `numericals7::set_partitions()`, the toolkit's single copy of the
#' unordered enumeration, by taking every permutation of each partition's
#' blocks. The count is the Fubini number: 1, 3, 13, 75 for `k` of 1 to 4.
#'
#' @param k A single positive integer, at most 4 in this package's use.
#'
#' @return A list of lists of integer vectors, each inner list one ordered
#'   partition.
#'
#' @seealso [inverse_of()], the only caller.
#'
#' @keywords internal
ordered_set_partitions <- function(k) {
  perms <- function(n) {
    if (n == 1L) return(list(1L))
    out <- list()
    for (i in seq_len(n)) {
      for (rest in perms(n - 1L)) {
        r <- rest + (rest >= i)
        out[[length(out) + 1L]] <- c(i, r)
      }
    }
    out
  }
  out <- list()
  for (pt in numericals7::set_partitions(k)) {
    q <- length(pt)
    for (o in perms(q)) out[[length(out) + 1L]] <- pt[o]
  }
  out
}


#' Derivatives of an Inverse Parameter of a Given Order
#'
#' @description
#' The ordered-block-partition sum of [inverse_of()], run over the enumeration
#' the class is keyed by.
#'
#' @details
#' The inner family's derivative arrays of every order up to `order` are
#' fetched once and keyed by their sorted index tuple, so a block appearing in
#' many partitions is read and not recomputed.
#'
#' @param s An [InverseParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of `s@dimension` square matrices, keyed and ordered as
#'   [param_tuple_names()] says.
#'
#' @seealso [inverse_of()] for the formula.
#'
#' @keywords internal
inverse_derivs <- function(s, eta, order) {
  inner <- .inv_inner(s)
  N <- unclass(param_value(s, eta))
  # every order up to `order`, keyed by sorted index tuple: a partition block
  # of size m asks for the inner family's order-m component
  tabs <- lapply(seq_len(order), function(m) block_derivs_by_tuple(inner, eta, m))
  parts <- ordered_set_partitions(order)
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(t) {
    acc <- matrix(0, s@dimension, s@dimension)
    for (pt in parts) {
      q <- length(pt)
      term <- N
      for (B in pt) {
        key <- paste(sort(t[B]), collapse = ",")
        term <- term %*% unclass(tabs[[length(B)]][[key]]) %*% N
      }
      acc <- acc + if (q %% 2L == 0L) term else -term
    }
    acc
  })
  stats::setNames(lapply(out, name_dims, s = s),
                  param_tuple_names(s, order))
}


#' @title Value and Inverse Map of an Inverse Parameter
#' @name param_value.InverseParam
#' @description
#' `param_value()` returns the inner family's inverse, read through
#' [param_solve()] so that a family with a closed inverse never reaches a
#' factorization. `param_free()` inverts the matrix given and hands it to the
#' inner family's own inverse map, so the free vector is the inner one and the
#' round trip closes on it.
#' @param s An [InverseParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param m A symmetric positive definite matrix of side `s@dimension`.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return `param_value()` a `s@dimension` square matrix, `param_free()` a
#'   numeric vector of length `s@n_free`.
#' @seealso [inverse_of()] for the family.
#' @keywords internal
S7::method(param_value, InverseParam) <- function(s, eta, ...) {
  name_dims(unclass(param_solve(.inv_inner(s), eta)), s)
}

#' @rdname param_value.InverseParam
#' @name param_free.InverseParam
#' @keywords internal
S7::method(param_free, InverseParam) <- function(s, m, ...) {
  param_free(.inv_inner(s), solve(unclass(m)))
}


#' @title Derivative Arrays of an Inverse Parameter
#' @name param_d1.InverseParam
#' @description
#' The four orders of \eqn{\partial N} for \eqn{N = S^{-1}}, each the ordered
#' set-partition sum [inverse_of()] writes out. Every factor comes from the
#' inner family's own arrays, so a closed form there stays closed here.
#' @param s An [InverseParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named list of `s@dimension` square matrices, keyed as
#'   [param_tuple_names()] says.
#' @seealso [inverse_of()] for the formula and [param_d1()] for the contract.
#' @keywords internal
S7::method(param_d1, InverseParam) <- function(s, eta, ...) {
  inverse_derivs(s, eta, 1L)
}

#' @rdname param_d1.InverseParam
#' @name param_d2.InverseParam
#' @keywords internal
S7::method(param_d2, InverseParam) <- function(s, eta, ...) {
  inverse_derivs(s, eta, 2L)
}

#' @rdname param_d1.InverseParam
#' @name param_d3.InverseParam
#' @keywords internal
S7::method(param_d3, InverseParam) <- function(s, eta, ...) {
  inverse_derivs(s, eta, 3L)
}

#' @rdname param_d1.InverseParam
#' @name param_d4.InverseParam
#' @keywords internal
S7::method(param_d4, InverseParam) <- function(s, eta, ...) {
  inverse_derivs(s, eta, 4L)
}


#' @title Log-Determinant of an Inverse Parameter
#' @name param_logdet.InverseParam
#' @description
#' \eqn{\log\lvert S^{-1}\rvert = -\log\lvert S\rvert}, so the value and every
#' derivative order are the inner family's negated. No determinant of the
#' inverted matrix is computed.
#' @param s An [InverseParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return `param_logdet()` a single number; the derivative methods a named
#'   list or vector shaped as the inner family's, with every entry negated.
#' @seealso [inverse_of()] for the family.
#' @keywords internal
S7::method(param_logdet, InverseParam) <- function(s, eta, ...) {
  -param_logdet(.inv_inner(s), eta)
}

#' @rdname param_logdet.InverseParam
#' @name param_dlogdet.InverseParam
#' @keywords internal
S7::method(param_dlogdet, InverseParam) <- function(s, eta, ...) {
  neg_like(param_dlogdet(.inv_inner(s), eta))
}

#' @rdname param_logdet.InverseParam
#' @name param_d2logdet.InverseParam
#' @keywords internal
S7::method(param_d2logdet, InverseParam) <- function(s, eta, ...) {
  neg_like(param_d2logdet(.inv_inner(s), eta))
}

#' @rdname param_logdet.InverseParam
#' @name param_d3logdet.InverseParam
#' @keywords internal
S7::method(param_d3logdet, InverseParam) <- function(s, eta, ...) {
  neg_like(param_d3logdet(.inv_inner(s), eta))
}

#' @rdname param_logdet.InverseParam
#' @name param_d4logdet.InverseParam
#' @keywords internal
S7::method(param_d4logdet, InverseParam) <- function(s, eta, ...) {
  neg_like(param_d4logdet(.inv_inner(s), eta))
}


#' Negate Without Changing the Container
#'
#' @description
#' Returns `-x` for a numeric vector and the elementwise negation for a list,
#' so that the log-determinant derivatives of [inverse_of()] come back in the
#' shape the inner family returned them in.
#'
#' @details
#' The container is part of the contract: `check_parameter()` subtracts the
#' result of [param_dlogdet()] from a reference, which a list does not
#' support, so an `lapply()` over a numeric vector turns a correct value into
#' an error several frames away.
#'
#' @param x A numeric vector or a list of numbers.
#'
#' @return `x` negated, of the same type.
#'
#' @seealso [inverse_of()], the only caller.
#'
#' @keywords internal
neg_like <- function(x) {
  if (is.list(x)) lapply(x, function(z) -z) else -x
}


#' @title Solve and Factor of an Inverse Parameter
#' @name param_solve.InverseParam
#' @description
#' `param_solve()` is the inner family's **value**, no inversion being needed:
#' the inverse of \eqn{S^{-1}} is \eqn{S}. `param_factor()` is the lower
#' triangular \eqn{L} with \eqn{LL^\top} the value, which is the package's
#' convention and the transpose of what [base::chol()] returns.
#' @param s An [InverseParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param b A numeric matrix or vector with `s@dimension` rows, or `NULL` for
#'   the whole inverse.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return `param_solve()` the inner value, or its product with `b`;
#'   `param_factor()` a lower triangular matrix `L` with `tcrossprod(L)` the
#'   value.
#' @seealso [inverse_of()] for the family.
#' @keywords internal
S7::method(param_solve, InverseParam) <- function(s, eta, b = NULL, ...) {
  m <- unclass(param_value(.inv_inner(s), eta))
  if (is.null(b)) name_dims(m, s) else m %*% b
}

#' @rdname param_solve.InverseParam
#' @name param_factor.InverseParam
#' @keywords internal
S7::method(param_factor, InverseParam) <- function(s, eta, ...) {
  t(chol(unclass(param_value(s, eta))))
}
