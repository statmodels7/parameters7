#' @include numerical_fallbacks.R
NULL

#' Distinct Blocks of a Matrix Parameter
#'
#' @description
#' The S7 class of a block-diagonal matrix built from several matrix parameters,
#' each block a family of its own and each carrying its own stretch of the free
#' vector. [block_diag()] builds one.
#'
#' It is the composition for independent groups of coefficients whose structures
#' differ, where [kron_identity()] is the one for identical blocks sharing a
#' single free vector.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `BlockDiagParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `blocks`, `labels`, `rows` and `free` (the ranges each block occupies in the
#'   matrix and in the free vector) and `owner` (the block each free value
#'   belongs to). `dimension`, `n_free` and `rank` are the sums of the blocks'.
#'
#' @seealso [block_diag()], the constructor, [kron_identity()] for identical
#'   blocks, and [matrix_parameter()] for the properties this inherits.
#'
#' @examples
#' # Everything the composite is, is the sum of what the blocks are.
#' s <- block_diag(subject = log_cholesky(2), time = ar1(3))
#' c(dimension = s@dimension, n_free = s@n_free, rank = s@rank)
#'
#' # And the free names carry the block's label, so two blocks of the same
#' # family stay distinguishable.
#' block_diag(a = ar1(3), b = ar1(4))@free_names
#'
#' @export
BlockDiagParam <- S7::new_class("BlockDiagParam", parent = matrix_parameter)


#' Construct a Block Diagonal of Matrix Parameters
#'
#' @description
#' Returns an object holding
#' \eqn{M(\eta) = \mathrm{diag}(S_1(\eta_1), \ldots, S_B(\eta_B))} for matrix
#' parameters \eqn{S_1, \ldots, S_B}: a block-diagonal matrix whose blocks are
#' **distinct families**, each governed by its own stretch of the free vector.
#'
#' It is the covariance (or precision) of several independent groups of
#' coefficients whose structures differ, as a model carrying more than one
#' random-effect term has.
#'
#' @details
#' # Nothing is rederived, and the reason
#'
#' The free values of one block do not enter another, so
#'
#' \deqn{\partial_k M = \mathrm{diag}(0, \ldots, \partial_k S_b, \ldots, 0),
#'       \qquad \log\lvert M \rvert_{+} = \sum_b \log\lvert S_b \rvert_{+},
#'       \qquad M^{-1} = \mathrm{diag}(S_1^{-1}, \ldots, S_B^{-1}),}
#'
#' and a derivative whose indices do not all belong to one block is
#' **identically zero**, at every order and for the log-determinant as well as
#' for the value. Measured on `block_diag(log_cholesky(2), ar1(3))`, whose five
#' free values split 3 and 2: 6 of the 15 second-order components are
#' cross-block, 21 of 35 at third order and 50 of 70 at fourth, and every one of
#' them is exactly 0. What is left is fetched from the block and placed in the
#' rows and columns that block occupies.
#'
#' The rank is the sum of the blocks' ranks and the null basis is their block
#' diagonal, both read from the components and never from an assembled matrix,
#' which is the rule this package follows everywhere: a rank is a property of the
#' family, and counting eigenvalues of the assembled matrix would make it a
#' property of the arithmetic. A deficient block is therefore admitted, and the
#' composite reports the deficiency: with a rank-one \eqn{P} of side 2,
#' `block_diag(scaled_matrix(P), ar1(3))` has rank 4 of 5 and a null basis of one
#' column, `param_logdet()` returns the log pseudo-determinant, and
#' `param_solve()` is refused by the generic.
#'
#' # Against kron_identity()
#'
#' [kron_identity()] repeats **one** block \eqn{m} times and they share a single
#' free vector, so its `n_free` does not grow with \eqn{m}. Here the blocks are
#' different objects and their free vectors are concatenated, so the composite has
#' \eqn{\sum_b d_b} free values. Use this one where the groups have different
#' structures, and that one where they have the same structure and the same
#' parameters.
#'
#' # Labels
#'
#' The free names are prefixed by the block's label, since two blocks of the same
#' family would otherwise report the same names and the class requires them to be
#' unique. Labels come from the names of the arguments where they are given, and
#' are `b1`, `b2`, ... otherwise. `param_name` records the blocks' own names, as
#' in `"blockdiag(log_cholesky, ar1)"`.
#'
#' @section Notation:
#' \eqn{B} is the number of blocks, \eqn{S_b} the \eqn{b}-th block's map,
#' \eqn{d_b} its number of free values, and \eqn{\eta_b} the stretch of the free
#' vector it owns.
#'
#' @param ... One or more objects inheriting from [matrix_parameter()], or a
#'   single list of them. Named arguments supply the block labels, which must be
#'   unique. A block that is not a `matrix_parameter` is rejected by position, so
#'   `block_diag(ar1(3), simplex(3))` reports that block 2 does not inherit from
#'   it.
#'
#' @return An object of class [BlockDiagParam()], with `dimension`, `n_free` and
#'   `rank` the sums of the blocks', `free_names` the blocks' own prefixed by the
#'   labels, and `null_basis` the block diagonal of the blocks'.
#'
#' @seealso [kron_identity()] for identical blocks, [dr_prod()] and
#'   [sum_struct()] for the other two compositions, and [matrix_parameter()] for
#'   the contract every block meets.
#'
#' @examples
#' s <- block_diag(subject = log_cholesky(2), time = ar1(3))
#' c(dimension = s@dimension, n_free = s@n_free)
#' s@free_names
#'
#' eta <- c(0.1, -0.2, 0.3, log(2), atanh(0.5))
#' M <- param_value(s, eta)
#'
#' # Each block is its own family's value, and the off-diagonal is exactly zero.
#' c(subject = max(abs(M[1:2, 1:2] - param_value(log_cholesky(2), eta[1:3]))),
#'   time = max(abs(M[3:5, 3:5] - param_value(ar1(3), eta[4:5]))),
#'   off = max(abs(M[1:2, 3:5])))
#'
#' # The log-determinant is the sum of the blocks'.
#' c(composite = param_logdet(s, eta),
#'   sum_of_blocks = param_logdet(log_cholesky(2), eta[1:3]) +
#'                   param_logdet(ar1(3), eta[4:5]))
#'
#' # And 50 of the 70 fourth-order components are exactly zero, their indices
#' # not all belonging to one block.
#' owner <- s@param_params$owner
#' cross <- vapply(param_tuple_indices(s, 4L),
#'                 function(t) length(unique(owner[t])) > 1L, logical(1))
#' c(cross_block = sum(cross), of = length(cross),
#'   largest = max(abs(unlist(param_d4(s, eta)[cross]))))
#'
#' # The round trip closes.
#' max(abs(param_free(s, M) - eta))
#'
#' @export
block_diag <- function(...) {
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
    param_params = list(blocks = blocks, labels = labels, rows = rows,
                        free = free, owner = owner)
  )
}


#' Consecutive Index Ranges of Given Widths
#'
#' @description
#' Turns the widths \eqn{n_1, \ldots, n_B} into the ranges they occupy when laid
#' end to end. It is called twice at construction, once for the blocks' rows in
#' the matrix and once for their stretches of the free vector.
#'
#' @details
#' A width of zero gives an **empty range**, `integer(0)`, instead of being
#' dropped, so the list stays aligned with the blocks and the \eqn{b}-th element
#' is always the \eqn{b}-th block's. That case is reachable: a fully known
#' [scaled_matrix()] has no free values at all, so
#' `split_ranges(c(2, 0, 3))` returns `1:2`, `integer(0)` and `3:5`.
#'
#' @param widths An integer vector of non-negative widths.
#'
#' @return A list of integer vectors, one per width, in order.
#'
#' @seealso [block_diag()], the only caller.
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
#' **sorted local index tuple**, as `"1"`, `"1,2"`, `"2,2"` and so on.
#'
#' @details
#' The re-keying lets a lookup from the composite's enumeration succeed without
#' assuming the two enumerations correspond. A block's own
#' `param_tuple_names()` are built from its own `free_names`, which the composite
#' has prefixed with a label, so the two name sets differ; the sorted index tuple
#' is the one key both sides can compute. Sorting matters because a derivative is
#' symmetric in its indices and the two enumerations need not order a tuple the
#' same way.
#'
#' @param block A [matrix_parameter()], one block of the composite.
#' @param eta The block's own stretch of the free vector, of length
#'   `block@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of the block's own derivative matrices of that order, named by
#'   the sorted local index tuples.
#'
#' @seealso [block_diag_derivs()], the only caller.
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
#' Places each block's own component in the rows and columns that block occupies,
#' and returns a zero matrix for a tuple whose indices are not all owned by one
#' block.
#'
#' @details
#' Each block's derivatives of the order are fetched at most once and held in a
#' list for the rest of the call, since the composite's enumeration visits a
#' block's tuples several times and a block's fourth-order array is the expensive
#' thing here. The zero matrix is one object shared by every cross-block
#' component: at order 4 with two blocks of 3 and 2 free values that is 50 of the
#' 70 components.
#'
#' @param s A [BlockDiagParam()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.
#'
#' @seealso [block_derivs_by_tuple()] for one block's components, and
#'   [param_d1.BlockDiagParam()], which calls this.
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
  # the dimnames convention every family's matrices carry; one point
  # here covers all four derivative orders, which route through this
  stats::setNames(lapply(out, name_dims, s = s),
                  param_tuple_names(s, order))
}


#' Assemble a Block Diagonal's Log-Determinant Derivatives
#'
#' @description
#' The log-determinant is the sum of the blocks', so it is separable across
#' blocks: a component is the owning block's own, and 0 for a tuple spanning two
#' blocks.
#'
#' @details
#' The structure is [block_diag_derivs()]'s with numbers in place of matrices,
#' including the per-block cache and the sorted-tuple key. It is a different
#' separability from the one [ar1()] or [compound_symmetry()] have: there the
#' log-determinant separates over the free values themselves, here it separates
#' over the blocks and a block's own mixed components survive.
#'
#' @param s A [BlockDiagParam()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
#'   as `param_tuple_names(s, order)` and in that order.
#'
#' @seealso [param_dlogdet.BlockDiagParam()], which calls this, and
#'   [block_diag_derivs()], its twin for the value.
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
#' @description
#' Evaluates each block at its own stretch of the free vector and writes it into
#' the rows and columns that block occupies, leaving the off-diagonal blocks at
#' the zeros the matrix was created with. Each block is exactly what its own
#' family returns, to the bit.
#'
#' The value is labeled `v1`, `v2`, ..., `vp` on both margins, the convention [name_dims()]
#' states and every family in the package follows.
#' @param s A [BlockDiagParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A symmetric `s@dimension` by `s@dimension` numeric matrix, block
#'   diagonal and labeled `v1`, `v2`, ..., `vp` on both margins.
#' @seealso [param_free.BlockDiagParam()] for the inverse, and [block_diag()] for
#'   the composition.
#' @keywords internal
S7::method(param_value, BlockDiagParam) <- function(s, eta, ...) {
  b <- .bd(s)
  out <- matrix(0, s@dimension, s@dimension)
  for (j in seq_along(b$blocks)) {
    out[b$rows[[j]], b$rows[[j]]] <- param_value(b$blocks[[j]], eta[b$free[[j]]])
  }
  name_dims(out, s)
}

#' @title Free Vector of a Block-Diagonal Parameter
#' @name param_free.BlockDiagParam
#' @description
#' Inverts each diagonal block through its own parameter and concatenates the
#' results, so the composite is exact wherever its blocks are: measured on
#' `block_diag(log_cholesky(2), ar1(3))`, the round trip closes to
#' \eqn{7 \times 10^{-17}}.
#' @details
#' A matrix whose off-diagonal blocks are not zero is rejected with
#' `'m' is not block diagonal in the blocks of this parameter.`, that being a
#' matrix this family cannot represent; the tolerance is \eqn{10^{-8}} relative
#' to `max(1, max(abs(m)))`. The blocks are inverted first, so a diagonal block
#' that is outside its own family is reported by that family's message, which is
#' the more specific of the two.
#' @param s A [BlockDiagParam()] object.
#' @param m A symmetric `s@dimension` by `s@dimension` matrix, block diagonal in
#'   this parameter's blocks, already checked for shape and symmetry by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, the blocks' free vectors
#'   concatenated.
#' @seealso [param_value.BlockDiagParam()], the map this inverts.
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
#' `param_d1()`, `param_d2()`, `param_d3()` and `param_d4()` for a
#' [block_diag()] parameter: each block's own derivatives, placed in the rows and
#' columns that block occupies, and **exactly zero** for a tuple spanning two
#' blocks. Nothing is rederived and nothing is differenced, so the composite is
#' as exact as its blocks are.
#' @param s A [BlockDiagParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return At order 1, a list of `s@n_free` symmetric matrices named by
#'   `s@free_names`; above it, `choose(s@n_free + k - 1, k)` of them keyed as
#'   `param_tuple_names(s, k)` and in that order. Each is `s@dimension` by
#'   `s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.
#' @seealso [block_diag_derivs()], which assembles them, and
#'   [param_dlogdet.BlockDiagParam()] for the log-determinant's own.
#' @details
#' The four share [block_diag_derivs()] and differ only in the order they pass.
#' The zeros are structural: they follow from the free values of one block not
#' entering another, and at order 4 with blocks of 3 and 2 free values they are
#' 50 of the 70 components. See [block_diag()] for the argument.
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
#' The sum of the blocks' log-(pseudo-)determinants, the determinant of a
#' block-diagonal matrix being the product of the blocks'. Each block computes
#' its own by whatever route it has, so a closed form in a block stays a closed
#' form in the composite, and no determinant of the assembled matrix is taken.
#'
#' Where a block is rank deficient its term is the log **pseudo**-determinant,
#' the sum over the non-zero eigenvalues, and the composite's is then a pseudo
#' one too.
#' @param s A [BlockDiagParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.BlockDiagParam()] for its derivatives.
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
#' `param_dlogdet()`, `param_d2logdet()`, `param_d3logdet()` and
#' `param_d4logdet()` for a [block_diag()] parameter: the blocks' own, placed by
#' owner, with every cross-block component exactly zero because the
#' log-determinant is a sum over the blocks.
#' @details
#' The four share [block_diag_logdet_derivs()] and differ only in the order they
#' pass. The first order names its result by `s@free_names`, one value per free
#' value; the orders above it are keyed by tuple, as the contract requires.
#'
#' A block's **own** mixed components are not zero, this separability being over
#' blocks, never over free values. Compare [ar1()], where it is over free values
#' and every mixed component vanishes.
#' @param s A [BlockDiagParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector: at order 1, `s@n_free` values named by
#'   `s@free_names`; above it, `choose(s@n_free + k - 1, k)` values keyed as
#'   `param_tuple_names(s, k)` and in that order.
#' @seealso [block_diag_logdet_derivs()], which assembles them, and
#'   [param_logdet.BlockDiagParam()] for the quantity differentiated.
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
#' `param_solve()` and `param_factor()` for a [block_diag()] parameter, both of
#' them blockwise: the inverse of a block-diagonal matrix is the block diagonal
#' of the inverses, and the same holds of a lower triangular factor. Each block
#' answers by whatever route it has, so `ar1()`'s tridiagonal inverse and
#' `log_cholesky()`'s free factor both survive into the composite.
#' @details
#' `param_solve()` takes the rows of `b` that belong to each block and hands them
#' to that block, so the whole matrix is never inverted. `param_factor()`
#' assembles the blocks' factors on the diagonal; the result is lower triangular
#' with \eqn{M = L L^\top}, which is the contract [param_factor()] states.
#' Measured against `solve()` and against the assembled matrix, the two agree to
#' \eqn{2 \times 10^{-16}} and \eqn{4 \times 10^{-16}}.
#'
#' Both are refused by the generic where any block is rank deficient, the
#' composite having no inverse and no Cholesky factor then.
#' @param s A [BlockDiagParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param b A numeric matrix with `s@dimension` rows, defaulted to the identity
#'   by the generic. `param_factor()` takes no `b`.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return `param_solve()` returns a numeric matrix with `s@dimension` rows and
#'   as many columns as `b`; `param_factor()` a lower triangular `s@dimension` by
#'   `s@dimension` matrix.
#' @seealso [param_solve()] and [param_factor()] for the two contracts.
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
