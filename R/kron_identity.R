#' @include numerical_fallbacks.R
NULL

#' Identical Blocks of a Matrix Parameter
#'
#' @description
#' The S7 class of a block-diagonal matrix built from \eqn{m} identical copies of
#' an inner matrix parameter, \eqn{I_m \otimes S(\eta)}. The blocks **share one
#' free vector**, so `n_free` is the inner parameter's however large \eqn{m} is,
#' and the free names are the inner ones unchanged.
#'
#' It is the first of the four composition wrappers, and the cheapest: every
#' quantity of the contract is a linear lift of the inner parameter's, so nothing
#' is rederived and nothing of size \eqn{(md)^2} is decomposed.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `KronIdentityParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `inner`, the per-block parameter, and `m`, the number of blocks.
#'   `param_name` is `kron(Im, <inner>)`.
#'
#' @seealso [kron_identity()], the constructor, [block_diag()] for blocks that
#'   are **not** identical and do not share a free vector, and
#'   [matrix_parameter()] for the properties this inherits.
#'
#' @examples
#' # Three copies of a 2 x 2 covariance: a 6 x 6 matrix with three free values.
#' s <- kron_identity(log_cholesky(2), 3)
#' S7::S7_inherits(s, KronIdentityParam)
#' c(dimension = s@dimension, n_free = s@n_free)
#' s@free_names
#'
#' # The free vector does not grow with the number of blocks.
#' vapply(c(1, 5, 100), function(m) kron_identity(log_cholesky(2), m)@n_free,
#'        integer(1))
#'
#' @export
KronIdentityParam <- S7::new_class("KronIdentityParam", parent = matrix_parameter)


#' Construct a Block Replication of a Matrix Parameter
#'
#' @description
#' Returns an object holding \eqn{M(\eta) = I_m \otimes S(\eta)} for an inner
#' matrix parameter \eqn{S} of side \eqn{d}: a block-diagonal matrix of side
#' \eqn{md} whose \eqn{m} diagonal blocks are the **same** matrix, sharing one
#' free vector.
#'
#' This is the covariance, or precision, of \eqn{m} independent groups whose
#' within-group structure is common, which is the shape a grouped random-effect
#' term needs: \eqn{S} is the per-group matrix over the coefficients of one group.
#' Because the blocks share their free values, `n_free` is the inner
#' parameter's however many groups there are, so a random slope over a thousand
#' subjects still estimates three numbers.
#'
#' @details
#' # Every quantity is a linear lift
#'
#' \deqn{\partial_k M = I_m \otimes \partial_k S, \qquad
#'       \log|M|_{+} = m \log|S|_{+}, \qquad
#'       M^{-1} = I_m \otimes S^{-1}.}
#'
#' Nothing is rederived and nothing of side \eqn{md} is decomposed: the rank is
#' \eqn{m} times the inner rank, the null basis is the inner one replicated
#' blockwise, and [param_solve()] loops over the blocks solving the inner
#' parameter \eqn{m} times.
#'
#' The multiplier on the log-determinant is worth expecting. At \eqn{m = 3} over a
#' 2 x 2 log-Cholesky covariance, `param_dlogdet()` answers `c(6, 6, 0)` where
#' the inner parameter answers `c(2, 2, 0)`: a step in a shared free value moves
#' every block.
#'
#' # A deficient inner parameter stays deficient
#'
#' The rank is exactly \eqn{m} times the inner rank, so replicating a
#' [scaled_matrix()] on a difference penalty over 2 blocks gives an 8 x 8 matrix
#' of rank 4, and [param_solve()] rejects it as the inner one would. The
#' replicated null basis is annihilated to \eqn{10^{-15}}.
#'
#' # What it does not do
#'
#' The value carries **no dimnames**, where every primitive family in the package
#' labels its rows and columns `v1`, `v2`, ...: `kronecker()` drops them. The
#' other three composition wrappers do the same. See [name_dims()].
#'
#' @section Notation:
#' \eqn{S} is the inner parameter, \eqn{d} its side, \eqn{m} the number of
#' blocks, \eqn{I_m} the \eqn{m \times m} identity and \eqn{\otimes} the
#' Kronecker product. \eqn{\eta} is the free vector, shared by every block.
#'
#' @param structure An object inheriting from [matrix_parameter()], the per-block
#'   parameter. A family that is not a matrix, such as [simplex()], throws
#'   `'structure' must inherit from 'matrix_parameter'.` The inner parameter may
#'   itself be rank deficient.
#' @param m The number of blocks, a single integer of at least 1. `0`, a
#'   fraction, `NA` and a vector all throw `'m' must be a single integer of at
#'   least 1.` At `m = 1` the result equals the inner parameter's value exactly,
#'   which makes it safe to call in a loop over group counts.
#'
#' @return An object of class [KronIdentityParam()], with `dimension` equal to
#'   `m * structure@dimension`, `n_free` and `free_names` the inner parameter's
#'   unchanged, `rank` equal to `m * structure@rank`, `null_basis` the inner one
#'   replicated blockwise, `role` copied from the inner parameter, and
#'   `param_params` holding `inner` and `m`.
#'
#' @seealso [block_diag()] for blocks that differ and each carry their own free
#'   values, [dr_prod()] and [sum_struct()] for the other two compositions, and
#'   [log_cholesky()] or [ar1()] for the usual inner parameters.
#'
#' @examples
#' # Three groups sharing one 2 x 2 covariance: 6 x 6, three free values.
#' inner <- log_cholesky(2)
#' s <- kron_identity(inner, 3)
#' c(dimension = s@dimension, n_free = s@n_free)
#' eta <- c(0.1, 0, -0.2)
#' round(param_value(s, eta), 3)
#'
#' # Every quantity is the inner one lifted, and the log-determinant is
#' # multiplied by the number of blocks.
#' c(outer = param_logdet(s, eta), m_times_inner = 3 * param_logdet(inner, eta))
#' rbind(outer = param_dlogdet(s, eta), inner = param_dlogdet(inner, eta))
#'
#' # The solve is blockwise and exact.
#' max(abs(param_solve(s, eta) - solve(param_value(s, eta))))
#'
#' # The round trip closes, and m = 1 is the inner parameter itself.
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#' max(abs(param_value(kron_identity(inner, 1), eta) - param_value(inner, eta)))
#'
#' # A deficient inner parameter stays deficient, block by block.
#' r <- kron_identity(scaled_matrix(crossprod(diff(diag(4), differences = 2))), 2)
#' c(dimension = r@dimension, rank = r@rank, null = ncol(r@null_basis))
#' max(abs(param_value(r, 0.3) %*% r@null_basis))
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

#' The Inner Parameter, the Block Count, and the Lift
#'
#' @description
#' The three one-line accessors every [KronIdentityParam()] method uses.
#' `.kron_inner()` returns the per-block parameter, `.kron_m()` the number of
#' blocks, and `.kron_lift()` places a `d` by `d` matrix into `m` identical
#' diagonal blocks as `kronecker(diag(m), x)`.
#'
#' @details
#' They exist so that the twelve methods below read as one line each and the
#' storage of `param_params` appears in one place. `.kron_lift()` is the whole
#' arithmetic of the composition: the value, every derivative component and the
#' factor are the inner quantity passed through it.
#'
#' Note that `kronecker()` drops dimnames, so a lifted matrix carries none where
#' the inner one carried `v1`, `v2`, ...
#'
#' @param s A [KronIdentityParam()] object.
#' @param m For `.kron_lift()`, the number of blocks.
#' @param x A numeric matrix to replicate.
#'
#' @return `.kron_inner()` a [matrix_parameter()], `.kron_m()` a single integer,
#'   `.kron_lift()` an `m * nrow(x)` by `m * ncol(x)` numeric matrix with `x` on
#'   the diagonal blocks and zeros elsewhere.
#'
#' @seealso [kron_identity()], which stores what these read.
#'
#' @name kron_accessors
#' @keywords internal
NULL

.kron_inner <- function(s) s@param_params$inner
.kron_m <- function(s) s@param_params$m
.kron_lift <- function(m, x) kronecker(diag(m), x)

#' @title Value of a Block Replication
#' @name param_value.KronIdentityParam
#' @description
#' #' Returns \eqn{I_m \otimes S(\eta)}, the inner parameter's matrix placed into
#' \eqn{m} identical diagonal blocks. One `kronecker()` call and no arithmetic of
#' its own: the value is the inner value lifted.
#'
#' The result carries **no dimnames**, `kronecker()` dropping the inner
#' parameter's `v1`, `v2`, ...
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` symmetric numeric matrix, block
#'   diagonal with `m` identical blocks, and with no dimnames. Positive definite
#'   exactly when the inner parameter is.
#' @seealso [param_free.KronIdentityParam()] for the inverse and [kron_identity()] for the construction.
#' @keywords internal
S7::method(param_value, KronIdentityParam) <- function(s, eta, ...) {
  .kron_lift(.kron_m(s), param_value(.kron_inner(s), eta))
}

#' @title Free Vector of a Block Replication
#' @name param_free.KronIdentityParam
#' @description
#' Reads the free vector from the **first** diagonal block, by inverting the
#' inner parameter there, and then checks the whole matrix against the
#' replication that free vector implies. A matrix whose blocks are not identical,
#' or whose first block is outside the inner family, is rejected.
#'
#' The rejection is complete because the whole rebuild is compared, not the
#' blocks pairwise: the first block might invert cleanly while the fourth is
#' something else entirely.
#' @param s A [KronIdentityParam()] object.
#' @param m A symmetric numeric matrix of side `s@dimension`, block diagonal
#'   with `s@param_params$m` identical blocks, already checked for shape and
#'   symmetry by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.KronIdentityParam()], the map this inverts.
#' @keywords internal
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

#' @title First Derivatives of a Block Replication
#' @name param_d1.KronIdentityParam
#' @description
#' Closed form: \eqn{\partial_k M = I_m \otimes \partial_k S}, the inner
#' parameter's first derivatives lifted one at a time. Exact whenever the inner
#' parameter's are, and numerical whenever they are not, so
#' [param_is_numerical()] on the composite reports the inner parameter's own
#' answer.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` symmetric matrices named by `s@free_names`, each
#'   `s@dimension` by `s@dimension` and block diagonal, with no dimnames.
#' @seealso [param_d2.KronIdentityParam()] for the order above.
#' @keywords internal
S7::method(param_d1, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d1(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

#' @title Second Derivatives of a Block Replication
#' @name param_d2.KronIdentityParam
#' @description
#' Closed form: \eqn{\partial_{kl} M = I_m \otimes \partial_{kl} S}. The blocks
#' share their free values, so there is no cross-block component to compute and
#' the list is exactly the inner parameter's, lifted.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
#'   `param_tuple_names(s)` and in that order, each block diagonal with no
#'   dimnames.
#' @seealso [param_d1.KronIdentityParam()] and [param_d3.KronIdentityParam()] for the neighbouring orders.
#' @keywords internal
S7::method(param_d2, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d2(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

#' @title Third Derivatives of a Block Replication
#' @name param_d3.KronIdentityParam
#' @description
#' Closed form: the inner parameter's third derivatives, each lifted into \eqn{m}
#' identical blocks.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
#'   `param_tuple_names(s, 3)` and in that order, each block diagonal with no
#'   dimnames.
#' @seealso [param_d4.KronIdentityParam()] for the order above.
#' @keywords internal
S7::method(param_d3, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d3(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

#' @title Fourth Derivatives of a Block Replication
#' @name param_d4.KronIdentityParam
#' @description
#' Closed form: the inner parameter's fourth derivatives, each lifted into \eqn{m}
#' identical blocks. The composition costs nothing at any order, which is why a
#' grouped random effect over many levels is affordable to fourth order.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
#'   `param_tuple_names(s, 4)` and in that order, each block diagonal with no
#'   dimnames.
#' @seealso [param_d3.KronIdentityParam()] for the order below.
#' @keywords internal
S7::method(param_d4, KronIdentityParam) <- function(s, eta, ...) {
  lapply(param_d4(.kron_inner(s), eta), .kron_lift, m = .kron_m(s))
}

#' @title Log-Determinant of a Block Replication
#' @name param_logdet.KronIdentityParam
#' @description
#' Closed form: \eqn{\log|M|_{+} = m \log|S|_{+}}, the inner parameter's
#' log-(pseudo-)determinant times the number of blocks. A block-diagonal
#' determinant is the product of the blocks' and they are identical, so nothing
#' of side \eqn{md} is decomposed.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.KronIdentityParam()] for its derivatives.
#' @keywords internal
S7::method(param_logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_logdet(.kron_inner(s), eta)
}

#' @title Log-Determinant Gradient of a Block Replication
#' @name param_dlogdet.KronIdentityParam
#' @description
#' Closed form: \eqn{m} times the inner parameter's gradient. A shared free value
#' moves every block, so the gradient is multiplied rather than copied: at
#' \eqn{m = 3} over a 2 x 2 log-Cholesky covariance it is `c(6, 6, 0)` where the
#' inner parameter answers `c(2, 2, 0)`.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_logdet.KronIdentityParam()] for the quantity differentiated.
#' @keywords internal
S7::method(param_dlogdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_dlogdet(.kron_inner(s), eta)
}

#' @title Log-Determinant Hessian of a Block Replication
#' @name param_d2logdet.KronIdentityParam
#' @description
#' Closed form: \eqn{m} times the inner parameter's Hessian, with the same keying.
#' Zero wherever the inner parameter's is, so replicating a [log_cholesky()]
#' gives exact zeros here and replicating an [ar1()] does not.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
#'   `param_tuple_names(s)`.
#' @seealso [param_dlogdet.KronIdentityParam()] for the order below.
#' @keywords internal
S7::method(param_d2logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_d2logdet(.kron_inner(s), eta)
}

#' @title Third Log-Determinant Derivatives of a Block Replication
#' @name param_d3logdet.KronIdentityParam
#' @description
#' Closed form: \eqn{m} times the inner parameter's third-order vector.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
#'   `param_tuple_names(s, 3)`.
#' @seealso [param_d4logdet.KronIdentityParam()] for the order above.
#' @keywords internal
S7::method(param_d3logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_d3logdet(.kron_inner(s), eta)
}

#' @title Fourth Log-Determinant Derivatives of a Block Replication
#' @name param_d4logdet.KronIdentityParam
#' @description
#' Closed form: \eqn{m} times the inner parameter's fourth-order vector. Exact
#' where the inner parameter's is, so the composition never introduces the
#' accuracy loss a fallback at this order would.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
#'   `param_tuple_names(s, 4)`.
#' @seealso [param_d3logdet.KronIdentityParam()] for the order below.
#' @keywords internal
S7::method(param_d4logdet, KronIdentityParam) <- function(s, eta, ...) {
  .kron_m(s) * param_d4logdet(.kron_inner(s), eta)
}

#' @title Solve of a Block Replication
#' @name param_solve.KronIdentityParam
#' @description
#' Solves blockwise: \eqn{M^{-1} = I_m \otimes S^{-1}}, so the rows of `b`
#' belonging to each block are handed to the inner parameter's own
#' [param_solve()] in turn. Nothing of side \eqn{md} is factorized, and the
#' method inherits whatever the inner parameter does, including a closed inverse
#' where it has one.
#'
#' Measured at \eqn{m = 3} over a 2 x 2 log-Cholesky covariance, the result
#' agrees with `base::solve()` on the assembled matrix to
#' \eqn{3 \times 10^{-17}}.
#'
#' The generic has already rejected a rank-deficient composite, which a deficient
#' inner parameter produces, and filled `b` with the identity when the caller
#' left it out.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param b A numeric matrix with `s@dimension` rows, already coerced from a
#'   vector and defaulted to the identity by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric matrix with `s@dimension` rows and as many columns as `b`.
#' @seealso [param_factor.KronIdentityParam()] for the lifted factor, and
#'   [kron_identity()] for why the deficiency is inherited.
#' @keywords internal
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

#' @title Factor of a Block Replication
#' @name param_factor.KronIdentityParam
#' @description
#' Closed form: \eqn{I_m \otimes L}, the inner parameter's Cholesky factor lifted.
#' A block-diagonal matrix of identical blocks factors blockwise, so one inner
#' factorization serves every block and nothing of side \eqn{md} is
#' decomposed.
#' @param s A [KronIdentityParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`,
#'   already checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` lower triangular numeric matrix,
#'   block diagonal with `m` identical blocks, satisfying
#'   `L %*% t(L) == param_value(s, eta)`, with no dimnames.
#' @seealso [param_solve.KronIdentityParam()] for the blockwise solve.
#' @keywords internal
S7::method(param_factor, KronIdentityParam) <- function(s, eta, ...) {
  .kron_lift(.kron_m(s), param_factor(.kron_inner(s), eta))
}
