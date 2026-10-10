#' @include numerical_fallbacks.R
NULL


#' Matrix Logarithm Parameter
#'
#' @description
#' The S7 class of unstructured symmetric positive definite matrices in the
#' matrix logarithm parametrization, \eqn{M = \exp(S)} with \eqn{S} symmetric and
#' its lower triangle read straight off the free vector. It is a second chart
#' onto the cone that [log_cholesky()] parametrizes, with different quantities
#' in closed form.
#'
#' [matrix_log()] builds one. Every quantity is exact, and every derivative is
#' computed through the eigendecomposition of \eqn{S}. The object stores only
#' the positions of the free values in `param_params$positions`.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `MatrixLogParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `positions`, the row and column of each free value. `n_free` is
#'   \eqn{p(p+1)/2} and `rank` is \eqn{p}.
#'
#' @seealso [matrix_log()], the constructor, [log_cholesky()] for the other
#'   unstructured chart, and [matrix_parameter()] for the properties that
#'   this class inherits.
#'
#' @examples
#' # The same cone as log_cholesky(), the same number of free values.
#' s <- matrix_log(3)
#' c(matrix_log = s@n_free, log_cholesky = log_cholesky(3)@n_free)
#'
#' # But different free names: no entry is transformed here.
#' rbind(matrix_log = s@free_names, log_cholesky = log_cholesky(3)@free_names)
#'
#' @export
MatrixLogParam <- S7::new_class("MatrixLogParam", parent = matrix_parameter)


#' Construct a Matrix Logarithm Parameter
#'
#' @description
#' Returns an object holding the matrix logarithm parametrization of a symmetric
#' positive definite matrix: \eqn{M = \exp(S)} with \eqn{S} symmetric and
#' **genuinely free**. No entry is transformed; the free values fill the lower
#' triangle of \eqn{S} directly, the diagonal first and then below the diagonal
#' column by column, and the exponential does the rest. Any vector in
#' \eqn{\mathbb{R}^{p(p+1)/2}} gives a positive definite matrix, the exponential
#' of a symmetric matrix having positive eigenvalues, as long as no eigenvalue
#' of \eqn{S} exceeds about 709.78, above which the exponential overflows in
#' double precision.
#'
#' It parametrizes the same cone as [log_cholesky()] with the same number of free
#' values. The choice between the two depends on which quantities should have a
#' closed form; see **Details**.
#'
#' @details
#' # The log-determinant and the inverse
#'
#' \deqn{\log|M| = \mathrm{tr}(S) = \sum_{i=1}^{p} \eta_i,}
#'
#' the sum of the diagonal free values, so the log-determinant is linear and its
#' second, third and fourth derivatives are exactly zero. And the inverse is
#' \eqn{M^{-1} = \exp(-S)}, evaluated through the same eigendecomposition, with
#' no factorization, so `param_solve(s, eta)` and `param_value(s, -eta)` agree
#' up to rounding.
#'
#' # The derivatives, and the reason for the Opitz route
#'
#' They are the Frechet derivatives of the matrix exponential, by the
#' Daleckii-Krein representation. With \eqn{S = Q \Lambda Q^\top} and the
#' directions rotated by \eqn{Q}, the \eqn{k}-th derivative contracts the
#' directions against divided differences of \eqn{e^x} of order \eqn{k+1} at the
#' eigenvalues, summed over the orderings of the directions.
#'
#' The divided differences come from the Opitz theorem, the exponential of a
#' small upper bidiagonal matrix read off its corner, in place of the recursive
#' quotient, which cancels catastrophically when two eigenvalues nearly
#' coincide. The Opitz value stays accurate to rounding there. Repeated
#' eigenvalues occur in ordinary use: an \eqn{S} proportional to the identity,
#' or any \eqn{S} with a symmetry, has them exactly.
#'
#' # Choosing between the two unstructured charts
#'
#' The derivatives of [log_cholesky()] are sparse products of single-entry
#' matrices and are far cheaper. Here a fourth-order component sums over up to
#' \eqn{4! = 24} orderings of its directions, each a contraction of
#' \eqn{O(p^5)} operations.
#'
#' This chart has a linear log-determinant and an inverse that is the map at
#' \eqn{-\eta}, while the inverse of [log_cholesky()] needs a triangular solve.
#' It suits a model that needs the covariance and the precision together, or
#' one whose free values should form an unconstrained symmetric matrix.
#' [log_cholesky()] suits a model that takes derivatives repeatedly.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = p(p+1)/2}, \eqn{S} the
#' symmetric matrix it fills, and \eqn{M = \exp(S)} the value. \eqn{Q} and
#' \eqn{\Lambda} are the eigenvectors and eigenvalues of \eqn{S}, and
#' \eqn{e[\lambda_1, \dots, \lambda_m]} a divided difference of the exponential.
#'
#' @param dimension The side \eqn{p} of the matrix. A single positive whole
#'   number, finite and at least 1; any other value signals the error
#'   `'dimension' must be a single positive integer.`
#'
#' @return An object of class [MatrixLogParam()], with `n_free` equal to
#'   \eqn{p(p+1)/2}, `free_names` `S1` ... `Sp` then `S2.1`, `S3.1`, ..., `rank`
#'   equal to `dimension`, an empty `null_basis`, `param_name` `"matrix_log"`,
#'   and `param_params` holding `positions`.
#'
#' @references
#' Daleckii, J. L. and Krein, S. G. (1965). Integration and differentiation of
#' functions of Hermitian operators. *American Mathematical Society
#' Translations* **47**, 1-30.
#'
#' Opitz, G. (1964). Steigungsmatrizen. *Zeitschrift fur Angewandte Mathematik
#' und Mechanik* **44**, T52-T54.
#'
#' @seealso [log_cholesky()] for the other unstructured chart onto the same cone,
#'   [dd_exp()] for the divided differences, and [param_solve()], which here is
#'   the map at \eqn{-\eta}.
#'
#' @examples
#' # A 3 x 3 unstructured covariance. The free names carry no transformation:
#' # the free values are the entries of S itself.
#' s <- matrix_log(3)
#' s@free_names
#' eta <- c(0.2, -0.3, 0.4, 0.1, -0.2, 0.15)
#' M <- param_value(s, eta)
#' round(M, 4)
#' eigen(M, only.values = TRUE)$values > 0
#'
#' # The log-determinant is the trace of S, so it is the sum of the first p
#' # free values, and it agrees with the eigenvalues.
#' c(closed = param_logdet(s, eta), trace = sum(eta[1:3]),
#'   from_eigen = sum(log(eigen(M, only.values = TRUE)$values)))
#'
#' # Which makes its gradient a vector of ones and zeros, and every higher
#' # order exactly zero.
#' param_dlogdet(s, eta)
#' c(second = max(abs(param_d2logdet(s, eta))),
#'   third = max(abs(param_d3logdet(s, eta))),
#'   fourth = max(abs(param_d4logdet(s, eta))))
#'
#' # The inverse is the map at -eta, computed with no factorization.
#' max(abs(param_solve(s, eta) - param_value(s, -eta)))
#' max(abs(param_solve(s, eta) - solve(M)))
#'
#' # The round trip closes.
#' max(abs(param_free(s, M) - eta))
#'
#' @export
matrix_log <- function(dimension) {
  p <- check_param_args(dimension)
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
    param_params = list(positions = pos)
  )
}


#' The Symmetric Matrix Behind a Free Vector
#'
#' @description
#' Fills the lower triangle of \eqn{S} with the free values, in the order
#' `param_params$positions` records, and mirrors it to the upper triangle. No
#' transformation is applied: the free values **are** the entries of \eqn{S}.
#' In [log_cholesky()], by contrast, the diagonal free values are logarithms.
#'
#' @param s A [MatrixLogParam()] object, whose `dimension` and
#'   `param_params$positions` are read.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A symmetric `s@dimension` by `s@dimension` numeric matrix, with no
#'   dimnames and no constraint on its entries.
#'
#' @seealso [mlog_basis()] for its derivative in one free value, and
#'   [param_value.MatrixLogParam()], which exponentiates it.
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
#' Returns \eqn{\partial S / \partial \eta_k}, which is a constant matrix,
#' \eqn{S} being linear in the free vector: a single 1 on the diagonal for a
#' diagonal free value, and a symmetric pair of 1s below and above the diagonal
#' for the rest. These are the directions against which the Frechet
#' derivatives contract.
#'
#' Because they do not depend on \eqn{\eta}, the whole nonlinearity of the family
#' sits in the exponential, and all four derivative orders are contractions of
#' the same fixed set of directions.
#'
#' @param s A [MatrixLogParam()] object, whose `dimension` and
#'   `param_params$positions` are read.
#' @param k The free-value index, a position in `1:s@n_free`.
#'
#' @return A symmetric `s@dimension` by `s@dimension` numeric matrix with one or
#'   two non-zero entries, all of them 1.
#'
#' @seealso [mlog_tables()], which rotates these by the eigenvectors of \eqn{S},
#'   and [mlog_contract()], which contracts them.
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
#' Scaling and squaring with a Taylor series, written out here for the tiny
#' bidiagonal matrices that [dd_exp()] builds. The matrices are at most 5 by 5,
#' so a general-purpose matrix exponential, and the dependency that it would
#' bring, is not needed.
#'
#' @param a A square numeric matrix, upper triangular in every call the package
#'   makes, and at most 5 by 5.
#'
#' @return Its exponential, a numeric matrix of the same shape.
#'
#' @seealso [dd_exp()], the only caller.
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
#' Returns \eqn{e[\lambda_1, \dots, \lambda_m]} by the **Opitz theorem**: build
#' the upper bidiagonal matrix with the arguments on the diagonal and ones above
#' it, exponentiate it, and read the divided difference off the top-right corner.
#' Exact under repeated arguments, where the recursive quotient
#' \eqn{(f[\lambda_2..] - f[\lambda_1..])/(\lambda_m - \lambda_1)} divides by
#' zero, and accurate under nearly repeated ones, where it cancels.
#'
#' @details
#' When two arguments nearly coincide, the quotient loses accuracy as the gap
#' shrinks, while the Opitz value stays accurate to rounding. Near-repeated
#' eigenvalues are ordinary here: any \eqn{S} with a symmetry has them, and a
#' scalar multiple of the identity has an exactly repeated pair.
#'
#' @param lams The arguments, in any order and with repeats allowed. A divided
#'   difference is symmetric in its arguments, so the order does not matter. At
#'   most 5 of them, the number fourth-order derivatives need.
#'
#' @return A single number.
#'
#' @seealso [mlog_expm_small()], which does the exponential, and
#'   [mlog_tables()], which builds whole tables of these.
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
#' of the requested orders, computed once per free vector.
#'
#' @details
#' Everything here depends on \eqn{\eta} alone, never on the index tuple, so it
#' is computed once and every component of the order reads it. The rotated
#' directions are \eqn{Q^\top E_k Q} for each basis direction \eqn{E_k}. The
#' table of \eqn{m} points is an array with one divided difference for each
#' ordered \eqn{m}-tuple of eigenvalue indices, \eqn{p^m} entries in all, for
#' \eqn{m} up to `order + 1`.
#'
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The highest derivative order wanted: 1, 2, 3 or 4.
#'
#' @return A list with `q` and `lam`, the eigenvectors and eigenvalues of
#'   \eqn{S}; `e`, a list of `s@n_free` rotated directions; and the
#'   divided-difference tables `dd2` up to `dd<order+1>`.
#'
#' @seealso [dd_exp()] for the entries of the tables, and [mlog_contract()],
#'   which reads them.
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
#' divided-difference table of the matching order and sums over the orderings
#' of the directions. The result is in the eigenbasis of \eqn{S}; the caller
#' rotates it back by \eqn{Q}.
#'
#' @details
#' The multilinear form sums over all orderings of the directions. The loop
#' runs over the distinct orderings, which [combinat_perms()] returns, and the
#' sum is multiplied by \eqn{\prod_j m_j!}, where the \eqn{m_j} are the
#' multiplicities of the repeated indices, because each distinct ordering
#' occurs that many times among all of them. Without that factor the result
#' would be too small by exactly that factor.
#'
#' The cost is the reason this family's higher orders are expensive: at order 4
#' a component with four different indices has 24 orderings, each a contraction
#' of \eqn{O(p^5)} operations.
#'
#' @param tb The tables of [mlog_tables()], carrying at least
#'   `length(dirs) + 1` points.
#' @param dirs The free-value indices of the tuple, possibly repeated. Its length
#'   is the derivative order.
#'
#' @return A symmetric numeric matrix with the side of `tb$q`, in the
#'   eigenbasis of \eqn{S} and with no dimnames.
#'
#' @seealso [mlog_tables()] for the input, [combinat_perms()] for the orderings,
#'   and [matrix_log()] for the representation.
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
#' Returns the distinct permutations of an index tuple as a list. When the tuple
#' has repeated entries, permutations that coincide appear once, so
#' `combinat_perms(c(1, 1))` gives one element. [mlog_contract()] restores the
#' multiplicity with the product of the factorials of the repeat counts.
#'
#' @param x An integer vector, of length 1 to 4 in every call the package makes.
#'
#' @return A list of the distinct orderings of `x`, each an integer vector:
#'   `factorial(length(x))` of them when the entries of `x` are all different,
#'   and fewer when some repeat.
#'
#' @seealso [mlog_contract()], the only caller.
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
#' against divided differences of \eqn{\exp} in the eigenvalues. The sum over
#' all orderings of the directions is computed as the sum over the distinct
#' orderings, weighted by the multiplicity factor of a tuple with repeated
#' indices (see [mlog_contract()]). The divided differences come from
#' the Opitz representation, an exponential of a small bidiagonal matrix read
#' off its corner, which stays exact where the quotient recursion cancels
#' catastrophically under near-repeated eigenvalues.
#'
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension`.
#'
#' @seealso [mlog_contract()], which computes one component, [mlog_tables()] for
#'   the shared tables, and [matrix_log()] for the cost this order carries.
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
#' @description
#' Returns \eqn{M = \exp(S)}, through the eigendecomposition of \eqn{S}:
#' \eqn{Q \exp(\Lambda) Q^\top}, the exponential applied to the eigenvalues.
#' Positive definiteness needs no checking, every \eqn{e^{\lambda}} being
#' positive whatever \eqn{S} is, up to the overflow of the exponential for an
#' eigenvalue above about 709.78. The cost is one \eqn{O(p^3)}
#' eigendecomposition.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A symmetric positive definite `s@dimension` by `s@dimension` numeric
#'   matrix with dimnames `v1`, `v2`, ...
#' @seealso [param_free.MatrixLogParam()] for the inverse, [mlog_s()] for the
#'   assembly of \eqn{S}, and [param_solve.MatrixLogParam()], which is this map
#'   at \eqn{-\eta}.
#' @keywords internal
S7::method(param_value, MatrixLogParam) <- function(s, eta, ...) {
  es <- eigen(mlog_s(s, eta), symmetric = TRUE)
  name_dims(es$vectors %*% (exp(es$values) * t(es$vectors)), s)
}


#' @title Free Vector of a Matrix Logarithm Parameter
#' @name param_free.MatrixLogParam
#' @description
#' Returns the matrix logarithm of `m`, read off its eigendecomposition as
#' \eqn{Q \log(\Lambda) Q^\top} and then read out of the lower triangle. Exact
#' and a true inverse of [param_value.MatrixLogParam()], the matrix logarithm of
#' a symmetric positive definite matrix being unique among symmetric matrices.
#' The round trip passes through two eigendecompositions, one in each
#' direction.
#' @details
#' A matrix that is not positive definite is rejected: the logarithm of a
#' non-positive eigenvalue is not a real number, so there is no free vector to
#' return. The test is on the eigenvalues of `m`.
#' @param s A [MatrixLogParam()] object.
#' @param m A symmetric positive definite `s@dimension` by `s@dimension` numeric
#'   matrix, already checked for shape and symmetry by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.MatrixLogParam()], the map this inverts.
#' @keywords internal
S7::method(param_free, MatrixLogParam) <- function(s, m, ...) {
  es <- eigen(m, symmetric = TRUE)
  if (any(es$values <= 0)) {
    stop(paste0(
      "'m' is not positive definite, so it is not in the set that\n",
      "  matrix_log() parametrizes."
    ), call. = FALSE)
  }
  sm <- es$vectors %*% (log(es$values) * t(es$vectors))
  pos <- s@param_params$positions
  stats::setNames(sm[cbind(pos$row, pos$col)], s@free_names)
}


#' @title First Derivatives of a Matrix Logarithm Parameter
#' @name param_d1.MatrixLogParam
#' @description
#' The Frechet derivative of the matrix exponential, by Daleckii-Krein: rotate
#' the basis direction \eqn{E_k} by the eigenvectors of \eqn{S}, weight it
#' entrywise by the two-point divided differences
#' \eqn{e[\lambda_i, \lambda_j]}, and rotate back. Exact, with no
#' differencing.
#'
#' The derivative \eqn{\partial_k M} differs from \eqn{E_k \exp(S)}, because the
#' exponential of a matrix does not commute with an arbitrary direction; the
#' weighting by divided differences accounts for this. Where the eigenvalues
#' are equal the weight is \eqn{e^{\lambda}} and the formula reduces to the
#' scalar one, which the Opitz route of [dd_exp()] delivers exactly.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` symmetric matrices named by `s@free_names`, each
#'   `s@dimension` by `s@dimension`.
#' @seealso [mlog_contract()], which does the contraction, [dd_exp()] for the
#'   divided differences, and [param_d2.MatrixLogParam()] for the order above.
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
#' The second Frechet derivative: chains of two rotated directions contracted
#' against three-point divided differences \eqn{e[\lambda_i, \lambda_j,
#' \lambda_k]}, summed over both orderings of the two directions.
#'
#' For the component in one free value twice the two orderings coincide;
#' [mlog_contract()] evaluates that ordering once and multiplies it by 2, which
#' gives the sum over both orderings.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
#'   `param_tuple_names(s)` and in that order.
#' @seealso [mlog_contract()], which does the contraction, and
#'   [param_d1.MatrixLogParam()] and [param_d3.MatrixLogParam()] for the
#'   neighboring orders.
#' @keywords internal
S7::method(param_d2, MatrixLogParam) <- function(s, eta, ...) {
  mlog_higher(s, eta, 2L)
}


#' @title Third Derivatives of a Matrix Logarithm Parameter
#' @name param_d3.MatrixLogParam
#' @description
#' Chains of three rotated directions contracted against four-point divided
#' differences, summed over the orderings of the three directions. Exact, with
#' no differencing.
#'
#' The cost grows quickly: up to six orderings per component, each an
#' \eqn{O(p^4)} contraction, and \eqn{\binom{d+2}{3}} components. See
#' [matrix_log()] for the comparison with [log_cholesky()], whose derivatives
#' are sparse products.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
#'   `param_tuple_names(s, 3)` and in that order.
#' @seealso [mlog_higher()], which assembles it, and
#'   [param_d4.MatrixLogParam()] for the order above.
#' @keywords internal
S7::method(param_d3, MatrixLogParam) <- function(s, eta, ...) {
  mlog_higher(s, eta, 3L)
}


#' @title Fourth Derivatives of a Matrix Logarithm Parameter
#' @name param_d4.MatrixLogParam
#' @description
#' Chains of four rotated directions contracted against five-point divided
#' differences, summed over the orderings of the four directions, up to
#' twenty-four per component. Exact, and much more expensive than the fourth
#' derivatives of [log_cholesky()] on a matrix of the same size.
#'
#' The cost is the price of a linear log-determinant and of an inverse equal to
#' the map at \eqn{-\eta}. When fourth derivatives are taken repeatedly and those
#' two properties are not needed, [log_cholesky()] parametrizes the same cone at
#' a much lower cost.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
#'   `param_tuple_names(s, 4)` and in that order.
#' @seealso [mlog_higher()], which assembles it, [param_d3.MatrixLogParam()] for
#'   the order below, and [log_cholesky()] for the cheaper chart.
#' @keywords internal
S7::method(param_d4, MatrixLogParam) <- function(s, eta, ...) {
  mlog_higher(s, eta, 4L)
}


#' @title Log-Determinant of a Matrix Logarithm Parameter
#' @name param_logdet.MatrixLogParam
#' @description
#' Closed form and linear. The eigenvalues of \eqn{M = \exp(S)} are
#' \eqn{e^{\lambda_i}}, so
#'
#' \deqn{\log|M| = \sum_i \lambda_i = \mathrm{tr}(S) = \sum_{i=1}^{p} \eta_i,}
#'
#' the sum of the diagonal free values. One `sum()` over \eqn{p} numbers: no
#' eigendecomposition and no determinant, with nothing growing in \eqn{p} beyond
#' the sum itself.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.MatrixLogParam()] for its gradient, and
#'   [param_logdet.LogCholeskyParam()], which is linear for the same kind of
#'   reason.
#' @keywords internal
S7::method(param_logdet, MatrixLogParam) <- function(s, eta, ...) {
  sum(eta[s@param_params$positions$on_diagonal])
}

#' @title Log-Determinant Gradient of a Matrix Logarithm Parameter
#' @name param_dlogdet.MatrixLogParam
#' @description
#' Closed form and constant: 1 in each of the \eqn{p} diagonal directions and 0
#' in the \eqn{p(p-1)/2} below-diagonal ones, the log-determinant being
#' \eqn{\sum_i \eta_i}. The value does not depend on `eta`, which is read only
#' for its length.
#'
#' The gradient is 1 where that of [log_cholesky()] is 2: there the diagonal
#' free value is \eqn{\log L_{ii}} and \eqn{|M| = |L|^2}, while here it is a
#' diagonal entry of \eqn{S}, and \eqn{\log|M| = \mathrm{tr}(S)} contains it
#' once.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`,
#'   holding 1 and 0.
#' @seealso [param_logdet.MatrixLogParam()] for the quantity differentiated, and
#'   [param_dlogdet.LogCholeskyParam()], which returns 2.
#' @keywords internal
S7::method(param_dlogdet, MatrixLogParam) <- function(s, eta, ...) {
  stats::setNames(
    ifelse(s@param_params$positions$on_diagonal, 1, 0), s@free_names
  )
}

#' @title Log-Determinant Hessian of a Matrix Logarithm Parameter
#' @name param_d2logdet.MatrixLogParam
#' @description
#' Closed form, and identically zero: \eqn{\log|M| = \mathrm{tr}(S)} is linear in
#' the free vector, so its second derivative vanishes at every \eqn{\eta}. The
#' zeros are exact.
#'
#' As with [log_cholesky()], the order is zero for this family, so a test of a
#' general log-determinant Hessian has nothing to compare here. [ar1()],
#' [compound_symmetry()] and [correlation_matrix()] are families where it is
#' not zero.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 1, 2)` zeros, keyed as
#'   `param_tuple_names(s)`.
#' @seealso [param_dlogdet.MatrixLogParam()] for the order below, and
#'   [param_d3logdet.MatrixLogParam()], zero for the same reason.
#' @keywords internal
S7::method(param_d2logdet, MatrixLogParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 2L)
  stats::setNames(rep(0, length(nm)), nm)
}

#' @title Third Log-Determinant Derivatives of a Matrix Logarithm Parameter
#' @name param_d3logdet.MatrixLogParam
#' @description
#' Closed form, and identically zero, the trace being linear in the free vector.
#' The zeros are exact, so a consumer can drop the term instead of carrying small
#' numbers through a contraction.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 2, 3)` zeros, keyed as
#'   `param_tuple_names(s, 3)`.
#' @seealso [param_d2logdet.MatrixLogParam()] for the order below, and
#'   [param_d4logdet.MatrixLogParam()] for the one above.
#' @keywords internal
S7::method(param_d3logdet, MatrixLogParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 3L)
  stats::setNames(rep(0, length(nm)), nm)
}

#' @title Fourth Log-Determinant Derivatives of a Matrix Logarithm Parameter
#' @name param_d4logdet.MatrixLogParam
#' @description
#' Closed form, and identically zero, for the reason it is zero at second and
#' third order: \eqn{\log|M| = \mathrm{tr}(S)} is linear in the free vector. The
#' zeros are exact.
#'
#' This is the order at which a numerical fallback is least accurate, so a
#' family whose log-determinant is not linear gains most from a closed form
#' here; see [param_d4logdet.matrix_parameter()] for the fallback.
#' @param s A [MatrixLogParam()] object.
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
#'   eigendecomposition as the value, applied to `b`.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values.
#' @param b A numeric matrix with `s@dimension` rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(param_solve, MatrixLogParam) <- function(s, eta, b = NULL, ...) {
  es <- eigen(mlog_s(s, eta), symmetric = TRUE)
  inv <- es$vectors %*% (exp(-es$values) * t(es$vectors))
  inv %*% b
}
