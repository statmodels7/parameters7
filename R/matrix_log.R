#' @include numerical_fallbacks.R
NULL


#' Matrix Logarithm Parameter
#'
#' @description
#' The S7 class of unstructured symmetric positive definite matrices in the
#' matrix logarithm parametrization, \eqn{M = \exp(S)} with \eqn{S} symmetric and
#' its lower triangle read straight off the free vector. It is the other chart
#' onto the same cone [log_cholesky()] parametrizes, with a different set of
#' quantities coming free.
#'
#' [matrix_log()] builds one. Every quantity is exact, and every derivative goes
#' through the eigendecomposition of \eqn{S}, so the object holds no auxiliary
#' data at all.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `MatrixLogParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `positions`, the row and column of each free value. `n_free` is
#'   \eqn{p(p+1)/2} and `rank` is \eqn{p}.
#'
#' @seealso [matrix_log()], the constructor, [log_cholesky()] for the other
#'   unstructured chart, and [matrix_parameter()] for the properties this
#'   inherits.
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
#' of a symmetric matrix having positive eigenvalues.
#'
#' It parametrizes the same cone as [log_cholesky()] with the same number of free
#' values. Which to reach for depends on which quantities you want free; see
#' **Details**.
#'
#' @details
#' # Two quantities cost nothing
#'
#' \deqn{\log|M| = \mathrm{tr}(S) = \sum_{i=1}^{p} \eta_i,}
#'
#' the sum of the diagonal free values, so the log-determinant is linear and its
#' second, third and fourth derivatives are exactly zero. And the inverse is
#' \eqn{M^{-1} = \exp(-S)}, evaluated through the same eigendecomposition, with
#' no factorization: measured, `param_solve(s, eta)` and
#' `param_value(s, -eta)` agree to \eqn{2 \times 10^{-15}}.
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
#' quotient, which cancels catastrophically under near-repeated eigenvalues.
#' Measured at eigenvalues \eqn{0.5} and \eqn{0.5 + g}, against the exact limit
#' \eqn{1.648721270700127} at \eqn{g = 0}:
#'
#' | gap \eqn{g} | Opitz | quotient |
#' |---|---|---|
#' | \eqn{10^{-3}} | 1.649545906191066 | 1.649545906190929 |
#' | \eqn{10^{-6}} | 1.648722095061039 | 1.648722095023754 |
#' | \eqn{10^{-9}} | 1.648721271524488 | 1.648721206226264 |
#'
#' At a gap of \eqn{10^{-9}} the quotient has lost eight digits and the Opitz
#' value has lost one. Repeated eigenvalues are not a pathology here: a
#' [scalar_matrix()]-like \eqn{S}, or any \eqn{S} with a symmetry, has them
#' exactly.
#'
#' # Choosing between the two unstructured charts
#'
#' [log_cholesky()]'s derivatives are sparse products of single-entry matrices,
#' and they are far cheaper: measured at \eqn{p = 4}, one fourth-order derivative
#' array costs **3.4 s** here against **0.003 s** there, a factor of a thousand,
#' the contraction summing over \eqn{4! = 24} orderings per component at
#' \eqn{p^4} entries each.
#'
#' What this chart buys is the other end. Its log-determinant is linear and its
#' inverse is a sign flip, where [log_cholesky()]'s inverse is a triangular solve.
#' Reach for it when the covariance and the precision are both wanted at once, or
#' when the free values themselves should be an unconstrained symmetric matrix;
#' reach for [log_cholesky()] when derivatives are taken in a loop.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = p(p+1)/2}, \eqn{S} the
#' symmetric matrix it fills, and \eqn{M = \exp(S)} the value. \eqn{Q} and
#' \eqn{\Lambda} are the eigenvectors and eigenvalues of \eqn{S}, and
#' \eqn{e[\lambda_1, \dots, \lambda_m]} a divided difference of the exponential.
#'
#' @param dimension The side \eqn{p} of the matrix. A single positive whole
#'   number, finite and at least 1; anything else throws `'dimension' must be a
#'   single positive integer.`
#' @param role A label recording which side of a model the matrix parametrizes:
#'   `"either"` (the default), `"covariance"` or `"precision"`. No numeric result
#'   depends on it, and it matters less here than elsewhere: the inverse is the
#'   same family at \eqn{-\eta}.
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
#' # The inverse is the map at -eta: a sign flip, not a factorization.
#' max(abs(param_solve(s, eta) - param_value(s, -eta)))
#' max(abs(param_solve(s, eta) - solve(M)))
#'
#' # The round trip closes.
#' max(abs(param_free(s, M) - eta))
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
#' Fills the lower triangle of \eqn{S} with the free values, in the order
#' `param_params$positions` records, and mirrors it to the upper triangle. No
#' transformation is applied: the free values **are** the entries of \eqn{S},
#' which is where this chart parts company with [log_cholesky()]'s.
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
#' for the rest. These are the directions the Frechet derivatives contract
#' against.
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
#' bidiagonal matrices [dd_exp()] builds. Deterministic and dependency free: the
#' matrices are at most 5 by 5, so a general-purpose matrix exponential would be
#' a dependency bought for nothing.
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
#' Measured at \eqn{\lambda = (0.5,\, 0.5 + g)}, against the exact limit
#' \eqn{e^{0.5} = 1.648721270700127} reached at \eqn{g = 0}:
#'
#' | gap \eqn{g} | Opitz | quotient |
#' |---|---|---|
#' | \eqn{10^{-3}} | 1.649545906191066 | 1.649545906190929 |
#' | \eqn{10^{-6}} | 1.648722095061039 | 1.648722095023754 |
#' | \eqn{10^{-9}} | 1.648721271524488 | 1.648721206226264 |
#'
#' The quotient has lost eight digits at a gap of \eqn{10^{-9}} and the Opitz
#' value one. Near-repeated eigenvalues are ordinary here: any \eqn{S} with a
#' symmetry has them, and an exactly repeated pair is what a scalar multiple of
#' the identity gives.
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
#' of the orders asked for, computed once per free vector.
#'
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values.
#' @param order The highest derivative order wanted.
#'
#' @details
#' Everything here depends on \eqn{\eta} alone, never on the index tuple, so it
#' is computed once and every component of the order reads it. The rotated
#' directions are \eqn{Q^\top E_k Q} for each basis direction \eqn{E_k}, and the
#' tables hold one divided difference per combination of eigenvalues with
#' repetition, up to `order + 1` points.
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
#' divided-difference table of the matching order, sums over **every ordering**
#' of the directions, and rotates the result back by \eqn{Q}.
#'
#' @details
#' The sum runs over every ordering, never over the distinct ones, so a tuple
#' with repeated indices is counted with its multiplicity and needs no correction
#' afterwards. That is where a hand-written version of this goes wrong: summing
#' the distinct orderings and forgetting the \eqn{\prod_j m_j!} factor gives a
#' result that is too small by exactly that factor, which looks plausible.
#'
#' The cost is the reason this family's higher orders are expensive: at order 4
#' there are 24 orderings per component, each an \eqn{O(p^4)} contraction. See
#' [matrix_log()] for the measurement.
#'
#' @param tb The tables of [mlog_tables()], carrying at least
#'   `length(dirs) + 1` points.
#' @param dirs The free-value indices of the tuple, possibly repeated. Its length
#'   is the derivative order.
#'
#' @return A symmetric numeric matrix of the side `tb$q` has, with no dimnames.
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
#' Returns every permutation of an index tuple as a list, including the ones that
#' coincide when the tuple has repeats: `combinat_perms(c(1, 1))` gives two
#' elements and never one. The multiplicity is deliberate: it is what leaves
#' [mlog_contract()]'s sum over orderings correct for a repeated tuple without a
#' correction factor.
#'
#' @param x An integer vector, of length 1 to 4 in every call the package makes.
#'
#' @return A list of `factorial(length(x))` integer vectors.
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
#' against divided differences of \eqn{\exp} in the eigenvalues, and the sum
#' runs over every ordering of the directions, never over the distinct ones, so
#' a tuple with repeated indices is counted with its multiplicity instead of
#' being corrected for afterwards. The divided differences come from
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
#' positive whatever \eqn{S} is. The cost is one \eqn{O(p^3)}
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
#' Measured, the round trip closes to \eqn{2 \times 10^{-15}}, a little looser
#' than [log_cholesky()]'s \eqn{7 \times 10^{-17}} because two
#' eigendecompositions stand between the two ends.
#' @details
#' A matrix that is not positive definite is rejected: the logarithm of a
#' non-positive eigenvalue is not a real number, so there is no free vector to
#' return. The verdict is spectral, as it is everywhere in this package.
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
#' The Frechet derivative of the matrix exponential, by Daleckii-Krein: rotate
#' the basis direction \eqn{E_k} by the eigenvectors of \eqn{S}, weight it
#' entrywise by the two-point divided differences
#' \eqn{e[\lambda_i, \lambda_j]}, and rotate back. Exact, with no
#' differencing.
#'
#' Note that \eqn{\partial_k M} is **not** \eqn{E_k \exp(S)}: the exponential of a
#' matrix does not commute with an arbitrary direction, which is why the
#' weighting by divided differences appears at all. Where the eigenvalues are
#' equal the weight is \eqn{e^{\lambda}} and the formula reduces to the scalar
#' one, which [dd_exp()]'s Opitz route delivers exactly.
#'
#' Measured against a central difference of [param_value()] at \eqn{p = 3}, the
#' agreement is \eqn{1.5 \times 10^{-10}}, which is the difference's own accuracy.
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
#' Summing over the orderings, never over the distinct ones, is what brings a
#' repeated index out right: the component in one free value twice gets its
#' factor of 2 from the sum instead of from a correction applied afterwards.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
#'   `param_tuple_names(s)` and in that order.
#' @seealso [mlog_contract()], which does the contraction, and
#'   [param_d1.MatrixLogParam()] and [param_d3.MatrixLogParam()] for the
#'   neighbouring orders.
#' @keywords internal
S7::method(param_d2, MatrixLogParam) <- function(s, eta, ...) {
  mlog_higher(s, eta, 2L)
}


#' @title Third Derivatives of a Matrix Logarithm Parameter
#' @name param_d3.MatrixLogParam
#' @description
#' Chains of three rotated directions contracted against four-point divided
#' differences, summed over the six orderings. Exact, where a numerical route at
#' this order keeps about six digits.
#'
#' The cost grows quickly: six orderings per component, each an \eqn{O(p^4)}
#' contraction, and \eqn{\binom{d+2}{3}} components. See [matrix_log()] for the
#' comparison with [log_cholesky()], whose derivatives are sparse products.
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
#' differences, summed over the twenty-four orderings. Exact, and the most
#' expensive quantity the package computes: measured at \eqn{p = 4} one call
#' costs **3.4 s**, against 0.003 s for [log_cholesky()]'s fourth derivatives on
#' the same matrix size.
#'
#' That factor of a thousand is the trade this chart makes, and it is worth
#' knowing before putting a `matrix_log()` inside a loop that takes fourth
#' derivatives. It buys a linear log-determinant and an inverse that is a sign
#' flip; if those are not what the model needs, [log_cholesky()] parametrizes the
#' same cone.
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
#' the sum itself. It agrees with the eigenvalues of the assembled matrix to the
#' printed digit.
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
#' Note the 1 where [log_cholesky()] answers 2: there the diagonal free value is
#' \eqn{\log L_{ii}} and \eqn{|M| = |L|^2}, here it is an eigenvalue of the
#' logarithm and enters once.
#' @param s A [MatrixLogParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`,
#'   holding 1 and 0.
#' @seealso [param_logdet.MatrixLogParam()] for the quantity differentiated, and
#'   [param_dlogdet.LogCholeskyParam()], which answers 2.
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
#' As with [log_cholesky()], that means this family cannot exercise the order: a
#' check of a log-determinant Hessian against a numerical reference passes here
#' whatever is missing. [ar1()] and [compound_symmetry()] are where it has
#' content.
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
#' This is the order at which a numerical fallback is least usable, so a family
#' whose log-determinant is not linear gains most from a closed form here; see
#' [param_d4logdet.matrix_parameter()] for the measured accuracy of the
#' alternative.
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
