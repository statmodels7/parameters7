#' @include numerical_fallbacks.R chain.R
NULL


#' Unstructured Positive Definite Parameter
#'
#' @description
#' The S7 class of unstructured symmetric positive definite matrices in the
#' log-Cholesky parametrization, \eqn{M = L L^\top} with \eqn{L} lower triangular
#' and positive on the diagonal. It is the class the methods of that family
#' dispatch on: given only [param_value()], everything else here is closed form,
#' including all four derivative orders and all four orders of the
#' log-determinant.
#'
#' Call [log_cholesky()] to build one. The class constructor takes the eight
#' properties directly and does no work; using it means computing `free_names`,
#' `rank` and `param_params` by hand.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `LogCholeskyParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. It carries
#'   `dimension`, `rank` (always \eqn{p}), `null_basis` (\eqn{p} by 0), `role`,
#'   `param_name` (`"log_cholesky"`), `n_free` (\eqn{p(p+1)/2}), `free_names` and
#'   `param_params`, whose only entry is `positions`.
#'
#' @seealso [log_cholesky()], the constructor, and [matrix_parameter()] for the
#'   properties this inherits.
#'
#' @examples
#' # The class is what dispatch keys on.
#' s <- log_cholesky(2)
#' c(S7::S7_inherits(s, LogCholeskyParam),
#'   S7::S7_inherits(s, matrix_parameter))
#'
#' # Belonging to it is what gives the family its closed forms: nothing here
#' # is obtained by finite differences.
#' any(param_is_numerical(s))
#'
#' @export
LogCholeskyParam <- S7::new_class("LogCholeskyParam", parent = matrix_parameter)


#' Positions of the Free Values in the Lower Triangle
#'
#' @description
#' Returns the row and column of each free value of a log-Cholesky parameter, in
#' the order the free vector uses: the \eqn{p} diagonal entries first, then the
#' strictly below-diagonal entries column by column. Everything in the family
#' reads it, so the ordering is decided once here and nowhere else.
#'
#' @details
#' The ordering is part of the interface, since `free_names` follows it and
#' consumers build their parameter tables from those names. At \eqn{p = 3} the
#' rows are `1 2 3 2 3 3`, the columns `1 2 3 1 1 2`, and `on_diagonal` is
#' `TRUE TRUE TRUE FALSE FALSE FALSE`, which labels the free values
#' `log_L1 log_L2 log_L3 L2.1 L3.1 L3.2`.
#'
#' @param p The side of the matrix, a single positive integer. At `p = 1` there
#'   are no below-diagonal entries and all three vectors have length 1.
#'
#' @return A list with three vectors of length \eqn{p(p+1)/2}: the integer `row`
#'   and `col` of each free value, and the logical `on_diagonal`.
#'
#' @seealso [log_cholesky()], which turns this into `free_names`, and
#'   [chol_assemble()], which uses it to place the values in \eqn{L}.
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
#' Estimating a covariance matrix directly is awkward. An optimizer moving
#' freely through \eqn{p(p+1)/2} numbers will eventually propose a matrix that is
#' not positive definite, and the likelihood is undefined there. This function
#' removes the problem instead of policing it: it returns an object holding the
#' log-Cholesky map, which sends an unconstrained vector to \eqn{M = L L^\top}
#' with \eqn{L} lower triangular and positive on the diagonal. Every vector in
#' \eqn{\mathbb{R}^{p(p+1)/2}} gives a valid matrix, so the optimizer never has
#' to be told about the constraint.
#'
#' Reach for it when nothing is known about the matrix. It is the parametrization
#' of Pinheiro and Bates (1996).
#'
#' @details
#' # Why the map is safe
#'
#' \deqn{M = L L^\top, \qquad L_{ii} = \exp(\eta_i) > 0}
#'
#' The exponential keeps the diagonal of \eqn{L} positive, and a triangular
#' matrix with a positive diagonal has full rank, so \eqn{L L^\top} is positive
#' definite at every finite \eqn{\eta}. The map is smooth in both directions and
#' one to one, the Cholesky factor with a positive diagonal being unique, so
#' there is no boundary on the free scale to run into.
#'
#' The logarithm on the diagonal belongs to the parametrization and is not a
#' swappable link, which is why it appears in the free names. Use
#' [diagonal_matrix()] if a choice of link is wanted.
#'
#' # Reading the free vector
#'
#' The \eqn{p} logarithms of the diagonal come first, then the strictly
#' below-diagonal entries column by column. For \eqn{p = 3}:
#'
#' \deqn{\eta = (\log L_{11}, \log L_{22}, \log L_{33},
#'              L_{21}, L_{31}, L_{32})}
#'
#' Its length is \eqn{p(p+1)/2}, which is 3, 6, 10 and 15 at \eqn{p} of 2 to 5.
#' The ordering is part of the interface: `free_names` follows it and downstream
#' parameter tables are built from those names. Treat it as fixed.
#'
#' # The log-determinant comes free
#'
#' A likelihood involving \eqn{M} almost always needs \eqn{\log|M|}, and here it
#' is twice the sum of the first \eqn{p} free values:
#'
#' \deqn{\log|M| = 2 \sum_{i=1}^{p} \log L_{ii} = 2 \sum_{i=1}^{p} \eta_i}
#'
#' Linear, so the gradient is the constant 2 in the diagonal directions and 0
#' elsewhere, and the second, third and fourth derivatives are exactly zero. All
#' four are returned without any factorization being taken.
#'
#' # When to use something else
#'
#' This is the parametrization for an unstructured matrix, where
#' \eqn{p(p+1)/2} free values is the price of assuming nothing. A structured
#' family estimates fewer and imposes its structure exactly: [diagonal_matrix()]
#' for independence, [compound_symmetry()] for exchangeable equicorrelation,
#' [ar1()] and [autoregressive()] for a time series, [correlation_matrix()] for a
#' correlation with a separate scale, [matrix_log()] for the same cone through
#' the matrix exponential.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, the point on the unconstrained scale that an
#' optimizer moves, of length \eqn{d = p(p+1)/2}. \eqn{p} is the side of the
#' matrix, \eqn{L} the lower triangular Cholesky factor and \eqn{M = L L^\top}
#' the matrix itself. \eqn{E_{ij}} is the matrix with 1 in position \eqn{(i, j)}
#' and 0 elsewhere.
#'
#' @param dimension The side \eqn{p} of the matrix. A single positive whole
#'   number, finite and at least 1. `0`, `2.5`, `c(1, 2)`, `"3"`, `Inf` and `NA`
#'   all throw `'dimension' must be a single positive integer.`
#' @param role A label recording which side of a model the matrix parametrizes:
#'   `"either"` (the default), `"covariance"` or `"precision"`. **No numeric
#'   result depends on it.** It is carried because the family name does not
#'   record it, the same object serving either side, and a consumer that prefixes
#'   a free name with the matrix it describes needs to know which. Matched with
#'   `match.arg()`, so an abbreviation such as `"cov"` is accepted and anything
#'   else throws.
#'
#' @return An object of class [LogCholeskyParam()], with properties
#'   \describe{
#'     \item{`dimension`}{integer, the side \eqn{p}.}
#'     \item{`n_free`}{integer, \eqn{p(p+1)/2}.}
#'     \item{`free_names`}{character, `log_L1` ... `log_Lp` then `L2.1`, `L3.1`,
#'       `L3.2`, ... in the order above.}
#'     \item{`rank`}{integer, always \eqn{p}: the map is onto the full-rank
#'       cone.}
#'     \item{`null_basis`}{a \eqn{p} by 0 matrix; there are no null directions.}
#'     \item{`role`}{character, as supplied.}
#'     \item{`param_name`}{`"log_cholesky"`.}
#'     \item{`param_params`}{a list with one entry, `positions`, holding the row,
#'       column and diagonal flag of each free value.}
#'   }
#'
#' @references
#' Pinheiro, J. C. and Bates, D. M. (1996). Unconstrained parametrizations for
#' variance-covariance matrices. *Statistics and Computing* **6**, 289-296.
#'
#' @seealso [param_value()] and [param_free()] for the map and its inverse,
#'   [param_factor()], which returns \eqn{L} itself at no cost,
#'   [param_logdet()] and [param_dlogdet()] for the determinant, and
#'   [check_parameter()] to verify a parametrization.
#'
#'   [diagonal_matrix()], [scaled_matrix()], [correlation_matrix()],
#'   [compound_symmetry()], [ar1()], [autoregressive()] and [matrix_log()] for
#'   the alternatives.
#'
#' @examples
#' # A 3 x 3 unstructured covariance: six free values, no constraints.
#' s <- log_cholesky(3)
#' s
#' s@n_free
#' s@free_names
#'
#' # Pick any six numbers at all; the result is positive definite.
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#' M <- param_value(s, eta)
#' round(M, 4)
#' eigen(M, only.values = TRUE)$values
#'
#' # Even absurd values stay in the cone, which is the point.
#' eigen(param_value(s, c(-8, 9, -7, 100, -100, 50)),
#'       only.values = TRUE)$values > 0
#'
#' # The inverse map recovers the free vector exactly.
#' max(abs(param_free(s, M) - eta))
#'
#' # log|M| is linear in eta: twice the sum of the first p entries.
#' all.equal(param_logdet(s, eta), 2 * sum(eta[1:3]))
#' all.equal(param_logdet(s, eta),
#'           sum(log(eigen(M, only.values = TRUE)$values)))
#' param_dlogdet(s, eta)
#'
#' # So every higher derivative of the log-determinant vanishes.
#' c(second = max(abs(param_d2logdet(s, eta))),
#'   third = max(abs(param_d3logdet(s, eta))),
#'   fourth = max(abs(param_d4logdet(s, eta))))
#'
#' # The factor is what the parametrization holds, so it costs nothing.
#' L <- param_factor(s, eta)
#' round(L, 4)
#' all.equal(diag(L), exp(eta[1:3]))
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
#' the first `dimension` values and the rest are placed below it, at the
#' positions [chol_positions()] recorded. Both [param_value()] and
#' [param_factor()] start here, and [param_factor()] returns it unchanged.
#'
#' @details
#' The exponential is applied by subsetting on `on_diagonal`, never through
#' `ifelse()`, which evaluates both branches over the whole vector. Here that
#' would only exponentiate values it then discards; the same shape in
#' [param_free()] would take a logarithm of below-diagonal entries that are free
#' to be negative, and warn about the `NaN`s it throws away.
#'
#' @param s A [LogCholeskyParam()] object, whose `dimension` and
#'   `param_params$positions` are read.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A `s@dimension` by `s@dimension` lower triangular numeric matrix with
#'   a strictly positive diagonal, and no dimnames.
#'
#' @seealso [chol_positions()] for the layout, and
#'   [param_factor.LogCholeskyParam()], which returns this directly.
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
#' @description
#' Returns \eqn{M = L L^\top}, with \eqn{L} assembled from the free vector by
#' [chol_assemble()] and the product taken as `tcrossprod(l)`. Positive
#' definiteness needs no checking: \eqn{L} is triangular with a diagonal that is
#' an exponential, so it has full rank at every finite \eqn{\eta}, and
#' \eqn{L L^\top} is therefore in the interior of the cone. The cost is one
#' \eqn{O(p^3)} product and no factorization.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A symmetric positive definite `s@dimension` by `s@dimension` numeric
#'   matrix, with dimnames `v1`, `v2`, ...
#' @seealso [param_free.LogCholeskyParam()] for the inverse,
#'   [param_factor.LogCholeskyParam()] for \eqn{L} itself, and [log_cholesky()]
#'   for the ordering of `eta`.
#' @keywords internal
S7::method(param_value, LogCholeskyParam) <- function(s, eta, ...) {
  l <- chol_assemble(s, eta)
  name_dims(tcrossprod(l), s)
}


#' @title Factor of a Log-Cholesky Parameter
#' @name param_factor.LogCholeskyParam
#' @description
#' Returns \eqn{L} by assembling it from the free vector, which is all this
#' family's parametrization is: no factorization is taken, because the factor is
#' what \eqn{\eta} holds. It is the only family here whose factor is free, and it
#' is \eqn{O(p^2)} against the base class's \eqn{O(p^3)} Cholesky.
#'
#' The diagonal is `exp(eta[1:p])` and the entries below it are the remaining
#' free values, so a caller who wants a standard deviation off the factor can
#' read it from `eta` directly.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` lower triangular numeric matrix with
#'   a positive diagonal, satisfying `L %*% t(L) == param_value(s, eta)` exactly,
#'   and carrying no dimnames.
#' @seealso [chol_assemble()], which does the assembly, and
#'   [param_factor.matrix_parameter()] for what every other family pays.
#' @keywords internal
S7::method(param_factor, LogCholeskyParam) <- function(s, eta, ...) {
  chol_assemble(s, eta)
}


#' @title Free Vector of a Log-Cholesky Parameter
#' @name param_free.LogCholeskyParam
#' @description
#' Returns the free vector behind a matrix: the Cholesky factor of `m`, read off
#' at the positions [chol_positions()] records, with the diagonal logged. Exact,
#' and a true inverse of [param_value.LogCholeskyParam()] because the triangular
#' factor with a positive diagonal is unique. Measured, the round trip closes to
#' \eqn{7 \times 10^{-17}}.
#' @details
#' `m` is rejected when it is not positive definite, with the verdict taken from
#' the eigenvalues through [chol_pd()] and the message saying so. That
#' distinction matters: a caught [base::chol()] error would be a statement about
#' the arithmetic and can differ between platforms on a matrix with an exactly
#' zero eigenvalue, while a test on the spectrum is a statement about the matrix.
#'
#' The logarithm is applied by subsetting on `on_diagonal`. Through `ifelse()` it
#' would be evaluated over the whole vector, including the below-diagonal entries
#' that are free to be negative, which produces `NaN`s and a warning about values
#' the function then discards.
#' @param s A [LogCholeskyParam()] object.
#' @param m A symmetric positive definite `s@dimension` by `s@dimension` numeric
#'   matrix, already checked for shape and symmetry by the generic. A matrix that
#'   is not positive definite throws.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.LogCholeskyParam()], the map this inverts, and
#'   [chol_pd()] for the definiteness test.
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
#' Closed form, by the Leibniz rule on \eqn{M = L L^\top}. Writing \eqn{L_k} for
#' the derivative of the factor in the \eqn{k}-th free value,
#'
#' \deqn{\partial_k M = L_k L^\top + L L_k^\top.}
#'
#' The factor's derivative is a single-entry matrix: \eqn{L_{ii} E_{ii}} for a
#' diagonal free value, because the parametrization holds its logarithm and
#' \eqn{\partial \exp(\eta_i)/\partial \eta_i = \exp(\eta_i) = L_{ii}}, and
#' \eqn{E_{ij}} for a value below the diagonal, which enters \eqn{L} linearly.
#' Each component therefore costs one row and one column of \eqn{L}, never a
#' matrix product.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` symmetric matrices named by `s@free_names`, each
#'   `s@dimension` by `s@dimension` with dimnames `v1`, `v2`, ...
#' @seealso [chol_leibniz()], which assembles it, [chol_dfactor()] for
#'   \eqn{\partial^S L}, and [param_d2.LogCholeskyParam()] for the order above.
#' @keywords internal
S7::method(param_d1, LogCholeskyParam) <- function(s, eta, ...) {
  chol_leibniz(s, eta, 1L)
}


#' @title Second Derivatives of a Log-Cholesky Parameter
#' @name param_d2.LogCholeskyParam
#' @description
#' Closed form. Differentiating \eqn{\partial_k M = L_k L^\top + L L_k^\top}
#' once more,
#'
#' \deqn{\partial_{kl} M = L_{kl} L^\top + L_k L_l^\top + L_l L_k^\top
#'                       + L L_{kl}^\top,}
#'
#' and the factor's second derivative \eqn{L_{kl}} is non-zero only where
#' \eqn{k = l} names a diagonal free value, where it is \eqn{L_{ii} E_{ii}}
#' again. A below-diagonal value enters \eqn{L} linearly, so any second
#' derivative touching one of those drops the outer terms and leaves the two
#' cross products.
#'
#' One consequence is worth knowing when reading a result: \eqn{M} is quadratic
#' in each below-diagonal free value, so the component keyed by that value twice
#' is a constant in \eqn{\eta}, and every third and fourth derivative repeating
#' it is exactly zero.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
#'   `param_tuple_names(s)` and in that order, each `s@dimension` by
#'   `s@dimension`.
#' @seealso [chol_leibniz()], which assembles it, and
#'   [param_d1.LogCholeskyParam()] for the order below.
#' @keywords internal
S7::method(param_d2, LogCholeskyParam) <- function(s, eta, ...) {
  chol_leibniz(s, eta, 2L)
}


#' @title Log-Determinant of a Log-Cholesky Parameter
#' @name param_logdet.LogCholeskyParam
#' @description
#' Closed form, and linear in the free vector. Since \eqn{|M| = |L|^2} and
#' \eqn{L} is triangular,
#'
#' \deqn{\log|M| = 2\sum_{i=1}^{p} \log L_{ii} = 2\sum_{i=1}^{p} \eta_i,}
#'
#' twice the sum of the free values on the diagonal. It is one `sum()` over
#' \eqn{p} numbers: no factorization, no determinant, no eigendecomposition, and
#' nothing that grows with \eqn{p} beyond the sum itself, against the base
#' class's \eqn{O(p^3)}.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number, always finite: `eta` is finite by the generic's check
#'   and the sum of \eqn{p} finite numbers cannot overflow at any usable \eqn{p}.
#' @seealso [param_dlogdet.LogCholeskyParam()] for its gradient, and
#'   [param_logdet.matrix_parameter()] for the route a family without a closed
#'   form takes.
#' @keywords internal
S7::method(param_logdet, LogCholeskyParam) <- function(s, eta, ...) {
  2 * sum(eta[s@param_params$positions$on_diagonal])
}


#' @title Log-Determinant Gradient of a Log-Cholesky Parameter
#' @name param_dlogdet.LogCholeskyParam
#' @description
#' Closed form and constant: 2 in each of the \eqn{p} diagonal directions and 0
#' in the \eqn{p(p-1)/2} below-diagonal ones, since \eqn{\log|M|} is
#' \eqn{2\sum_i \eta_i} and does not involve the rest of \eqn{\eta} at all. The
#' value does not depend on `eta`, which is read only for its length.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`,
#'   holding 2 and 0.
#' @seealso [param_logdet.LogCholeskyParam()] for the quantity differentiated,
#'   and [param_d2logdet.LogCholeskyParam()], which is zero for the same reason.
#' @keywords internal
S7::method(param_dlogdet, LogCholeskyParam) <- function(s, eta, ...) {
  stats::setNames(
    ifelse(s@param_params$positions$on_diagonal, 2, 0),
    s@free_names
  )
}


#' @title Log-Determinant Hessian of a Log-Cholesky Parameter
#' @name param_d2logdet.LogCholeskyParam
#' @description
#' Closed form, and identically zero: \eqn{\log|M| = 2\sum_i \eta_i} is linear in
#' the free vector, so its second derivative vanishes at every \eqn{\eta}. The
#' zeros are exact, so a consumer can drop the term entirely instead of
#' carrying small numbers through a contraction.
#'
#' This is worth knowing when checking another family: a comparison of a
#' log-determinant Hessian against a numerical reference can pass on
#' `log_cholesky()` while a term is missing, both sides being zero. [ar1()] and
#' [compound_symmetry()] are the families where this order has something to get
#' wrong.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 1, 2)` zeros, keyed as
#'   `param_tuple_names(s)`.
#' @seealso [param_dlogdet.LogCholeskyParam()] for the order below, and
#'   [param_d3logdet.LogCholeskyParam()], zero for the same reason.
#' @keywords internal
S7::method(param_d2logdet, LogCholeskyParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s)
  stats::setNames(rep(0, length(nm)), nm)
}


#' The Factor of a Log-Cholesky Parameter, and Its Derivatives
#'
#' @description
#' Returns \eqn{\partial^S L} for a multiset \eqn{S} of free-value indices, or
#' `NULL` where that derivative is identically zero. Almost all of them are, and
#' that is what keeps the Leibniz sum of [chol_leibniz()] short.
#'
#' @details
#' Three cases exhaust it. An empty \eqn{S} gives \eqn{L} itself. A single index
#' gives a **single-entry matrix**: \eqn{L_{ii} E_{ii}} for a diagonal free
#' value, whose parametrization is a logarithm, and \eqn{E_{ij}} for one below
#' the diagonal, which enters \eqn{L} linearly. A repeated index gives
#' \eqn{L_{ii} E_{ii}} again where every repetition names the same diagonal
#' value, every derivative of \eqn{\exp(\eta_i)} being itself, and `NULL`
#' otherwise.
#'
#' So the only surviving higher derivatives are the pure diagonal ones, and the
#' single-entry shape is what the compiled kernel exploits: each Leibniz term
#' becomes one row, one column or one cell of a product, never a matrix
#' multiplication.
#'
#' @param s A [LogCholeskyParam()] object, whose `dimension` and
#'   `param_params$positions` are read.
#' @param l The factor at the point, from [chol_assemble()].
#' @param ks A multiset of free-value indices, possibly empty. Positions in
#'   `1:s@n_free`; the empty multiset gives \eqn{L} itself.
#'
#' @return A `s@dimension` by `s@dimension` numeric matrix with at most one
#'   non-zero entry, or `l` unchanged for an empty `ks`, or `NULL` where the
#'   derivative is identically zero.
#'
#' @seealso [chol_leibniz()], the caller, and [leibniz_gram()], the generic
#'   Leibniz sum the R twin uses.
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
#' Assembles a whole derivative order of \eqn{M = L L^\top} by the Leibniz rule,
#' in compiled code. Each component of \eqn{\partial^S M} distributes the
#' differentiations of \eqn{S} over the two factors, and since every surviving
#' \eqn{\partial^T L} is a single-entry matrix (see [chol_dfactor()]), each term
#' of the sum is one row, one column or one cell of a product. The kernel walks
#' those entries directly and never multiplies two matrices.
#'
#' The four methods [param_d1()], [param_d2()], [param_d3()] and [param_d4()] of
#' this family are one call each to this function.
#'
#' @details
#' Measured against [.chol_leibniz_r()], the dense R twin, on the same free
#' vector: **30x** at \eqn{p = 8} and order 4 (0.295 s against 8.89 s over 82251
#' components), 30x at \eqn{p = 8} order 2, 35x at \eqn{p = 5} order 4 and 19x at
#' \eqn{p = 3} order 2. The two agree **exactly**, to 0, at \eqn{p = 5} order 4.
#'
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension` with dimnames `v1`, `v2`, ...
#'
#' @seealso [.chol_leibniz_r()], the R twin the tests hold it against,
#'   [chol_dfactor()] for the factor's derivatives, and [leibniz_gram()].
#'
#' @keywords internal
chol_leibniz <- function(s, eta, order) {
  l <- chol_assemble(s, eta)
  idx <- param_tuple_indices(s, order)
  tup <- matrix(unlist(idx, use.names = FALSE),
                ncol = order, byrow = TRUE)
  pos <- s@param_params$positions
  chol_leibniz_cpp(s@dimension, tup, pos$row, pos$col, pos$on_diagonal, l,
                   s@free_names, paste0("v", seq_len(s@dimension)))
}

#' The R Twin of the Compiled Leibniz Assembly
#'
#' @description
#' Computes the same components as `chol_leibniz_cpp`, through the dense matrix
#' products of [leibniz_gram()]. It exists as the independent reference the tests
#' hold that kernel against, so a change to one side that is not a change to both
#' shows up as a disagreement. Not called on any production path;
#' [chol_leibniz()] is.
#'
#' @details
#' The two routes agree **exactly**: measured at \eqn{p = 5} and order 4, over
#' 3060 components, the largest absolute difference is 0. That identity is
#' licensed rather than hoped for, both routes summing the same Leibniz terms,
#' and the compiled one skipping only additions of structural zeros.
#'
#' The compiled route is the production one because every derivative of the
#' factor is a single-entry matrix, so a Leibniz term is one row, one column or
#' one cell where this twin forms a full \eqn{p \times p} product. Measured:
#' 0.295 s against 8.89 s at \eqn{p = 8} and order 4, a factor of 30, and between
#' 19x and 35x over \eqn{p} of 3 to 8 and orders 2 and 4.
#'
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)`, identical to [chol_leibniz()]'s.
#'
#' @seealso [chol_leibniz()], the compiled route this mirrors, and
#'   [leibniz_gram()], the dense Leibniz sum it uses.
#'
#' @keywords internal
.chol_leibniz_r <- function(s, eta, order) {
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
#' Closed form, by the Leibniz rule on \eqn{M = L L^\top}. Each component
#' distributes its three differentiations over the two factors, and a factor
#' differentiated more than once survives only where the repetitions name the
#' same diagonal free value, every derivative of \eqn{e^{\eta_k}} being itself.
#'
#' Two consequences a reader will meet in the output. A component repeating a
#' below-diagonal free value three times is exactly zero, \eqn{M} being quadratic
#' in it. A component repeating a diagonal value three times is not, that value
#' entering through an exponential.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
#'   `param_tuple_names(s, 3)` and in that order.
#' @seealso [chol_leibniz()], which assembles it, and
#'   [param_d4.LogCholeskyParam()] for the order above.
#' @keywords internal
S7::method(param_d3, LogCholeskyParam) <- function(s, eta, ...) {
  chol_leibniz(s, eta, 3L)
}

#' @title Fourth Derivatives of a Log-Cholesky Parameter
#' @name param_d4.LogCholeskyParam
#' @description
#' Closed form, by the same Leibniz rule on \eqn{M = L L^\top} as at third order,
#' with the four differentiations distributed over the two factors. The
#' derivative of \eqn{L} in a multiset of indices is non-zero only for an empty
#' multiset, a single index, or a repetition of one diagonal free value, so of
#' the \eqn{2^4} ways to split a quadruple between the factors almost all
#' contribute nothing.
#'
#' Exact at this order, which is where the contract stops: a fourth-order chain
#' rule through a link needs this much. Nothing is differenced, so the accuracy
#' is machine precision, against the \eqn{10^{-4}} a stencil would give.
#'
#' A component repeating a below-diagonal free value three or more times is
#' exactly zero, \eqn{M} being quadratic in each of those.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
#'   `param_tuple_names(s, 4)` and in that order. At \eqn{p = 3} that is 126
#'   components, at \eqn{p = 8} it is 82251.
#' @seealso [chol_leibniz()], which assembles it, [param_d3.LogCholeskyParam()]
#'   for the order below, and [numerical_d4()] for what a family without a closed
#'   form gets instead.
#' @keywords internal
S7::method(param_d4, LogCholeskyParam) <- function(s, eta, ...) {
  chol_leibniz(s, eta, 4L)
}

#' @title Third Log-Determinant Derivatives of a Log-Cholesky Parameter
#' @name param_d3logdet.LogCholeskyParam
#' @description
#' Closed form, and identically zero. \eqn{\log|M| = 2\sum_i \eta_i} is linear in
#' the free vector, so every derivative above the first vanishes at every
#' \eqn{\eta}. The zeros are exact, so a consumer can drop the term instead of
#' carrying small numbers through a contraction.
#'
#' It also means this family cannot exercise the order: a check of a
#' log-determinant's third derivative against a numerical reference passes here
#' whatever is missing, both sides being zero. [ar1()] and [compound_symmetry()]
#' are where this order has content.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 2, 3)` zeros, keyed as
#'   `param_tuple_names(s, 3)`.
#' @seealso [param_d2logdet.LogCholeskyParam()] for the order below, and
#'   [param_d4logdet.LogCholeskyParam()] for the one above.
#' @keywords internal
S7::method(param_d3logdet, LogCholeskyParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 3L)
  stats::setNames(rep(0, length(nm)), nm)
}

#' @title Fourth Log-Determinant Derivatives of a Log-Cholesky Parameter
#' @name param_d4logdet.LogCholeskyParam
#' @description
#' Closed form, and identically zero, for the reason it is zero at second and
#' third order: \eqn{\log|M| = 2\sum_i \eta_i} is linear in the free vector. The
#' zeros are exact.
#'
#' This is the order at which a numerical fallback is least usable, so a family
#' whose log-determinant is not linear gains most from a closed form here. See
#' [param_d4logdet.matrix_parameter()] for the measured accuracy of the
#' alternative.
#' @param s A [LogCholeskyParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. Its values do not enter the result.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 3, 4)` zeros, keyed as
#'   `param_tuple_names(s, 4)`.
#' @seealso [param_d3logdet.LogCholeskyParam()] for the order below, and
#'   [param_d4logdet.matrix_parameter()] for the numerical route.
#' @keywords internal
S7::method(param_d4logdet, LogCholeskyParam) <- function(s, eta, ...) {
  nm <- param_tuple_names(s, 4L)
  stats::setNames(rep(0, length(nm)), nm)
}
