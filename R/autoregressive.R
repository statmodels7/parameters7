#' @include compound_symmetry.R
NULL


#' Autoregressive Parameter
#'
#' @description
#' The S7 class of covariance matrices of \eqn{p} consecutive observations of a
#' stationary autoregression of order \eqn{q}, parametrized by the marginal
#' variance and the \eqn{q} partial autocorrelations. [autoregressive()] builds
#' one.
#'
#' It is the one family here whose free vector does not grow with the dimension:
#' `n_free` is \eqn{q + 1} at \eqn{p = 6} and \eqn{q + 1} at \eqn{p = 200}. What
#' grows with \eqn{p} is the matrix, not the parametrization.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `AutoregressiveParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `order`, `link_scale` and `link_pacf`. `n_free` is \eqn{q + 1} and `rank` is
#'   `dimension`.
#'
#' @seealso [autoregressive()], the constructor, [ar1()] for the case
#'   \eqn{q = 1} written out, and [matrix_parameter()] for the properties this
#'   inherits.
#'
#' @examples
#' # The free vector does not grow with the dimension.
#' vapply(c(6, 20, 200), function(p) autoregressive(p, order = 2)@n_free,
#'        numeric(1))
#'
#' # One scale and q partial autocorrelations, each on its own chart.
#' autoregressive(20, order = 3)@free_names
#'
#' @export
AutoregressiveParam <- S7::new_class("AutoregressiveParam",
  parent = matrix_parameter)


#' Construct an Autoregressive Parameter
#'
#' @description
#' Returns an object holding the covariance of \eqn{p} consecutive observations
#' of a stationary autoregression of order \eqn{q},
#'
#' \deqn{y_t = \phi_1 y_{t-1} + \cdots + \phi_q y_{t-q} + \varepsilon_t,}
#'
#' parametrized by its marginal variance \eqn{\gamma_0} and its \eqn{q} partial
#' autocorrelations \eqn{r_1, \dots, r_q}. That is \eqn{q + 1} free values
#' **whatever the dimension**: measured, `n_free` is 3 for an order-2 process
#' observed 6 times and 3 for one observed 200 times.
#'
#' @details
#' # Why the partial autocorrelations carry the parametrization
#'
#' The coefficients \eqn{\phi} are stationary exactly when the roots of
#' \eqn{1 - \phi_1 z - \cdots - \phi_q z^{q}} lie outside the unit circle, and
#' that set is not a box. At \eqn{q = 2} it is the open triangle with vertices
#' \eqn{(-2, -1)}, \eqn{(2, -1)} and \eqn{(0, 1)}, of area 4 inside a bounding
#' box \eqn{(-2, 2) \times (-1, 1)} of area 8: **exactly half the box is
#' non-stationary**, and \eqn{\phi = (1.5, 0.6)}, which sits comfortably inside
#' the box, has a root of modulus 0.547. So no collection of scalar links onto
#' intervals can cover the region, whatever intervals are chosen.
#'
#' The partial autocorrelations do not have this problem. Each lies in
#' \eqn{(-1, 1)} independently of the others, and the Levinson-Durbin recursion
#' carries them onto the stationary coefficients bijectively, which is the
#' transformation of Barndorff-Nielsen and Schou (1973) and Monahan (1984). Each
#' takes a [linkfunctions7::rhobit_link()], and every free vector then gives a
#' stationary positive definite matrix.
#'
#' In double precision that last statement has a boundary, and it belongs to the
#' chart, not to this family. At \eqn{\lvert \eta_k \rvert = 6} the matrix is
#' strictly positive definite with an eigenvalue ratio of \eqn{5 \times
#' 10^{-12}}; by \eqn{\lvert \eta_k \rvert = 10} the partial autocorrelation is
#' \eqn{1 - 4 \times 10^{-9}} and the smallest eigenvalue has crossed zero at the
#' rounding floor, \eqn{-3 \times 10^{-17}} of the largest.
#'
#' # The map is polynomial, so nothing is differenced
#'
#' The autocorrelations follow from the same recursion. Writing
#' \eqn{\phi^{(k)}} for the coefficients of the order-\eqn{k} predictor,
#'
#' \deqn{\phi^{(k)}_k = r_k, \qquad
#'   \phi^{(k)}_j = \phi^{(k-1)}_j - r_k \phi^{(k-1)}_{k-j},}
#'
#' and the Yule-Walker equations at order \eqn{k} give
#' \eqn{\rho_k = \sum_{j<k} \phi^{(k)}_j \rho_{k-j} + r_k}, after which
#' \eqn{\rho_h = \sum_j \phi_j \rho_{h-j}} for every lag beyond the order. The
#' whole map from the partial autocorrelations to the matrix is therefore
#' **polynomial**, built from sums and products alone, and its derivatives to
#' fourth order come from propagating the derivative arrays through the recursion
#' in compiled code, the product rule written out per order. Measured against one
#' central difference of [param_value()], the first derivatives agree to
#' \eqn{5 \times 10^{-11}}, which is the difference's own accuracy.
#'
#' # Two quantities are closed form
#'
#' The innovation variances of the Levinson-Durbin recursion give
#'
#' \deqn{\log\lvert M \rvert = p\log\gamma_0
#'   + \sum_{k=1}^{q} (p - k)\log(1 - r_k^{2}),}
#'
#' one term per free value, so the log-determinant is **separable** and every
#' mixed derivative of it is exactly zero: at order 4 and \eqn{q = 2}, 12 of the
#' 15 components are 0 by construction. It agrees with `determinant()` to
#' \eqn{2 \times 10^{-15}} at \eqn{p = 5} and \eqn{4 \times 10^{-14}} at
#' \eqn{p = 100}, and costs a sum of \eqn{q + 1} terms at either size.
#'
#' The inverse is **banded of bandwidth \eqn{q}**: an autoregression of order
#' \eqn{q} is Markov of that order, so its precision carries no entry beyond the
#' \eqn{q}-th off-diagonal. Measured at \eqn{q = 3} and \eqn{p = 9}, every one of
#' the 30 entries outside the band is exactly 0, not merely small. It is
#' assembled from the prediction form \eqn{M^{-1} = U^\top D^{-1} U}, with
#' \eqn{U} unit lower triangular holding the predictor coefficients and \eqn{D}
#' the innovation variances, so no factorization is taken.
#'
#' # What a call costs
#'
#' [ar_taylor()] runs the recursion once and packs **every** order up to the
#' fourth, so `param_d1()` pays most of what `param_d4()` pays and the difference
#' between them is the R-level assembly. Seconds per call, over repetition loops
#' sized by elapsed time:
#'
#' | \eqn{q} | \eqn{p} | `param_value` | `param_d1` | `param_d4` |
#' |---|---|---|---|---|
#' | 1 | 10 | 0.00018 | 0.00024 | 0.00033 |
#' | 1 | 200 | 0.00143 | 0.00191 | 0.00375 |
#' | 2 | 200 | 0.00398 | 0.00516 | 0.01203 |
#' | 4 | 200 | 0.03313 | 0.03563 | 0.08313 |
#'
#' # Against ar1()
#'
#' [ar1()] is the case \eqn{q = 1} written out: there the autocorrelation is
#' \eqn{\rho^{h}}, the determinant is \eqn{(1-\rho^2)^{p-1}} and the inverse is
#' tridiagonal in three lines, so it keeps its own closed forms and does not go
#' through the recursion. The two agree: at \eqn{p = 8} and \eqn{\rho = 0.6} the
#' values differ by \eqn{1.4 \times 10^{-17}}, the log-determinants by 0, the
#' inverses by \eqn{2.2 \times 10^{-16}} and the four derivative orders by
#' \eqn{1.1 \times 10^{-16}} to \eqn{1.1 \times 10^{-14}}. Cost is not the reason
#' to prefer one: they measure 0.00033 s and 0.00050 s for one fourth-order
#' array. `ar1()`'s free name is `z_rho` where this one's is `z_pacf1`, and at
#' \eqn{q = 1} the two are the same number.
#'
#' # The name
#'
#' The name is `autoregressive()` rather than `ar()` because [stats::ar()] is a
#' function of \pkg{stats}, and a package meant to be attached alongside others
#' should not mask one.
#'
#' @section Notation:
#' \eqn{p} is the dimension, \eqn{q} the order, \eqn{\gamma_0} the marginal
#' variance, \eqn{r_k} the \eqn{k}-th partial autocorrelation, \eqn{\phi_j} the
#' autoregressive coefficients, \eqn{\rho_h} the autocorrelation at lag \eqn{h},
#' and \eqn{\eta = (\log \gamma_0, \mathrm{atanh}\, r_1, \dots)} the free vector.
#'
#' @param dimension The side \eqn{p} of the matrix: the number of consecutive
#'   observations. Must exceed `order`, since a stretch of \eqn{q} observations
#'   does not identify \eqn{q} partial autocorrelations; \eqn{p = q + 1} is legal
#'   and is the smallest legal dimension.
#' @param order The order \eqn{q} of the autoregression: a single positive whole
#'   number, no smaller than 1.
#' @param link_scale A \pkg{linkfunctions7} link carrying the marginal variance,
#'   `linkfunctions7::log_link()` by default. It must map onto the positive half
#'   line, so `identity_link()` is rejected, and from the whole real line, which
#'   rules out `sqrt_link()` and its relatives; see [diagonal_matrix()] for the
#'   two conditions.
#' @param role A label recording which side of a model the matrix parametrizes:
#'   `"either"` (the default), `"covariance"` or `"precision"`. No numeric result
#'   depends on it.
#'
#' @return An object of class [AutoregressiveParam()], with `n_free` equal to
#'   \eqn{q + 1}, `free_names` the tagged `log_scale`, `z_pacf1`, ..., `z_pacfq`,
#'   `rank` equal to `dimension`, an empty `null_basis`, and `param_name`
#'   `"ar(q)"`.
#'
#' @references
#' Barndorff-Nielsen, O. and Schou, G. (1973). On the parametrization of
#' autoregressive models by partial autocorrelations. *Journal of Multivariate
#' Analysis* **3**, 408-419.
#'
#' Monahan, J. F. (1984). A note on enforcing stationarity in autoregressive
#' moving average models. *Biometrika* **71**, 403-404.
#'
#' @seealso [ar1()] for \eqn{q = 1} written out, [compound_symmetry()] for the
#'   other two-value family, and [param_readable()], which reports the
#'   autoregressive coefficients this parametrization does not carry.
#'
#' @examples
#' s <- autoregressive(6, order = 2)
#' s@free_names
#'
#' eta <- c(log(2), atanh(0.7), atanh(-0.3))
#' M <- param_value(s, eta)
#' round(M, 4)
#'
#' # The first partial autocorrelation is the lag-one correlation.
#' c(from_matrix = M[1, 2] / M[1, 1], from_eta = tanh(eta[2]))
#'
#' # The precision is banded of bandwidth two, the process being Markov of
#' # order two: everything beyond the second off-diagonal is exactly zero.
#' lag <- abs(outer(1:6, 1:6, "-"))
#' max(abs(param_solve(s, eta)[lag > 2]))
#'
#' # The log-determinant is a sum of three terms, whatever p is.
#' c(closed = param_logdet(s, eta),
#'   from_determinant = as.numeric(determinant(M, logarithm = TRUE)$modulus))
#'
#' # And it is separable, so every mixed derivative of it is exactly zero.
#' tuples <- param_tuple_indices(s, 2L)
#' mixed <- vapply(tuples, function(t) t[1] != t[2], logical(1))
#' max(abs(param_d2logdet(s, eta)[mixed]))
#'
#' # The round trip closes.
#' max(abs(param_free(s, M) - eta))
#'
#' # The coefficients are not free values, and param_readable() reports them
#' # with the Jacobian a delta method needs.
#' param_readable(s, eta)$value
#'
#' @export
autoregressive <- function(dimension, order,
                           link_scale = linkfunctions7::log_link(),
                           role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
    order < 1 || order != round(order)) {
    stop("'order' must be a single positive integer.", call. = FALSE)
  }
  q <- as.integer(order)
  if (p <= q) {
    stop(sprintf(paste0(
      "'dimension' is %d and 'order' is %d. A stretch of %d observations does\n",
      "  not identify %d partial autocorrelations: the dimension must exceed\n",
      "  the order."
    ), p, q, p, q), call. = FALSE)
  }
  check_positive_link(link_scale)
  link_pacf <- linkfunctions7::rhobit_link()

  AutoregressiveParam(
    param_name = sprintf("ar(%d)", q),
    dimension = p,
    n_free = q + 1L,
    free_names = c(tagged_name(link_scale, "scale"),
                   tagged_name(link_pacf, paste0("pacf", seq_len(q)))),
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(
      order = q,
      link_scale = link_scale,
      link_pacf = link_pacf
    )
  )
}


#' The Levinson-Durbin Recursion With Its Derivatives
#'
#' @description
#' Runs the compiled recursion `ar_taylor_cpp`. The scale and the partial
#' autocorrelations enter as their link inverses carrying four derivatives each,
#' and the autocorrelations, the coefficients and every partial derivative to
#' fourth order come back as packed arrays.
#'
#' @details
#' The recursion is sums and products only, so the propagation rules are the
#' product rule written out per order: every derivative is exact and nothing is
#' differenced.
#'
#' It always fills **all four orders**, whatever order the caller went on to
#' want. Measured at \eqn{q = 2}, `gamma` comes back 6 by 121 for \eqn{p = 6},
#' and \eqn{121 = 1 + 3 + 3^2 + 3^3 + 3^4}. That is why `param_d1()` costs nearly
#' what `param_d4()` costs (0.00516 s against 0.01203 s at \eqn{q = 2},
#' \eqn{p = 200}): the difference between them is the R-level assembly, not the
#' recursion.
#'
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list with `n`, the number of free values \eqn{q + 1}; `gamma`, a
#'   matrix with one row per lag \eqn{0, \dots, p-1}; and `phi`, one row per
#'   autoregressive coefficient. Each row packs the value first and then the full
#'   derivative tensors of orders one to four in row-major order, so it has
#'   \eqn{1 + n + n^2 + n^3 + n^4} columns.
#'
#' @seealso [ar_pack_col()] for the column a component sits in, [ar_assemble()]
#'   for the matrix built from one column, and [autoregressive()] for the
#'   recursion itself.
#'
#' @keywords internal
ar_taylor <- function(s, eta) {
  q <- s@param_params$order
  grab <- function(link, e) {
    c(
      linkfunctions7::linkinv(link, e),
      linkfunctions7::dlinkinv(link, e),
      linkfunctions7::d2linkinv(link, e),
      linkfunctions7::d3linkinv(link, e),
      linkfunctions7::d4linkinv(link, e)
    )
  }
  seeds <- rbind(
    grab(s@param_params$link_scale, eta[1L]),
    t(vapply(seq_len(q), function(k) {
      grab(s@param_params$link_pacf, eta[k + 1L])
    }, numeric(5)))
  )
  out <- ar_taylor_cpp(s@dimension, q, seeds)
  out$n <- q + 1L
  out
}


#' The Matrix and Its Derivatives, From the Packed Arrays
#'
#' @description
#' Fills the Toeplitz matrix \eqn{M_{ij} = \gamma_{\lvert i - j \rvert}} from one
#' column of [ar_taylor()]'s packed rows, taking either the value column or one
#' derivative component. Every derivative of the matrix is Toeplitz too, the
#' Toeplitz structure being a property of the family, fixed as the point moves,
#' so one indexing operation serves all five cases.
#'
#' @param s An [AutoregressiveParam()] object, whose `dimension` is read.
#' @param tay The arrays of [ar_taylor()].
#' @param order The derivative order 1 to 4, or 0 for the value.
#' @param tuple The index tuple of that order, ignored at order 0.
#'
#' @return A symmetric `s@dimension` by `s@dimension` numeric matrix, with no
#'   dimnames.
#'
#' @seealso [ar_pack_col()], which locates the column, and [ar_derivative()],
#'   which loops this over the tuples of an order.
#'
#' @keywords internal
ar_assemble <- function(s, tay, order = 0L, tuple = NULL) {
  p <- s@dimension
  lag <- abs(outer(seq_len(p), seq_len(p), "-"))
  vals <- tay$gamma[, ar_pack_col(tay$n, order, tuple)]
  matrix(vals[lag + 1L], p, p)
}


#' The Column of a Packed Derivative Record
#'
#' @description
#' Locates a derivative component in a row of [ar_taylor()]'s output. The value
#' sits in column 1 and the full tensors of orders one to four follow it in
#' row-major order, so order \eqn{k} begins at column
#' \eqn{1 + \sum_{j<k} n^{j}} and the tuple \eqn{(t_1, \dots, t_k)} sits
#' \eqn{\sum_i (t_i - 1) n^{k-i}} further on.
#'
#' @details
#' At \eqn{n = 3} the value is column 1, the three first-order components are
#' columns 2 to 4, and the nine second-order ones are columns 5 to 13:
#' \eqn{(1,1)} is 5, \eqn{(1,2)} is 6 and \eqn{(3,3)} is 13. The tensor is stored
#' **full**, with one column per ordered tuple, so \eqn{(1,2)} and \eqn{(2,1)} are two
#' columns holding the same number; the caller reads whichever one
#' [param_tuple_indices()] hands it.
#'
#' @param n The number of free values, \eqn{q + 1}.
#' @param order The derivative order 1 to 4, or 0 for the value.
#' @param tuple The index tuple, 1-based, of length `order`. Ignored at order 0.
#'
#' @return A single column index.
#'
#' @seealso [ar_taylor()] for the layout this describes.
#'
#' @keywords internal
ar_pack_col <- function(n, order = 0L, tuple = NULL) {
  if (order == 0L) return(1L)
  off <- 1L + cumsum(c(0L, n^(1:3)))[order]
  idx <- 0L
  for (t in tuple) idx <- idx * n + (t - 1L)
  off + idx + 1L
}


#' @title Value of an Autoregressive Parameter
#' @name param_value.AutoregressiveParam
#' @description
#' Returns the Toeplitz matrix \eqn{\gamma_0 \rho_{\lvert i-j \rvert}}, the
#' autocorrelations coming from the Levinson-Durbin recursion of [ar_taylor()]
#' and the marginal variance from the scale link. Positive definite at every free
#' vector, the partial autocorrelations being inside \eqn{(-1, 1)} by
#' construction.
#'
#' The cost is the recursion, which is linear in \eqn{p}: 0.00018 s at
#' \eqn{q = 1, p = 10} and 0.00143 s at \eqn{q = 1, p = 200}.
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A positive definite Toeplitz `s@dimension` by `s@dimension` numeric
#'   matrix with dimnames `v1`, `v2`, ...
#' @seealso [param_free.AutoregressiveParam()] for the inverse, and
#'   [ar_assemble()], which fills it.
#' @keywords internal
S7::method(param_value, AutoregressiveParam) <- function(s, eta, ...) {
  name_dims(ar_assemble(s, ar_taylor(s, eta)), s)
}


#' Derivative Components of an Autoregressive Parameter
#'
#' @description
#' Assembles one derivative order by running [ar_taylor()] once and reading the
#' matching column out of the packed arrays for each tuple of the order. The four
#' methods differ only in the order they pass.
#'
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension` with dimnames.
#'
#' @seealso [ar_taylor()] for the recursion, [ar_assemble()] for one component,
#'   and [param_d1.AutoregressiveParam()], which calls this.
#'
#' @keywords internal
ar_derivative <- function(s, eta, order) {
  tay <- ar_taylor(s, eta)
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(t) {
    name_dims(ar_assemble(s, tay, order, t), s)
  })
  stats::setNames(out, param_tuple_names(s, order))
}


#' @title Derivatives of an Autoregressive Parameter
#' @name param_d1.AutoregressiveParam
#' @description
#' `param_d1()`, `param_d2()`, `param_d3()` and `param_d4()` for an
#' [autoregressive()] parameter, closed form at every order. The map from the
#' partial autocorrelations to the matrix is polynomial, so the derivative arrays
#' propagated through the Levinson-Durbin recursion give each derivative exactly
#' and nothing is differenced. Measured against one central difference of
#' [param_value()], the first derivatives agree to \eqn{5 \times 10^{-11}}, which
#' is the difference's own accuracy.
#'
#' Every component is Toeplitz, the structure being a property of the family and
#' fixed as the point moves.
#' @details
#' The four share [ar_derivative()] and differ only in the order they pass.
#' [ar_taylor()] fills all four orders in one pass whatever is asked, so the
#' first order costs nearly what the fourth costs; see [autoregressive()] for the
#' table.
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of symmetric `s@dimension` by `s@dimension` matrices with
#'   dimnames, keyed as `param_tuple_names(s, order)` and in that order:
#'   `s@n_free` of them at order 1, then `choose(s@n_free + k - 1, k)` at order
#'   \eqn{k}.
#' @seealso [ar_derivative()], which assembles them, and
#'   [param_dlogdet.AutoregressiveParam()] for the log-determinant's own
#'   derivatives.
#' @keywords internal
S7::method(param_d1, AutoregressiveParam) <- function(s, eta, ...) {
  ar_derivative(s, eta, 1L)
}

#' @rdname param_d1.AutoregressiveParam
#' @name param_d2.AutoregressiveParam
#' @keywords internal
S7::method(param_d2, AutoregressiveParam) <- function(s, eta, ...) {
  ar_derivative(s, eta, 2L)
}

#' @rdname param_d1.AutoregressiveParam
#' @name param_d3.AutoregressiveParam
#' @keywords internal
S7::method(param_d3, AutoregressiveParam) <- function(s, eta, ...) {
  ar_derivative(s, eta, 3L)
}

#' @rdname param_d1.AutoregressiveParam
#' @name param_d4.AutoregressiveParam
#' @keywords internal
S7::method(param_d4, AutoregressiveParam) <- function(s, eta, ...) {
  ar_derivative(s, eta, 4L)
}


#' @title Free Vector of an Autoregressive Parameter
#' @name param_free.AutoregressiveParam
#' @description
#' Reads the marginal variance off the common diagonal entry and the partial
#' autocorrelations off the Levinson-Durbin recursion run forwards on the
#' autocorrelations. Exact where `m` is in the set: measured at \eqn{p = 6},
#' \eqn{q = 2}, the round trip closes to \eqn{1.1 \times 10^{-16}}.
#' @details
#' Four things are rejected rather than approximated, each with a message of its
#' own, because the set this family describes is much smaller than the set of
#' symmetric positive definite matrices and a caller who lands outside it has a
#' model error, not a rounding one:
#'
#' - a diagonal that is not constant, which is not stationary;
#' - a matrix that is not Toeplitz, which is rejected instead of being averaged
#'   along its diagonals;
#' - an implied partial autocorrelation at or beyond \eqn{\pm 1}, or a vanishing
#'   prediction variance, neither of which is stationary;
#' - a Toeplitz matrix whose autocorrelations stop following the Yule-Walker
#'   recursion beyond lag \eqn{q}. This last one is what separates an AR(2) from
#'   an AR(3) at \eqn{p \ge 4}, and it is checked by rebuilding the matrix from
#'   the recovered free vector and comparing.
#'
#' The tolerances are relative to `max(1, max(abs(m)))`: \eqn{10^{-8}} for the
#' structure and \eqn{10^{-7}} for the rebuild.
#' @param s An [AutoregressiveParam()] object.
#' @param m The covariance matrix of a stationary autoregression of order
#'   `s@param_params$order`, already checked for shape and symmetry by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.AutoregressiveParam()], the map this inverts.
#' @keywords internal
S7::method(param_free, AutoregressiveParam) <- function(s, m, ...) {
  p <- s@dimension
  q <- s@param_params$order
  scl <- max(1, max(abs(m)))
  d <- diag(m)
  if (diff(range(d)) > 1e-8 * scl) {
    stop("'m' does not have a constant diagonal, so it is not stationary.",
      call. = FALSE)
  }
  if (d[1L] <= 0) stop("'m' has a non-positive diagonal.", call. = FALSE)
  lag <- abs(outer(seq_len(p), seq_len(p), "-"))
  rho <- m[1L, ] / d[1L]
  if (max(abs(m - d[1L] * matrix(rho[lag + 1L], p, p))) > 1e-8 * scl) {
    stop(paste0(
      "'m' is not Toeplitz, so it is not the covariance of a stationary\n",
      "  process. It is rejected rather than averaged along its diagonals."
    ), call. = FALSE)
  }

  r <- numeric(q)
  phi <- numeric(0)
  for (k in seq_len(q)) {
    if (k == 1L) {
      r[1L] <- rho[2L]
    } else {
      num <- rho[k + 1L] - sum(phi * rho[k:2L])
      den <- 1 - sum(phi * rho[2L:k])
      if (abs(den) < 1e-12) {
        stop("'m' is degenerate: a prediction variance vanishes.", call. = FALSE)
      }
      r[k] <- num / den
    }
    if (abs(r[k]) >= 1) {
      stop(paste0(
        "'m' implies a partial autocorrelation outside (-1, 1), so it is not\n",
        "  the covariance of a stationary autoregression."
      ), call. = FALSE)
    }
    phi <- if (k == 1L) r[1L] else c(phi - r[k] * rev(phi), r[k])
  }

  eta <- c(
    linkfunctions7::linkfun(s@param_params$link_scale, d[1L]),
    linkfunctions7::linkfun(s@param_params$link_pacf, r)
  )
  want <- param_value(s, eta)
  if (max(abs(unname(want) - unname(m))) > 1e-7 * scl) {
    stop(sprintf(paste0(
      "'m' is Toeplitz but does not follow the Yule-Walker recursion of an\n",
      "  order-%d autoregression beyond lag %d, so it is not in the set this\n",
      "  parameter describes."
    ), q, q), call. = FALSE)
  }
  stats::setNames(eta, s@free_names)
}


#' The Prediction Form of an Autoregressive Parameter
#'
#' @description
#' Returns the unit lower triangular matrix \eqn{U} of one-step predictor
#' coefficients and the vector \eqn{v} of innovation variances, which factor the
#' matrix as \eqn{M = U^{-1} D U^{-\top}} with \eqn{D = \mathrm{diag}(v)}. It is
#' where [param_solve.AutoregressiveParam()]'s banded precision comes from.
#'
#' @details
#' Row \eqn{t} holds the coefficients of the best linear predictor of \eqn{y_t}
#' from its predecessors, which for \eqn{t} beyond the order are the
#' autoregression's own coefficients, so \eqn{U} has bandwidth \eqn{q}. The
#' innovation variances fall by a factor \eqn{1 - r_k^2} at each of the first
#' \eqn{q} steps and are constant thereafter. At \eqn{p = 6}, \eqn{q = 2},
#' \eqn{\gamma_0 = 2} and partial autocorrelations \eqn{(0.7, -0.3)}, \eqn{v} is
#' \eqn{(2, 1.02, 0.9282, 0.9282, 0.9282, 0.9282)} and every row of \eqn{U} from
#' the third on repeats \eqn{(-0.91, 0.30)}, which are the coefficients
#' \eqn{\phi = (0.91, -0.30)} with a sign.
#'
#' Both come from the same recursion the derivatives use, run here without
#' derivatives, since a solve needs no arrays.
#'
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list with `u`, a unit lower triangular `s@dimension` by
#'   `s@dimension` numeric matrix of bandwidth \eqn{q}, and `v`, a positive
#'   numeric vector of length `s@dimension`.
#'
#' @seealso [param_solve.AutoregressiveParam()], the only caller, and
#'   [autoregressive()] for the Levinson-Durbin recursion.
#'
#' @keywords internal
ar_prediction <- function(s, eta) {
  p <- s@dimension
  q <- s@param_params$order
  r <- vapply(seq_len(q), function(k) {
    linkfunctions7::linkinv(s@param_params$link_pacf, eta[k + 1L])
  }, numeric(1))
  v0 <- linkfunctions7::linkinv(s@param_params$link_scale, eta[1L])

  # the predictor coefficients of every order, from the same recursion
  coef <- vector("list", q)
  phi <- numeric(0)
  for (k in seq_len(q)) {
    phi <- if (k == 1L) r[1L] else c(phi - r[k] * rev(phi), r[k])
    coef[[k]] <- phi
  }

  u <- diag(p)
  for (t in seq_len(p)) {
    k <- min(t - 1L, q)
    if (k >= 1L) u[t, seq.int(t - 1L, t - k)] <- -coef[[k]]
  }
  v <- v0 * cumprod(c(1, (1 - r^2)))[pmin(seq_len(p), q + 1L)]
  list(u = u, v = v)
}


#' @title Solve of an Autoregressive Parameter
#' @name param_solve.AutoregressiveParam
#' @description
#' Returns \eqn{M^{-1} b}, exactly and without a factorization: the precision is
#' \eqn{U^\top D^{-1} U} in the prediction form of [ar_prediction()], assembled
#' and applied.
#' @details
#' The precision is **banded of bandwidth \eqn{q}**, which is the content of an
#' order-\eqn{q} Markov property: no partial correlation beyond the lag.
#' Measured at \eqn{q = 3} and \eqn{p = 9}, all 30 entries outside the band are
#' exactly 0 and 51 of the 81 entries are non-zero; against `solve()` on the
#' assembled matrix the agreement is \eqn{7 \times 10^{-16}}.
#'
#' The bandedness is not exploited for speed here, the whole \eqn{p} by \eqn{p}
#' precision being formed and multiplied. What it buys a consumer is the
#' structure: a model assembling a penalty from this parameter gets a banded
#' block, and a sparse solver can use it.
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param b A numeric matrix with `s@dimension` rows, defaulted to the identity
#'   by the generic, so `param_solve(s, eta)` is the whole precision.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric matrix with `s@dimension` rows and as many columns as `b`.
#' @seealso [ar_prediction()] for the factor, and
#'   [param_logdet.AutoregressiveParam()], which reads the same innovation
#'   variances.
#' @keywords internal
S7::method(param_solve, AutoregressiveParam) <- function(s, eta, b = NULL, ...) {
  pf <- ar_prediction(s, eta)
  (crossprod(pf$u, pf$u / pf$v)) %*% b
}


#' @title Log-Determinant of an Autoregressive Parameter
#' @name param_logdet.AutoregressiveParam
#' @description
#' Closed form, from the innovation variances of the Levinson-Durbin recursion:
#'
#' \deqn{\log\lvert M \rvert = p\log\gamma_0
#'   + \sum_{k=1}^{q} (p - k)\log(1 - r_k^{2}).}
#'
#' A sum of \eqn{q + 1} terms at any dimension, with no factorization and no
#' determinant taken. Measured against `determinant()`, the gap is
#' \eqn{2 \times 10^{-15}} at \eqn{p = 5} and \eqn{4 \times 10^{-14}} at
#' \eqn{p = 100}.
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.AutoregressiveParam()] for its derivatives, and
#'   [ar_prediction()], which returns the same innovation variances.
#' @keywords internal
S7::method(param_logdet, AutoregressiveParam) <- function(s, eta, ...) {
  p <- s@dimension
  q <- s@param_params$order
  r <- vapply(seq_len(q), function(k) {
    linkfunctions7::linkinv(s@param_params$link_pacf, eta[k + 1L])
  }, numeric(1))
  v0 <- linkfunctions7::linkinv(s@param_params$link_scale, eta[1L])
  p * log(v0) + sum((p - seq_len(q)) * log(1 - r^2))
}


#' Log-Determinant Components of an Autoregressive Parameter
#'
#' @description
#' Assembles one derivative order of the log-determinant. It is a sum with one
#' term per free value, so a mixed component is exactly zero and is returned
#' without arithmetic; only the \eqn{q + 1} pure components are computed.
#'
#' @details
#' The scale term is \eqn{p \log \gamma_0}, whose derivatives in
#' \eqn{\gamma_0} are \eqn{(-1)^{k-1}(k-1)!/\gamma_0^{k}}, and each partial
#' autocorrelation contributes \eqn{(p-k)\log(1 - r_k^2)}, whose derivatives come
#' from [log_affine_derivs()] applied to the two factors \eqn{1 - r} and
#' \eqn{1 + r}. Both are then carried onto the free scale by [compose4()], the
#' Faa di Bruno chain with a one-dimensional inner map.
#'
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
#'   as `param_tuple_names(s, order)` and in that order.
#'
#' @seealso [compose4()] and [log_affine_derivs()] for the two pieces, and
#'   [param_dlogdet.AutoregressiveParam()], which calls this.
#'
#' @keywords internal
ar_logdet_derivative <- function(s, eta, order) {
  p <- s@dimension
  q <- s@param_params$order
  ls <- s@param_params$link_scale
  lr <- s@param_params$link_pacf
  chain <- function(link, e, outer_fun) {
    gd <- list(
      linkfunctions7::dlinkinv(link, e),
      linkfunctions7::d2linkinv(link, e),
      linkfunctions7::d3linkinv(link, e),
      linkfunctions7::d4linkinv(link, e)
    )
    compose4(outer_fun(linkfunctions7::linkinv(link, e)), gd)
  }
  d_scale <- chain(ls, eta[1L], function(x) {
    lapply(1:4, function(k) (-1)^(k - 1L) * factorial(k - 1L) / x^k)
  })
  d_pacf <- lapply(seq_len(q), function(k) {
    chain(lr, eta[k + 1L], function(x) {
      log_affine_derivs(x, list(c(p - k, 1, -1), c(p - k, 1, 1)))
    })
  })

  idx <- param_tuple_indices(s, order)
  out <- vapply(idx, function(t) {
    if (any(t != t[1L])) return(0)
    if (t[1L] == 1L) p * d_scale[[order]] else d_pacf[[t[1L] - 1L]][[order]]
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, order))
}


#' @title Log-Determinant Derivatives of an Autoregressive Parameter
#' @name param_dlogdet.AutoregressiveParam
#' @description
#' `param_dlogdet()`, `param_d2logdet()`, `param_d3logdet()` and
#' `param_d4logdet()` for an [autoregressive()] parameter, closed form at every
#' order and **separable**: the log-determinant is a sum with one term per free
#' value, so every mixed component is exactly zero. At \eqn{q = 2} that is 3 of
#' the 6 components at order 2, 7 of 10 at order 3 and 12 of 15 at order 4, and
#' the zeros are exact, not merely small.
#' @details
#' The four share [ar_logdet_derivative()] and differ only in the order they
#' pass. Compare [compound_symmetry()] and [ar1()], whose log-determinants are
#' separable for the same reason, against [correlation_matrix()], where a mixed
#' component has content.
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector keyed as `param_tuple_names(s, order)` and in that
#'   order: `s@n_free` values at order 1, then
#'   `choose(s@n_free + k - 1, k)` at order \eqn{k}.
#' @seealso [ar_logdet_derivative()], which assembles them, and
#'   [param_logdet.AutoregressiveParam()] for the quantity differentiated.
#' @keywords internal
S7::method(param_dlogdet, AutoregressiveParam) <- function(s, eta, ...) {
  ar_logdet_derivative(s, eta, 1L)
}

#' @rdname param_dlogdet.AutoregressiveParam
#' @name param_d2logdet.AutoregressiveParam
#' @keywords internal
S7::method(param_d2logdet, AutoregressiveParam) <- function(s, eta, ...) {
  ar_logdet_derivative(s, eta, 2L)
}

#' @rdname param_dlogdet.AutoregressiveParam
#' @name param_d3logdet.AutoregressiveParam
#' @keywords internal
S7::method(param_d3logdet, AutoregressiveParam) <- function(s, eta, ...) {
  ar_logdet_derivative(s, eta, 3L)
}

#' @rdname param_dlogdet.AutoregressiveParam
#' @name param_d4logdet.AutoregressiveParam
#' @keywords internal
S7::method(param_d4logdet, AutoregressiveParam) <- function(s, eta, ...) {
  ar_logdet_derivative(s, eta, 4L)
}
