#' @include compound_symmetry.R
NULL


#' AR(1) Parameter
#'
#' @description
#' The S7 class of first-order autoregressive covariance matrices: equal
#' variances and a correlation falling geometrically with the lag,
#' \eqn{M_{ij} = \sigma^2 \rho^{|i-j|}}. Two free values at every dimension, and
#' the correlation is unrestricted in \eqn{(-1, 1)}, unlike
#' [compound_symmetry()]'s.
#'
#' [ar1()] builds one. It shares its whole derivative skeleton with
#' [compound_symmetry()] through [econ_scalars()], [econ_derivative()] and
#' [econ_logdet_derivative()]; what differs is the pattern and the
#' log-determinant terms.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `Ar1Param`, a subclass of [matrix_parameter()]
#'   adding no properties of its own. `param_params` holds `link_scale` and
#'   `link_rho`, the second a `linkfunctions7::rhobit_link()` at every dimension.
#'   `n_free` is 2 and `rank` is \eqn{p}.
#'
#' @seealso [ar1()], the constructor, [compound_symmetry()] for the sibling
#'   family, [autoregressive()] for higher orders, and [matrix_parameter()] for
#'   the properties that this class inherits.
#'
#' @examples
#' # Two free values whatever the dimension, and a rhobit correlation.
#' s <- ar1(4)
#' S7::S7_inherits(s, Ar1Param)
#' c(p4 = ar1(4)@n_free, p20 = ar1(20)@n_free)
#' s@free_names
#'
#' # Unlike compound symmetry, the bound does not move with the dimension.
#' t(vapply(c(2, 4, 20), function(p) ar1(p)@param_params$link_rho@link_bounds,
#'          numeric(2)))
#'
#' @export
Ar1Param <- S7::new_class("Ar1Param", parent = matrix_parameter)


#' Construct an AR(1) Parameter
#'
#' @description
#' Returns an object holding the first-order autoregressive covariance
#'
#' \deqn{M(\eta)_{ij} = \sigma^2 \rho^{|i - j|},}
#'
#' with \eqn{\sigma^2} positive and \eqn{\rho \in (-1, 1)}: two free values at
#' every dimension. Use it when the measurements are ordered and the correlation
#' should fall with the distance between them, which is a time series, a spatial
#' transect, or repeated measures where an exchangeable correlation is too
#' strong an assumption.
#'
#' @details
#' # The correlation is unrestricted
#'
#' Unlike [compound_symmetry()], the bound is \eqn{|\rho| < 1} at every
#' dimension, so the correlation is mapped by `linkfunctions7::rhobit_link()`, the
#' inverse hyperbolic tangent, and every free value gives a positive definite
#' matrix in exact arithmetic. In double precision the correlation rounds to
#' \eqn{-1} or 1 from \eqn{|\eta_2|} of about 19, and the matrix is then
#' singular to rounding; the log-determinant, its derivatives and the inverse
#' are evaluated from \eqn{\eta_2} and stay accurate there. A free value of 0 is a correlation of 0 here, where in
#' [compound_symmetry()] it is the midpoint of a dimension-dependent interval.
#'
#' # The inverse is tridiagonal
#'
#' An AR(1) process is Markov, so its precision has no entries beyond the first
#' off-diagonal:
#'
#' \deqn{M^{-1} = \frac{1}{\sigma^2(1-\rho^2)}
#'   \begin{pmatrix}
#'     1 & -\rho & & \\
#'     -\rho & 1+\rho^2 & -\rho & \\
#'      & \ddots & \ddots & \ddots \\
#'      & & -\rho & 1
#'   \end{pmatrix}.}
#'
#' [param_solve()] builds it from that form, so no factorization is performed;
#' the matrix is assembled densely and multiplied by `b`, at \eqn{O(p^2)} per
#' column of `b`.
#'
#' The precision is tridiagonal but is not itself an AR(1) matrix. Its lag-1
#' correlation is not constant along the diagonal (at \eqn{p = 4} and
#' \eqn{\rho = 0.6} it is \eqn{-0.514} at the ends and \eqn{-0.441} in the
#' middle), and its lag-2 correlation is 0, where an AR(1) matrix has
#' \eqn{\rho^2}. An AR(1) covariance and an AR(1) precision therefore define
#' different models. This family is the AR(1) pattern itself, whichever side a
#' consumer puts it on; [ar1_inv()] is the family whose value is the matrix
#' above, so that the AR(1) process is the one written on the other side.
#'
#' # The log-determinant
#'
#' The determinant of the correlation pattern is \eqn{(1-\rho^2)^{p-1}}, so
#'
#' \deqn{\log|M| = p\log\sigma^2 + (p-1)\log(1-\rho^2),}
#'
#' a sum of a function of one free value and a function of the other. Every
#' mixed derivative of it is exactly zero, and no determinant is computed.
#'
#' # Derivatives of the pattern
#'
#' Unlike the pattern of [compound_symmetry()], this pattern is not linear in
#' the correlation. An entry is \eqn{\rho^{m}}
#' for the lag \eqn{m}, so its derivatives in the free value are a power composed
#' with the link, taken to fourth order by [compose_order()]. Each distinct lag is
#' composed once and written into every entry that carries it, the matrix having
#' only \eqn{p} distinct values.
#'
#' @section Notation:
#' \eqn{\eta = (\eta_1, \eta_2)} is the free vector, \eqn{\sigma^2} the common
#' variance and \eqn{\rho} the lag-one correlation. \eqn{p} is the side of the
#' matrix and \eqn{m = |i-j|} the lag of an entry.
#'
#' @param dimension The side \eqn{p} of the matrix, **at least 2**. A one by one
#'   matrix has no correlation, so `ar1(1)` signals an error stating that the
#'   family would carry a free value with no effect.
#' @param link_scale A \pkg{linkfunctions7} link carrying the first free value
#'   onto the positive variance, `linkfunctions7::log_link()` by default. It must
#'   map onto the positive half line, so `identity_link()` is rejected, and from
#'   the whole real line, which rules out `sqrt_link()` and its relatives; see
#'   [diagonal_matrix()] for the two conditions.
#'
#' @return An object of class [Ar1Param()], with `n_free` 2, `free_names`
#'   `log_scale` and `z_rho` under the defaults, `rank` equal to `dimension`, an
#'   empty `null_basis`, `param_name` `"ar1"`, and `param_params` holding
#'   `link_scale` and `link_rho`.
#'
#' @seealso [compound_symmetry()] for an exchangeable correlation,
#'   [autoregressive()] for orders above one, [correlation_matrix()] for an
#'   unrestricted correlation, and [param_solve()] for the tridiagonal inverse.
#'
#' @examples
#' # A 4 x 4 AR(1) covariance at a variance of 2 and a correlation of 0.6.
#' s <- ar1(4)
#' eta <- c(log(2), atanh(0.6))
#' m <- param_value(s, eta)
#' round(m, 4)
#'
#' # The entries really are sigma^2 rho^|i-j|.
#' all.equal(m[1, 4], 2 * 0.6^3, check.attributes = FALSE)
#'
#' # A free value of 0 is a correlation of 0, the link being rhobit.
#' param_value(s, c(0, 0))[1, 2]
#'
#' # The inverse is tridiagonal, and exact: an AR(1) process is Markov.
#' om <- param_solve(s, eta)
#' round(om, 4)
#' max(abs(om[abs(row(om) - col(om)) > 1]))
#' max(abs(om - solve(m)))
#'
#' # But it is not itself AR(1): its lag-1 correlation varies along the
#' # diagonal, and its lag-2 correlation is 0 where rho^2 would be 0.36.
#' d <- sqrt(diag(om))
#' round(om / outer(d, d), 4)
#'
#' # The log-determinant is closed form and agrees with the eigenvalues.
#' c(closed = param_logdet(s, eta),
#'   from_eigen = sum(log(eigen(m, only.values = TRUE)$values)))
#'
#' # The round trip closes exactly, and every mixed second derivative of the
#' # log-determinant is zero.
#' max(abs(param_free(s, m) - eta))
#' param_d2logdet(s, eta)
#'
#' @export
ar1 <- function(dimension,
                link_scale = linkfunctions7::log_link()) {
  p <- check_param_args(dimension)
  if (p < 2L) {
    stop(paste0(
      "'dimension' must be at least 2: a one by one matrix has no\n",
      "  correlation, and the family would carry a free value with no effect."
    ), call. = FALSE)
  }
  check_positive_link(link_scale)
  link_rho <- linkfunctions7::rhobit_link()

  Ar1Param(
    param_name = "ar1",
    dimension = p,
    n_free = 2L,
    free_names = c(tagged_name(link_scale, "scale"),
                   tagged_name(link_rho, "rho")),
    rank = p,
    null_basis = empty_null_basis(p),
    param_params = list(
      link_scale = link_scale,
      link_rho = link_rho
    )
  )
}


#' The Pattern of an AR(1) Parameter
#'
#' @description
#' Returns \eqn{P(\rho)_{ij} = \rho^{|i-j|}}, the correlation pattern of an AR(1)
#' matrix, together with its derivatives in the second free value to the
#' order that the scalars carry. Each
#' entry is a **power** of the correlation, so unlike [cs_pattern()]'s the
#' pattern is not linear and each order needs a genuine chain: [compose_order()]
#' composes \eqn{\rho \mapsto \rho^m} with the link's own derivatives, one
#' order at a time.
#'
#' @details
#' Only \eqn{p} distinct lags occur, so each is composed once and the result
#' written into every entry that carries it. The lag-0 entries are the constant 1,
#' so the diagonal of every derivative is exactly zero.
#'
#' @param s An [Ar1Param()] object, whose `dimension` supplies the lags.
#' @param sc The scalars of [econ_scalars()], whose `rho` component supplies the
#'   correlation and its derivatives in the free value.
#'
#' @return A list of `s@dimension` by `s@dimension` matrices, one more than the
#'   number of derivatives that `sc` carries: the pattern at index 1 and its derivatives at
#'   the following indices, each derivative with a zero diagonal.
#'
#' @seealso [cs_pattern()], the compound-symmetric counterpart, which is linear in
#'   the correlation, [compose_order()] for the chain, and [econ_derivative()], the
#'   caller.
#'
#' @keywords internal
ar1_pattern <- function(s, sc) {
  p <- s@dimension
  lag <- abs(outer(seq_len(p), seq_len(p), "-"))
  r <- sc$rho
  K <- length(r) - 1L
  by_lag <- lapply(0:(p - 1L), function(m) {
    pw <- power_derivs(r[[1L]], m, K)
    c(list(r[[1L]]^m), lapply(seq_len(K), function(o) {
      compose_order(pw, r[-1L], o)
    }))
  })
  lapply(seq_len(K + 1L), function(k) {
    matrix(vapply(as.vector(lag), function(m) by_lag[[m + 1L]][[k]], numeric(1)),
           p, p)
  })
}


#' @title Value of an AR(1) Parameter
#' @name param_value.Ar1Param
#' @description
#' Returns \eqn{M_{ij} = \sigma^2 \rho^{|i-j|}}: the common variance on the
#' diagonal and a covariance falling geometrically with the lag. Positive
#' definiteness holds at every free vector in exact arithmetic, the rhobit link
#' keeping \eqn{|\rho| < 1}, so nothing is tested here; in double precision
#' \eqn{\rho} rounds to \eqn{-1} or 1 from \eqn{|\eta_2|} of about 19.
#' @param s An [Ar1Param()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` symmetric positive definite Toeplitz
#'   matrix with dimnames `v1`, `v2`, ...
#' @seealso [param_free.Ar1Param()] for the inverse map,
#'   [param_solve.Ar1Param()] for the tridiagonal matrix inverse, and
#'   `ar1_pattern()` for the pattern.
#' @keywords internal
S7::method(param_value, Ar1Param) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta, 0L)
  name_dims(sc$scale[[1L]] * ar1_pattern(s, sc)[[1L]], s)
}


#' @title Free Vector of an AR(1) Parameter
#' @name param_free.Ar1Param
#' @description
#' Returns the two free values behind an AR(1) matrix. The variance is the common
#' diagonal entry and the correlation is the first off-diagonal entry divided by
#' it, both read exactly; the remaining entries of the matrix are then compared
#' with the pattern that those two values imply.
#' @details
#' A matrix that does not match the implied pattern is rejected and is never
#' projected onto the nearest AR(1) matrix. A Toeplitz matrix whose lag-2 entry
#' is not \eqn{\sigma^2\rho^2} is not in the family; a caller who wants the
#' nearest AR(1) matrix computes that projection separately and passes the
#' result to this method.
#' @param s An [Ar1Param()] object.
#' @param m An AR(1) numeric matrix of side `s@dimension`, already checked for
#'   shape and symmetry by the generic. Its diagonal must be constant, and every
#'   entry must equal \eqn{\sigma^2\rho^{|i-j|}} for the two values read off the
#'   diagonal and the first off-diagonal.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length 2, named by `s@free_names`.
#' @seealso [param_value.Ar1Param()], the map this inverts.
#' @keywords internal
S7::method(param_free, Ar1Param) <- function(s, m, ...) {
  p <- s@dimension
  d <- diag(m)
  scl <- max(1, max(abs(m)))
  if (diff(range(d)) > 1e-8 * scl) {
    stop("'m' does not have a constant diagonal, so it is not AR(1).",
      call. = FALSE)
  }
  if (d[1L] <= 0) stop("'m' has a non-positive diagonal.", call. = FALSE)
  r <- m[2L, 1L] / d[1L]
  if (abs(r) >= 1) {
    stop("'m' implies a correlation outside (-1, 1), so it is not AR(1).",
      call. = FALSE)
  }
  want <- d[1L] * r^abs(outer(seq_len(p), seq_len(p), "-"))
  if (max(abs(m - want)) > 1e-8 * scl) {
    stop(paste0(
      "'m' does not follow the geometric pattern that its first two entries\n",
      "  imply, so it is not AR(1)."
    ), call. = FALSE)
  }
  stats::setNames(
    c(linkfunctions7::linkfun(s@param_params$link_scale, d[1L]),
      linkfunctions7::linkfun(s@param_params$link_rho, r)),
    s@free_names
  )
}


#' @title Solve of an AR(1) Parameter
#' @name param_solve.Ar1Param
#' @description
#' Exact, and tridiagonal. The precision of an AR(1) covariance is
#' \eqn{\{\sigma^2(1-\rho^2)\}^{-1}} times the matrix with 1 at the two corners
#' of the diagonal, \eqn{1+\rho^2} elsewhere on it, and \eqn{-\rho} on the two
#' first off-diagonals. Every other entry is exactly zero, an AR(1) process being
#' Markov, so the inverse is written down and no factorization is performed. The
#' matrix is assembled densely and multiplied by `b`, at \eqn{O(p^2)} per column
#' of `b`. The factor \eqn{1 - \rho^2} is evaluated from the free value \eqn{z}
#' as \eqn{\mathrm{sech}^2 z = 4e^{-2|z|}/(1 + e^{-2|z|})^2}, so the inverse stays
#' accurate where \eqn{\rho} rounds to \eqn{-1} or 1.
#' @param s An [Ar1Param()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param b A numeric matrix with `s@dimension` rows, already coerced from a
#'   vector and defaulted to the identity by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric matrix with `s@dimension` rows and as many columns as `b`.
#' @seealso [ar1()], which explains why the tridiagonal precision is not itself
#'   AR(1), and [param_solve.matrix_parameter()] for the factorization this
#'   avoids.
#' @keywords internal
S7::method(param_solve, Ar1Param) <- function(s, eta, b = NULL, ...) {
  p <- s@dimension
  sc <- econ_scalars(s, eta, 0L)
  v <- sc$scale[[1L]]
  r <- sc$rho[[1L]]
  # 1 - rho^2 = sech(z)^2, written in exp(-2|z|) so that it does not cancel
  # where rho rounds to -1 or 1
  a <- exp(-2 * abs(eta[2L]))
  one_minus <- 4 * a / (1 + a)^2
  t_mat <- diag(c(1, rep(1 + r^2, p - 2L), 1), nrow = p)
  off <- cbind(seq_len(p - 1L), seq.int(2L, p))
  t_mat[off] <- -r
  # drop = FALSE, because at p = 2 the index matrix has one row and would
  # otherwise collapse to a length-two vector, which R reads as two LINEAR
  # positions rather than as one row-column pair -- writing into [1, 1].
  t_mat[off[, c(2L, 1L), drop = FALSE]] <- -r
  (t_mat / (v * one_minus)) %*% b
}


#' @title Derivatives of an AR(1) Parameter
#' @name param_d1.Ar1Param
#' @description
#' Closed form at all four orders. One page covers [param_d1()], [param_d2()],
#' [param_d3()] and [param_d4()] for this family because the four are the same
#' product at different orders: \eqn{M = \sigma^2 P(\rho)}, so a component with
#' \eqn{a} scale indices and \eqn{b} correlation indices is the \eqn{a}-th
#' derivative of the scale times the \eqn{b}-th derivative of the pattern.
#'
#' The pattern's derivatives are the ones that need work here. An entry is
#' \eqn{\rho^{|i-j|}}, a power where [compound_symmetry()]'s pattern is linear,
#' so each is composed with the rhobit link to the requested order, through
#' [compose_order()]. Each distinct lag is composed once.
#' @details
#' The four methods return lists of `order + 1` matrices, keyed by the tuple names
#' of their own order:
#'
#' - [param_d1()]: 2 entries, keyed by `free_names`.
#' - [param_d2()]: 3 entries, `param_tuple_names(s, 2)`.
#' - [param_d3()]: 4 entries, `param_tuple_names(s, 3)`.
#' - [param_d4()]: 5 entries, `param_tuple_names(s, 4)`.
#'
#' Every component with at least one correlation index has a zero diagonal, the
#' lag-0 entries of the pattern being the constant 1.
#' @param s An [Ar1Param()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named list of `order + 1` symmetric `s@dimension` by `s@dimension`
#'   matrices; see **Details** for the keying of each order.
#' @seealso [econ_derivative()], which assembles all four, `ar1_pattern()` for
#'   the pattern and its derivatives, and
#'   [param_d1.CompoundSymmetryParam()] for the sibling family, whose pattern is
#'   linear.
#' @keywords internal
S7::method(param_d1, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 1L, ar1_pattern)
}

#' @rdname param_d1.Ar1Param
#' @name param_d2.Ar1Param
#' @keywords internal
S7::method(param_d2, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 2L, ar1_pattern)
}

#' @rdname param_d1.Ar1Param
#' @name param_d3.Ar1Param
#' @keywords internal
S7::method(param_d3, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 3L, ar1_pattern)
}

#' @rdname param_d1.Ar1Param
#' @name param_d4.Ar1Param
#' @keywords internal
S7::method(param_d4, Ar1Param) <- function(s, eta, ...) {
  econ_derivative(s, eta, 4L, ar1_pattern)
}


#' Log-Determinant Chain of an AR(1) Parameter
#'
#' @description
#' Returns a function of the second free value \eqn{z} and an order that gives
#' the order-th derivative of \eqn{q = (p-1)\log(1-\rho^2)} in \eqn{z}, with
#' \eqn{\rho = \tanh z} under the rhobit link.
#'
#' @details
#' With \eqn{t = \tanh z} and \eqn{u = 1 - t^2 = \mathrm{sech}^2 z},
#' \eqn{q = -2(p-1)\log\cosh z}, and from \eqn{t' = u} and \eqn{u' = -2tu} the
#' four derivatives are \eqn{-2(p-1)t}, \eqn{-2(p-1)u}, \eqn{4(p-1)tu} and
#' \eqn{4(p-1)u(u - 2t^2)}. The factor \eqn{u} is evaluated as
#' \eqn{4e^{-2|z|}/(1 + e^{-2|z|})^2}, never as \eqn{1 - t^2}, which loses all
#' accuracy once \eqn{t} rounds to \eqn{-1} or 1.
#'
#' @param s An [Ar1Param()] object, whose `dimension` supplies \eqn{p}.
#'
#' @return A function of `(e, order)`, with `e` a single number and `order` an
#'   integer from 1 to 4, returning a single number.
#'
#' @seealso [econ_logdet_derivative()], which calls it, and [cs_logdet_chain()]
#'   for the compound-symmetric counterpart.
#'
#' @keywords internal
ar1_logdet_chain <- function(s) {
  p <- s@dimension
  function(e, order) {
    t <- tanh(e)
    a <- exp(-2 * abs(e))
    u <- 4 * a / (1 + a)^2
    switch(order,
      -2 * (p - 1) * t,
      -2 * (p - 1) * u,
      4 * (p - 1) * t * u,
      4 * (p - 1) * u * (u - 2 * t^2))
  }
}


#' @title Log-Determinant of an AR(1) Parameter
#' @name param_logdet.Ar1Param
#' @description
#' Closed form. The determinant of the correlation pattern is
#' \eqn{(1-\rho^2)^{p-1}}, so
#'
#' \deqn{\log|M| = p\log\sigma^2 + (p-1)\log(1-\rho^2).}
#'
#' Two logarithms and no factorization, whatever \eqn{p} is, where the base
#' method takes an \eqn{O(p^3)} eigendecomposition. With \eqn{\rho = \tanh z},
#' the second term is evaluated as
#' \eqn{(p-1)\{2\log 2 - 2|z| - 2\log(1 + e^{-2|z|})\}}, which does not lose
#' accuracy where \eqn{\rho} rounds to \eqn{-1} or 1.
#' @param s An [Ar1Param()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number, finite at every free vector.
#' @seealso [param_dlogdet.Ar1Param()] for its four derivative orders, and
#'   [ar1_logdet_chain()] for the derivatives of the correlation's term.
#' @keywords internal
S7::method(param_logdet, Ar1Param) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta, 0L)
  p <- s@dimension
  z <- abs(eta[2L])
  p * log(sc$scale[[1L]]) + (p - 1) * (2 * log(2) - 2 * z - 2 * log1p(exp(-2 * z)))
}

#' @title Log-Determinant Derivatives of an AR(1) Parameter
#' @name param_dlogdet.Ar1Param
#' @description
#' Closed form at all four orders. One page covers [param_dlogdet()],
#' [param_d2logdet()], [param_d3logdet()] and [param_d4logdet()] because the four
#' differ only in the order taken:
#' \eqn{\log|M| = p\log\sigma^2 + (p-1)\log(1-\rho^2)} is a **sum** of a function
#' of one free value and a function of the other, so every mixed component is
#' exactly zero and the pure ones are the two chains taken separately.
#' @details
#' The scale's chain is \eqn{p} times the derivatives of \eqn{\log h}, which under
#' the default log link is \eqn{p} at first order and 0 above it. The
#' correlation's is [ar1_logdet_chain()], the derivatives in the free value
#' written in closed form in \eqn{\tanh z} and \eqn{\mathrm{sech}^2 z}, so that
#' they stay accurate where \eqn{\rho} rounds to \eqn{-1} or 1.
#'
#' The four methods return vectors of `order + 1` entries, keyed by the tuple
#' names of their own order. At \eqn{p = 4} and \eqn{\rho = 0.6} the second order
#' is \eqn{(0, -3.84, 0)} over `log_scale:log_scale`, `z_rho:z_rho` and
#' `log_scale:z_rho`.
#'
#' The higher log-determinant orders are not zero here, as for
#' [compound_symmetry()], [autoregressive()] and [correlation_matrix()], so a
#' check of them compares non-zero numbers, where on [log_cholesky()] it would
#' compare two zeros.
#' @param s An [Ar1Param()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named numeric vector of `order + 1` entries, keyed as
#'   `param_tuple_names(s, order)`, with every mixed entry exactly zero.
#' @seealso [param_logdet.Ar1Param()] for the quantity differentiated,
#'   [econ_logdet_derivative()], which assembles all four, and
#'   [param_dlogdet.CompoundSymmetryParam()] for the sibling family.
#' @keywords internal
S7::method(param_dlogdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 1L, ar1_logdet_chain(s))
}

#' @rdname param_dlogdet.Ar1Param
#' @name param_d2logdet.Ar1Param
#' @keywords internal
S7::method(param_d2logdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 2L, ar1_logdet_chain(s))
}

#' @rdname param_dlogdet.Ar1Param
#' @name param_d3logdet.Ar1Param
#' @keywords internal
S7::method(param_d3logdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 3L, ar1_logdet_chain(s))
}

#' @rdname param_dlogdet.Ar1Param
#' @name param_d4logdet.Ar1Param
#' @keywords internal
S7::method(param_d4logdet, Ar1Param) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 4L, ar1_logdet_chain(s))
}
