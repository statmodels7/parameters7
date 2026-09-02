#' @include inverse_of.R
NULL

#' The Precision of an AR(1) Process
#'
#' @description
#' The S7 class of the family whose value is the **inverse** of an [ar1()]
#' matrix: tridiagonal, the precision of a first-order autoregression.
#' [ar1_inv()] builds one.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `Ar1InvParam`, a subclass of [InverseParam()]
#'   adding no properties of its own. `param_params` holds `inner`, the
#'   [ar1()] family being inverted.
#'
#' @seealso [ar1_inv()], the constructor, [ar1()] for the family it inverts,
#'   and [inverse_of()] for the general composition it specializes.
#'
#' @examples
#' s <- ar1_inv(5)
#' s@free_names
#' round(param_value(s, c(0, atanh(0.6))), 4)
#'
#' @export
Ar1InvParam <- S7::new_class("Ar1InvParam", parent = InverseParam)


#' Construct the Precision of an AR(1) Process
#'
#' @description
#' Returns the family whose value is \eqn{\Sigma(\eta)^{-1}} for
#' \eqn{\Sigma} an [ar1()] matrix: the **tridiagonal** precision of a
#' first-order autoregression. It carries [ar1()]'s free vector unchanged, so
#' its two coordinates are the scale and the correlation of the process whose
#' precision this is.
#'
#' @details
#' # Which side is which
#'
#' [ar1()] is the AR(1) **pattern**, \eqn{\sigma^2\rho^{|i-j|}}. A family that
#' takes a matrix parameter decides what it does with it, so
#' `mvgaussian1_distrib(p, ar1(p))` is the process whose covariance is AR(1)
#' and `mvgaussian2_distrib(p, ar1(p))` the model whose precision has those
#' entries, which is a different and legitimate model. This family is the third
#' object in that picture: the matrix whose inverse is an AR(1), so
#' `mvgaussian2_distrib(p, ar1_inv(p))` is again the AR(1) process, written on
#' the precision side.
#'
#' # The value
#'
#' A first-order autoregression is Markov, so its precision carries no entry
#' beyond the first off-diagonal. With \eqn{\tau = \sigma^{-2}} and
#' \eqn{w = (1-\rho^2)^{-1}},
#'
#' \deqn{\Omega = \tau\,G(\rho), \qquad
#'   G_{11} = G_{pp} = w, \quad G_{ii} = (1+\rho^2)w, \quad
#'   G_{i,i\pm 1} = -\rho w,}
#'
#' every other entry exactly zero. The value is read from
#' [param_solve.Ar1Param()], which writes that matrix down and performs no
#' factorization, so there is one copy of it in the package.
#'
#' # The derivatives
#'
#' The value is a **product** of a function of the first free value and a
#' matrix function of the second, so a component with \eqn{a} scale indices and
#' \eqn{b} correlation indices is \eqn{(\partial^a\tau)(\partial^b G)} and no
#' mixed expansion is needed. The scale factor is the reciprocal of the link's
#' inverse, chained through \eqn{x\mapsto x^{-1}}; the pattern's three distinct
#' entries are written in \eqn{w} alone,
#'
#' \deqn{w^{(k)} = \frac{k!}{2}\left\{(1-\rho)^{-(k+1)}
#'   + (-1)^k (1+\rho)^{-(k+1)}\right\},}
#'
#' with \eqn{G_{11}^{(k)} = w^{(k)}}, \eqn{G_{ii} = 2w - 1} so
#' \eqn{G_{ii}^{(k)} = 2w^{(k)}} above order zero, and
#' \eqn{G_{i,i\pm1}^{(k)} = -(\rho\,w^{(k)} + k\,w^{(k-1)})} by the Leibniz
#' rule. Each is then chained onto the free value through the correlation's
#' link.
#'
#' This is what the family adds over `inverse_of(ar1(p))`, which reaches the
#' same numbers through the ordered-block-partition sum: there a fourth-order
#' component is 75 products of five matrices, here it is one elementwise
#' product. The two agree to machine precision, which is what a test asserts.
#'
#' # The log-determinant
#'
#' \eqn{\log\lvert\Omega\rvert = -\log\lvert\Sigma\rvert}, so the value and all
#' four orders are [ar1()]'s negated, and the closed form
#' \eqn{p\log\sigma^2 + (p-1)\log(1-\rho^2)} is what is negated.
#'
#' @section Notation:
#' \eqn{p} is the matrix side, \eqn{\sigma^2} the AR(1) process's variance,
#' \eqn{\rho} its lag-one correlation, \eqn{\tau = \sigma^{-2}} and
#' \eqn{w = (1-\rho^2)^{-1}}.
#'
#' @param dimension The side of the matrix, at least 2.
#' @param link_scale The link carrying the AR(1) process's variance, passed to
#'   [ar1()]. Must map onto the positive half line.
#'
#' @return An object of class [Ar1InvParam()], with `n_free` 2, `free_names`
#'   [ar1()]'s and `rank` equal to `dimension`.
#'
#' @seealso [ar1()] for the family this inverts, [autoregressive_inv()] for the
#'   order-\eqn{q} case, and [inverse_of()] for any other family's inverse.
#'
#' @examples
#' s <- ar1_inv(5)
#' eta <- c(log(2), atanh(0.6))
#'
#' # The value is tridiagonal: an AR(1) is Markov, so the precision carries no
#' # entry beyond the first off-diagonal.
#' round(param_value(s, eta), 4)
#'
#' # And it is exactly the inverse of the AR(1) it names.
#' max(abs(param_value(s, eta) %*% param_value(ar1(5), eta) - diag(5)))
#'
#' # The written-out derivatives agree with the general composition's.
#' g <- inverse_of(ar1(5))
#' max(abs(unlist(param_d4(s, eta)) - unlist(param_d4(g, eta))))
#'
#' @export
ar1_inv <- function(dimension, link_scale = linkfunctions7::log_link()) {
  inner <- ar1(dimension, link_scale = link_scale)
  p <- inner@dimension
  Ar1InvParam(
    param_name = "ar1_inv",
    n_free = inner@n_free,
    free_names = inner@free_names,
    dimension = p,
    rank = p,
    null_basis = empty_null_basis(p),
    param_params = list(inner = inner)
  )
}


#' Derivatives of a Reciprocal, for Composition
#'
#' @description
#' Returns the value and four derivatives of \eqn{\eta \mapsto 1/h(\eta)} from
#' the value and four derivatives of \eqn{h}, by composing \eqn{x^{-1}} onto
#' them.
#'
#' @details
#' [power_derivs()] cannot serve here: it truncates beyond the exponent, which
#' is right for a non-negative integer power and wrong for \eqn{-1}, where
#' every order is non-zero. The outer derivatives written out are
#' \eqn{-x^{-2}}, \eqn{2x^{-3}}, \eqn{-6x^{-4}} and \eqn{24x^{-5}}.
#'
#' @param v A list of five numbers: the value of \eqn{h} and its four
#'   derivatives, as `econ_scalars()` returns for one link.
#'
#' @return A list of five numbers in the same shape, for \eqn{1/h}.
#'
#' @seealso [ar1_inv()], the caller, and [compose4()] for the chain.
#'
#' @keywords internal
recip_derivs <- function(v) {
  x <- v[[1L]]
  outer_d <- list(-x^-2, 2 * x^-3, -6 * x^-4, 24 * x^-5)
  c(list(1 / x), compose4(outer_d, v[-1L]))
}


#' Derivatives of the Reciprocal of One Minus a Squared Correlation
#'
#' @description
#' Returns the value and four derivatives of
#' \eqn{\eta \mapsto w = (1 - \rho(\eta)^2)^{-1}}, the quantity every inverse
#' autoregression is written in: it is the reciprocal of the innovation
#' variance's share at one lag, and it multiplies every entry of an inverse
#' AR(1) and every diagonal of an inverse AR(\eqn{q}).
#'
#' @details
#' Written by partial fractions,
#' \eqn{w = \tfrac{1}{2}\{(1-\rho)^{-1} + (1+\rho)^{-1}\}}, every order is an
#' exact expression rather than a repeated quotient rule,
#'
#' \deqn{w^{(k)} = \frac{k!}{2}\left\{(1-\rho)^{-(k+1)}
#'   + (-1)^k (1+\rho)^{-(k+1)}\right\},}
#'
#' and it is finite throughout \eqn{|\rho| < 1}, which the correlation's link
#' guarantees. The result is then chained onto the free value.
#'
#' @param v A list of five numbers: \eqn{\rho} and its four derivatives in the
#'   free value, as `econ_scalars()` returns for one link.
#'
#' @return A list of five numbers: \eqn{w} and its four derivatives in the free
#'   value.
#'
#' @seealso [ar1_inv()] and [autoregressive_inv()], which share it.
#'
#' @keywords internal
w_derivs <- function(v) {
  r <- v[[1L]]
  d <- vapply(0:4, function(k) {
    factorial(k) / 2 * ((1 - r)^-(k + 1L) + (-1)^k * (1 + r)^-(k + 1L))
  }, numeric(1))
  c(list(d[1L]), compose4(as.list(d[-1L]), v[-1L]))
}


#' The Pattern of an Inverse AR(1) Parameter
#'
#' @description
#' Returns \eqn{G(\rho)}, the tridiagonal correlation pattern of the precision
#' of an AR(1), together with its four derivatives in the second free value.
#'
#' @details
#' The three distinct entries are \eqn{w}, \eqn{2w-1} and \eqn{-\rho w} for
#' \eqn{w = (1-\rho^2)^{-1}}, so one sequence of derivatives of \eqn{w} in
#' \eqn{\rho} serves all three. Writing \eqn{w} by partial fractions as
#' \eqn{\tfrac{1}{2}\{(1-\rho)^{-1} + (1+\rho)^{-1}\}} makes every order an
#' exact expression rather than a repeated quotient rule, and it is finite
#' throughout \eqn{|\rho| < 1}, which the rhobit link guarantees.
#'
#' @param s An [Ar1InvParam()] object, whose `dimension` supplies \eqn{p}.
#' @param sc The scalars of [econ_scalars()] read on the inner [ar1()], whose
#'   `rho` entry supplies the correlation and its link's derivatives.
#'
#' @return A list of five `s@dimension` square matrices: the pattern and its
#'   four derivatives in the second free value.
#'
#' @seealso [ar1_inv()] for the formulas and `ar1_pattern()` for the
#'   counterpart on the covariance side.
#'
#' @keywords internal
ar1_inv_pattern <- function(s, sc) {
  p <- s@dimension
  r <- sc$rho[[1L]]
  # w and its four derivatives in RHO, by partial fractions. The chained
  # version is w_derivs(); here the derivatives in rho are wanted first,
  # because the other two entries are built from them before the chain.
  w <- vapply(0:4, function(k) {
    factorial(k) / 2 * ((1 - r)^-(k + 1L) + (-1)^k * (1 + r)^-(k + 1L))
  }, numeric(1))

  # the three distinct entries, value then four derivatives in rho
  corner <- w
  interior <- c(2 * w[1L] - 1, 2 * w[-1L])
  off <- vapply(0:4, function(k) {
    -(r * w[k + 1L] + if (k >= 1L) k * w[k] else 0)
  }, numeric(1))

  # chain each onto the free value through the correlation's link
  chain <- function(d) c(list(d[1L]), compose4(as.list(d[-1L]), sc$rho[-1L]))
  cn <- chain(corner)
  it <- chain(interior)
  of <- chain(off)

  lag <- abs(outer(seq_len(p), seq_len(p), "-"))
  is_corner <- outer(seq_len(p), seq_len(p), function(i, j) {
    i == j & (i == 1L | i == p)
  })
  lapply(1:5, function(k) {
    m <- matrix(0, p, p)
    m[lag == 1L] <- of[[k]]
    m[lag == 0L] <- it[[k]]
    m[is_corner] <- cn[[k]]
    m
  })
}


#' Derivative Arrays of an Inverse AR(1) Parameter
#'
#' @description
#' Assembles one order from the product structure \eqn{\Omega = \tau G(\rho)}:
#' a component with \eqn{a} scale indices and \eqn{b} correlation indices is
#' the \eqn{a}-th derivative of the reciprocal scale times the \eqn{b}-th
#' derivative of the pattern.
#'
#' @param s An [Ar1InvParam()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A named list of `s@dimension` square matrices, keyed and ordered as
#'   [param_tuple_names()] says.
#'
#' @seealso [ar1_inv()] for the formulas and `econ_derivative()`, the same
#'   assembly on the covariance side.
#'
#' @keywords internal
ar1_inv_derivative <- function(s, eta, order) {
  inner <- .inv_inner(s)
  sc <- econ_scalars(inner, eta)
  tau <- recip_derivs(sc$scale)
  pt <- ar1_inv_pattern(s, sc)
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(t) {
    a <- sum(t == 1L)
    b <- sum(t == 2L)
    name_dims(tau[[a + 1L]] * pt[[b + 1L]], s)
  })
  stats::setNames(out, param_tuple_names(s, order))
}


#' @title Derivative Arrays of an Inverse AR(1) Parameter
#' @name param_d1.Ar1InvParam
#' @description
#' The four orders of \eqn{\partial\Omega} for \eqn{\Omega} the precision of an
#' AR(1), written out from the product structure rather than assembled by the
#' ordered-block-partition sum [inverse_of()] uses. The two routes agree to
#' machine precision.
#' @param s An [Ar1InvParam()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named list of `s@dimension` square matrices, keyed as
#'   [param_tuple_names()] says.
#' @seealso [ar1_inv()] for the formulas.
#' @keywords internal
S7::method(param_d1, Ar1InvParam) <- function(s, eta, ...) {
  ar1_inv_derivative(s, eta, 1L)
}

#' @rdname param_d1.Ar1InvParam
#' @name param_d2.Ar1InvParam
#' @keywords internal
S7::method(param_d2, Ar1InvParam) <- function(s, eta, ...) {
  ar1_inv_derivative(s, eta, 2L)
}

#' @rdname param_d1.Ar1InvParam
#' @name param_d3.Ar1InvParam
#' @keywords internal
S7::method(param_d3, Ar1InvParam) <- function(s, eta, ...) {
  ar1_inv_derivative(s, eta, 3L)
}

#' @rdname param_d1.Ar1InvParam
#' @name param_d4.Ar1InvParam
#' @keywords internal
S7::method(param_d4, Ar1InvParam) <- function(s, eta, ...) {
  ar1_inv_derivative(s, eta, 4L)
}


#' The Precision of an Autoregression of Order q
#'
#' @description
#' The S7 class of the family whose value is the **inverse** of an
#' [autoregressive()] matrix: banded of bandwidth \eqn{q}, the precision of an
#' autoregression of that order. [autoregressive_inv()] builds one.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `AutoregressiveInvParam`, a subclass of
#'   [InverseParam()] adding no properties of its own. `param_params` holds
#'   `inner`, the [autoregressive()] family being inverted.
#'
#' @seealso [autoregressive_inv()], the constructor, and [ar1_inv()] for the
#'   order-one case.
#'
#' @examples
#' s <- autoregressive_inv(6, order = 2)
#' s@free_names
#'
#' @export
AutoregressiveInvParam <- S7::new_class("AutoregressiveInvParam",
                                        parent = InverseParam)


#' Construct the Precision of an Autoregression of Order q
#'
#' @description
#' Returns the family whose value is \eqn{\Sigma(\eta)^{-1}} for \eqn{\Sigma}
#' an [autoregressive()] matrix: the precision of an autoregression of order
#' \eqn{q}, **banded of bandwidth \eqn{q}**. It carries
#' [autoregressive()]'s free vector unchanged, so its coordinates are the scale
#' and the partial autocorrelations of the process whose precision this is.
#'
#' @details
#' # The value and the log-determinant are closed
#'
#' An autoregression of order \eqn{q} is Markov of that order, so its precision
#' carries no entry beyond the \eqn{q}-th diagonal. The value comes from
#' [param_solve.AutoregressiveParam()], which reads it off the prediction form
#' of the process rather than factorizing, and the log-determinant is
#' [autoregressive()]'s negated. Both are therefore \eqn{O(p)} in the entries
#' that matter and exact.
#'
#' # The derivative arrays
#'
#' Written out, from the same prediction form the value comes from:
#'
#' \deqn{\Omega = U^\top \mathrm{diag}(\tau)\, U,}
#'
#' with \eqn{U} unit lower triangular of bandwidth \eqn{q}, row \eqn{t}
#' carrying the coefficients of the best linear predictor of \eqn{y_t} from its
#' predecessors, and \eqn{\tau_t = 1/v_t} the reciprocal innovation variances.
#' Two facts make every order exact without a new recursion.
#'
#' The lower-order rows of \eqn{U} are the coefficients of the SAME family at
#' that order, measured to 0, so their derivative arrays come from the compiled
#' Levinson-Durbin recursion of [ar_taylor()] run once per order, and a
#' component differentiating in a partial autocorrelation an order does not
#' reach is exactly zero. And \eqn{\tau_t} is a PRODUCT of one factor per free
#' value, \eqn{1/v_0} from the scale and \eqn{(1-r_j^2)^{-1}} from each
#' correlation the prediction has reached, so a mixed derivative of it is a
#' product of univariate derivatives and is exactly zero where it
#' differentiates in a factor a row does not carry.
#'
#' What is left is the Leibniz rule over three factors, taken twice so that a
#' component costs \eqn{2^m} matrix products rather than \eqn{3^m}.
#'
#' Measured against `inverse_of(autoregressive(p, q))`, which reaches the same
#' numbers through the ordered-block-partition sum, at order four: **11.2x** at
#' \eqn{p = 6, q = 2}, 11.4x at \eqn{p = 20, q = 3} and **22.8x** at
#' \eqn{p = 100}, with the two agreeing to 7e-11 over six shapes and two free
#' vectors each. That agreement is what licenses the written-out route, the two
#' sharing no arithmetic.
#'
#' @param dimension The side of the matrix.
#' @param order The autoregressive order \eqn{q}, passed to
#'   [autoregressive()].
#' @param ... Further arguments passed to [autoregressive()], such as
#'   `link_scale`.
#'
#' @return An object of class [AutoregressiveInvParam()], with `n_free`,
#'   `free_names` and `dimension` [autoregressive()]'s.
#'
#' @seealso [autoregressive()] for the family this inverts, [ar1_inv()] for the
#'   order-one case, and [inverse_of()] for the general composition both of
#'   them specialize.
#'
#' @examples
#' s <- autoregressive_inv(6, order = 2)
#' eta <- c(log(1.5), atanh(0.5), atanh(-0.2))
#'
#' # Banded of bandwidth two: the process is Markov of order two.
#' lag <- abs(outer(1:6, 1:6, "-"))
#' max(abs(param_value(s, eta)[lag > 2]))
#'
#' # And it is the inverse of the autoregression it names.
#' max(abs(param_value(s, eta) %*% param_value(autoregressive(6, 2), eta) -
#'         diag(6)))
#'
#' @export
autoregressive_inv <- function(dimension, order, ...) {
  inner <- autoregressive(dimension, order, ...)
  p <- inner@dimension
  q <- inner@param_params$order
  # The lower-order families whose coefficients fill U's first q rows. They
  # depend on the dimension, the order and the link and not on the free
  # vector, so they are built once here: measured, rebuilding them at every
  # call was 2.33 ms of a 14.17 ms fourth derivative.
  lower <- lapply(seq_len(q), function(k) {
    if (k == q) inner else autoregressive(p, order = k, ...)
  })
  AutoregressiveInvParam(
    param_name = "autoregressive_inv",
    n_free = inner@n_free,
    free_names = inner@free_names,
    dimension = p,
    rank = p,
    null_basis = empty_null_basis(p),
    param_params = list(inner = inner, lower = lower,
                        cache = new.env(parent = emptyenv()))
  )
}


#' The Sub-Multiset Codes of a Derivative Order
#'
#' @description
#' Enumerates, once per order, every sub-multiset of every index tuple the
#' order carries, together with the integer code each is stored under. A
#' sub-multiset is held as its COUNT VECTOR, one entry per free value, and its
#' code is that vector read as a base-five numeral.
#'
#' @details
#' The strings this replaces were 65 per cent of the cost of a fourth-order
#' component, measured: `sort()` and `paste()` inside the key, not the
#' arithmetic. Counts need neither, `tabulate()` being one C call, and the code
#' indexes a list directly. Five is a safe base because a derivative order here
#' is at most four, so no count can reach it.
#'
#' @param idx The index tuples of the order, as [param_tuple_indices()]
#'   returns them.
#' @param n The number of free values.
#'
#' @return A list with `pow`, the base-five weights; `codes`, one integer
#'   vector per tuple holding the code of `t[take]` for every subset of
#'   positions in the order the bits enumerate them; `comp`, the matching codes
#'   of the complements; and `all`, every distinct code with `cnt`, the count
#'   vector each stands for.
#'
#' @seealso [ar_inv_derivative()], the only caller.
#'
#' @keywords internal
ar_inv_codes <- function(idx, n) {
  pow <- 5^(seq_len(n) - 1L)
  top <- as.integer(sum(4 * pow)) + 1L
  seen <- vector("list", top)
  codes <- comp <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    t <- idx[[i]]
    m <- length(t)
    bits <- 2^(seq_len(m) - 1L)
    cc <- dd <- integer(2^m)
    for (b in seq_len(2^m) - 1L) {
      take <- as.logical(bitwAnd(b, bits))
      a <- tabulate(t[take], nbins = n)
      z <- tabulate(t[!take], nbins = n)
      ca <- as.integer(sum(a * pow)) + 1L
      cz <- as.integer(sum(z * pow)) + 1L
      cc[b + 1L] <- ca
      dd[b + 1L] <- cz
      seen[[ca]] <- a
      seen[[cz]] <- z
    }
    codes[[i]] <- cc
    comp[[i]] <- dd
  }
  seen[[1L]] <- integer(n)
  have <- which(!vapply(seen, is.null, logical(1)))
  list(pow = pow, codes = codes, comp = comp,
       all = have, cnt = seen[have])
}


#' The Prediction Factors of an Inverse Autoregression, With Derivatives
#'
#' @description
#' Returns the two factors of \eqn{\Omega = U^\top \mathrm{diag}(\tau) U},
#' each with its derivative arrays to fourth order, stored under the integer
#' code of the sub-multiset differentiated in.
#'
#' @details
#' \eqn{U} is unit lower triangular of bandwidth \eqn{q}, row \eqn{t} holding
#' the coefficients of the best linear predictor of \eqn{y_t} from its
#' predecessors, which for \eqn{t} beyond the order are the autoregression's
#' own. Those lower-order coefficients are the coefficients of the SAME family
#' at that order, measured exactly: `ar_prediction()`'s row \eqn{k+1} and
#' `autoregressive(p, order = k)`'s `phi` agree to 0. So the intermediate
#' derivative arrays come from the compiled Levinson-Durbin recursion run once
#' per order, and nothing is rederived here. A component differentiating in a
#' partial autocorrelation the order does not reach is exactly zero.
#'
#' \eqn{\tau_t = 1/v_t} is a PRODUCT of one factor per free value:
#' \eqn{1/v_0} from the scale, and \eqn{(1-r_j^2)^{-1}} from each partial
#' autocorrelation the prediction at \eqn{t} has reached. A mixed derivative of
#' a product of univariate factors is the product of their own derivatives, and
#' it is exactly zero whenever it differentiates in a factor that row does not
#' carry.
#'
#' @param s An [AutoregressiveInvParam()] object.
#' @param eta A numeric vector of length `s@n_free`.
#' @param cd The codes of [ar_inv_codes()].
#'
#' @return A list with `u` and `tau`, each a list indexed by code: `u` of
#'   `s@dimension` square matrices, `tau` of numeric vectors of that length.
#'
#' @seealso [autoregressive_inv()] for the formula and [ar_taylor()] for the
#'   recursion the coefficients come from.
#'
#' @keywords internal
ar_inv_factors <- function(s, eta, cd) {
  inner <- .inv_inner(s)
  p <- inner@dimension
  q <- inner@param_params$order
  n <- q + 1L

  # the coefficients of every order, with their derivative arrays, from the
  # lower-order families the constructor built once
  tay <- lapply(seq_len(q), function(k) {
    list(n = k + 1L,
         tay = ar_taylor(s@param_params$lower[[k]], eta[seq_len(k + 1L)]))
  })

  # the univariate factors of tau: 1 / v0, then 1 / (1 - r_j^2)
  fac <- vector("list", n)
  fac[[1L]] <- recip_derivs(link_derivs(inner@param_params$link_scale, eta[1L]))
  for (j in seq_len(q)) {
    fac[[j + 1L]] <- w_derivs(link_derivs(inner@param_params$link_pacf,
                                          eta[j + 1L]))
  }
  # which factors the prediction at row t has reached: the scale always, and
  # the first min(t - 1, q) correlations
  reach <- pmin(seq_len(p) - 1L, q)
  active <- cbind(TRUE, outer(reach, seq_len(q), ">="))

  # the band positions, laid out once
  rows_of <- lapply(seq_len(q), function(k) if (k < q) k + 1L else
    seq.int(q + 1L, p))

  u_list <- tau_list <- vector("list", max(cd$all))
  for (i in seq_along(cd$all)) {
    code <- cd$all[i]
    cnt <- cd$cnt[[i]]
    ord <- sum(cnt)
    tup <- rep(seq_len(n), cnt)

    m <- matrix(0, p, p)
    if (ord == 0L) diag(m) <- 1
    for (k in seq_len(q)) {
      tk <- tay[[k]]
      cf <- if (ord && any(tup > tk$n)) rep(0, k) else
        tk$tay$phi[, ar_pack_col(tk$n, ord, tup)]
      for (t in rows_of[[k]]) m[t, seq.int(t - 1L, t - k)] <- -cf
    }
    u_list[[code]] <- m

    # tau, vectorized over rows: a row contributes zero when it is
    # differentiated in a factor it does not carry
    keep <- !apply(active[, cnt > 0L, drop = FALSE], 1L, function(z) any(!z))
    val <- rep(1, p)
    for (v in seq_len(n)) {
      d <- fac[[v]][[cnt[v] + 1L]]
      val <- val * ifelse(active[, v], d, if (cnt[v] > 0L) 0 else 1)
    }
    tau_list[[code]] <- val * keep
  }

  list(u = u_list, tau = tau_list)
}


#' A Link's Inverse and Its Four Derivatives
#'
#' @description
#' The five-element list every chain in this package consumes: the link's
#' inverse at a free value and its first four derivatives.
#'
#' @param link A \pkg{linkfunctions7} link.
#' @param e A single free value.
#'
#' @return A list of five numbers.
#'
#' @seealso `econ_scalars()`, which builds the same list for two links at once.
#'
#' @keywords internal
link_derivs <- function(link, e) {
  list(
    linkfunctions7::linkinv(link, e),
    linkfunctions7::dlinkinv(link, e),
    linkfunctions7::d2linkinv(link, e),
    linkfunctions7::d3linkinv(link, e),
    linkfunctions7::d4linkinv(link, e)
  )
}


#' Derivative Arrays of an Inverse Autoregressive Parameter
#'
#' @description
#' Assembles one order of \eqn{\partial\Omega} for
#' \eqn{\Omega = U^\top \mathrm{diag}(\tau) U} by the Leibniz rule, taken twice
#' so that the three factors cost \eqn{2^m} products per component rather than
#' \eqn{3^m}.
#'
#' @details
#' Writing \eqn{N = \mathrm{diag}(\tau)U}, the inner Leibniz gives
#' \eqn{\partial_T N = \sum_{R \subseteq T} \mathrm{diag}(\partial_R \tau)
#' \partial_{T \setminus R} U}, and the outer one
#' \eqn{\partial_I \Omega = \sum_{S \subseteq I} (\partial_S U)^\top
#' \partial_{I \setminus S} N}. Both sums run over subsets of the index
#' POSITIONS, which is what makes a repeated index count with its multiplicity.
#' Every \eqn{\partial_S N} is built once and read by every component that
#' needs it.
#'
#' @param s An [AutoregressiveInvParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A named list of `s@dimension` square matrices, keyed and ordered as
#'   [param_tuple_names()] says.
#'
#' @seealso [autoregressive_inv()] for the formula.
#'
#' @keywords internal
ar_inv_derivative <- function(s, eta, order) {
  p <- s@dimension
  st <- ar_inv_structure(s, order)
  cd <- st$cd
  f <- ar_inv_factors(s, eta, cd)

  # N = diag(tau) U and its arrays, by the inner Leibniz, one per distinct code
  nmat <- vector("list", length(f$u))
  for (i in seq_along(cd$all)) {
    acc <- matrix(0, p, p)
    for (b in seq_along(st$sub$codes[[i]])) {
      acc <- acc + f$tau[[st$sub$codes[[i]][b]]] *
        f$u[[st$sub$comp[[i]][b]]]
    }
    nmat[[cd$all[i]]] <- acc
  }

  out <- lapply(seq_along(cd$codes), function(i) {
    acc <- matrix(0, p, p)
    for (b in seq_along(cd$codes[[i]])) {
      acc <- acc + crossprod(f$u[[cd$codes[[i]][b]]], nmat[[cd$comp[[i]][b]]])
    }
    name_dims(acc, s)
  })
  stats::setNames(out, st$names)
}


#' The Subset Structure of One Derivative Order, Memoized
#'
#' @description
#' The index tuples of an order, the codes of every sub-multiset they carry,
#' and the codes of the sub-multisets of THOSE, which the inner Leibniz sums
#' over. None of it depends on the free vector.
#'
#' @details
#' Memoized in an environment held on the object, because it is a function of
#' the order and the number of free values alone and was 59 per cent of a
#' fourth derivative when rebuilt at every call. The cache is pure -- the same
#' order always gives the same structure -- so sharing it across the copies S7
#' makes of the object is safe.
#'
#' @param s An [AutoregressiveInvParam()] object.
#' @param order The derivative order, 1 to 4.
#'
#' @return A list with `cd` and `sub`, both as [ar_inv_codes()] returns, and
#'   `names`, the component names of the order.
#'
#' @seealso [ar_inv_derivative()], the only caller.
#'
#' @keywords internal
ar_inv_structure <- function(s, order) {
  cache <- s@param_params$cache
  key <- as.character(order)
  hit <- cache[[key]]
  if (!is.null(hit)) return(hit)
  n <- s@n_free
  cd <- ar_inv_codes(param_tuple_indices(s, order), n)
  sub <- ar_inv_codes(lapply(cd$cnt, function(cnt) rep(seq_len(n), cnt)), n)
  out <- list(cd = cd, sub = sub, names = param_tuple_names(s, order))
  assign(key, out, envir = cache)
  out
}


#' @title Derivative Arrays of an Inverse Autoregressive Parameter
#' @name param_d1.AutoregressiveInvParam
#' @description
#' The four orders of \eqn{\partial\Omega} for \eqn{\Omega} the precision of an
#' autoregression of order \eqn{q}, written out from the prediction
#' factorization rather than assembled by the ordered-block-partition sum
#' [inverse_of()] uses. The two routes agree to machine precision.
#' @param s An [AutoregressiveInvParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named list of `s@dimension` square matrices, keyed as
#'   [param_tuple_names()] says.
#' @seealso [autoregressive_inv()] for the formulas.
#' @keywords internal
S7::method(param_d1, AutoregressiveInvParam) <- function(s, eta, ...) {
  ar_inv_derivative(s, eta, 1L)
}

#' @rdname param_d1.AutoregressiveInvParam
#' @name param_d2.AutoregressiveInvParam
#' @keywords internal
S7::method(param_d2, AutoregressiveInvParam) <- function(s, eta, ...) {
  ar_inv_derivative(s, eta, 2L)
}

#' @rdname param_d1.AutoregressiveInvParam
#' @name param_d3.AutoregressiveInvParam
#' @keywords internal
S7::method(param_d3, AutoregressiveInvParam) <- function(s, eta, ...) {
  ar_inv_derivative(s, eta, 3L)
}

#' @rdname param_d1.AutoregressiveInvParam
#' @name param_d4.AutoregressiveInvParam
#' @keywords internal
S7::method(param_d4, AutoregressiveInvParam) <- function(s, eta, ...) {
  ar_inv_derivative(s, eta, 4L)
}
