#' @include compound_symmetry.R
NULL


#' Autoregressive Parameter
#'
#' @description
#' The S7 class of covariance matrices of a stationary autoregression of any
#' order. Constructed by \code{\link{autoregressive}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{AutoregressiveParam}.
#'
#' @seealso \code{\link{autoregressive}}
#'
#' @examples
#' S7::S7_inherits(autoregressive(5, order = 2), AutoregressiveParam)
#'
#' @export
AutoregressiveParam <- S7::new_class("AutoregressiveParam",
  parent = matrix_parameter)


#' Construct an Autoregressive Parameter
#'
#' @description
#' The covariance of \eqn{p} consecutive observations of a stationary
#' autoregression of order \eqn{q},
#' \deqn{y_t = \phi_1 y_{t-1} + \cdots + \phi_q y_{t-q} + \varepsilon_t,}
#' parametrized by its marginal variance and its \eqn{q} partial
#' autocorrelations, so \eqn{q + 1} free values whatever the dimension.
#'
#' @details
#' The parametrization is forced by the shape of the stationary region. The
#' coefficients \eqn{\phi} are stationary exactly when the roots of
#' \eqn{1 - \phi_1 z - \cdots - \phi_q z^{q}} lie outside the unit circle, and
#' that set is not a box: already at \eqn{q = 2} it is a triangle, so no
#' collection of scalar links onto intervals can cover it. The partial
#' autocorrelations do not have this problem. Each lies in \eqn{(-1, 1)}
#' independently of the others, and the Levinson-Durbin recursion carries them
#' onto the stationary coefficients bijectively, which is the transformation
#' of Barndorff-Nielsen and Schou (1973) and Monahan (1984). Each partial
#' autocorrelation therefore takes an \code{\link[linkfunctions7]{rhobit_link}}
#' and every free vector gives a stationary, positive definite matrix.
#'
#' The autocorrelations follow from the same recursion. Writing
#' \eqn{\phi^{(k)}} for the coefficients of the order-\eqn{k} predictor,
#' \deqn{\phi^{(k)}_k = r_k, \qquad
#'   \phi^{(k)}_j = \phi^{(k-1)}_j - r_k \phi^{(k-1)}_{k-j},}
#' and the Yule-Walker equations at order \eqn{k} give
#' \eqn{\rho_k = \sum_{j<k} \phi^{(k)}_j \rho_{k-j} + r_k}, after which
#' \eqn{\rho_h = \sum_j \phi_j \rho_{h-j}} for every lag beyond the order. The
#' whole map from the partial autocorrelations to the matrix is therefore
#' \strong{polynomial}, built from sums and products alone, and its
#' derivatives to fourth order are obtained by propagating the derivative
#' arrays through the recursion in compiled code, the product rule written
#' out per order, rather than by expanding it.
#'
#' Two quantities are closed form. The innovation variances of the
#' Levinson-Durbin recursion give
#' \deqn{\log\lvert M \rvert = p\log\gamma_0
#'   + \sum_{k=1}^{q} (p - k)\log(1 - r_k^{2}),}
#' one term per free value, so the log-determinant is separable and every
#' mixed derivative of it is exactly zero. And the inverse is
#' \strong{banded of bandwidth \eqn{q}}: an autoregression of order \eqn{q} is
#' Markov of that order, so its precision carries no entry beyond the
#' \eqn{q}-th off-diagonal. It is assembled from the prediction form
#' \eqn{M^{-1} = U^\top D^{-1} U}, with \eqn{U} unit lower triangular holding
#' the predictor coefficients and \eqn{D} the innovation variances, rather
#' than by a factorization.
#'
#' \code{\link{ar1}} is the case \eqn{q = 1} written out: there the
#' autocorrelation is simply \eqn{\rho^{h}}, the determinant is
#' \eqn{(1-\rho^2)^{p-1}} and the inverse is tridiagonal in three lines, so it
#' keeps its own closed forms and does not go through the recursion.
#'
#' The name is \code{autoregressive()} rather than \code{ar()} because
#' \code{\link[stats]{ar}} is a function of \pkg{stats}, and a package meant
#' to be attached alongside others should not mask one.
#'
#' @param dimension The side \eqn{p} of the matrix: the number of consecutive
#'   observations. Must exceed the order, since a shorter stretch does not
#'   identify the last partial autocorrelation.
#' @param order The order \eqn{q} of the autoregression, at least 1.
#' @param link_scale A \pkg{linkfunctions7} link onto the positive scale,
#'   carrying the marginal variance. Defaults to
#'   \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{AutoregressiveParam}}.
#'
#' @references
#' Barndorff-Nielsen, O. and Schou, G. (1973). On the parametrization of
#' autoregressive models by partial autocorrelations. \emph{Journal of
#' Multivariate Analysis} 3, 408-419.
#'
#' Monahan, J. F. (1984). A note on enforcing stationarity in autoregressive
#' moving average models. \emph{Biometrika} 71, 403-404.
#'
#' @seealso \code{\link{ar1}}, \code{\link{compound_symmetry}}
#'
#' @examples
#' s <- autoregressive(6, order = 2)
#' s@free_names
#'
#' eta <- c(log(2), atanh(0.7), atanh(-0.3))
#' round(param_value(s, eta), 4)
#'
#' # the precision is banded of bandwidth two, the process being Markov of
#' # order two
#' round(param_solve(s, eta), 4)
#'
#' # the round trip closes exactly
#' max(abs(param_free(s, param_value(s, eta)) - eta))
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
#' Runs the compiled recursion of \code{ar_taylor_cpp}: the scale and
#' the partial autocorrelations enter as their link inverses with four
#' derivatives each, and the autocorrelations, the coefficients and every
#' partial derivative to fourth order come out as packed arrays.
#'
#' @details
#' The recursion is sums and products only, so the propagation rules are the
#' product rule written out per order; every derivative is exact and nothing
#' is differenced.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with \code{n}, the number of free values; \code{gamma}, a
#'   matrix with one row per lag; and \code{phi}, one row per coefficient.
#'   Each row packs the value, then the full derivative tensors of orders one
#'   to four, in row-major order.
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
#' Fills the Toeplitz matrix of the scaled autocorrelations, taking either the
#' value column or one derivative component out of the packed rows of
#' \code{\link{ar_taylor}}.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param tay The arrays of \code{\link{ar_taylor}}.
#' @param order The derivative order, or 0 for the value.
#' @param tuple The index tuple of that order, ignored at order 0.
#'
#' @return A symmetric numeric matrix.
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
#' Where a derivative component sits in a row of \code{\link{ar_taylor}}'s
#' output: the value first, then the tensors of orders one to four in
#' row-major order.
#'
#' @param n The number of free values.
#' @param order The derivative order, or 0 for the value.
#' @param tuple The index tuple, 1-based.
#'
#' @return A single column index.
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
#' The Toeplitz matrix \eqn{\gamma_0 \rho_{\lvert i-j \rvert}}, the
#' autocorrelations coming from the Levinson-Durbin recursion.
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A positive definite Toeplitz matrix.
#' @keywords internal
S7::method(param_value, AutoregressiveParam) <- function(s, eta, ...) {
  name_dims(ar_assemble(s, ar_taylor(s, eta)), s)
}


#' Derivative Components of an Autoregressive Parameter
#'
#' @description
#' Assembles one derivative order by reading the matching component out of
#' the packed arrays.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of symmetric matrices.
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
#' Closed form at every order. The map from the partial autocorrelations to
#' the matrix is polynomial, so the derivative arrays propagated through the
#' Levinson-Durbin recursion give each derivative exactly; nothing is
#' differenced.
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
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
#' The marginal variance is the common diagonal entry and the partial
#' autocorrelations come from the Levinson-Durbin recursion run forwards on
#' the autocorrelations. The matrix is then checked against the pattern those
#' values imply, so a Toeplitz matrix that is not the covariance of an
#' autoregression of this order is refused rather than fitted.
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param m A covariance matrix of a stationary autoregression.
#' @param ... Unused.
#' @return A named numeric vector of free values.
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
      "  process. It is refused rather than averaged along its diagonals."
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
#' The unit lower triangular matrix of one-step predictor coefficients and the
#' innovation variances, which factor the matrix as
#' \eqn{M = U^{-1} D U^{-\top}}.
#'
#' @details
#' Row \eqn{t} holds the coefficients of the best linear predictor of
#' \eqn{y_t} from its predecessors, which for \eqn{t} beyond the order are the
#' autoregression's own coefficients, so \eqn{U} has bandwidth \eqn{q}. The
#' innovation variances fall by a factor \eqn{1 - r_k^2} at each of the first
#' \eqn{q} steps and are constant thereafter.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with the matrix \code{u} and the vector \code{v}.
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
#' Exact and banded of bandwidth \eqn{q}. The precision is
#' \eqn{U^\top D^{-1} U} in the prediction form of
#' \code{\link{ar_prediction}}, which is what an order-\eqn{q} Markov property
#' means: no partial correlation beyond the lag.
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param b A numeric matrix with \code{s@dimension} rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(param_solve, AutoregressiveParam) <- function(s, eta, b = NULL, ...) {
  pf <- ar_prediction(s, eta)
  (crossprod(pf$u, pf$u / pf$v)) %*% b
}


#' @title Log-Determinant of an Autoregressive Parameter
#' @name param_logdet.AutoregressiveParam
#' @description
#' Closed form from the innovation variances:
#' \eqn{p\log\gamma_0 + \sum_k (p-k)\log(1 - r_k^2)}. No factorization and no
#' determinant is computed.
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
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
#' term per free value, so every mixed component is exactly zero.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named numeric vector.
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
#' Closed form at every order, with every mixed component exactly zero: the
#' log-determinant is a sum with one term per free value.
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
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
