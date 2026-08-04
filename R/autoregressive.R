#' @include jet.R compound_symmetry.R
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
#' parametrised by its marginal variance and its \eqn{q} partial
#' autocorrelations, so \eqn{q + 1} free values whatever the dimension.
#'
#' @details
#' The parametrisation is forced by the shape of the stationary region. The
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
#' derivatives to fourth order are obtained by carrying jets through the
#' recursion rather than by expanding it.
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
#' than by a factorisation.
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

  AutoregressiveParam(
    param_name = sprintf("ar(%d)", q),
    dimension = p,
    n_free = q + 1L,
    free_names = c("scale", paste0("pacf", seq_len(q))),
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(
      order = q,
      link_scale = link_scale,
      link_pacf = linkfunctions7::rhobit_link()
    )
  )
}


#' The Levinson-Durbin Recursion, Carried in Jets
#'
#' @description
#' From the partial autocorrelations as jets, the predictor coefficients of
#' every order and the autocorrelations out to the dimension, each a jet in
#' the free values.
#'
#' @details
#' The recursion is sums and products only, so carrying jets through it gives
#' every derivative to fourth order exactly. The coefficients of the last
#' order are the autoregression's own, and the autocorrelations beyond the
#' order come from the Yule-Walker equations.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param r A list of \eqn{q} jets, the partial autocorrelations.
#' @param lay A layout from \code{\link{jet_layout}}.
#'
#' @return A list with \code{phi}, the coefficient jets of the last order, and
#'   \code{rho}, the autocorrelation jets from lag 0 to \code{dimension - 1}.
#'
#' @keywords internal
ar_levinson <- function(s, r, lay) {
  q <- s@param_params$order
  p <- s@dimension
  one <- jet_const(1, lay)

  phi <- list()
  rho <- list(one)
  for (k in seq_len(q)) {
    new <- vector("list", k)
    new[[k]] <- r[[k]]
    if (k > 1L) {
      for (j in seq_len(k - 1L)) {
        new[[j]] <- jet_add(phi[[j]],
                            jet_cmul(jet_mul(r[[k]], phi[[k - j]], lay), -1))
      }
    }
    phi <- new
    acc <- r[[k]]
    if (k > 1L) {
      for (j in seq_len(k - 1L)) {
        acc <- jet_add(acc, jet_mul(phi[[j]], rho[[k - j + 1L]], lay))
      }
    }
    rho[[k + 1L]] <- acc
  }
  # beyond the order, the Yule-Walker equations carry the recursion on
  if (p - 1L > q) {
    for (h in seq.int(q + 1L, p - 1L)) {
      acc <- jet_const(0, lay)
      for (j in seq_len(q)) {
        acc <- jet_add(acc, jet_mul(phi[[j]], rho[[h - j + 1L]], lay))
      }
      rho[[h + 1L]] <- acc
    }
  }
  list(phi = phi, rho = rho)
}


#' The Jets an Autoregressive Parameter Is Built From
#'
#' @description
#' The scale and the partial autocorrelations as jets, and the recursion's
#' output, computed once for a free vector.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with \code{lay}, \code{scale}, \code{r}, \code{phi} and
#'   \code{rho}.
#'
#' @keywords internal
ar_jets <- function(s, eta) {
  q <- s@param_params$order
  lay <- jet_layout(s@n_free)
  grab <- function(link, e) {
    list(
      linkfunctions7::linkinv(link, e),
      linkfunctions7::dlinkinv(link, e),
      linkfunctions7::d2linkinv(link, e),
      linkfunctions7::d3linkinv(link, e),
      linkfunctions7::d4linkinv(link, e)
    )
  }
  scale <- jet_var(1L, grab(s@param_params$link_scale, eta[1L]), lay)
  r <- lapply(seq_len(q), function(k) {
    jet_var(k + 1L, grab(s@param_params$link_pacf, eta[k + 1L]), lay)
  })
  ld <- ar_levinson(s, r, lay)
  list(lay = lay, scale = scale, r = r, phi = ld$phi, rho = ld$rho)
}


#' The Matrix and Its Derivatives, From the Jets
#'
#' @description
#' Fills the Toeplitz matrix of the scaled autocorrelations, taking either the
#' value or one derivative component from each jet.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param j The jets of \code{\link{ar_jets}}.
#' @param order The derivative order, or 0 for the value.
#' @param position Which component of that order.
#'
#' @return A symmetric numeric matrix.
#'
#' @keywords internal
ar_assemble <- function(s, j, order = 0L, position = 1L) {
  p <- s@dimension
  lag <- abs(outer(seq_len(p), seq_len(p), "-"))
  pick <- function(x) if (order == 0L) x$v else x$d[[order]][position]
  vals <- vapply(seq_len(p), function(h) {
    pick(jet_mul(j$scale, j$rho[[h]], j$lay))
  }, numeric(1))
  matrix(vals[lag + 1L], p, p)
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
  name_dims(ar_assemble(s, ar_jets(s, eta)), s)
}


#' Derivative Components of an Autoregressive Parameter
#'
#' @description
#' Assembles one derivative order by reading the matching component out of
#' every jet.
#'
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of symmetric matrices.
#'
#' @keywords internal
ar_derivative <- function(s, eta, order) {
  j <- ar_jets(s, eta)
  nm <- param_tuple_names(s, order)
  out <- lapply(seq_along(nm), function(i) {
    name_dims(ar_assemble(s, j, order, i), s)
  })
  stats::setNames(out, nm)
}


#' @title Derivatives of an Autoregressive Parameter
#' @name param_d1.AutoregressiveParam
#' @description
#' Closed form at every order. The map from the partial autocorrelations to
#' the matrix is polynomial, so jets carried through the Levinson-Durbin
#' recursion give each derivative exactly; nothing is differenced.
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
  j <- ar_jets(s, eta)
  r <- vapply(j$r, function(x) x$v, numeric(1))
  v0 <- j$scale$v

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
#' \eqn{p\log\gamma_0 + \sum_k (p-k)\log(1 - r_k^2)}. No factorisation and no
#' determinant is computed.
#' @param s An \code{\link{AutoregressiveParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, AutoregressiveParam) <- function(s, eta, ...) {
  p <- s@dimension
  q <- s@param_params$order
  j <- ar_jets(s, eta)
  r <- vapply(j$r, function(x) x$v, numeric(1))
  p * log(j$scale$v) + sum((p - seq_len(q)) * log(1 - r^2))
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
