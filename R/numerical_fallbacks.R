#' @include generics.R
NULL

# What a parameter gets for free when it implements only param_value().
# Every method here is registered on the base class, so a subclass that
# supplies a closed form takes over through dispatch with no registration step.


#' Is This the Package's Own Base Class?
#'
#' @description
#' Asks whether an S7 class is the abstract \code{\link{parameter}} class,
#' which is how a method registered on it is told apart from one a subclass
#' supplied.
#'
#' @details
#' Identity is tried first because it is the usual case and costs nothing, then
#' the name and the package, because identity is not preserved when a package's
#' code is re-evaluated rather than loaded. Coverage tools do exactly that, so
#' an identity-only test passes every ordinary check and fails there.
#'
#' @param cls An S7 class.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
is_base_param_class <- function(cls) {
  for (base in list(parameter, matrix_parameter)) {
    if (identical(cls, base)) return(TRUE)
    if (identical(attr(cls, "name"), attr(base, "name")) &&
        identical(attr(cls, "package"), attr(base, "package"))) {
      return(TRUE)
    }
  }
  FALSE
}


#' Which of a Parameter's Quantities Come From the Base Class
#'
#' @description
#' Reports, for each of the five derivative quantities, whether the parameter
#' supplies its own method or falls back to the one registered on
#' \code{\link{parameter}}.
#'
#' @details
#' The distinction that matters is whether an independent check exists. A
#' derivative computed by finite differences cannot be checked against a finite
#' difference, and a log-determinant read off an eigendecomposition cannot be
#' checked against an eigendecomposition; the comparison is the same arithmetic
#' twice, and it agrees however wrong the parameter is.
#' \code{\link{check_parameter}} uses this to report such a quantity as not
#' checked rather than as passed.
#'
#' \code{param_solve()} and \code{param_factor()} are deliberately absent.
#' Their base-class versions are a Cholesky factorization, which is exact
#' whoever performs it, and the validator compares them with
#' \code{base::solve} either way. Calling them numerical would suggest an
#' approximation that is not there.
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#'
#' @return A named logical vector over \code{param_d1},
#'   \code{param_d2}, \code{param_logdet}, \code{param_dlogdet} and
#'   \code{param_d2logdet}, \code{TRUE} where the base-class method is in
#'   force.
#'
#' @examples
#' # everything closed form except the second-order log-determinant of a
#' # parameter that does not supply one
#' param_is_numerical(log_cholesky(3))
#' param_is_numerical(diagonal_matrix(3))
#'
#' @export
param_is_numerical <- function(s) {
  cls <- S7::S7_class(s)
  gens <- list(
    param_d1 = param_d1,
    param_d2 = param_d2,
    param_d3 = param_d3,
    param_d4 = param_d4,
    param_logdet = param_logdet,
    param_dlogdet = param_dlogdet,
    param_d2logdet = param_d2logdet,
    param_d3logdet = param_d3logdet,
    param_d4logdet = param_d4logdet
  )
  if (!S7::S7_inherits(s, matrix_parameter)) {
    gens <- gens[c("param_d1", "param_d2", "param_d3", "param_d4")]
  }
  vapply(gens, function(g) {
    m <- tryCatch(S7::method(g, cls), error = function(e) NULL)
    if (is.null(m)) return(TRUE)
    is_base_param_class(attr(m, "signature")[[1L]])
  }, logical(1))
}


#' The Spectral Decomposition a Parameter's Quantities Are Read From
#'
#' @description
#' The eigenvalues and eigenvectors of \eqn{M(\eta)}, with the eigenvalues the
#' declared rank says are zero identified by position rather than by size.
#'
#' @details
#' The rank is taken from the object and never re-derived here. Counting
#' eigenvalues above a relative tolerance is not scale invariant, so a
#' parameter whose components differ by many orders of magnitude would be
#' assigned a different rank at different \eqn{\eta} -- and a fitted model with
#' smoothing parameters that far apart is ordinary. The object settled the
#' question once, at construction, from the components.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with \code{values}, \code{vectors}, and \code{keep}, a
#'   logical vector marking the \code{s@rank} directions that carry the matrix.
#'
#' @keywords internal
param_spectrum <- function(s, eta) {
  m <- param_value(s, eta)
  e <- eigen((m + t(m)) / 2, symmetric = TRUE)
  keep <- rep(FALSE, s@dimension)
  if (s@rank > 0L) keep[seq_len(s@rank)] <- TRUE
  list(values = e$values, vectors = e$vectors, keep = keep, matrix = m)
}


#' Moore-Penrose Inverse From a Parameter's Spectrum
#'
#' @description
#' The pseudo-inverse of \eqn{M(\eta)}, formed from the directions the declared
#' rank keeps.
#'
#' @param sp The result of \code{\link{param_spectrum}}.
#'
#' @return A symmetric numeric matrix.
#'
#' @keywords internal
spectrum_pinv <- function(sp) {
  v <- sp$vectors[, sp$keep, drop = FALSE]
  d <- sp$values[sp$keep]
  if (!length(d)) return(matrix(0, nrow(sp$vectors), nrow(sp$vectors)))
  v %*% (t(v) / d)
}


#' Cholesky Factorization, With the Rank Decided Before It
#'
#' @description
#' The lower triangular Cholesky factor of a symmetric matrix, or \code{NULL}
#' when the matrix is not positive definite to the given relative tolerance.
#'
#' @details
#' The verdict comes from the eigenvalues rather than from whether
#' \code{\link[base]{chol}} raises. On a matrix with an exactly zero eigenvalue
#' the pivot that should be zero comes out positive or negative according to
#' rounding, so \code{chol()} succeeds on some platforms and fails on others,
#' and a construction that asks it whether a matrix is usable gets a different
#' answer on different machines.
#'
#' @param m A symmetric numeric matrix.
#' @param tol The relative tolerance below which the smallest eigenvalue counts
#'   as zero.
#'
#' @return The lower triangular factor \eqn{L} with \eqn{M = L L^\top}, or
#'   \code{NULL}.
#'
#' @keywords internal
chol_pd <- function(m, tol = 1e-12) {
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  if (!length(ev) || anyNA(ev)) return(NULL)
  if (max(ev) <= 0 || min(ev) <= tol * max(ev)) return(NULL)
  r <- tryCatch(chol(m), error = function(e) NULL)
  if (is.null(r)) NULL else t(r)
}


#' Finite-Difference Step for a Free Value
#'
#' @description
#' The step for a central difference in one component of \eqn{\eta}, scaled by
#' the size of that component.
#'
#' @details
#' The free scale is unbounded, so unlike the response and parameter steps of
#' \pkg{distributions7} this one has no boundary to be clamped away from: the
#' whole point of the unconstrained scale is that there is nowhere to fall off.
#'
#' @param eta_k The value of the component.
#' @param order The derivative order the step is for.
#'
#' @return A single positive number.
#'
#' @keywords internal
fd_step <- function(eta_k, order = 1L) {
  .Machine$double.eps^(1 / (order + 2)) * max(1, abs(eta_k))
}


#' Numerical First Derivatives of a Parameter's Matrix
#'
#' @description
#' One central difference of \code{\link{param_value}} in each component of
#' \eqn{\eta}.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A named list of symmetric matrices.
#'
#' @seealso \code{\link{param_d1}}
#' @examples
#' numerical_d1(scalar_matrix(2), 0.3)
#' @export
numerical_d1 <- function(s, eta) {
  out <- vector("list", s@n_free)
  names(out) <- s@free_names
  for (k in seq_len(s@n_free)) {
    h <- fd_step(eta[k], 1L)
    up <- eta
    dn <- eta
    up[k] <- eta[k] + h
    dn[k] <- eta[k] - h
    out[[k]] <- (param_value(s, up) - param_value(s, dn)) / (2 * h)
  }
  out
}


#' Numerical Second Derivatives of a Parameter's Matrix
#'
#' @description
#' One central difference of the analytic first derivatives where a parameter
#' supplies them, and a single mixed stencil on \code{\link{param_value}}
#' where it does not.
#'
#' @details
#' Differentiating the analytic first derivative costs one finite-difference
#' layer rather than two, which is the rule the whole toolkit follows: never
#' compose differences in the same variable. When the first derivatives are
#' themselves numerical the two differences act on different components for an
#' off-diagonal pair, so they commute into a four-point mixed stencil rather
#' than compounding; the diagonal pairs use the three-point second-difference
#' stencil directly, which is again one layer and not two.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A named list of symmetric matrices, keyed as
#'   \code{param_tuple_names(s)}.
#'
#' @seealso \code{\link{param_d2}}
#' @examples
#' numerical_d2(scalar_matrix(2), 0.3)
#' @export
numerical_d2 <- function(s, eta) {
  idx <- param_tuple_indices(s)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s)
  if (!length(idx)) return(out)

  analytic <- !param_is_numerical(s)[["param_d1"]]

  for (i in seq_along(idx)) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    if (analytic) {
      # One layer on the analytic first derivative in the other component.
      h <- fd_step(eta[l], 1L)
      up <- eta
      dn <- eta
      up[l] <- eta[l] + h
      dn[l] <- eta[l] - h
      out[[i]] <- (param_d1(s, up)[[k]] - param_d1(s, dn)[[k]]) /
        (2 * h)
    } else if (k == l) {
      h <- fd_step(eta[k], 2L)
      up <- eta
      dn <- eta
      up[k] <- eta[k] + h
      dn[k] <- eta[k] - h
      out[[i]] <- (param_value(s, up) - 2 * param_value(s, eta) +
        param_value(s, dn)) / h^2
    } else {
      hk <- fd_step(eta[k], 2L)
      hl <- fd_step(eta[l], 2L)
      pp <- pm <- mp <- mm <- eta
      pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
      pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
      mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
      mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
      out[[i]] <- (param_value(s, pp) - param_value(s, pm) -
        param_value(s, mp) + param_value(s, mm)) / (4 * hk * hl)
    }
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


#' @title Default First Derivatives
#' @name param_d1.parameter
#' @description Fallback: one central difference of
#'   \code{\link{param_value}} per component (see
#'   \code{\link{numerical_d1}}).
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, parameter) <- function(s, eta, ...) {
  numerical_d1(s, eta)
}

#' @title Default Second Derivatives
#' @name param_d2.parameter
#' @description Fallback: see \code{\link{numerical_d2}}.
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d2, parameter) <- function(s, eta, ...) {
  numerical_d2(s, eta)
}

#' @title Default Log-Determinant
#' @name param_logdet.parameter
#' @description
#' Fallback: the sum of the logs of the eigenvalues the declared rank keeps,
#' which is the log-determinant for a full-rank family and the log
#' pseudo-determinant otherwise.
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, matrix_parameter) <- function(s, eta, ...) {
  sp <- param_spectrum(s, eta)
  v <- sp$values[sp$keep]
  if (any(v <= 0)) {
    stop(sprintf(paste0(
      "'%s' produced %d non-positive eigenvalue(s) among the %d its rank\n",
      "  keeps, so it is outside the set it claims to parametrize."
    ), s@param_name, sum(v <= 0), s@rank), call. = FALSE)
  }
  sum(log(v))
}

#' @title Default Log-Determinant Gradient
#' @name param_dlogdet.parameter
#' @description
#' Fallback: \eqn{\mathrm{tr}(M^{+} \partial_k M)}, with the pseudo-inverse
#' formed from the directions the declared rank keeps, which is the ordinary
#' inverse when the family is of full rank.
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, matrix_parameter) <- function(s, eta, ...) {
  sp <- param_spectrum(s, eta)
  mi <- spectrum_pinv(sp)
  d <- param_d1(s, eta)
  stats::setNames(
    vapply(d, function(dk) sum(mi * dk), numeric(1)),
    s@free_names
  )
}

#' @title Default Log-Determinant Hessian
#' @name param_d2logdet.parameter
#' @description
#' Fallback: \eqn{\mathrm{tr}(M^{+} \partial_{kl} M) -
#' \mathrm{tr}(M^{+} \partial_k M\, M^{+} \partial_l M)}, the derivative of the
#' identity behind \code{\link{param_dlogdet}}.
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d2logdet, matrix_parameter) <- function(s, eta, ...) {
  sp <- param_spectrum(s, eta)
  mi <- spectrum_pinv(sp)
  d <- param_d1(s, eta)
  d2 <- param_d2(s, eta)
  idx <- param_tuple_indices(s)
  out <- vapply(seq_along(idx), function(i) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    sum(mi * d2[[i]]) - sum(t(mi %*% d[[k]]) * (mi %*% d[[l]]))
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s))
}

#' @title Default Solve
#' @name param_solve.parameter
#' @description Fallback: a Cholesky of \code{\link{param_value}}, with the
#'   definiteness verdict taken spectrally (see \code{\link{chol_pd}}).
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param b A numeric matrix with \code{s@dimension} rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(param_solve, matrix_parameter) <- function(s, eta, b = NULL, ...) {
  l <- param_factor(s, eta)
  backsolve(t(l), forwardsolve(l, b))
}

#' @title Default Factor
#' @name param_factor.parameter
#' @description Fallback: the lower Cholesky factor of
#'   \code{\link{param_value}}, rejected when the matrix is not positive
#'   definite spectrally.
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A lower triangular numeric matrix.
#' @keywords internal
S7::method(param_factor, matrix_parameter) <- function(s, eta, ...) {
  l <- chol_pd(param_value(s, eta))
  if (is.null(l)) {
    stop(sprintf(paste0(
      "'%s' is not positive definite at this eta, although it declares full\n",
      "  rank. The verdict is spectral, so this is a statement about the\n",
      "  matrix rather than about a factorization that happened to fail."
    ), s@param_name), call. = FALSE)
  }
  l
}

#' @title Rejection to Invert Without a Closed Form
#' @name param_free.parameter
#' @description
#' The base class rejects rather than inverting the map numerically: an
#' optimization-based inverse would return a plausible \eqn{\eta} for a matrix
#' outside the set the family parametrizes.
#' @param s A \code{\link{parameter}} object.
#' @param m A symmetric numeric matrix.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @keywords internal
S7::method(param_free, parameter) <- function(s, m, ...) {
  stop(sprintf(paste0(
    "'%s' does not implement param_free(). The inverse map is exact or\n",
    "  rejected, never obtained by optimization, because a numerical inverse\n",
    "  would return a plausible eta for a matrix outside the set."
  ), s@param_name), call. = FALSE)
}


#' Numerical Third Derivatives of a Parameter's Value
#'
#' @description
#' One product stencil applied directly to \code{\link{param_value}} for each
#' distinct index tuple: a single stencil on the map itself, never a stencil
#' on a lower-order numerical derivative, so the no-nesting rule holds
#' whatever the family implements. A repeated component uses the matching
#' higher-order one-dimensional factor; distinct components each contribute a
#' central two-point factor, and the product is evaluated in one pass.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A named list keyed as \code{\link{param_tuple_names}(s, 3)}.
#'
#' @seealso \code{\link{param_d3}}
#' @examples
#' numerical_d3(scalar_matrix(2), 0.3)
#' @export
numerical_d3 <- function(s, eta) {
  idx <- param_tuple_indices(s, 3L)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, 3L)
  f <- function(e) param_value(s, e)
  for (i in seq_along(idx)) {
    out[[i]] <- mixed_stencil(f, eta, idx[[i]])
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


#' Numerical Fourth Derivatives of a Parameter's Value
#'
#' @description
#' The order-four analogue of \code{\link{numerical_d3}}: one stencil per
#' component, applied directly to the map. Rounding is amplified by the
#' fourth power of the step, so this is accurate to roughly four significant
#' digits -- a starting point, and the reason every shipped family carries
#' closed forms instead.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A named list keyed as \code{\link{param_tuple_names}(s, 4)}.
#'
#' @seealso \code{\link{param_d4}}
#' @examples
#' numerical_d4(scalar_matrix(2), 0.3)
#' @export
numerical_d4 <- function(s, eta) {
  idx <- param_tuple_indices(s, 4L)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, 4L)
  f <- function(e) param_value(s, e)
  for (i in seq_along(idx)) {
    out[[i]] <- mixed_stencil(f, eta, idx[[i]])
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


#' One Product Stencil for a Mixed Partial Derivative
#'
#' @description
#' Differentiates \code{f} once in each component the index tuple names, a
#' central factor per distinct component of the order its multiplicity asks:
#' two points at order one, three at two, the five-point forms at three and
#' four. The result is a single product stencil, not a composition of
#' lower-order numerical derivatives.
#'
#' @param f A function of the free vector.
#' @param eta The point.
#' @param tuple An integer vector of component indices, possibly repeated.
#'
#' @return The stencil's value, shaped like \code{f(eta)}.
#'
#' @keywords internal
mixed_stencil <- function(f, eta, tuple) {
  counts <- table(tuple)
  ks <- as.integer(names(counts))
  ms <- as.integer(counts)
  hs <- vapply(seq_along(ks), function(j) {
    fd_step(eta[ks[j]], sum(ms))
  }, numeric(1))

  one_dim <- function(m) {
    switch(m,
      list(offsets = c(-1, 1), weights = c(-0.5, 0.5), power = 1),
      list(offsets = c(-1, 0, 1), weights = c(1, -2, 1), power = 2),
      list(offsets = c(-2, -1, 1, 2), weights = c(-0.5, 1, -1, 0.5), power = 3),
      list(offsets = c(-2, -1, 0, 1, 2), weights = c(1, -4, 6, -4, 1), power = 4)
    )
  }
  facs <- lapply(ms, one_dim)

  combos <- expand.grid(lapply(facs, function(fc) seq_along(fc$offsets)))
  acc <- NULL
  for (r in seq_len(nrow(combos))) {
    e <- eta
    w <- 1
    for (j in seq_along(ks)) {
      pick <- combos[r, j]
      e[ks[j]] <- e[ks[j]] + facs[[j]]$offsets[pick] * hs[j]
      w <- w * facs[[j]]$weights[pick]
    }
    term <- w * f(e)
    acc <- if (is.null(acc)) term else acc + term
  }
  denom <- prod(vapply(seq_along(ks), function(j) {
    hs[j]^facs[[j]]$power
  }, numeric(1)))
  acc / denom
}


#' @title Default Third Derivatives
#' @name param_d3.parameter
#' @description Fallback: see \code{\link{numerical_d3}}.
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list.
#' @keywords internal
S7::method(param_d3, parameter) <- function(s, eta, ...) {
  numerical_d3(s, eta)
}

#' @title Default Fourth Derivatives
#' @name param_d4.parameter
#' @description Fallback: see \code{\link{numerical_d4}}.
#' @param s A \code{\link{parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list.
#' @keywords internal
S7::method(param_d4, parameter) <- function(s, eta, ...) {
  numerical_d4(s, eta)
}

#' @title Default Higher Log-Determinant Derivatives
#' @name param_d3logdet.matrix_parameter
#' @description
#' Fallback: one central stencil on \code{\link{param_d2logdet}}, which is
#' the exact trace identity given the matrix derivatives -- a single layer on
#' an analytic quantity, per the toolkit's rule.
#' @param s A \code{\link{matrix_parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d3logdet, matrix_parameter) <- function(s, eta, ...) {
  idx3 <- param_tuple_indices(s, 3L)
  idx2 <- param_tuple_indices(s, 2L)
  key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
  out <- vapply(seq_along(idx3), function(i) {
    t3 <- idx3[[i]]
    k <- t3[3L]
    pos <- match(paste(sort(t3[1:2]), collapse = ","), key2)
    h <- fd_step(eta[k], 1L)
    up <- eta
    dn <- eta
    up[k] <- eta[k] + h
    dn[k] <- eta[k] - h
    (param_d2logdet(s, up)[[pos]] - param_d2logdet(s, dn)[[pos]]) / (2 * h)
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, 3L))
}

#' @title Default Fourth Log-Determinant Derivatives
#' @name param_d4logdet.matrix_parameter
#' @description
#' Fallback: one second-order stencil on \code{\link{param_d2logdet}} in the
#' last two components of the tuple -- one layer, mixed across components
#' where they differ.
#' @param s A \code{\link{matrix_parameter}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d4logdet, matrix_parameter) <- function(s, eta, ...) {
  idx4 <- param_tuple_indices(s, 4L)
  idx2 <- param_tuple_indices(s, 2L)
  key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
  out <- vapply(seq_along(idx4), function(i) {
    t4 <- idx4[[i]]
    pos <- match(paste(sort(t4[1:2]), collapse = ","), key2)
    k <- t4[3L]
    l <- t4[4L]
    g <- function(e) param_d2logdet(s, e)[[pos]]
    if (k == l) {
      h <- fd_step(eta[k], 2L)
      up <- eta
      dn <- eta
      up[k] <- eta[k] + h
      dn[k] <- eta[k] - h
      (g(up) - 2 * g(eta) + g(dn)) / h^2
    } else {
      hk <- fd_step(eta[k], 2L)
      hl <- fd_step(eta[l], 2L)
      pp <- pm <- mp <- mm <- eta
      pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
      pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
      mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
      mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
      (g(pp) - g(pm) - g(mp) + g(mm)) / (4 * hk * hl)
    }
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, 4L))
}
