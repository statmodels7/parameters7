#' @include generics.R
NULL

# What a structure gets for free when it implements only struct_matrix().
# Every method here is registered on the base class, so a subclass that
# supplies a closed form takes over through dispatch with no registration step.


#' Is This the Package's Own Base Class?
#'
#' @description
#' Asks whether an S7 class is the abstract \code{\link{covstruct}} class,
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
is_base_struct_class <- function(cls) {
  if (identical(cls, covstruct)) return(TRUE)
  identical(attr(cls, "name"), attr(covstruct, "name")) &&
    identical(attr(cls, "package"), attr(covstruct, "package"))
}


#' Which of a Structure's Quantities Come From the Base Class
#'
#' @description
#' Reports, for each of the five derivative quantities, whether the structure
#' supplies its own method or falls back to the one registered on
#' \code{\link{covstruct}}.
#'
#' @details
#' The distinction that matters is whether an independent check exists. A
#' derivative computed by finite differences cannot be checked against a finite
#' difference, and a log-determinant read off an eigendecomposition cannot be
#' checked against an eigendecomposition; the comparison is the same arithmetic
#' twice, and it agrees however wrong the structure is.
#' \code{\link{check_covstruct}} uses this to report such a quantity as not
#' checked rather than as passed.
#'
#' \code{struct_solve()} and \code{struct_factor()} are deliberately absent.
#' Their base-class versions are a Cholesky factorisation, which is exact
#' whoever performs it, and the validator compares them with
#' \code{base::solve} either way. Calling them numerical would suggest an
#' approximation that is not there.
#'
#' @param s An object inheriting from class \code{\link{covstruct}}.
#'
#' @return A named logical vector over \code{struct_dmatrix},
#'   \code{struct_d2matrix}, \code{struct_logdet}, \code{struct_dlogdet} and
#'   \code{struct_d2logdet}, \code{TRUE} where the base-class method is in
#'   force.
#'
#' @examples
#' # everything closed form except the second-order log-determinant of a
#' # structure that does not supply one
#' struct_is_numerical(log_cholesky(3))
#' struct_is_numerical(diag_struct(3))
#'
#' @export
struct_is_numerical <- function(s) {
  cls <- S7::S7_class(s)
  gens <- list(
    struct_dmatrix = struct_dmatrix,
    struct_d2matrix = struct_d2matrix,
    struct_logdet = struct_logdet,
    struct_dlogdet = struct_dlogdet,
    struct_d2logdet = struct_d2logdet
  )
  vapply(gens, function(g) {
    m <- tryCatch(S7::method(g, cls), error = function(e) NULL)
    if (is.null(m)) return(TRUE)
    is_base_struct_class(attr(m, "signature")[[1L]])
  }, logical(1))
}


#' The Spectral Decomposition a Structure's Quantities Are Read From
#'
#' @description
#' The eigenvalues and eigenvectors of \eqn{M(\eta)}, with the eigenvalues the
#' declared rank says are zero identified by position rather than by size.
#'
#' @details
#' The rank is taken from the object and never re-derived here. Counting
#' eigenvalues above a relative tolerance is not scale invariant, so a
#' structure whose components differ by many orders of magnitude would be
#' assigned a different rank at different \eqn{\eta} -- and a fitted model with
#' smoothing parameters that far apart is ordinary. The object settled the
#' question once, at construction, from the components.
#'
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with \code{values}, \code{vectors}, and \code{keep}, a
#'   logical vector marking the \code{s@rank} directions that carry the matrix.
#'
#' @keywords internal
struct_spectrum <- function(s, eta) {
  m <- struct_matrix(s, eta)
  e <- eigen((m + t(m)) / 2, symmetric = TRUE)
  keep <- rep(FALSE, s@dimension)
  if (s@rank > 0L) keep[seq_len(s@rank)] <- TRUE
  list(values = e$values, vectors = e$vectors, keep = keep, matrix = m)
}


#' Moore-Penrose Inverse From a Structure's Spectrum
#'
#' @description
#' The pseudo-inverse of \eqn{M(\eta)}, formed from the directions the declared
#' rank keeps.
#'
#' @param sp The result of \code{\link{struct_spectrum}}.
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


#' Cholesky Factorisation, With the Rank Decided Before It
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


#' Numerical First Derivatives of a Structure's Matrix
#'
#' @description
#' One central difference of \code{\link{struct_matrix}} in each component of
#' \eqn{\eta}.
#'
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A named list of symmetric matrices.
#'
#' @seealso \code{\link{struct_dmatrix}}
#' @export
numerical_dmatrix <- function(s, eta) {
  out <- vector("list", s@n_free)
  names(out) <- s@free_names
  for (k in seq_len(s@n_free)) {
    h <- fd_step(eta[k], 1L)
    up <- eta
    dn <- eta
    up[k] <- eta[k] + h
    dn[k] <- eta[k] - h
    out[[k]] <- (struct_matrix(s, up) - struct_matrix(s, dn)) / (2 * h)
  }
  out
}


#' Numerical Second Derivatives of a Structure's Matrix
#'
#' @description
#' One central difference of the analytic first derivatives where a structure
#' supplies them, and a single mixed stencil on \code{\link{struct_matrix}}
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
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A named list of symmetric matrices, keyed as
#'   \code{struct_pair_names(s)}.
#'
#' @seealso \code{\link{struct_d2matrix}}
#' @export
numerical_d2matrix <- function(s, eta) {
  idx <- struct_pair_indices(s)
  out <- vector("list", length(idx))
  names(out) <- struct_pair_names(s)
  if (!length(idx)) return(out)

  analytic <- !struct_is_numerical(s)[["struct_dmatrix"]]

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
      out[[i]] <- (struct_dmatrix(s, up)[[k]] - struct_dmatrix(s, dn)[[k]]) /
        (2 * h)
    } else if (k == l) {
      h <- fd_step(eta[k], 2L)
      up <- eta
      dn <- eta
      up[k] <- eta[k] + h
      dn[k] <- eta[k] - h
      out[[i]] <- (struct_matrix(s, up) - 2 * struct_matrix(s, eta) +
        struct_matrix(s, dn)) / h^2
    } else {
      hk <- fd_step(eta[k], 2L)
      hl <- fd_step(eta[l], 2L)
      pp <- pm <- mp <- mm <- eta
      pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
      pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
      mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
      mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
      out[[i]] <- (struct_matrix(s, pp) - struct_matrix(s, pm) -
        struct_matrix(s, mp) + struct_matrix(s, mm)) / (4 * hk * hl)
    }
    out[[i]] <- (out[[i]] + t(out[[i]])) / 2
  }
  out
}


#' @title Default First Derivatives
#' @name struct_dmatrix.covstruct
#' @description Fallback: one central difference of
#'   \code{\link{struct_matrix}} per component (see
#'   \code{\link{numerical_dmatrix}}).
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(struct_dmatrix, covstruct) <- function(s, eta, ...) {
  numerical_dmatrix(s, eta)
}

#' @title Default Second Derivatives
#' @name struct_d2matrix.covstruct
#' @description Fallback: see \code{\link{numerical_d2matrix}}.
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(struct_d2matrix, covstruct) <- function(s, eta, ...) {
  numerical_d2matrix(s, eta)
}

#' @title Default Log-Determinant
#' @name struct_logdet.covstruct
#' @description
#' Fallback: the sum of the logs of the eigenvalues the declared rank keeps,
#' which is the log-determinant for a full-rank family and the log
#' pseudo-determinant otherwise.
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(struct_logdet, covstruct) <- function(s, eta, ...) {
  sp <- struct_spectrum(s, eta)
  v <- sp$values[sp$keep]
  if (any(v <= 0)) {
    stop(sprintf(paste0(
      "'%s' produced %d non-positive eigenvalue(s) among the %d its rank\n",
      "  keeps, so it is outside the set it claims to parametrise."
    ), s@struct_name, sum(v <= 0), s@rank), call. = FALSE)
  }
  sum(log(v))
}

#' @title Default Log-Determinant Gradient
#' @name struct_dlogdet.covstruct
#' @description
#' Fallback: \eqn{\mathrm{tr}(M^{+} \partial_k M)}, with the pseudo-inverse
#' formed from the directions the declared rank keeps, which is the ordinary
#' inverse when the family is of full rank.
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(struct_dlogdet, covstruct) <- function(s, eta, ...) {
  sp <- struct_spectrum(s, eta)
  mi <- spectrum_pinv(sp)
  d <- struct_dmatrix(s, eta)
  stats::setNames(
    vapply(d, function(dk) sum(mi * dk), numeric(1)),
    s@free_names
  )
}

#' @title Default Log-Determinant Hessian
#' @name struct_d2logdet.covstruct
#' @description
#' Fallback: \eqn{\mathrm{tr}(M^{+} \partial_{kl} M) -
#' \mathrm{tr}(M^{+} \partial_k M\, M^{+} \partial_l M)}, the derivative of the
#' identity behind \code{\link{struct_dlogdet}}.
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(struct_d2logdet, covstruct) <- function(s, eta, ...) {
  sp <- struct_spectrum(s, eta)
  mi <- spectrum_pinv(sp)
  d <- struct_dmatrix(s, eta)
  d2 <- struct_d2matrix(s, eta)
  idx <- struct_pair_indices(s)
  out <- vapply(seq_along(idx), function(i) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    sum(mi * d2[[i]]) - sum(t(mi %*% d[[k]]) * (mi %*% d[[l]]))
  }, numeric(1))
  stats::setNames(out, struct_pair_names(s))
}

#' @title Default Solve
#' @name struct_solve.covstruct
#' @description Fallback: a Cholesky of \code{\link{struct_matrix}}, with the
#'   definiteness verdict taken spectrally (see \code{\link{chol_pd}}).
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#' @param b A numeric matrix with \code{s@dimension} rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(struct_solve, covstruct) <- function(s, eta, b = NULL, ...) {
  l <- struct_factor(s, eta)
  backsolve(t(l), forwardsolve(l, b))
}

#' @title Default Factor
#' @name struct_factor.covstruct
#' @description Fallback: the lower Cholesky factor of
#'   \code{\link{struct_matrix}}, refused when the matrix is not positive
#'   definite spectrally.
#' @param s A \code{\link{covstruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A lower triangular numeric matrix.
#' @keywords internal
S7::method(struct_factor, covstruct) <- function(s, eta, ...) {
  l <- chol_pd(struct_matrix(s, eta))
  if (is.null(l)) {
    stop(sprintf(paste0(
      "'%s' is not positive definite at this eta, although it declares full\n",
      "  rank. The verdict is spectral, so this is a statement about the\n",
      "  matrix rather than about a factorisation that happened to fail."
    ), s@struct_name), call. = FALSE)
  }
  l
}

#' @title Refusal to Invert Without a Closed Form
#' @name struct_free.covstruct
#' @description
#' The base class refuses rather than inverting the map numerically: an
#' optimisation-based inverse would return a plausible \eqn{\eta} for a matrix
#' outside the set the family parametrises.
#' @param s A \code{\link{covstruct}} object.
#' @param m A symmetric numeric matrix.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @keywords internal
S7::method(struct_free, covstruct) <- function(s, m, ...) {
  stop(sprintf(paste0(
    "'%s' does not implement struct_free(). The inverse map is exact or\n",
    "  refused, never obtained by optimisation, because a numerical inverse\n",
    "  would return a plausible eta for a matrix outside the set."
  ), s@struct_name), call. = FALSE)
}
