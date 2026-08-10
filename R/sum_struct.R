#' @include numerical_fallbacks.R
NULL

#' A Non-Negative Combination of Fixed Matrices
#'
#' @description
#' The S7 class of a matrix parameter that is a sum of fixed symmetric positive
#' semidefinite matrices, each carried by one positive free value. Constructed
#' by \code{\link{sum_struct}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{SumStructParam}.
#'
#' @seealso \code{\link{sum_struct}}
#'
#' @examples
#' S7::S7_inherits(sum_struct(list(diag(3), matrix(1, 3, 3))), SumStructParam)
#'
#' @export
SumStructParam <- S7::new_class("SumStructParam", parent = matrix_parameter)


#' Construct a Sum of Fixed Matrices
#'
#' @description
#' \eqn{M(\eta) = \sum_k c_k(\eta_k) P_k} for fixed symmetric positive
#' semidefinite \eqn{P_1, \ldots, P_K} and positive weights carried through a
#' link.
#'
#' @details
#' This is the variance-components covariance, \eqn{\sum_k \sigma_k^2 Z_kZ_k'},
#' and it is also the matrix a penalty with one smoothing parameter per
#' component assembles. \pkg{penalties7}'s \code{additive_penalty()} builds the
#' same sum for its own purposes; the difference is that a penalty is a
#' function of the coefficients while this is a matrix map, so a distribution
#' can take it as a covariance.
#'
#' The value is linear in the weights, so a derivative component is zero unless
#' every index of the tuple names the same free value, and it then carries the
#' corresponding derivative of the inverse link times that component:
#' \deqn{\partial^{m}_{\eta_k} M = c_k^{(m)}(\eta_k)\, P_k.}
#' The log-determinant is not separable, and its derivatives in the weights
#' come from the expansion
#' \deqn{\frac{\partial^{n}\log\lvert M\rvert}
#'            {\partial c_{k_1}\cdots\partial c_{k_n}}
#'       = (-1)^{n-1}\sum_{\sigma}
#'         \operatorname{tr}\bigl(M^{-1}P_{\sigma(1)}\cdots
#'                                M^{-1}P_{\sigma(n)}\bigr),}
#' the sum running over the \eqn{(n-1)!} cyclic orderings counted with
#' multiplicity, and are then carried onto the free scale by a chain rule whose
#' Jacobian is diagonal.
#'
#' \strong{Rank.} The null space of a sum of positive semidefinite matrices is
#' the intersection of theirs, so it does not move with the weights and the
#' rank is fixed at construction. It is read from the components stacked and
#' individually normalized, never from an assembled matrix: a count of small
#' eigenvalues of \eqn{M(\eta)} falls as the weights spread apart, which is an
#' ordinary fitted model and not a pathology.
#'
#' @param components A list of symmetric matrices of the same side, each
#'   positive semidefinite. Named entries supply the free-value labels.
#' @param link The positive link carrying each weight onto the free scale.
#'   Defaults to \code{\link[linkfunctions7]{log_link}()}.
#' @param role One of \code{"covariance"}, \code{"precision"} or
#'   \code{"either"}.
#'
#' @return An object of class \code{\link{SumStructParam}}.
#'
#' @seealso \code{\link{scaled_matrix}}, \code{\link{block_diag}},
#'   \code{\link[penalties7]{additive_penalty}}
#'
#' @examples
#' # two variance components on three coefficients
#' s <- sum_struct(list(between = matrix(1, 3, 3), within = diag(3)))
#' s@free_names
#' param_value(s, c(log(0.5), log(2)))
#' param_logdet(s, c(log(0.5), log(2)))
#'
#' @export
sum_struct <- function(components, link = linkfunctions7::log_link(),
                       role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  if (!is.list(components) || !length(components)) {
    stop("'components' must be a non-empty list of matrices.", call. = FALSE)
  }
  components <- lapply(components, function(m) {
    if (!is.matrix(m) || !is.numeric(m)) {
      stop("Every component must be a numeric matrix.", call. = FALSE)
    }
    unname(as.matrix(m))
  })
  p <- nrow(components[[1L]])
  for (m in components) {
    if (!identical(dim(m), c(p, p))) {
      stop("Every component must be square and of the same side.", call. = FALSE)
    }
    if (max(abs(m - t(m))) > 1e-10 * max(1, max(abs(m)))) {
      stop("Every component must be symmetric.", call. = FALSE)
    }
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8 * max(1, max(ev))) {
      stop("Every component must be positive semidefinite.", call. = FALSE)
    }
  }
  check_positive_link(link)
  p <- as.integer(p)
  K <- length(components)

  labels <- names(components)
  if (is.null(labels)) labels <- rep("", K)
  blank <- !nzchar(labels)
  labels[blank] <- paste0("w", seq_len(K))[blank]
  if (anyDuplicated(labels)) {
    stop("Component labels must be unique.", call. = FALSE)
  }

  nb <- sum_struct_null_basis(components)
  SumStructParam(
    param_name = sprintf("sum_struct(%d)", K),
    dimension = p,
    n_free = K,
    free_names = tagged_name(link, labels),
    rank = as.integer(p - ncol(nb)),
    null_basis = nb,
    role = role,
    param_params = list(components = components, link = link, labels = labels)
  )
}


#' The Shared Null Space of a Set of Positive Semidefinite Matrices
#'
#' @description
#' Returns an orthonormal basis of the intersection of the components' null
#' spaces, which is the null space of every non-negative combination of them
#' and so a property of the family rather than of a point.
#'
#' @details
#' The components are stacked after each is divided by its own largest entry.
#' Without that normalization a component whose scale is many orders below
#' another's sinks below the tolerance and is read as absent, which is the
#' failure a rank taken from an assembled matrix already shows.
#'
#' @param components A list of symmetric positive semidefinite matrices.
#'
#' @return A matrix whose columns are an orthonormal basis, with zero columns
#'   when the sum has full rank.
#'
#' @keywords internal
sum_struct_null_basis <- function(components) {
  p <- nrow(components[[1L]])
  stacked <- do.call(rbind, lapply(components, function(m) {
    s <- max(abs(m))
    if (s > 0) m / s else m
  }))
  sv <- svd(stacked)
  tol <- max(dim(stacked)) * .Machine$double.eps * max(sv$d)
  keep <- sv$d <= tol
  if (length(sv$d) < p) keep <- c(keep, rep(TRUE, p - length(sv$d)))
  sv$v[, keep, drop = FALSE]
}

.ss <- function(s) s@param_params
.ss_weights <- function(s, eta) linkfunctions7::linkinv(.ss(s)$link, eta)

#' Derivatives of the Weights of a Sum of Fixed Matrices
#'
#' @description
#' Returns each weight and its first four derivatives in the free value that
#' carries it, as a matrix with one row per order.
#'
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#'
#' @return A 5 by \code{K} numeric matrix, rows being orders 0 to 4.
#'
#' @keywords internal
sum_struct_weight_derivs <- function(s, eta) {
  lk <- .ss(s)$link
  rbind(linkfunctions7::linkinv(lk, eta),
        linkfunctions7::dlinkinv(lk, eta),
        linkfunctions7::d2linkinv(lk, eta),
        linkfunctions7::d3linkinv(lk, eta),
        linkfunctions7::d4linkinv(lk, eta))
}

#' Assemble a Sum of Fixed Matrices' Derivatives of a Given Order
#'
#' @description
#' The value is linear in the weights, so a component is zero unless every
#' index of the tuple names the same free value.
#'
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of matrices keyed as \code{param_tuple_names(s, order)}.
#'
#' @keywords internal
sum_struct_derivs <- function(s, eta, order) {
  cd <- sum_struct_weight_derivs(s, eta)
  comp <- .ss(s)$components
  p <- s@dimension
  zero <- matrix(0, p, p)
  out <- lapply(param_tuple_indices(s, order), function(t) {
    k <- unique(t)
    if (length(k) > 1L) return(zero)
    cd[order + 1L, k] * comp[[k]]
  })
  stats::setNames(out, param_tuple_names(s, order))
}

#' Orderings of a Multiset, Counted With Multiplicity
#'
#' @description
#' All \eqn{n!} orderings of a vector of indices, without deduplicating the
#' ones that coincide because an index repeats.
#'
#' @details
#' The distinction is load bearing. The cyclic sum behind the log-determinant
#' expansion runs over \eqn{(n-1)!} orderings, and two that happen to be equal
#' still count twice; deduplicating them makes the third derivative of a
#' component in one weight too small by a factor of two and the fourth by six.
#'
#' @param v An integer vector.
#'
#' @return A list of integer vectors.
#'
#' @keywords internal
multiset_orderings <- function(v) {
  if (!length(v)) return(list(integer(0)))
  if (length(v) == 1L) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    for (rest in multiset_orderings(v[-i])) {
      out[[length(out) + 1L]] <- c(v[i], rest)
    }
  }
  out
}

#' Log-Determinant Derivatives in the Weights
#'
#' @description
#' Evaluates the cyclic trace expansion for one tuple of weight indices.
#'
#' @param minv The inverse of the assembled matrix.
#' @param comp The list of components.
#' @param t An integer vector of weight indices, with repeats.
#'
#' @return A single number.
#'
#' @keywords internal
sum_struct_trace_term <- function(minv, comp, t) {
  n <- length(t)
  acc <- 0
  for (rest in multiset_orderings(t[-1L])) {
    ord <- c(t[1L], rest)
    term <- minv %*% comp[[ord[1L]]]
    for (j in ord[-1L]) term <- term %*% minv %*% comp[[j]]
    acc <- acc + sum(diag(term))
  }
  (-1)^(n - 1L) * acc
}

#' Assemble a Sum of Fixed Matrices' Log-Determinant Derivatives
#'
#' @description
#' Carries the trace expansion in the weights onto the free scale. The map is
#' diagonal, so the chain rule groups the tuple by index and takes one set
#' partition per group, each block contributing a derivative of the inverse
#' link and one differentiation in that weight.
#'
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named numeric vector keyed as \code{param_tuple_names(s, order)}.
#'
#' @keywords internal
sum_struct_logdet_derivs <- function(s, eta, order) {
  comp <- .ss(s)$components
  cd <- sum_struct_weight_derivs(s, eta)
  minv <- solve(param_value(s, eta))
  out <- vapply(param_tuple_indices(s, order), function(t) {
    idx <- sort(unique(t))
    mult <- tabulate(match(t, idx), length(idx))
    # one set partition per distinct index; the product over the choices is
    # the diagonal chain rule
    parts <- lapply(mult, numericals7::set_partitions)
    grid <- expand.grid(lapply(parts, seq_along))
    acc <- 0
    for (r in seq_len(nrow(grid))) {
      coef <- 1
      inner <- integer(0)
      for (g in seq_along(idx)) {
        pk <- parts[[g]][[grid[r, g]]]
        for (b in pk) coef <- coef * cd[length(b) + 1L, idx[g]]
        inner <- c(inner, rep(idx[g], length(pk)))
      }
      acc <- acc + coef * sum_struct_trace_term(minv, comp, inner)
    }
    acc
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, order))
}


#' @title Value of a Sum of Fixed Matrices
#' @name param_value.SumStructParam
#' @description The weighted sum of the components.
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A symmetric matrix of side \code{s@dimension}.
#' @keywords internal
S7::method(param_value, SumStructParam) <- function(s, eta, ...) {
  w <- .ss_weights(s, eta)
  comp <- .ss(s)$components
  out <- w[[1L]] * comp[[1L]]
  for (k in seq_along(comp)[-1L]) out <- out + w[[k]] * comp[[k]]
  out
}

#' @title Free Vector of a Sum of Fixed Matrices
#' @name param_free.SumStructParam
#' @description
#' Recovers the weights by least squares on the components' entries and rejects
#' a matrix the combination cannot reproduce. Non-positive weights are rejected
#' too, the family carrying them through a positive link.
#' @param s A \code{\link{SumStructParam}} object.
#' @param m A symmetric matrix of side \code{s@dimension}.
#' @param ... Unused.
#' @return A numeric vector of length \code{s@n_free}.
#' @keywords internal
S7::method(param_free, SumStructParam) <- function(s, m, ...) {
  comp <- .ss(s)$components
  X <- do.call(cbind, lapply(comp, as.numeric))
  w <- tryCatch(qr.solve(X, as.numeric(m)), error = function(e) NULL)
  if (is.null(w) ||
      max(abs(X %*% w - as.numeric(m))) > 1e-8 * max(1, max(abs(m)))) {
    stop("'m' is not a combination of this parameter's components.",
         call. = FALSE)
  }
  if (any(w <= 0)) {
    stop("'m' needs a non-positive weight, which the link cannot carry.",
         call. = FALSE)
  }
  linkfunctions7::linkfun(.ss(s)$link, w)
}

#' @title Derivatives of a Sum of Fixed Matrices
#' @name param_d1.SumStructParam
#' @description
#' The value being linear in the weights, a component is zero unless every
#' index names the same free value, and is then that weight's derivative times
#' its component.
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, SumStructParam) <- function(s, eta, ...) {
  stats::setNames(sum_struct_derivs(s, eta, 1L), s@free_names)
}

#' @rdname param_d1.SumStructParam
#' @name param_d2.SumStructParam
#' @keywords internal
S7::method(param_d2, SumStructParam) <- function(s, eta, ...) {
  sum_struct_derivs(s, eta, 2L)
}

#' @rdname param_d1.SumStructParam
#' @name param_d3.SumStructParam
#' @keywords internal
S7::method(param_d3, SumStructParam) <- function(s, eta, ...) {
  sum_struct_derivs(s, eta, 3L)
}

#' @rdname param_d1.SumStructParam
#' @name param_d4.SumStructParam
#' @keywords internal
S7::method(param_d4, SumStructParam) <- function(s, eta, ...) {
  sum_struct_derivs(s, eta, 4L)
}

#' @title Log-Determinant of a Sum of Fixed Matrices
#' @name param_logdet.SumStructParam
#' @description
#' The log-determinant of the assembled matrix, or its log pseudo-determinant
#' over the non-zero eigenvalues when the family is rank deficient.
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, SumStructParam) <- function(s, eta, ...) {
  m <- param_value(s, eta)
  if (s@rank == s@dimension) {
    return(as.numeric(determinant(m, logarithm = TRUE)$modulus))
  }
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  sum(log(sort(ev, decreasing = TRUE)[seq_len(s@rank)]))
}

#' @title Log-Determinant Derivatives of a Sum of Fixed Matrices
#' @name param_dlogdet.SumStructParam
#' @description
#' The cyclic trace expansion in the weights, carried onto the free scale by a
#' chain rule with a diagonal Jacobian.
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, SumStructParam) <- function(s, eta, ...) {
  stats::setNames(sum_struct_logdet_derivs(s, eta, 1L), s@free_names)
}

#' @rdname param_dlogdet.SumStructParam
#' @name param_d2logdet.SumStructParam
#' @keywords internal
S7::method(param_d2logdet, SumStructParam) <- function(s, eta, ...) {
  sum_struct_logdet_derivs(s, eta, 2L)
}

#' @rdname param_dlogdet.SumStructParam
#' @name param_d3logdet.SumStructParam
#' @keywords internal
S7::method(param_d3logdet, SumStructParam) <- function(s, eta, ...) {
  sum_struct_logdet_derivs(s, eta, 3L)
}

#' @rdname param_dlogdet.SumStructParam
#' @name param_d4logdet.SumStructParam
#' @keywords internal
S7::method(param_d4logdet, SumStructParam) <- function(s, eta, ...) {
  sum_struct_logdet_derivs(s, eta, 4L)
}

#' @title Solve and Factor of a Sum of Fixed Matrices
#' @name param_solve.SumStructParam
#' @description
#' Both come from the assembled matrix: a sum of fixed matrices has no
#' structure a solve could exploit, unlike the families whose factor is written
#' out.
#' @param s A \code{\link{SumStructParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param b A matrix with \code{s@dimension} rows, or \code{NULL}.
#' @param ... Unused.
#' @return A matrix.
#' @keywords internal
S7::method(param_solve, SumStructParam) <- function(s, eta, b = NULL, ...) {
  solve(param_value(s, eta), b)
}

#' @rdname param_solve.SumStructParam
#' @name param_factor.SumStructParam
#' @keywords internal
S7::method(param_factor, SumStructParam) <- function(s, eta, ...) {
  t(chol(param_value(s, eta)))
}
