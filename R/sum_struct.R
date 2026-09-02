#' @include numerical_fallbacks.R
NULL

#' A Non-Negative Combination of Fixed Matrices
#'
#' @description
#' The S7 class of a matrix parameter that is a sum of **fixed** symmetric
#' positive semidefinite matrices, each carried by one positive free value.
#' [sum_struct()] builds one.
#'
#' It is the variance-components covariance, and the one family here whose free
#' values are weights rather than entries: `n_free` is the number of components,
#' whatever the dimension.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `SumStructParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `components`, `link` and `labels`. `n_free` is the number of components and
#'   `rank` is fixed at construction from their shared null space.
#'
#' @seealso [sum_struct()], the constructor, [scaled_matrix()] for the
#'   one-component case with the matrix fixed, and [matrix_parameter()] for the
#'   properties this inherits.
#'
#' @examples
#' # One free value per component, whatever the side of the matrices.
#' s <- sum_struct(list(between = matrix(1, 3, 3), within = diag(3)))
#' c(dimension = s@dimension, n_free = s@n_free)
#' s@free_names
#'
#' @export
SumStructParam <- S7::new_class("SumStructParam", parent = matrix_parameter)


#' Construct a Sum of Fixed Matrices
#'
#' @description
#' Returns an object holding \eqn{M(\eta) = \sum_k c_k(\eta_k) P_k}, for fixed
#' symmetric positive semidefinite \eqn{P_1, \ldots, P_K} and positive weights
#' carried through a link. The matrices are given once and never move; what the
#' free vector carries is the \eqn{K} weights.
#'
#' @details
#' # What it is for
#'
#' This is the variance-components covariance \eqn{\sum_k \sigma_k^2 Z_kZ_k^\top},
#' and it is also the matrix a penalty with one smoothing parameter per component
#' assembles. \pkg{penalties7}'s `additive_penalty()` builds the same sum for its
#' own purposes; the difference is that a penalty is a function of the
#' coefficients while this is a matrix map, so a distribution can take it as a
#' covariance.
#'
#' # The value is linear, so most derivative components vanish
#'
#' \deqn{\partial^{m}_{\eta_k} M = c_k^{(m)}(\eta_k)\, P_k,}
#'
#' and a component whose indices name two different weights is **exactly zero**:
#' at \eqn{K = 2} that is 1 of the 3 second-order components, 2 of 4 at third
#' order and 3 of 5 at fourth. What survives is one number times a fixed matrix.
#'
#' # The log-determinant is not separable
#'
#' Its derivatives in the weights come from the cyclic trace expansion
#'
#' \deqn{\frac{\partial^{n}\log\lvert M\rvert}
#'            {\partial c_{k_1}\cdots\partial c_{k_n}}
#'       = (-1)^{n-1}\sum_{\sigma}
#'         \operatorname{tr}\bigl(M^{-1}P_{\sigma(1)}\cdots
#'                                M^{-1}P_{\sigma(n)}\bigr),}
#'
#' the sum running over the \eqn{(n-1)!} cyclic orderings **counted with
#' multiplicity**, and are then carried onto the free scale by a chain rule whose
#' Jacobian is diagonal. Counting with multiplicity is load bearing and the cost
#' of getting it wrong is measured: deduplicating the orderings that coincide
#' when an index repeats leaves a third-order component too small by exactly 2
#' and a fourth-order one by exactly 6.
#'
#' Against one stencil on the analytic order below, the four orders agree to
#' \eqn{7 \times 10^{-11}}, \eqn{1 \times 10^{-11}}, \eqn{3 \times 10^{-12}} and
#' \eqn{3 \times 10^{-13}}.
#'
#' # The rank is fixed at construction
#'
#' The null space of a sum of positive semidefinite matrices is the
#' **intersection** of theirs, so it does not move with the weights, and it is
#' read from the components stacked and individually normalized, never from an
#' assembled matrix. The distinction is measurable, not merely conceptual. Take
#' \eqn{P_1 = \mathbf{1}\mathbf{1}^\top} of side 4, of rank 1, and \eqn{P_2} the
#' first-difference penalty, of rank 3: their null spaces meet only at the
#' origin, so the family has rank 4, and counting eigenvalues of
#' \eqn{M(\eta)} above \eqn{10^{-10}} of the largest gives
#'
#' | weight ratio | \eqn{10^{0}} | \eqn{10^{3.5}} | \eqn{10^{6.9}} | \eqn{10^{10.4}} | \eqn{10^{13.9}} |
#' |---|---|---|---|---|---|
#' | eigenvalue count | 4 | 4 | 4 | 1 | 1 |
#'
#' while the family reports 4 throughout. Weights ten orders of magnitude apart
#' are an ordinary fitted model. Where there **is** a shared null space it is
#' held exactly: on the first- and second-difference penalties over five points,
#' whose null spaces meet in the constants, the residual
#' \eqn{\lVert M N\rVert / \lVert M\rVert} stays at \eqn{3 \times 10^{-16}} over
#' fourteen decades of weight ratio.
#'
#' # What a call costs
#'
#' The log-determinant's fourth derivatives are the expensive quantity, the
#' expansion evaluating \eqn{(n-1)!} matrix chains per component and one set
#' partition per repeated index. Seconds per call at side 6, over repetition
#' loops sized by elapsed time:
#'
#' | \eqn{K} | `param_value` | `param_d4` | `param_dlogdet` | `param_d4logdet` |
#' |---|---|---|---|---|
#' | 2 | 0.00003 | 0.00017 | 0.00029 | 0.00172 |
#' | 3 | 0.00002 | 0.00022 | 0.00037 | 0.00344 |
#' | 5 | 0.00003 | 0.00068 | 0.00039 | 0.01328 |
#'
#' @section Notation:
#' \eqn{K} is the number of components, \eqn{P_k} the \eqn{k}-th fixed matrix,
#' \eqn{c_k = h(\eta_k)} its weight, \eqn{p} the side, and \eqn{N} a basis of the
#' shared null space.
#'
#' @param components A non-empty list of symmetric positive semidefinite numeric
#'   matrices of the same side. Each is checked for all three properties, the
#'   semidefiniteness spectrally at a relative tolerance of \eqn{10^{-8}}. Named
#'   entries supply the free-value labels, which must be unique; unnamed ones are
#'   `w1`, `w2`, ...
#' @param link The positive link carrying each weight onto the free scale,
#'   [linkfunctions7::log_link()] by default. It must map onto the positive half
#'   line and from the whole real line; see [diagonal_matrix()] for the two
#'   conditions.
#'
#' @return An object of class [SumStructParam()], with `n_free` equal to the
#'   number of components, `free_names` the tagged labels, `null_basis` an
#'   orthonormal basis of the components' shared null space, and `rank` the side
#'   less its width.
#'
#' @seealso [scaled_matrix()] for one fixed matrix and one weight,
#'   [block_diag()], [kron_identity()] and [dr_prod()] for the other
#'   compositions, and `penalties7::additive_penalty()`, which assembles the same
#'   sum as a penalty.
#'
#' @examples
#' # Two variance components on three coefficients.
#' s <- sum_struct(list(between = matrix(1, 3, 3), within = diag(3)))
#' s@free_names
#' eta <- c(log(0.5), log(2))
#' param_value(s, eta)
#'
#' # The value is linear in the weights, so a component naming two of them is
#' # exactly zero, and a pure one is the weight's derivative times its matrix.
#' c(mixed = max(abs(param_d2(s, eta)[["log_between:log_within"]])),
#'   pure = max(abs(param_d1(s, eta)[[1]] - 0.5 * matrix(1, 3, 3))))
#'
#' # The log-determinant is not separable: the mixed second derivative is as
#' # large as the pure ones.
#' param_d2logdet(s, eta)
#'
#' # A single component of rank 1 gives a deficient family, and the
#' # log-determinant is then the log pseudo-determinant.
#' d <- sum_struct(list(matrix(1, 3, 3)))
#' c(rank = d@rank, logdet = param_logdet(d, log(2)), check = log(3 * 2))
#'
#' # The round trip closes.
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' @export
sum_struct <- function(components, link = linkfunctions7::log_link()) {
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
    param_params = list(components = components, link = link, labels = labels)
  )
}


#' The Shared Null Space of a Set of Positive Semidefinite Matrices
#'
#' @description
#' Returns an orthonormal basis of the **intersection** of the components' null
#' spaces, which is the null space of every non-negative combination of them: a
#' quadratic form \eqn{v^\top M v = \sum_k c_k\, v^\top P_k v} is a sum of
#' non-negative terms, so it vanishes only where every term does. It is therefore
#' a property of the family and can be computed once, at construction.
#'
#' @details
#' The components are stacked into one tall matrix and its right singular vectors
#' at negligible singular values are returned, which is the intersection of the
#' row spaces' orthogonal complements.
#'
#' Each component is **divided by its own largest entry** before stacking.
#' Without that normalization a component whose scale is many orders below
#' another's sinks below the tolerance and is read as absent, which is exactly
#' the failure a rank taken from an assembled matrix shows: see [sum_struct()]
#' for the measured table, where the eigenvalue count of \eqn{M(\eta)} falls from
#' 4 to 1 as the weights spread while the family's rank stays 4.
#'
#' The tolerance is LAPACK's usual one, `max(dim) * eps * max(d)`. Where the
#' stacked matrix has fewer singular values than the side, the missing directions
#' are null by construction and are appended.
#'
#' @param components A list of symmetric positive semidefinite matrices of the
#'   same side, already checked by [sum_struct()].
#'
#' @return A numeric matrix with `nrow(components[[1]])` rows whose columns are
#'   an orthonormal basis of the shared null space, with **zero columns** where
#'   the family has full rank.
#'
#' @seealso [sum_struct()], the only caller, and [param_null_basis()] for what
#'   the basis is used for.
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
#' carries it, for every component at once, as a matrix with one row per order.
#'
#' @details
#' Each weight depends on one free value only, so the table is complete: there
#' are no cross-derivatives between weights to record, and that is the whole
#' reason a derivative of the value naming two weights is zero. It is computed
#' once per call and read for every component of the order.
#'
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`.
#'
#' @return A 5 by \eqn{K} numeric matrix, row \eqn{m+1} holding the \eqn{m}-th
#'   derivative of the inverse link at each free value, so row 1 is the weights
#'   themselves.
#'
#' @seealso [sum_struct_derivs()] and [sum_struct_logdet_derivs()], the two
#'   callers.
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
#' The value is linear in the weights, so a component is zero unless every index
#' of the tuple names the same free value, and it is then that weight's
#' \eqn{m}-th derivative times its own fixed matrix.
#'
#' @details
#' No arithmetic is done on the matrices at all: a surviving component is one
#' scalar times a component the object has held since construction, and the rest
#' share a single zero matrix. At \eqn{K = 2} that is 1 of the 3 second-order
#' components, 2 of 4 at third order and 3 of 5 at fourth.
#'
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.
#'
#' @seealso [sum_struct_weight_derivs()] for the scalars, and
#'   [param_d1.SumStructParam()], which calls this.
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
  # the dimnames convention every family's matrices carry; one point
  # here covers all four derivative orders, which route through this
  stats::setNames(lapply(out, name_dims, s = s),
                  param_tuple_names(s, order))
}

#' Orderings of a Multiset, Counted With Multiplicity
#'
#' @description
#' All \eqn{n!} orderings of a vector of indices, **without** deduplicating the
#' ones that coincide because an index repeats: `multiset_orderings(c(1, 1))`
#' returns two elements.
#'
#' @details
#' The distinction is load bearing, and its cost is measured.
#' The cyclic sum behind the log-determinant expansion runs over \eqn{(n-1)!}
#' orderings, and two that happen to be equal still count twice; deduplicating
#' them leaves a third derivative in one weight too small by exactly 2 and a
#' fourth by exactly 6, which is \eqn{2!} and \eqn{3!}, the orderings of the tail.
#' Both are numbers a reader would accept without noticing.
#'
#' @param v An integer vector, of length 0 to 3 in every call the package makes,
#'   the expansion at order \eqn{n} permuting the \eqn{n-1} indices after the
#'   first.
#'
#' @return A list of `factorial(length(v))` integer vectors.
#'
#' @seealso [sum_struct_trace_term()], the only caller, and [combinat_perms()],
#'   which does the same for the matrix logarithm's contraction.
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
#' Evaluates the cyclic trace expansion
#' \eqn{(-1)^{n-1}\sum_{\sigma} \operatorname{tr}(M^{-1}P_{\sigma(1)}\cdots
#' M^{-1}P_{\sigma(n)})} for one tuple of weight indices, the sum running over
#' the \eqn{(n-1)!} orderings of the tail with the first index held.
#'
#' @details
#' Holding the first index is what leaves the orderings cyclic: a cyclic
#' permutation leaves the trace unchanged, so fixing one position counts each
#' distinct cycle once. The orderings come from [multiset_orderings()] and are
#' counted with multiplicity, which see.
#'
#' This is a derivative in the **weights**, not in the free values;
#' [sum_struct_logdet_derivs()] carries it onto the free scale.
#'
#' @param minv The inverse of the assembled matrix, computed once by the caller
#'   and reused for every component of the order.
#' @param comp The list of fixed components.
#' @param t An integer vector of weight indices, with repeats, of length 1 to 4.
#'
#' @return A single number.
#'
#' @seealso [sum_struct_logdet_derivs()], the only caller, and [sum_struct()] for
#'   the expansion.
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
#' Carries [sum_struct_trace_term()]'s expansion in the weights onto the free
#' scale. The map from free values to weights is diagonal, so the chain rule
#' groups the tuple by index and takes one set partition per group, each block
#' contributing a derivative of the inverse link and one differentiation in that
#' weight.
#'
#' @details
#' This is Faa di Bruno with a diagonal inner map, written out here instead of
#' going through [compose4()], because the outer function is a derivative in
#' several weights at once, never a univariate composition. The product over the groups'
#' partitions is the grid the loop walks; the partitions come from
#' [numericals7::set_partitions()], the one enumeration the toolkit keeps.
#'
#' It is the dearest quantity this family computes: at side 6 and \eqn{K = 5} the
#' fourth order costs 0.013 s against 0.0007 s for the value's own fourth
#' derivatives, the expansion evaluating a chain of matrix products per ordering.
#' Against one stencil on the analytic order below, the four orders agree to
#' \eqn{7 \times 10^{-11}}, \eqn{1 \times 10^{-11}}, \eqn{3 \times 10^{-12}} and
#' \eqn{3 \times 10^{-13}}.
#'
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
#'   as `param_tuple_names(s, order)` and in that order.
#'
#' @seealso [sum_struct_trace_term()] for the inner expansion, and
#'   [param_dlogdet.SumStructParam()], which calls this.
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
#' @description
#' The weighted sum \eqn{\sum_k c_k(\eta_k) P_k}, accumulated in place over the
#' components. Positive semidefinite at every free vector, a non-negative
#' combination of positive semidefinite matrices being one; positive **definite**
#' only where the components' null spaces meet at the origin, which is the
#' condition `rank` records.
#'
#' The value is labeled `v1`, `v2`, ..., `vp` on both margins, the convention [name_dims()]
#' states and every family in the package follows.
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A symmetric positive semidefinite `s@dimension` by `s@dimension`
#'   numeric matrix, labeled `v1`, `v2`, ..., `vp` on both margins.
#' @seealso [param_free.SumStructParam()] for the inverse, and [sum_struct()] for
#'   the parametrization.
#' @keywords internal
S7::method(param_value, SumStructParam) <- function(s, eta, ...) {
  w <- .ss_weights(s, eta)
  comp <- .ss(s)$components
  out <- w[[1L]] * comp[[1L]]
  for (k in seq_along(comp)[-1L]) out <- out + w[[k]] * comp[[k]]
  name_dims(out, s)
}

#' @title Free Vector of a Sum of Fixed Matrices
#' @name param_free.SumStructParam
#' @description
#' Recovers the weights by least squares on the components' entries, stacking the
#' \eqn{P_k} as columns and solving against the entries of `m`. Exact where `m` is
#' in the span: measured on two variance components over three coefficients, the
#' round trip closes to \eqn{3 \times 10^{-16}}.
#' @details
#' Two refusals, and they are different failures:
#'
#' - a matrix the combination cannot reproduce, checked by residual at
#'   \eqn{10^{-8}} relative to `max(1, max(abs(m)))`. The span is
#'   \eqn{K}-dimensional inside the \eqn{p(p+1)/2}-dimensional space of symmetric
#'   matrices, so almost every matrix is outside it and this is the ordinary
#'   failure.
#' - a matrix in the span needing a **non-positive weight**, which the link
#'   cannot carry. Such a matrix may still be positive semidefinite, so the
#'   refusal is about the parametrization, never about the matrix.
#'
#' The least-squares solve is exact wherever it succeeds, the residual check
#' being what separates a solution from a projection.
#' @param s A [SumStructParam()] object.
#' @param m A symmetric `s@dimension` by `s@dimension` matrix lying in the
#'   non-negative span of the components, already checked for shape and symmetry
#'   by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, the linked weights.
#' @seealso [param_value.SumStructParam()], the map this inverts.
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
#' `param_d1()`, `param_d2()`, `param_d3()` and `param_d4()` for a
#' [sum_struct()] parameter. The value being linear in the weights, a component
#' is **exactly zero** unless every index names the same free value, and is then
#' that weight's \eqn{m}-th derivative times its own fixed matrix. Nothing is
#' differenced and no matrix arithmetic is done.
#' @details
#' The four share [sum_struct_derivs()] and differ only in the order they pass.
#' At \eqn{K = 2} the zeros are 1 of the 3 second-order components, 2 of 4 at
#' third order and 3 of 5 at fourth. Compare
#' [param_dlogdet.SumStructParam()], where nothing vanishes: the value is linear
#' in the weights and its log-determinant is not.
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return At order 1, a list of `s@n_free` symmetric matrices named by
#'   `s@free_names`; above it, `choose(s@n_free + k - 1, k)` of them keyed as
#'   `param_tuple_names(s, k)` and in that order. Each is `s@dimension` by
#'   `s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.
#' @seealso [sum_struct_derivs()], which assembles them.
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
#' The log-determinant of the assembled matrix where the family has full rank,
#' and its log **pseudo**-determinant, the sum over the `s@rank` largest
#' eigenvalues, where it does not. A sum of fixed matrices has no structure a
#' closed form could exploit, so this is the one family here that assembles and
#' decomposes.
#' @details
#' Which branch is taken is decided by `s@rank`, fixed at construction from the
#' components' shared null space, and never by counting eigenvalues at the point:
#' that count falls as the weights spread apart, which [sum_struct()] measures.
#' So the number of eigenvalues summed does not move as a fit walks the free
#' vector.
#'
#' Measured on `sum_struct(list(matrix(1, 3, 3)))`, of rank 1, at a weight of 2:
#' the pseudo-determinant is \eqn{\log 6}, its one non-zero eigenvalue, to the
#' printed digit.
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.SumStructParam()] for its derivatives, and
#'   [param_null_basis()] for the directions the pseudo-determinant leaves out.
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
#' `param_dlogdet()`, `param_d2logdet()`, `param_d3logdet()` and
#' `param_d4logdet()` for a [sum_struct()] parameter: the cyclic trace expansion
#' in the weights, carried onto the free scale by a chain rule with a diagonal
#' Jacobian.
#' @details
#' The four share [sum_struct_logdet_derivs()] and differ only in the order they
#' pass. **Nothing vanishes here**, unlike the value's own derivatives: at two
#' variance components the mixed second derivative is \eqn{-0.245} against
#' \eqn{0.245} for each pure one, the same size. The log-determinant of a sum is
#' not a sum, which is why this family needs an expansion where [block_diag()]
#' and [dr_prod()] need none.
#'
#' It is the dearest quantity the family computes; [sum_struct()] carries the
#' timings and the accuracy against a stencil.
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector: at order 1, `s@n_free` values named by
#'   `s@free_names`; above it, `choose(s@n_free + k - 1, k)` values keyed as
#'   `param_tuple_names(s, k)` and in that order.
#' @seealso [sum_struct_logdet_derivs()], which assembles them, and
#'   [param_logdet.SumStructParam()] for the quantity differentiated.
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
#' `param_solve()` and `param_factor()` for a [sum_struct()] parameter. Both come
#' from the assembled matrix, by `solve()` and by a transposed `chol()`: a sum of
#' fixed matrices has no structure a solve could exploit, where [ar1()] has a
#' tridiagonal inverse and [log_cholesky()] holds its factor already.
#' @details
#' `param_solve()` therefore agrees with `solve()` on the assembled matrix
#' exactly, being the same call; the factor is `t(chol(M))`, lower triangular
#' with \eqn{M = L L^\top} as [param_factor()] requires, and reproduces the matrix
#' to \eqn{4 \times 10^{-16}}.
#'
#' Both are refused by the generic where the family is rank deficient, a singular
#' matrix having neither an inverse nor a Cholesky factor.
#' @param s A [SumStructParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param b A numeric matrix with `s@dimension` rows, defaulted to the identity
#'   by the generic. `param_factor()` takes no `b`.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return `param_solve()` returns a numeric matrix with `s@dimension` rows and
#'   as many columns as `b`; `param_factor()` a lower triangular `s@dimension` by
#'   `s@dimension` matrix. Both are bare, where the value and the derivative
#'   arrays carry `v1`, `v2`, ...: this is the one family whose two answers are
#'   built from its own value rather than from its structure, so they are
#'   unnamed explicitly to match the twelve families that never label them.
#' @seealso [param_solve()] and [param_factor()] for the two contracts.
#' @keywords internal
S7::method(param_solve, SumStructParam) <- function(s, eta, b = NULL, ...) {
  # unname because this is the one family whose solve delegates to the value:
  # every other computes it from its own structure, so all of them return a
  # bare matrix, and the dimnames convention covers the value and the four
  # derivative orders rather than everything derived from them.
  solve(unname(param_value(s, eta)), b)
}

#' @rdname param_solve.SumStructParam
#' @name param_factor.SumStructParam
#' @keywords internal
S7::method(param_factor, SumStructParam) <- function(s, eta, ...) {
  # unname for the reason the solve does: the primitives disagree about whether
  # a factor carries the labels, and this commit leaves that where it found it.
  t(chol(unname(param_value(s, eta))))
}
