#' @include numerical_fallbacks.R
NULL

#' Scales Times a Correlation Matrix
#'
#' @description
#' The S7 class of a covariance written as \eqn{D R D}, for a diagonal matrix
#' \eqn{D} of positive scales and a correlation matrix \eqn{R}. [dr_prod()]
#' builds one.
#'
#' Its free values are the standard deviations and the correlation's own
#' coordinates, so the quantities a reader takes off a fitted covariance are the
#' coordinates themselves.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `DrProdParam`, a subclass of [matrix_parameter()]
#'   adding no properties of its own. `param_params` holds `cor` (the correlation
#'   parameter), `link` and `p`. `n_free` is \eqn{p} plus the correlation's, and
#'   `rank` is \eqn{p}, this family admitting no deficiency.
#'
#' @seealso [dr_prod()], the constructor, [correlation_matrix()] for the default
#'   block, and [matrix_parameter()] for the properties this inherits.
#'
#' @examples
#' # The scales come first and the correlation's coordinates follow.
#' s <- dr_prod(3)
#' s@free_names
#'
#' # p standard deviations plus the correlation's p(p-1)/2 angles.
#' c(p = 4, n_free = dr_prod(4)@n_free)
#'
#' @export
DrProdParam <- S7::new_class("DrProdParam", parent = matrix_parameter)


#' Construct a Covariance as Scales Times a Correlation
#'
#' @description
#' Returns an object holding
#' \eqn{\Sigma(\eta) = D(\eta_D)\, R(\eta_R)\, D(\eta_D)}, with
#' \eqn{D = \mathrm{diag}(d_1, \ldots, d_p)} carrying the standard deviations
#' through a positive link and \eqn{R} a correlation matrix parameter.
#'
#' @details
#' # What the separation buys
#'
#' The quantities a reader takes off a fitted covariance are the standard
#' deviations and the correlations, and here they **are** the coordinates. A
#' [log_cholesky()] factor produces the same set of matrices with no coordinate
#' meaning anything on its own. Measured at
#' \eqn{\eta = (0, \log 2, \log 0.5, 1, 1.2, 0.9)}: `sqrt(diag(M))` is
#' \eqn{(1, 2, 0.5)}, which is `exp(eta[1:3])`, and dividing them out returns the
#' correlation block's own value exactly.
#'
#' # Every derivative factorizes
#'
#' Because \eqn{\Sigma_{ij} = d_i d_j R_{ij}} and the two groups of free values
#' are disjoint,
#'
#' \deqn{\partial^{S}\Sigma_{ij} =
#'       \bigl[\partial^{S_D}(d_i d_j)\bigr]\,
#'       \bigl[\partial^{S_R} R_{ij}\bigr],}
#'
#' where \eqn{S_D} and \eqn{S_R} are the parts of the multiset \eqn{S} falling in
#' each group. Nothing of the correlation family is rederived: its own components
#' are fetched and multiplied entrywise.
#'
#' Both factors are **sparse**, and where their supports miss each other the
#' component is exactly zero. The scale factor is supported on the rows and
#' columns the indices of \eqn{S_D} name, and vanishes altogether once \eqn{S_D}
#' names three distinct scales; a correlation derivative is supported on the two
#' entries its angle governs. So `log_sd1:z3.2` is zero, the first factor living
#' in row and column 1 and the second on entries \eqn{(3,2)} and \eqn{(2,3)}.
#' Measured at \eqn{p = 3}, the disjoint-support rule predicts every exact zero:
#' 1 of the 21 second-order components, 10 of 56 at third order and 37 of 126 at
#' fourth.
#'
#' # The log-determinant
#'
#' \deqn{\log\lvert\Sigma\rvert = 2\sum_j \log d_j + \log\lvert R\rvert,}
#'
#' separable in the scales and separable from the correlation, so a component
#' mixing two scales, or a scale with a correlation, is exactly zero: 12 of the 21
#' second-order components, 43 of 56 and 108 of 126. The first derivative in a
#' scale is 2 whatever the point, the scales entering on both sides.
#'
#' # The correlation block must have full rank
#'
#' A rank-deficient \eqn{R} would give \eqn{\Sigma} the null space
#' \eqn{D^{-1}\ker R}, which **moves with the free vector**, while the class
#' records the rank and the null space as properties of the family. The
#' constructor rejects one instead of recording a rank that is not stable.
#'
#' @section Notation:
#' \eqn{p} is the dimension, \eqn{d_j} the \eqn{j}-th standard deviation,
#' \eqn{D = \mathrm{diag}(d)}, \eqn{R} the correlation matrix, and \eqn{\eta_D}
#' and \eqn{\eta_R} the two stretches of the free vector.
#'
#' @param dimension The side \eqn{p} of the matrix, at least 2: a 1 by 1
#'   correlation carries nothing, and the constructor says so.
#' @param correlation A [matrix_parameter()] of side `dimension` **producing
#'   correlation matrices**, that is with a unit diagonal at every free vector.
#'   Defaults to [`correlation_matrix(dimension)`][correlation_matrix], which is
#'   the only shipped family with that property. The requirement is not checked;
#'   a block carrying a scale of its own, such as [ar1()], makes the composite
#'   unidentified, the same matrix arising from a whole ray of free vectors, and
#'   `check_parameter()`'s round trip is what reports it.
#' @param link The positive link carrying each standard deviation onto the free
#'   scale, [linkfunctions7::log_link()] by default. It must map onto the
#'   positive half line and from the whole real line; see [diagonal_matrix()] for
#'   the two conditions.
#' @param role One of `"covariance"` (the default), `"precision"` or `"either"`.
#'   No numeric result depends on it.
#'
#' @return An object of class [DrProdParam()], with `n_free` equal to \eqn{p}
#'   plus the correlation's, `free_names` the tagged `log_sd1` ... `log_sdp`
#'   followed by the correlation's own, `rank` equal to `dimension`, and an empty
#'   `null_basis`.
#'
#' @seealso [correlation_matrix()] for the default block, [log_cholesky()] for
#'   the same set of matrices in coordinates that mean nothing separately, and
#'   [block_diag()], [kron_identity()] and [sum_struct()] for the other
#'   compositions.
#'
#' @examples
#' s <- dr_prod(3)
#' s@free_names
#'
#' eta <- c(log(1), log(2), log(0.5), 1.0, 1.2, 0.9)
#' M <- param_value(s, eta)
#'
#' # The standard deviations are the coordinates, read back off the diagonal.
#' rbind(from_matrix = sqrt(diag(M)), from_eta = exp(eta[1:3]))
#'
#' # And dividing them out leaves the correlation block's own value.
#' R <- M / outer(sqrt(diag(M)), sqrt(diag(M)))
#' max(abs(R - param_value(correlation_matrix(3), eta[4:6])))
#'
#' # A second derivative is exactly zero where the two factors' supports miss
#' # each other: the scale factor lives in row and column 1, the correlation
#' # derivative on the (3, 2) entry.
#' max(abs(param_d2(s, eta)[["log_sd1:z3.2"]]))
#'
#' # The log-determinant separates, so the first derivative in a scale is 2.
#' param_dlogdet(s, eta)[1:3]
#'
#' # The round trip closes.
#' max(abs(param_free(s, M) - eta))
#'
#' @export
dr_prod <- function(dimension, correlation = NULL,
                    link = linkfunctions7::log_link(),
                    role = c("covariance", "precision", "either")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)
  if (p < 2L) {
    stop("'dimension' must be at least 2: a 1 by 1 correlation carries nothing.",
         call. = FALSE)
  }
  check_positive_link(link)
  if (is.null(correlation)) correlation <- correlation_matrix(p)
  if (!S7::S7_inherits(correlation, matrix_parameter)) {
    stop("'correlation' must inherit from 'matrix_parameter'.", call. = FALSE)
  }
  if (correlation@dimension != p) {
    stop(sprintf("'correlation' has side %d, not %d.",
                 correlation@dimension, p), call. = FALSE)
  }
  if (correlation@rank < p) {
    stop(paste0("'correlation' must have full rank: the null space of D R D is\n",
                "  D^-1 ker R, which moves with the free vector, while the rank\n",
                "  and the null space are recorded as properties of the family."),
         call. = FALSE)
  }

  DrProdParam(
    param_name = sprintf("dr_prod(%s)", correlation@param_name),
    dimension = p,
    n_free = as.integer(p + correlation@n_free),
    free_names = c(tagged_name(link, paste0("sd", seq_len(p))),
                   correlation@free_names),
    rank = p,
    null_basis = matrix(0, p, 0),
    role = role,
    param_params = list(cor = correlation, link = link, p = p)
  )
}

.dr <- function(s) s@param_params
.dr_scales <- function(s, eta) {
  linkfunctions7::linkinv(.dr(s)$link, eta[seq_len(.dr(s)$p)])
}
.dr_eta_cor <- function(s, eta) eta[-seq_len(.dr(s)$p)]

#' Derivatives of the Inverse Link at Every Scale Coordinate
#'
#' @description
#' Returns \eqn{d_j} and its first four derivatives in the free value that
#' carries it, for every \eqn{j} at once, as a matrix with one row per order.
#'
#' @details
#' Each scale depends on one free value only, so the table is complete: there are
#' no cross-derivatives between scales to record. It is computed once per call to
#' [dr_prod_derivs()] and read by [dr_scale_factor()] for every component of the
#' order.
#'
#' The correlation's part of `eta` is ignored, so one table serves every
#' component whatever the correlation indices are.
#'
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`; only its first \eqn{p}
#'   entries are read.
#'
#' @return A 5 by \eqn{p} numeric matrix, row \eqn{k+1} holding the \eqn{k}-th
#'   derivative of the inverse link at each scale coordinate, so row 1 is the
#'   scales themselves.
#'
#' @seealso [dr_scale_factor()], which reads it, and [diag_dlog()], the
#'   log-determinant's own chain.
#'
#' @keywords internal
dr_scale_derivs <- function(s, eta) {
  lk <- .dr(s)$link
  e <- eta[seq_len(.dr(s)$p)]
  rbind(linkfunctions7::linkinv(lk, e),
        linkfunctions7::dlinkinv(lk, e),
        linkfunctions7::d2linkinv(lk, e),
        linkfunctions7::d3linkinv(lk, e),
        linkfunctions7::d4linkinv(lk, e))
}

#' The Scale Factor of a Derivative Component
#'
#' @description
#' Evaluates \eqn{\partial^{S_D}(d_i d_j)} for every pair \eqn{(i, j)} at once,
#' given the multiset \eqn{S_D} of scale indices. It is zero wherever \eqn{S_D}
#' contains an index naming neither \eqn{i} nor \eqn{j}, so the result is
#' **supported on the rows and columns those indices name**, and that sparsity is
#' half of why the composition costs so little.
#'
#' @details
#' Three cases, and the arithmetic differs in each:
#'
#' - **no indices.** The factor is \eqn{d_i d_j}, an outer product, and that is
#'   the case a component differentiating the correlation alone falls into.
#' - **one distinct index \eqn{k}, with multiplicity \eqn{m}.** Off the diagonal
#'   the entry carries one factor of \eqn{d_k}, so it is
#'   \eqn{d_k^{(m)} d_j}; on the diagonal it carries two, so it is the Leibniz
#'   expansion \eqn{\sum_r \binom{m}{r} d_k^{(r)} d_k^{(m-r)}} of
#'   \eqn{\partial^m d_k^2}.
#' - **two distinct indices \eqn{i} and \eqn{j}.** Only the entries \eqn{(i,j)}
#'   and \eqn{(j,i)} survive, each \eqn{d_i^{(m_i)} d_j^{(m_j)}}.
#'
#' Three or more distinct indices give the zero matrix: an entry of
#' \eqn{\Sigma} carries at most two scales, so a third differentiation in a new
#' scale annihilates it.
#'
#' @param sd A 5 by \eqn{p} matrix of inverse-link derivatives, as returned by
#'   [dr_scale_derivs()].
#' @param tuple An integer vector of scale indices, possibly empty and possibly
#'   with repeats. Its length is the number of scale indices in the component,
#'   which is at most the derivative order.
#'
#' @return A symmetric \eqn{p} by \eqn{p} numeric matrix.
#'
#' @seealso [dr_prod_derivs()], the only caller, and [dr_prod()] for the
#'   factorization this is half of.
#'
#' @keywords internal
dr_scale_factor <- function(sd, tuple) {
  p <- ncol(sd)
  d <- sd[1L, ]
  if (!length(tuple)) return(outer(d, d))
  used <- sort(unique(tuple))
  if (length(used) > 2L) return(matrix(0, p, p))
  mult <- tabulate(match(tuple, used), length(used))
  out <- matrix(0, p, p)
  if (length(used) == 1L) {
    k <- used
    m <- mult
    # the diagonal entry is d_k^2 differentiated m times, and an off-diagonal
    # entry in row or column k carries one factor only
    out[k, ] <- sd[m + 1L, k] * d
    out[, k] <- sd[m + 1L, k] * d
    out[k, k] <- sum(choose(m, 0:m) * sd[0:m + 1L, k] * sd[m - (0:m) + 1L, k])
  } else {
    i <- used[1L]
    j <- used[2L]
    v <- sd[mult[1L] + 1L, i] * sd[mult[2L] + 1L, j]
    out[i, j] <- v
    out[j, i] <- v
  }
  out
}

#' Assemble a Scales-Times-Correlation Derivative of a Given Order
#'
#' @description
#' Splits each tuple of the composite's enumeration into its scale indices and
#' its correlation indices, and multiplies [dr_scale_factor()] by the
#' correlation's own component of the matching order, elementwise.
#'
#' @details
#' The correlation's arrays are fetched at most once per order and re-keyed by
#' the sorted local index tuple, the same device [block_derivs_by_tuple()] uses
#' and for the same reason: the composite's names are not the block's, but the
#' sorted index tuple is a key both sides can compute.
#'
#' A component whose correlation part is empty reads the correlation's **value**,
#' which is why `param_value()` appears in a derivative routine.
#'
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension` and without dimnames.
#'
#' @seealso [dr_scale_factor()] for one factor, and [param_d1.DrProdParam()],
#'   which calls this.
#'
#' @keywords internal
dr_prod_derivs <- function(s, eta, order) {
  pp <- .dr(s)
  p <- pp$p
  sd <- dr_scale_derivs(s, eta)
  ec <- .dr_eta_cor(s, eta)
  cor_cache <- list()
  cor_component <- function(t) {
    k <- length(t)
    if (k == 0L) return(param_value(pp$cor, ec))
    key <- as.character(k)
    if (is.null(cor_cache[[key]])) {
      d <- switch(k,
        param_d1(pp$cor, ec), param_d2(pp$cor, ec),
        param_d3(pp$cor, ec), param_d4(pp$cor, ec)
      )
      cor_cache[[key]] <<- stats::setNames(
        d, vapply(param_tuple_indices(pp$cor, k),
                  function(u) paste(sort(u), collapse = ","), character(1))
      )
    }
    cor_cache[[key]][[paste(sort(t), collapse = ",")]]
  }
  out <- lapply(param_tuple_indices(s, order), function(t) {
    sdi <- t[t <= p]
    cri <- t[t > p] - p
    unname(dr_scale_factor(sd, sdi) * cor_component(cri))
  })
  stats::setNames(out, param_tuple_names(s, order))
}

#' Assemble a Scales-Times-Correlation Log-Determinant Derivative
#'
#' @description
#' The log-determinant is \eqn{2\sum_j \log d_j + \log\lvert R\rvert}, a sum of
#' one term per scale and one term for the correlation. A component is therefore
#' the scale's own where every index names **one** scale coordinate, the
#' correlation's own where every index is a correlation coordinate, and 0
#' otherwise.
#'
#' @details
#' Two separabilities at once, and it is worth keeping them apart: the scales
#' separate from each other, so a component naming two different scales is zero,
#' and the correlation separates from the scales, so a mixed component is zero.
#' What survives is \eqn{2\,\partial^m \log d_k} on the diagonal of the scale
#' block, through [diag_dlog()], and whatever the correlation family answers.
#'
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
#'   as `param_tuple_names(s, order)` and in that order.
#'
#' @seealso [diag_dlog()] for the scale terms, and
#'   [param_dlogdet.DrProdParam()], which calls this.
#'
#' @keywords internal
dr_prod_logdet_derivs <- function(s, eta, order) {
  pp <- .dr(s)
  p <- pp$p
  ec <- .dr_eta_cor(s, eta)
  cor_v <- switch(order,
    param_dlogdet(pp$cor, ec), param_d2logdet(pp$cor, ec),
    param_d3logdet(pp$cor, ec), param_d4logdet(pp$cor, ec)
  )
  cor_v <- stats::setNames(
    cor_v, vapply(param_tuple_indices(pp$cor, order),
                  function(u) paste(sort(u), collapse = ","), character(1))
  )
  es <- eta[seq_len(p)]
  out <- vapply(param_tuple_indices(s, order), function(t) {
    if (all(t <= p)) {
      k <- unique(t)
      if (length(k) > 1L) return(0)
      return(2 * diag_dlog(pp$link, es[k], length(t)))
    }
    if (all(t > p)) {
      return(cor_v[[paste(sort(t - p), collapse = ",")]])
    }
    0
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, order))
}


#' @title Value of a Scales-Times-Correlation Parameter
#' @name param_value.DrProdParam
#' @description
#' Forms \eqn{D R D} as `outer(d, d) * R`, the elementwise product being what two
#' diagonal multiplications amount to. The correlation block is asked for its own
#' value at its own stretch of the free vector, so the composite is exactly what
#' that family returns, rescaled.
#'
#' The value carries **no dimnames**, where the primitive families label their
#' margins `v1`, `v2`, ...; the four compositions share that.
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A symmetric positive definite `s@dimension` by `s@dimension` numeric
#'   matrix, without dimnames.
#' @seealso [param_free.DrProdParam()] for the inverse, and [dr_prod()] for the
#'   parametrization.
#' @keywords internal
S7::method(param_value, DrProdParam) <- function(s, eta, ...) {
  d <- .dr_scales(s, eta)
  unname(outer(d, d) * param_value(.dr(s)$cor, .dr_eta_cor(s, eta)))
}

#' @title Free Vector of a Scales-Times-Correlation Parameter
#' @name param_free.DrProdParam
#' @description
#' Reads the standard deviations off the diagonal as \eqn{\sqrt{\Sigma_{jj}}},
#' divides them out, and hands the resulting correlation matrix to the
#' correlation block. Exact wherever the block is: measured at \eqn{p = 3} with
#' the default block, the round trip closes to \eqn{7 \times 10^{-16}}.
#' @details
#' A matrix with a non-positive diagonal entry is rejected with
#' `'m' must have a positive diagonal.`; anything else outside the family is
#' reported by the correlation block's own message, which is the more specific
#' of the two.
#'
#' This is where a `correlation` block that carries a scale of its own shows: the
#' division always leaves a unit diagonal, so such a block is handed a matrix its
#' own scale cannot be recovered from, and the round trip fails without any
#' error being signaled. See [dr_prod()] on that requirement.
#' @param s A [DrProdParam()] object.
#' @param m A symmetric positive definite `s@dimension` by `s@dimension` matrix,
#'   already checked for shape and symmetry by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`: the linked standard deviations
#'   followed by the correlation block's own free vector.
#' @seealso [param_value.DrProdParam()], the map this inverts.
#' @keywords internal
S7::method(param_free, DrProdParam) <- function(s, m, ...) {
  dg <- diag(m)
  if (any(dg <= 0)) {
    stop("'m' must have a positive diagonal.", call. = FALSE)
  }
  d <- sqrt(dg)
  r <- m / outer(d, d)
  c(linkfunctions7::linkfun(.dr(s)$link, d),
    param_free(.dr(s)$cor, r))
}

#' @title Derivatives of a Scales-Times-Correlation Parameter
#' @name param_d1.DrProdParam
#' @description
#' `param_d1()`, `param_d2()`, `param_d3()` and `param_d4()` for a [dr_prod()]
#' parameter. Each component is the scale factor times the correlation's own
#' component, elementwise, the two groups of free values being disjoint, so
#' nothing of the correlation family is rederived and nothing is differenced.
#' @details
#' The four share [dr_prod_derivs()] and differ only in the order they pass. Both
#' factors are sparse, and where their supports miss each other the component is
#' exactly zero: measured at \eqn{p = 3} that is 1 of the 21 second-order
#' components, 10 of 56 at third order and 37 of 126 at fourth. See [dr_prod()]
#' for the rule.
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return At order 1, a list of `s@n_free` symmetric matrices named by
#'   `s@free_names`; above it, `choose(s@n_free + k - 1, k)` of them keyed as
#'   `param_tuple_names(s, k)` and in that order. Each is `s@dimension` by
#'   `s@dimension` and carries no dimnames.
#' @seealso [dr_prod_derivs()], which assembles them, and
#'   [param_dlogdet.DrProdParam()] for the log-determinant's own.
#' @keywords internal
S7::method(param_d1, DrProdParam) <- function(s, eta, ...) {
  stats::setNames(dr_prod_derivs(s, eta, 1L), s@free_names)
}

#' @rdname param_d1.DrProdParam
#' @name param_d2.DrProdParam
#' @keywords internal
S7::method(param_d2, DrProdParam) <- function(s, eta, ...) {
  dr_prod_derivs(s, eta, 2L)
}

#' @rdname param_d1.DrProdParam
#' @name param_d3.DrProdParam
#' @keywords internal
S7::method(param_d3, DrProdParam) <- function(s, eta, ...) {
  dr_prod_derivs(s, eta, 3L)
}

#' @rdname param_d1.DrProdParam
#' @name param_d4.DrProdParam
#' @keywords internal
S7::method(param_d4, DrProdParam) <- function(s, eta, ...) {
  dr_prod_derivs(s, eta, 4L)
}

#' @title Log-Determinant of a Scales-Times-Correlation Parameter
#' @name param_logdet.DrProdParam
#' @description
#' Closed form, \eqn{2\sum_j \log d_j + \log\lvert R \rvert}, the scales
#' contributing twice because they multiply on both sides. The correlation's term
#' comes from that family by whatever route it has, so a closed form there stays a
#' closed form here and no determinant of the assembled matrix is taken.
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.DrProdParam()] for its derivatives, and
#'   [param_logdet.CorrelationParam()] for the term the default block supplies.
#' @keywords internal
S7::method(param_logdet, DrProdParam) <- function(s, eta, ...) {
  2 * sum(log(.dr_scales(s, eta))) +
    param_logdet(.dr(s)$cor, .dr_eta_cor(s, eta))
}

#' @title Log-Determinant Derivatives of a Scales-Times-Correlation Parameter
#' @name param_dlogdet.DrProdParam
#' @description
#' `param_dlogdet()`, `param_d2logdet()`, `param_d3logdet()` and
#' `param_d4logdet()` for a [dr_prod()] parameter. The log-determinant is
#' separable in the scales and separable from the correlation, so a component
#' mixing two scales, or a scale with a correlation, is exactly zero: measured at
#' \eqn{p = 3}, 12 of the 21 second-order components, 43 of 56 at third order and
#' 108 of 126 at fourth.
#' @details
#' The four share [dr_prod_logdet_derivs()] and differ only in the order they
#' pass. At first order the scale entries are 2 whatever the point, the scales
#' entering \eqn{\log\lvert\Sigma\rvert} through \eqn{2\log d_j} and the log link
#' canceling its own derivative.
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector: at order 1, `s@n_free` values named by
#'   `s@free_names`; above it, `choose(s@n_free + k - 1, k)` values keyed as
#'   `param_tuple_names(s, k)` and in that order.
#' @seealso [dr_prod_logdet_derivs()], which assembles them, and
#'   [param_logdet.DrProdParam()] for the quantity differentiated.
#' @keywords internal
S7::method(param_dlogdet, DrProdParam) <- function(s, eta, ...) {
  stats::setNames(dr_prod_logdet_derivs(s, eta, 1L), s@free_names)
}

#' @rdname param_dlogdet.DrProdParam
#' @name param_d2logdet.DrProdParam
#' @keywords internal
S7::method(param_d2logdet, DrProdParam) <- function(s, eta, ...) {
  dr_prod_logdet_derivs(s, eta, 2L)
}

#' @rdname param_dlogdet.DrProdParam
#' @name param_d3logdet.DrProdParam
#' @keywords internal
S7::method(param_d3logdet, DrProdParam) <- function(s, eta, ...) {
  dr_prod_logdet_derivs(s, eta, 3L)
}

#' @rdname param_dlogdet.DrProdParam
#' @name param_d4logdet.DrProdParam
#' @keywords internal
S7::method(param_d4logdet, DrProdParam) <- function(s, eta, ...) {
  dr_prod_logdet_derivs(s, eta, 4L)
}

#' @title Solve and Factor of a Scales-Times-Correlation Parameter
#' @name param_solve.DrProdParam
#' @description
#' `param_solve()` and `param_factor()` for a [dr_prod()] parameter. Both come
#' from the correlation block with a scaling on either side:
#' \eqn{\Sigma^{-1} = D^{-1} R^{-1} D^{-1}}, and
#' \eqn{\Sigma = (DL)(DL)^\top} for \eqn{L} the correlation's own factor. Neither
#' forms \eqn{\Sigma} and neither factorizes it.
#' @details
#' The scaling is done by dividing and multiplying rows, which for a diagonal
#' matrix is what the products amount to. Measured at \eqn{p = 3} with the
#' default block, the solve agrees with `solve()` on the assembled matrix to
#' \eqn{7 \times 10^{-15}} and `tcrossprod(param_factor(s, eta))` with the matrix
#' to \eqn{4 \times 10^{-16}}. The factor is lower triangular, as
#' [param_factor()] requires, because scaling the rows of a lower triangular
#' matrix leaves it lower triangular.
#'
#' Neither is ever refused for rank: this family admits no deficient correlation
#' block, so it is always of full rank.
#' @param s A [DrProdParam()] object.
#' @param eta A numeric vector of length `s@n_free`, already checked by the
#'   generic.
#' @param b A numeric matrix with `s@dimension` rows, defaulted to the identity
#'   by the generic. `param_factor()` takes no `b`.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return `param_solve()` returns a numeric matrix with `s@dimension` rows and
#'   as many columns as `b`; `param_factor()` a lower triangular `s@dimension` by
#'   `s@dimension` matrix.
#' @seealso [param_solve()] and [param_factor()] for the two contracts.
#' @keywords internal
S7::method(param_solve, DrProdParam) <- function(s, eta, b = NULL, ...) {
  d <- .dr_scales(s, eta)
  inner <- param_solve(.dr(s)$cor, .dr_eta_cor(s, eta), b / d)
  unname(inner / d)
}

#' @rdname param_solve.DrProdParam
#' @name param_factor.DrProdParam
#' @keywords internal
S7::method(param_factor, DrProdParam) <- function(s, eta, ...) {
  d <- .dr_scales(s, eta)
  unname(d * param_factor(.dr(s)$cor, .dr_eta_cor(s, eta)))
}
