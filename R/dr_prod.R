#' @include numerical_fallbacks.R
NULL

#' Scales Times a Correlation Matrix
#'
#' @description
#' The S7 class of a covariance written as \eqn{D R D} for a diagonal matrix
#' of positive scales and a correlation matrix. Constructed by
#' \code{\link{dr_prod}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{DrProdParam}.
#'
#' @seealso \code{\link{dr_prod}}
#'
#' @examples
#' S7::S7_inherits(dr_prod(3), DrProdParam)
#'
#' @export
DrProdParam <- S7::new_class("DrProdParam", parent = matrix_parameter)


#' Construct a Covariance as Scales Times a Correlation
#'
#' @description
#' \eqn{\Sigma(\eta) = D(\eta_D)\, R(\eta_R)\, D(\eta_D)} with
#' \eqn{D = \mathrm{diag}(d_1, \ldots, d_p)} carrying the standard deviations
#' through a positive link, and \eqn{R} a correlation matrix parameter.
#'
#' @details
#' The separation is what makes the parametrization worth having: the
#' quantities a reader reads off a fitted covariance are the standard
#' deviations and the correlations, and here they are the coordinates rather
#' than a function of them. A log-Cholesky factor produces the same set of
#' matrices with no coordinate meaning anything on its own.
#'
#' Every derivative factorizes, because \eqn{\Sigma_{ij} = d_i d_j R_{ij}} and
#' the two groups of free values are disjoint:
#' \deqn{\partial^{S}\Sigma_{ij} =
#'       \bigl[\partial^{S_D}(d_i d_j)\bigr]\,
#'       \bigl[\partial^{S_R} R_{ij}\bigr],}
#' where \eqn{S_D} and \eqn{S_R} are the parts of the multiset \eqn{S} falling
#' in each group. The scale factor is zero unless every index of \eqn{S_D}
#' names \eqn{i} or \eqn{j}, so most components vanish; when \eqn{i = j} it is
#' the Leibniz expansion of \eqn{\partial^m d_i^2}. Nothing of the correlation
#' family is rederived.
#'
#' The log-determinant is \eqn{2\sum_j \log d_j + \log\lvert R\rvert}, which is
#' separable in the scales and separable from the correlation, so a
#' log-determinant component mixing two scales, or a scale with a correlation,
#' is exactly zero.
#'
#' The correlation block must have full rank. A rank-deficient \eqn{R} would
#' give \eqn{\Sigma} the null space \eqn{D^{-1}\ker R}, which moves with the
#' free vector, and the class records the rank and the null space as
#' properties of the family rather than of a point.
#'
#' @param dimension The side \eqn{p} of the matrix, at least 2.
#' @param correlation A \code{\link{matrix_parameter}} of side
#'   \code{dimension} producing correlation matrices. Defaults to
#'   \code{\link{correlation_matrix}(dimension)}.
#' @param link The positive link carrying each standard deviation onto the
#'   free scale. Defaults to \code{\link[linkfunctions7]{log_link}()}.
#' @param role One of \code{"covariance"}, \code{"precision"} or
#'   \code{"either"}.
#'
#' @return An object of class \code{\link{DrProdParam}}.
#'
#' @seealso \code{\link{correlation_matrix}}, \code{\link{log_cholesky}},
#'   \code{\link{block_diag}}
#'
#' @examples
#' s <- dr_prod(3)
#' s@free_names
#' round(param_value(s, c(0, log(2), 0, 0.4, 0.4, 0.4)), 3)
#'
#' # p standard deviations and the p(p-1)/2 angles of the correlation
#' c(p = 4, n_free = dr_prod(4)@n_free)
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
#' carries it, as a matrix with one row per order.
#'
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#'
#' @return A 5 by \code{p} numeric matrix, rows being orders 0 to 4.
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
#' Evaluates \eqn{\partial^{S_D}(d_i d_j)} for every pair \eqn{(i, j)}, given
#' the multiset \eqn{S_D} of scale indices. It is zero wherever \eqn{S_D}
#' contains an index other than \eqn{i} or \eqn{j}, so the result is supported
#' on the rows and columns those indices name.
#'
#' @param sd A 5 by \code{p} matrix of inverse-link derivatives, as returned by
#'   \code{\link{dr_scale_derivs}}.
#' @param tuple An integer vector of scale indices, possibly empty and possibly
#'   with repeats.
#'
#' @return A \code{p} by \code{p} numeric matrix.
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
#' Multiplies the scale factor by the correlation's own component, elementwise,
#' for every tuple of the composite's enumeration.
#'
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of matrices keyed as \code{param_tuple_names(s, order)}.
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
#' The log-determinant is \eqn{2\sum_j \log d_j + \log\lvert R\rvert}, so a
#' component is the scale's own when every index names one scale coordinate,
#' the correlation's own when every index is a correlation coordinate, and zero
#' otherwise.
#'
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named numeric vector keyed as \code{param_tuple_names(s, order)}.
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
#' @description Forms \eqn{D R D} from the scales and the correlation block.
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A symmetric matrix of side \code{s@dimension}.
#' @keywords internal
S7::method(param_value, DrProdParam) <- function(s, eta, ...) {
  d <- .dr_scales(s, eta)
  unname(outer(d, d) * param_value(.dr(s)$cor, .dr_eta_cor(s, eta)))
}

#' @title Free Vector of a Scales-Times-Correlation Parameter
#' @name param_free.DrProdParam
#' @description
#' Reads the standard deviations off the diagonal, divides them out, and hands
#' the resulting correlation matrix to the correlation block.
#' @param s A \code{\link{DrProdParam}} object.
#' @param m A symmetric positive definite matrix of side \code{s@dimension}.
#' @param ... Unused.
#' @return A numeric vector of length \code{s@n_free}.
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
#' Each component is the scale factor times the correlation's own component,
#' elementwise, the two groups of free values being disjoint.
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
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
#' \eqn{2\sum_j \log d_j + \log\lvert R \rvert}, the scales contributing twice
#' because they multiply on both sides.
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, DrProdParam) <- function(s, eta, ...) {
  2 * sum(log(.dr_scales(s, eta))) +
    param_logdet(.dr(s)$cor, .dr_eta_cor(s, eta))
}

#' @title Log-Determinant Derivatives of a Scales-Times-Correlation Parameter
#' @name param_dlogdet.DrProdParam
#' @description
#' Separable in the scales and separable from the correlation, so a component
#' mixing two scales, or a scale with a correlation, is exactly zero.
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param ... Unused.
#' @return A named numeric vector.
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
#' \eqn{\Sigma^{-1} = D^{-1} R^{-1} D^{-1}} and \eqn{\Sigma = (DL)(DL)'} for
#' \eqn{L} the correlation's factor, so both come from the correlation block
#' with a scaling on either side.
#' @param s A \code{\link{DrProdParam}} object.
#' @param eta A numeric vector of length \code{s@n_free}.
#' @param b A matrix with \code{s@dimension} rows, or \code{NULL}.
#' @param ... Unused.
#' @return A matrix.
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
