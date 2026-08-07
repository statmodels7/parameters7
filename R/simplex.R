#' @include numerical_fallbacks.R
NULL


#' Simplex Parameter
#'
#' @description
#' The S7 class of probability vectors on the open simplex, in the additive
#' log-ratio parametrization. Constructed by \code{\link{simplex}}.
#'
#' @inheritParams parameter
#'
#' @return An object of class \code{SimplexParam}.
#'
#' @seealso \code{\link{simplex}}
#'
#' @examples
#' S7::S7_inherits(simplex(3), SimplexParam)
#'
#' @export
SimplexParam <- S7::new_class("SimplexParam", parent = parameter)


#' Construct a Simplex Parameter
#'
#' @description
#' A probability vector of length \eqn{K} on the open simplex, carried by the
#' additive log-ratio map with the last category as reference:
#' \deqn{\pi_a = \dfrac{e^{\eta_a}}{1 + \sum_b e^{\eta_b}}, \qquad
#'       \pi_K = \dfrac{1}{1 + \sum_b e^{\eta_b}},}
#' with \eqn{\eta \in \mathbb{R}^{K-1}} unconstrained.
#'
#' @details
#' The derivatives of the map close over the value itself. With
#' \eqn{Z = 1 + \sum_b e^{\eta_b}}, each \eqn{\pi_a} with \eqn{a < K} is a
#' derivative of \eqn{\log Z}, so the derivative tensors of \eqn{\pi} are the
#' cumulants of a categorical indicator and follow from the one rule
#' \eqn{\partial_b \pi_a = \pi_a(\delta_{ab} - \pi_b)} applied repeatedly.
#' All four orders are closed form.
#'
#' The inverse is \eqn{\eta_a = \log(\pi_a / \pi_K)}, exact; a vector outside
#' the open simplex, or one that does not sum to one, is rejected rather than
#' repaired, because a silent renormalization would mask the caller's defect.
#'
#' Stick-breaking is the other common chart and was not chosen: its
#' derivatives chain through \eqn{K - 1} nested logistic maps and are not
#' symmetric in the categories, while the additive log-ratio's close over
#' \eqn{\pi} in one rule.
#'
#' @param n_cat The number of categories \eqn{K}, at least 2.
#'
#' @return An object of class \code{\link{SimplexParam}}.
#'
#' @seealso \code{\link{transition_matrix}}, \code{\link{param_value}}
#'
#' @examples
#' s <- simplex(3)
#' round(param_value(s, c(0.5, -0.2)), 4)
#'
#' # the round trip closes exactly
#' eta <- c(0.5, -0.2)
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' @export
simplex <- function(n_cat) {
  if (!is.numeric(n_cat) || length(n_cat) != 1L || !is.finite(n_cat) ||
    n_cat < 2 || n_cat != round(n_cat)) {
    stop("'n_cat' must be a single integer of at least 2.", call. = FALSE)
  }
  k <- as.integer(n_cat)
  SimplexParam(
    param_name = "simplex",
    n_free = k - 1L,
    free_names = paste0("alr", seq_len(k - 1L)),
    param_params = list(n_cat = k)
  )
}


#' The Softmax Point Behind a Free Vector
#'
#' @description
#' The value \eqn{\pi}, computed through the log-sum-exp shift so that a large
#' free value saturates gracefully instead of overflowing.
#'
#' @param eta A numeric vector of length \eqn{K - 1}.
#'
#' @return A numeric vector of length \eqn{K}.
#'
#' @keywords internal
simplex_point <- function(eta) {
  z <- c(eta, 0)
  z <- z - max(z)
  e <- exp(z)
  e / sum(e)
}


#' Derivative Tensors of the Softmax Map
#'
#' @description
#' The arrays \eqn{D_1[a, b]} to \eqn{D_4[a, b, c, d, e]} of derivatives of
#' \eqn{\pi} in the free values, built by applying the product rule to
#' \eqn{\partial_b \pi_a = \pi_a(\delta_{ab} - \pi_b)} as many times as the
#' order asks. The first index runs over the \eqn{K} categories of the value,
#' the rest over the \eqn{K - 1} free values.
#'
#' @param pi_full The value, a vector of length \eqn{K}.
#' @param order The highest order wanted, 1 to 4.
#'
#' @return A list with elements \code{d1} to \code{d\{order\}}.
#'
#' @keywords internal
simplex_tensors <- function(pi_full, order = 4L) {
  k <- length(pi_full)
  d <- k - 1L
  delta <- diag(1, k, d)

  d1 <- array(0, c(k, d))
  for (a in seq_len(k)) {
    for (b in seq_len(d)) {
      d1[a, b] <- pi_full[a] * (delta[a, b] - pi_full[b])
    }
  }
  out <- list(d1 = d1)
  if (order < 2L) return(out)

  # d2[a,b,c] = d_c d1[a,b] = d1[a,c](delta_ab - pi_b) - pi_a d1[b,c]
  d2 <- array(0, c(k, d, d))
  for (a in seq_len(k)) {
    for (b in seq_len(d)) {
      for (cc in seq_len(d)) {
        d2[a, b, cc] <- d1[a, cc] * (delta[a, b] - pi_full[b]) -
          pi_full[a] * d1[b, cc]
      }
    }
  }
  out$d2 <- d2
  if (order < 3L) return(out)

  # d3[a,b,c,e] = d_e d2[a,b,c]
  #   = d2[a,c,e](delta_ab - pi_b) - d1[a,c] d1[b,e]
  #     - d1[a,e] d1[b,c] - pi_a d2[b,c,e]
  d3 <- array(0, c(k, d, d, d))
  for (a in seq_len(k)) {
    for (b in seq_len(d)) {
      for (cc in seq_len(d)) {
        for (e in seq_len(d)) {
          d3[a, b, cc, e] <- d2[a, cc, e] * (delta[a, b] - pi_full[b]) -
            d1[a, cc] * d1[b, e] -
            d1[a, e] * d1[b, cc] -
            pi_full[a] * d2[b, cc, e]
        }
      }
    }
  }
  out$d3 <- d3
  if (order < 4L) return(out)

  # d4[a,b,c,e,f] = d_f d3[a,b,c,e], one product rule per term of d3
  d4 <- array(0, c(k, d, d, d, d))
  for (a in seq_len(k)) {
    for (b in seq_len(d)) {
      for (cc in seq_len(d)) {
        for (e in seq_len(d)) {
          for (f in seq_len(d)) {
            d4[a, b, cc, e, f] <-
              d3[a, cc, e, f] * (delta[a, b] - pi_full[b]) -
              d2[a, cc, e] * d1[b, f] -
              d2[a, cc, f] * d1[b, e] -
              d1[a, cc] * d2[b, e, f] -
              d2[a, e, f] * d1[b, cc] -
              d1[a, e] * d2[b, cc, f] -
              d1[a, f] * d2[b, cc, e] -
              pi_full[a] * d3[b, cc, e, f]
          }
        }
      }
    }
  }
  out$d4 <- d4
  out
}


#' Extract Named Components From Softmax Tensors
#'
#' @description
#' Slices the tensors of \code{\link{simplex_tensors}} into the named lists
#' the derivative generics return.
#'
#' @param s The parameter the tuples belong to.
#' @param tens The tensor list.
#' @param order The order to extract.
#' @param wrap A function applied to each raw slice, for the transition
#'   matrix to embed a row; the identity here.
#'
#' @return A named list keyed as \code{param_tuple_names(s, order)}.
#'
#' @keywords internal
simplex_components <- function(s, tens, order, wrap = identity) {
  idx <- param_tuple_indices(s, order)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, order)
  arr <- tens[[paste0("d", order)]]
  for (i in seq_along(idx)) {
    t <- idx[[i]]
    slice <- switch(order,
      arr[, t[1L]],
      arr[, t[1L], t[2L]],
      arr[, t[1L], t[2L], t[3L]],
      arr[, t[1L], t[2L], t[3L], t[4L]]
    )
    out[[i]] <- wrap(slice)
  }
  out
}


#' @title Value of a Simplex Parameter
#' @name param_value.SimplexParam
#' @description The softmax point, named \code{p1..pK}.
#' @param s A \code{\link{SimplexParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A numeric vector of length \eqn{K} summing to one.
#' @keywords internal
S7::method(param_value, SimplexParam) <- function(s, eta, ...) {
  stats::setNames(simplex_point(eta), paste0("p", seq_len(s@param_params$n_cat)))
}


#' @title Free Vector of a Simplex Parameter
#' @name param_free.SimplexParam
#' @description
#' The additive log-ratio, exact: \eqn{\eta_a = \log(\pi_a/\pi_K)}. Rejected
#' outside the open simplex or when the vector does not sum to one.
#' @param s A \code{\link{SimplexParam}} object.
#' @param m A probability vector of length \eqn{K}.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(param_free, SimplexParam) <- function(s, m, ...) {
  k <- s@param_params$n_cat
  if (!is.numeric(m) || length(m) != k) {
    stop(sprintf("'m' must be a numeric vector of length %d.", k), call. = FALSE)
  }
  if (anyNA(m) || any(m <= 0)) {
    stop("'m' must be strictly inside the simplex: every entry positive.",
      call. = FALSE
    )
  }
  if (abs(sum(m) - 1) > 1e-8) {
    stop(paste0(
      "'m' does not sum to one, so it is not on the simplex. It is rejected\n",
      "  rather than renormalized, because a silent repair would mask the\n",
      "  caller's defect."
    ), call. = FALSE)
  }
  stats::setNames(log(m[seq_len(k - 1L)] / m[k]), s@free_names)
}


#' @title First Derivatives of a Simplex Parameter
#' @name param_d1.SimplexParam
#' @description
#' Closed form: \eqn{\partial_b \pi_a = \pi_a(\delta_{ab} - \pi_b)}, the
#' covariance structure of a categorical indicator.
#' @param s A \code{\link{SimplexParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of vectors of length \eqn{K}.
#' @keywords internal
S7::method(param_d1, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 1L)
  out <- lapply(seq_len(s@n_free), function(b) tens$d1[, b])
  stats::setNames(out, s@free_names)
}


#' @title Second Derivatives of a Simplex Parameter
#' @name param_d2.SimplexParam
#' @description Closed form, one further product rule.
#' @param s A \code{\link{SimplexParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of vectors of length \eqn{K}.
#' @keywords internal
S7::method(param_d2, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 2L)
  simplex_components(s, tens, 2L)
}


#' @title Third Derivatives of a Simplex Parameter
#' @name param_d3.SimplexParam
#' @description Closed form, the cumulant recursion at third order.
#' @param s A \code{\link{SimplexParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of vectors of length \eqn{K}.
#' @keywords internal
S7::method(param_d3, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 3L)
  simplex_components(s, tens, 3L)
}


#' @title Fourth Derivatives of a Simplex Parameter
#' @name param_d4.SimplexParam
#' @description Closed form, the cumulant recursion at fourth order.
#' @param s A \code{\link{SimplexParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of vectors of length \eqn{K}.
#' @keywords internal
S7::method(param_d4, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 4L)
  simplex_components(s, tens, 4L)
}
