#' @include numerical_fallbacks.R
NULL


#' Simplex Parameter
#'
#' @description
#' The S7 class of probability vectors on the open simplex, in the additive
#' log-ratio parametrization. It inherits [parameter()] **directly**, never
#' [matrix_parameter()], its value being a vector: there is no `dimension`, no
#' `rank`, no `null_basis` and no `role`, and [param_logdet()], [param_solve()]
#' and [param_factor()] have no method for it, so asking for the
#' log-determinant of a probability vector fails at dispatch.
#'
#' [simplex()] builds one. The four derivative orders are all closed form.
#'
#' @inheritParams parameter
#'
#' @return An object of class `SimplexParam`, a subclass of [parameter()] adding
#'   no properties of its own: `param_name` is `"simplex"`, `n_free` is
#'   \eqn{K - 1}, `free_names` is `alr1` ... `alr(K-1)`, and `param_params` holds
#'   `n_cat`.
#'
#' @seealso [simplex()], the constructor, [transition_matrix()], which is one of
#'   these per row, and [parameter()] for the properties this inherits and the
#'   four generics it does not get.
#'
#' @examples
#' # A probability vector is not a matrix, so it inherits parameter() alone.
#' s <- simplex(3)
#' c(parameter = S7::S7_inherits(s, parameter),
#'   matrix_parameter = S7::S7_inherits(s, matrix_parameter))
#'
#' # K - 1 free values for K categories.
#' vapply(2:5, function(k) simplex(k)@n_free, integer(1))
#'
#' # And no log-determinant to ask for.
#' try(param_logdet(s, c(0, 0)))
#'
#' @export
SimplexParam <- S7::new_class("SimplexParam", parent = parameter)


#' Construct a Simplex Parameter
#'
#' @description
#' Returns an object holding the additive log-ratio map onto the open simplex: a
#' probability vector of length \eqn{K} from \eqn{K - 1} unconstrained free
#' values, with the last category as reference,
#'
#' \deqn{\pi_a = \frac{e^{\eta_a}}{1 + \sum_b e^{\eta_b}}, \qquad
#'       \pi_K = \frac{1}{1 + \sum_b e^{\eta_b}}.}
#'
#' Every free vector gives positive entries summing to exactly 1, so a mixture
#' weight, a categorical probability or a latent-state distribution can be
#' estimated without a constraint. This is the first family here whose value is
#' **not a matrix**; see [SimplexParam()] for what that costs.
#'
#' @details
#' # The derivatives close over the value
#'
#' With \eqn{Z = 1 + \sum_b e^{\eta_b}}, each \eqn{\pi_a} for \eqn{a < K} is a
#' derivative of \eqn{\log Z}, so the derivative tensors of \eqn{\pi} are the
#' cumulants of a categorical indicator and every order follows from one rule
#' applied repeatedly:
#'
#' \deqn{\partial_b \pi_a = \pi_a(\delta_{ab} - \pi_b).}
#'
#' At first order that is the covariance matrix of a categorical indicator, which
#' is why the tensors are the same objects a multinomial score already carries.
#' All four orders are closed form and no stencil is used.
#'
#' # An identity worth checking against
#'
#' \eqn{\sum_a \pi_a = 1} at every \eqn{\eta}, so differentiating it gives
#' \eqn{\sum_a \partial \pi_a = 0}, and the same at every higher order. Every
#' derivative component therefore sums to zero over the value index. Measured
#' over all four orders the worst sum is \eqn{1.7 \times 10^{-17}}, and
#' [check_parameter()] runs exactly this check: a derivative array that does not
#' sum to zero is wrong whatever else it agrees with.
#'
#' # Large free values saturate instead of overflowing
#'
#' [simplex_point()] applies the log-sum-exp shift, so a free value of 800 gives
#' \eqn{\pi_1 = 1} and \eqn{\pi_K = 0} with the vector still summing to exactly 1,
#' where the naive expression would divide `Inf` by `Inf`. Note that the value
#' then sits on the **boundary** of the simplex, so [param_free()] cannot invert
#' it: \eqn{\log(0)} is not finite.
#'
#' # Why not stick-breaking
#'
#' The other common chart chains through \eqn{K - 1} nested logistic maps, so its
#' derivatives compose that many times and are not symmetric in the categories.
#' The additive log-ratio's close over \eqn{\pi} in one rule, and every category
#' but the reference enters the same way.
#'
#' @section Notation:
#' \eqn{K} is the number of categories, \eqn{\eta \in \mathbb{R}^{K-1}} the free
#' vector, \eqn{\pi} the probability vector and \eqn{\delta_{ab}} the Kronecker
#' delta. The reference category is the last, \eqn{K}.
#'
#' @param n_cat The number of categories \eqn{K}, **at least 2**. A single
#'   integer; `1`, a fraction, `NA` and a vector all throw `'n_cat' must be a
#'   single integer of at least 2.` A one-category simplex is the constant 1 and
#'   has nothing to estimate.
#'
#' @return An object of class [SimplexParam()], with `n_free` equal to
#'   `n_cat - 1`, `free_names` `alr1` ... `alr(K-1)`, `param_name` `"simplex"`,
#'   and `param_params` holding `n_cat`. It carries no `dimension`, `rank`,
#'   `null_basis` or `role`, being a [parameter()] and never a
#'   [matrix_parameter()].
#'
#' @references
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
#' Chapman and Hall, London. The additive log-ratio chart is chapter 6.
#'
#' @seealso [transition_matrix()], which is one of these per row,
#'   [param_value()] and [param_free()] for the map and its inverse, and
#'   [check_parameter()], whose seven-check battery is what a non-matrix family
#'   gets.
#'
#' @examples
#' # Four categories, three free values.
#' s <- simplex(4)
#' s@free_names
#' eta <- c(0.5, -0.2, 1.1)
#' pi <- param_value(s, eta)
#' round(pi, 5)
#' sum(pi)
#'
#' # It is the closed form, with the last category as reference.
#' all.equal(unname(pi),
#'           c(exp(eta), 1) / (1 + sum(exp(eta))))
#'
#' # The round trip closes exactly.
#' max(abs(param_free(s, pi) - eta))
#'
#' # The first derivative is the covariance of a categorical indicator.
#' d1 <- param_d1(s, eta)
#' round(d1[["alr1"]], 6)
#' round(pi * ((seq_along(pi) == 1) - pi[1]), 6)
#'
#' # Every derivative component sums to zero over the value index, at every
#' # order, because sum(pi) is the constant 1.
#' vapply(1:4, function(k) {
#'   dk <- do.call(paste0("param_d", k), list(s, eta))
#'   max(abs(vapply(dk, sum, numeric(1))))
#' }, numeric(1))
#'
#' # A large free value saturates instead of overflowing, and the vector still
#' # sums to exactly 1.
#' v <- param_value(s, c(800, 0, 0))
#' c(p1 = v[[1]], pK = v[[4]], total = sum(v))
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
#' Computes the value \eqn{\pi} from a free vector, appending the reference
#' category's implicit 0 and taking the softmax through the **log-sum-exp
#' shift**: the largest entry is subtracted before exponentiating, so a large
#' free value saturates instead of overflowing.
#'
#' @details
#' The naive expression \eqn{e^{\eta_a}/(1 + \sum_b e^{\eta_b})} is `Inf/Inf` from
#' about \eqn{\eta = 710}. With the shift, a free value of 800 returns
#' \eqn{\pi_1 = 1} and \eqn{\pi_K = 0} with the vector summing to exactly 1, and
#' at 500 it returns \eqn{\pi_K = 7 \times 10^{-218}}, still representable. The
#' value is then on the boundary of the simplex, which [param_free()] cannot
#' invert.
#'
#' @param eta A numeric vector of length \eqn{K - 1}. Not checked; the callers
#'   have been through the generic.
#'
#' @return A numeric vector of length \eqn{K} with non-negative entries summing to
#'   1, and no names. The names `p1` ... `pK` are applied by
#'   [param_value.SimplexParam()].
#'
#' @seealso [param_value.SimplexParam()] and
#'   [param_value.TransitionMatrixParam()], the two callers, and
#'   [simplex_tensors()] for the derivatives of the same map.
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
#' @details
#' Each array is one more application of the product rule to the one below, so the
#' cost of the fourth is \eqn{K(K-1)^4} entries and the whole set is built in one
#' pass. Nothing is differenced.
#'
#' Every slice sums to zero over its first index, \eqn{\sum_a \pi_a} being the
#' constant 1, which is the identity [check_parameter()] tests.
#'
#' @param pi_full The value, a vector of length \eqn{K} summing to 1, as
#'   [simplex_point()] returns it.
#' @param order The highest order wanted: 1, 2, 3 or 4.
#'
#' @return A list with elements `d1` to `d<order>`: `d1` is `K` by `K-1`, `d2` is
#'   `K` by `K-1` by `K-1`, and so on, the first index running over the
#'   categories of the value and the rest over the free values.
#'
#' @seealso [simplex_components()], which slices these into the named lists the
#'   generics return, and [simplex_point()] for the value they are built from.
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
#' Slices the tensors of [simplex_tensors()] into the named list a derivative
#' generic returns, one entry per distinct index tuple, keyed as
#' `param_tuple_names(s, order)`.
#'
#' @details
#' The `wrap` argument is how [transition_matrix()] reuses this. A simplex
#' component is a vector of length \eqn{K} and is returned as it stands; a
#' transition matrix's is that vector placed into one row of a \eqn{K \times K}
#' matrix of zeros, and the caller passes the function that does the placing.
#'
#' @param s The parameter the tuples belong to, a [SimplexParam()] or, through
#'   `wrap`, one row's worth of a [TransitionMatrixParam()].
#' @param tens The tensor list from [simplex_tensors()].
#' @param order The order to extract: 1, 2, 3 or 4. `tens` must carry at least
#'   this order.
#' @param wrap A function applied to each raw slice before it is stored.
#'   `identity` by default, which is the simplex's own case;
#'   [transition_matrix()] passes a function embedding the slice in a row.
#'
#' @return A named list of `choose(s@n_free + order - 1, order)` entries, keyed
#'   as `param_tuple_names(s, order)` and in that order, each entry `wrap()`'s
#'   result.
#'
#' @seealso [simplex_tensors()] for the arrays sliced, and [tm_derivative()] for
#'   the `wrap` that embeds a row.
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
#' @description
#' Returns the probability vector \eqn{\pi}, the softmax of the free vector with
#' the reference category's implicit 0 appended, named `p1` ... `pK`. The entries
#' are positive and sum to exactly 1 at every free vector, so nothing is tested
#' and nothing is renormalized. Computed through [simplex_point()]'s log-sum-exp
#' shift, so a large free value saturates instead of overflowing.
#' @param s A [SimplexParam()] object, whose `param_params$n_cat` is read.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length \eqn{K} named `p1` ... `pK`, with
#'   non-negative entries summing to 1.
#' @seealso [param_free.SimplexParam()] for the inverse, and [simplex_point()]
#'   for the arithmetic.
#' @keywords internal
S7::method(param_value, SimplexParam) <- function(s, eta, ...) {
  stats::setNames(simplex_point(eta), paste0("p", seq_len(s@param_params$n_cat)))
}


#' @title Free Vector of a Simplex Parameter
#' @name param_free.SimplexParam
#' @description
#' Returns the additive log-ratio, \eqn{\eta_a = \log(\pi_a/\pi_K)}, exact and a
#' true inverse of [param_value.SimplexParam()]: the round trip closes to
#' \eqn{2 \times 10^{-16}}.
#' @details
#' Two rejections, both with their own message. A vector with a non-positive
#' entry is outside the **open** simplex, and \eqn{\log 0} is not finite; a
#' vector that does not sum to 1 is not a probability vector, and it is **not
#' renormalized**, a silent repair being the kind of thing that hides a caller's
#' defect for a long time.
#'
#' The first rejection is the one a fit runs into. [simplex_point()] saturates a
#' large free value to an exact 0 in the reference category, so a value produced
#' at the boundary cannot be inverted back.
#' @param s A [SimplexParam()] object.
#' @param m A probability vector of length \eqn{K}: strictly positive and summing
#'   to 1.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.SimplexParam()], the map this inverts.
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
#' Closed form:
#'
#' \deqn{\partial_b \pi_a = \pi_a(\delta_{ab} - \pi_b),}
#'
#' which is the covariance structure of a categorical indicator: the same array a
#' multinomial score already carries. Every component sums to zero over the
#' category index, \eqn{\sum_a \pi_a} being the constant 1.
#' @param s A [SimplexParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` numeric vectors of length \eqn{K}, named by
#'   `s@free_names`, each summing to 0.
#' @seealso [simplex_tensors()], which builds the array, and
#'   [param_d2.SimplexParam()] for the order above.
#' @keywords internal
S7::method(param_d1, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 1L)
  out <- lapply(seq_len(s@n_free), function(b) tens$d1[, b])
  stats::setNames(out, s@free_names)
}


#' @title Second Derivatives of a Simplex Parameter
#' @name param_d2.SimplexParam
#' @description
#' Closed form, one further application of the product rule to
#' \eqn{\partial_b \pi_a = \pi_a(\delta_{ab} - \pi_b)}. The array is the
#' third cumulant of a categorical indicator, and every component sums to zero
#' over the category index.
#'
#' Unlike a diagonal family, nothing here is separable: the entries of \eqn{\pi}
#' are coupled through the shared normalizing constant, so a component in two
#' different free values is generally non-zero.
#' @param s A [SimplexParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` numeric vectors of length \eqn{K},
#'   keyed as `param_tuple_names(s)` and in that order, each summing to 0.
#' @seealso [param_d1.SimplexParam()] and [param_d3.SimplexParam()] for the
#'   neighboring orders.
#' @keywords internal
S7::method(param_d2, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 2L)
  simplex_components(s, tens, 2L)
}


#' @title Third Derivatives of a Simplex Parameter
#' @name param_d3.SimplexParam
#' @description
#' Closed form, the same product rule applied a third time, which is the cumulant
#' recursion of a categorical indicator at third order. Exact, where a family
#' without a closed form would get a product stencil good to about six digits.
#' Every component sums to zero over the category index.
#' @param s A [SimplexParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` numeric vectors of length \eqn{K},
#'   keyed as `param_tuple_names(s, 3)` and in that order, each summing to 0.
#' @seealso [simplex_tensors()], which builds the array, and
#'   [param_d4.SimplexParam()] for the order above.
#' @keywords internal
S7::method(param_d3, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 3L)
  simplex_components(s, tens, 3L)
}


#' @title Fourth Derivatives of a Simplex Parameter
#' @name param_d4.SimplexParam
#' @description
#' Closed form, the cumulant recursion at fourth order, which is where the
#' contract stops. Exact, and that matters most here: a product stencil at fourth
#' order keeps about five digits, and this array is what a fourth-order chain rule
#' through a link reads. Every component sums to zero over the category index,
#' measured at \eqn{1.7 \times 10^{-17}}.
#'
#' The array holds \eqn{K(K-1)^4} entries, so it is the largest object the family
#' builds; at \eqn{K = 4} that is 324.
#' @param s A [SimplexParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` numeric vectors of length \eqn{K},
#'   keyed as `param_tuple_names(s, 4)` and in that order, each summing to 0.
#' @seealso [simplex_tensors()], which builds the array,
#'   [param_d3.SimplexParam()] for the order below, and [numerical_d4()] for the
#'   alternative.
#' @keywords internal
S7::method(param_d4, SimplexParam) <- function(s, eta, ...) {
  tens <- simplex_tensors(simplex_point(eta), 4L)
  simplex_components(s, tens, 4L)
}
