#' @include parameter_class.R
NULL

# The two pieces of arithmetic every family in this package ends up needing:
# composing a scalar map with another one, and assembling a Gram product from
# the derivatives of its factor. Both are written once here, because a family
# that transcribes either of them can get a coefficient wrong in a way only a
# fourth-order check would notice.


#' Compose Two Scalar Maps, One Order
#'
#' @description
#' The derivative of order `k` of \eqn{f(g(x))}, from the derivatives of
#' \eqn{f} at \eqn{g(x)} and of \eqn{g} at \eqn{x} up to that order. Each order
#' is its own formula and forms nothing above it; a caller that needs several
#' orders calls it once per order.
#'
#' @details
#' Faa di Bruno's formula, whose coefficients are the numbers of set
#' partitions of a given shape:
#' \deqn{(f \circ g)' = f' g',}
#' \deqn{(f \circ g)'' = f'' g'^2 + f' g'',}
#' \deqn{(f \circ g)''' = f''' g'^3 + 3 f'' g' g'' + f' g''',}
#' \deqn{(f \circ g)'''' = f'''' g'^4 + 6 f''' g'^2 g'' + 3 f'' g''^2
#'   + 4 f'' g' g''' + f' g''''.}
#' The four coefficients of the last line count the partitions of four
#' elements into four singletons, a pair and two singletons, two pairs, a
#' triple and a singleton, and one block.
#'
#' Every argument may be a vector, in which case the composition is applied
#' elementwise and the result has the same length. That is how the families use
#' it: one call composes a whole table of angles or lags at once.
#'
#' @param fd A list or vector of the derivatives of the outer map at
#'   \eqn{g(x)}, in order, at least `k` of them. Read with `[[`, so a list, a
#'   numeric vector or a list of numeric vectors all serve. The **value** of
#'   \eqn{f} is not needed and is not read.
#' @param gd A list or vector of the derivatives of the inner map at \eqn{x},
#'   in order, at least `k` of them, read the same way.
#' @param k The order, an integer from 1 to 4.
#'
#' @return The derivative of order `k` of the composition, the shape the
#'   arithmetic on `fd` and `gd` produces: a single number where both are
#'   scalar, a vector where either is, a matrix where either is.
#'
#' @seealso [leibniz_gram()] for the other piece of shared arithmetic,
#'   [power_derivs()] for the commonest outer map here, and
#'   [numericals7::set_partitions()], which enumerates the partitions the
#'   coefficients count.
#'
#' @keywords internal
compose_order <- function(fd, gd, k) {
  switch(k,
    fd[[1L]] * gd[[1L]],
    fd[[2L]] * gd[[1L]]^2 + fd[[1L]] * gd[[2L]],
    fd[[3L]] * gd[[1L]]^3 + 3 * fd[[2L]] * gd[[1L]] * gd[[2L]] +
      fd[[1L]] * gd[[3L]],
    fd[[4L]] * gd[[1L]]^4 + 6 * fd[[3L]] * gd[[1L]]^2 * gd[[2L]] +
      3 * fd[[2L]] * gd[[2L]]^2 + 4 * fd[[2L]] * gd[[1L]] * gd[[3L]] +
      fd[[1L]] * gd[[4L]]
  )
}


#' Derivatives of a Power, for Composition
#'
#' @description
#' Returns the derivatives of orders 1 to `order` of \eqn{r \mapsto r^{m}} at
#' \eqn{r}, ready to be the outer map of a [compose_order()] call. They are
#' \eqn{m(m-1)\cdots(m-k+1)\,r^{m-k}} and vanish beyond order \eqn{m}, the power
#' being a polynomial. `ar1_pattern()` calls it once per distinct lag.
#'
#' @param r The point, a single number or a vector.
#' @param m The exponent, a non-negative integer. At `m = 0` all are 0, the
#'   power being the constant 1; at `m = 2` the first four are `c(2r, 2, 0, 0)`.
#' @param order The highest order wanted, an integer from 0 to 4.
#'
#' @return A list of `order` elements, the first to the `order`-th
#'   derivative, each the shape of `r`.
#'
#' @seealso [compose_order()], which chains these onto a link's derivatives,
#'   and `ar1_pattern()`, the caller.
#'
#' @keywords internal
power_derivs <- function(r, m, order) {
  lapply(seq_len(order), function(k) {
    if (k > m) return(0)
    prod(m - seq_len(k) + 1L) * r^(m - k)
  })
}


#' A Link Inverse and Its Derivatives, to a Given Order
#'
#' @description
#' Returns \eqn{g^{-1}(\eta)} and its derivatives of orders 1 to `order`, the
#' seed every family composes its map from. Each order is one call of the
#' link's own generic (`dlinkinv()` to `d4linkinv()`), and no order above
#' `order` is evaluated.
#'
#' @param link A \pkg{linkfunctions7} link.
#' @param e A numeric vector of free values.
#' @param order The highest order wanted, an integer from 0 to 4.
#'
#' @return A list of `order + 1` elements: the value, then the derivatives in
#'   order, each the shape of `e`.
#'
#' @seealso [compose_order()], which composes onto these.
#'
#' @keywords internal
linkinv_upto <- function(link, e, order) {
  out <- list(linkfunctions7::linkinv(link, e))
  if (order >= 1L) out[[2L]] <- linkfunctions7::dlinkinv(link, e)
  if (order >= 2L) out[[3L]] <- linkfunctions7::d2linkinv(link, e)
  if (order >= 3L) out[[4L]] <- linkfunctions7::d3linkinv(link, e)
  if (order >= 4L) out[[5L]] <- linkfunctions7::d4linkinv(link, e)
  out
}


#' A Gram Product's Derivatives From Its Factor's
#'
#' @description
#' The derivative of \eqn{M = L L^\top} for one index tuple, from the
#' derivatives of \eqn{L}.
#'
#' @details
#' The Leibniz rule distributes the differentiations of a product over its
#' two factors in every way, so
#' \deqn{\partial^T (L L^\top) = \sum_{S \subseteq T}
#'   (\partial^S L)(\partial^{T \setminus S} L)^\top,}
#' the sum running over subsets of *positions* in the tuple, which
#' handles a repeated index correctly without a multiplicity bookkeeping of
#' its own.
#'
#' The result is symmetric by construction, with no symmetrizing step: the term
#' for a subset \eqn{S} and the term for its complement are transposes of each
#' other, so the sum pairs off. Measured at orders 1 to 3 on a random factor, the
#' asymmetry is exactly 0.
#'
#' The loop walks the \eqn{2^{|T|}} subsets through a bit mask, so the cost is
#' \eqn{2^{\text{order}}} matrix products at worst, and far fewer in practice: a
#' term whose `dfactor` returns `NULL` is skipped before the product is formed,
#' and for both families that use this most terms do.
#'
#' @param dfactor A function of a (possibly empty, possibly repeating) integer
#'   vector of free-value indices, returning the corresponding derivative of
#'   \eqn{L}, or `NULL` where that derivative is identically zero. The empty
#'   vector must give \eqn{L} itself.
#' @param tuple The index tuple, an integer vector whose length is the derivative
#'   order. Repeated entries are handled correctly, the sum running over subsets
#'   of positions.
#' @param p The side of the matrix.
#'
#' @return A symmetric `p` by `p` numeric matrix, with no dimnames.
#'
#' @seealso [compose_order()] for the other piece of shared arithmetic,
#'   [chol_dfactor()] and [corr_dfactor()] for the two `dfactor` arguments the
#'   package supplies, and [chol_leibniz()], the compiled route that replaces
#'   this for the log-Cholesky family.
#'
#' @keywords internal
leibniz_gram <- function(dfactor, tuple, p) {
  order <- length(tuple)
  positions <- seq_len(order)
  acc <- matrix(0, p, p)
  for (b in 0:(2^order - 1L)) {
    take <- as.logical(bitwAnd(b, bitwShiftL(1L, positions - 1L)) > 0L)
    a <- dfactor(tuple[take])
    if (is.null(a)) next
    d <- dfactor(tuple[!take])
    if (is.null(d)) next
    acc <- acc + a %*% t(d)
  }
  acc
}
