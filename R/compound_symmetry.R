#' @include chain.R numerical_fallbacks.R
NULL

# Compound symmetry and AR(1) are the economical families: whatever p is, two
# free values, a scale and a correlation. They share the whole skeleton -- a
# positive scale through one link, a bounded correlation through another, and
# a value that is the scale times a pattern depending on the correlation alone
# -- so the skeleton is written here once and ar1.R adds the second family to
# it, the way transition_matrix.R builds on simplex.R.
#
# The economy is the point rather than a convenience. A random effect over
# twenty occasions is 210 free values unstructured and two here, and the two
# are what a model can actually estimate from a handful of subjects.


#' Compound Symmetry Parameter
#'
#' @description
#' The S7 class of compound-symmetric covariance matrices: equal variances and
#' one common correlation, so **two** free values at every dimension. It is the
#' covariance of an exchangeable set of measurements, and the one that a random
#' intercept induces.
#'
#' [compound_symmetry()] builds one. `param_params` records the two links, and
#' the correlation's is a `bounded_link(-1/(p-1), 1)` whose lower end depends on
#' `dimension`: two objects of different sizes carry different links.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `CompoundSymmetryParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `link_scale` and `link_rho`. `n_free` is 2 at every `dimension`, and `rank`
#'   is \eqn{p}.
#'
#' @seealso [compound_symmetry()], the constructor, [ar1()] for the other
#'   two-value family, and [matrix_parameter()] for the properties that
#'   this class inherits.
#'
#' @examples
#' # Two free values whatever the dimension.
#' s <- compound_symmetry(3)
#' S7::S7_inherits(s, CompoundSymmetryParam)
#' vapply(2:6, function(p) compound_symmetry(p)@n_free, integer(1))
#'
#' # The correlation's link depends on the dimension, its lower bound being
#' # -1/(p-1): below that the matrix would not be positive definite.
#' t(vapply(c(2, 4, 10),
#'          function(p) compound_symmetry(p)@param_params$link_rho@link_bounds,
#'          numeric(2)))
#'
#' @export
CompoundSymmetryParam <- S7::new_class("CompoundSymmetryParam",
  parent = matrix_parameter)


#' Construct a Compound Symmetry Parameter
#'
#' @description
#' Returns an object holding the exchangeable covariance
#'
#' \deqn{M(\eta) = \sigma^2\{(1-\rho)I + \rho J\},}
#'
#' with \eqn{\sigma^2} positive and \eqn{\rho} the correlation that every pair
#' shares. Two free values at every dimension, against the \eqn{p(p+1)/2} of
#' [log_cholesky()], so it is the parametrization to use when the measurements
#' are exchangeable: repeated measures on a subject, items in a block, the
#' covariance that a random intercept induces.
#'
#' @details
#' # Bounds of the correlation
#'
#' The eigenvalues are \eqn{\sigma^2\{1 + (p-1)\rho\}} once and
#' \eqn{\sigma^2(1-\rho)} with multiplicity \eqn{p-1}, so the matrix is positive
#' definite exactly when
#'
#' \deqn{-\frac{1}{p-1} < \rho < 1.}
#'
#' The correlation is therefore mapped by
#' `linkfunctions7::bounded_link(-1/(p-1), 1)`. A `rhobit_link()` onto
#' \eqn{(-1, 1)} would give an indefinite matrix at an ordinary free value: at \eqn{p = 4} a correlation of
#' \eqn{-0.5} is inside \eqn{(-1, 1)} and outside the cone. The bound depends on
#' the dimension, so two objects of different sizes carry different links.
#'
#' A free value of 0 is the **midpoint** of that interval, so at \eqn{p = 4} it
#' gives a correlation of \eqn{1/3}.
#'
#' # The log-determinant and the inverse
#'
#' From the eigenvalues,
#'
#' \deqn{\log|M| = p\log\sigma^2 + \log\{1 + (p-1)\rho\} + (p-1)\log(1-\rho),}
#'
#' a sum of a function of one free value and a function of the other, so every
#' mixed derivative of the log-determinant is exactly zero and no determinant is
#' computed.
#'
#' The inverse is compound symmetric again, by Sherman-Morrison, so
#' [param_solve()] returns it in closed form and factorizes nothing. An
#' exchangeable covariance has an exchangeable precision, which leaves the family
#' closed under the choice of side.
#'
#' # The derivatives
#'
#' The value is the scale times a pattern **linear** in the correlation,
#' \eqn{I + \rho(J - I)}, so a component with \eqn{a} scale indices and \eqn{b}
#' correlation indices is the \eqn{a}-th derivative of the scale times the
#' \eqn{b}-th derivative of the pattern. All four orders follow from the two
#' links' own derivatives with no further algebra: a component with \eqn{b}
#' correlation indices carries the \eqn{b}-th derivative of the correlation's
#' link times the constant matrix \eqn{J - I}.
#'
#' @section Notation:
#' \eqn{\eta = (\eta_1, \eta_2)} is the free vector, \eqn{\sigma^2} the common
#' variance and \eqn{\rho} the common correlation. \eqn{p} is the side of the
#' matrix, \eqn{I} the identity and \eqn{J} the matrix of ones.
#'
#' @param dimension The side \eqn{p} of the matrix, **at least 2**. A one by one
#'   matrix has no correlation, so `compound_symmetry(1)` signals an error
#'   stating that the family would carry a free value with no effect.
#' @param link_scale A \pkg{linkfunctions7} link carrying the first free value
#'   onto the positive variance, `linkfunctions7::log_link()` by default. It must
#'   map onto the positive half line, so `identity_link()` is rejected, and from
#'   the whole real line, which rules out `sqrt_link()` and its relatives; see
#'   [diagonal_matrix()] for the two conditions.
#'
#' @return An object of class [CompoundSymmetryParam()], with `n_free` 2,
#'   `free_names` `log_scale` and `logit_rho` under the defaults, `rank` equal to
#'   `dimension`, an empty `null_basis`, `param_name` `"compound_symmetry"`, and
#'   `param_params` holding `link_scale` and `link_rho`.
#'
#' @seealso [ar1()], the other two-value family, where the correlation decays
#'   with distance instead of being shared; [correlation_matrix()] for an
#'   unrestricted correlation; [kron_identity()] to replicate this over
#'   independent groups; and [param_solve()] for the closed inverse.
#'
#' @examples
#' # Two free values at p = 4: a variance and one correlation.
#' s <- compound_symmetry(4)
#' c(n_free = s@n_free)
#' s@free_names
#'
#' eta <- c(log(2), 0.8)
#' m <- param_value(s, eta)
#' round(m, 4)
#'
#' # The eigenvalues are the two the bound comes from.
#' sig2 <- m[1, 1]
#' rho <- m[1, 2] / sig2
#' round(eigen(m, only.values = TRUE)$values, 6)
#' round(c(sig2 * (1 + 3 * rho), rep(sig2 * (1 - rho), 3)), 6)
#'
#' # The log-determinant is closed form and agrees with them.
#' c(closed = param_logdet(s, eta),
#'   from_eigen = sum(log(eigen(m, only.values = TRUE)$values)))
#'
#' # The inverse is compound symmetric too, and is not factorized.
#' round(param_solve(s, eta), 4)
#' max(abs(param_solve(s, eta) - solve(m)))
#'
#' # The round trip closes up to rounding.
#' max(abs(param_free(s, m) - eta))
#'
#' # A free value of 0 is the midpoint of (-1/(p-1), 1), not a correlation of 0.
#' param_value(s, c(0, 0))[1, 2]
#'
#' # The bound tightens with the dimension.
#' t(vapply(c(2, 4, 10),
#'          function(p) compound_symmetry(p)@param_params$link_rho@link_bounds,
#'          numeric(2)))
#'
#' @export
compound_symmetry <- function(dimension,
                              link_scale = linkfunctions7::log_link()) {
  p <- check_param_args(dimension)
  if (p < 2L) {
    stop(paste0(
      "'dimension' must be at least 2: a one by one matrix has no\n",
      "  correlation, and the family would carry a free value with no effect."
    ), call. = FALSE)
  }
  check_positive_link(link_scale)
  link_rho <- linkfunctions7::bounded_link(lwr = -1 / (p - 1), upr = 1)

  CompoundSymmetryParam(
    param_name = "compound_symmetry",
    dimension = p,
    n_free = 2L,
    free_names = c(tagged_name(link_scale, "scale"),
                   tagged_name(link_rho, "rho")),
    rank = p,
    null_basis = empty_null_basis(p),
    param_params = list(
      link_scale = link_scale,
      link_rho = link_rho
    )
  )
}


#' The Scale and the Correlation of an Economical Parameter
#'
#' @description
#' Returns the two scalars from which a [compound_symmetry()] or [ar1()]
#' parameter is built, each with its value and its derivatives to order `order` in its
#' **own** free value. Each scalar depends on one free value alone, so the two
#' families are separable and their derivative assembly is a product of two
#' chains. A value alone is `order = 0`, and no derivative is evaluated then.
#'
#' @param s A [CompoundSymmetryParam()] or [Ar1Param()] object, whose
#'   `param_params$link_scale` and `param_params$link_rho` are read.
#' @param eta A numeric vector of two free values.
#' @param order The highest derivative order wanted, an integer from 0 to 4.
#'
#' @return A list with two components, `scale` and `rho`, each a list of
#'   `order + 1` numbers: the value at index 1 and the derivatives in that free
#'   value at the following indices.
#'
#' @seealso [econ_derivative()] and [cs_pattern()], which consume it, and
#'   [param_value.CompoundSymmetryParam()].
#'
#' @keywords internal
econ_scalars <- function(s, eta, order) {
  list(scale = linkinv_upto(s@param_params$link_scale, eta[1L], order),
       rho = linkinv_upto(s@param_params$link_rho, eta[2L], order))
}


#' Derivative Components of an Economical Parameter
#'
#' @description
#' Assembles one derivative order of \eqn{M = \sigma^2 P(\rho)} for a
#' [compound_symmetry()] or [ar1()] parameter, from the scale's derivatives and
#' the pattern's. The four derivative methods of both families are one call each
#' to this function, differing only in the `pattern` passed.
#'
#' @details
#' The value is a **product** of a function of the first free value and a
#' function of the second, so a component with \eqn{a} scale indices and \eqn{b}
#' correlation indices is the \eqn{a}-th derivative of the scale times the
#' \eqn{b}-th derivative of the pattern. Every derivative is exact, and a tuple
#' enters only through the counts \eqn{a} and \eqn{b}.
#'
#' @param s A [CompoundSymmetryParam()] or [Ar1Param()] object.
#' @param eta A numeric vector of two free values.
#' @param order The derivative order: 1, 2, 3 or 4.
#' @param pattern A function of the parameter and the scalars of
#'   [econ_scalars()], returning a list of matrices: the pattern
#'   \eqn{P(\rho)} and its derivatives in the second free value, to the order
#'   that the scalars carry. See [cs_pattern()] and `ar1_pattern()`.
#'
#' @return A list of `choose(order + 1, order)`, that is `order + 1`, symmetric
#'   matrices keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension`.
#'
#' @seealso [econ_scalars()] for the two chains, [cs_pattern()] for one of the
#'   two patterns, and [econ_logdet_derivative()], its log-determinant
#'   counterpart.
#'
#' @keywords internal
econ_derivative <- function(s, eta, order, pattern) {
  sc <- econ_scalars(s, eta, order)
  pt <- pattern(s, sc)
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(t) {
    a <- sum(t == 1L)
    b <- sum(t == 2L)
    name_dims(sc$scale[[a + 1L]] * pt[[b + 1L]], s)
  })
  stats::setNames(out, param_tuple_names(s, order))
}


#' Log-Determinant Components of an Economical Parameter
#'
#' @description
#' Assembles one derivative order of
#' \eqn{\log|M| = p\log\sigma^2 + q(\rho)} for a [compound_symmetry()] or
#' [ar1()] parameter. The four log-determinant derivative methods of both
#' families are one call each to this function, differing only in the `terms`
#' passed.
#'
#' @details
#' The log-determinant is a **sum** of a function of one free value and a
#' function of the other, where the value is a product, so the structure is
#' stronger: every mixed component is exactly zero, and the pure ones are the two
#' chains taken separately. Under the default log link the scale's own
#' contribution is \eqn{p\eta_1}, linear, so its second, third and fourth
#' derivatives vanish too and only the correlation's chain survives above first
#' order.
#'
#' @param s A [CompoundSymmetryParam()] or [Ar1Param()] object.
#' @param eta A numeric vector of two free values.
#' @param order The derivative order: 1, 2, 3 or 4.
#' @param terms A function of the second free value and the order that returns
#'   the order-th derivative of \eqn{q} in that free value, as
#'   [cs_logdet_chain()] and [ar1_logdet_chain()] do.
#'
#' @return A numeric vector of `order + 1` entries, keyed as
#'   `param_tuple_names(s, order)` and in that order, with every mixed entry
#'   exactly zero.
#'
#' @seealso [cs_logdet_chain()] and [ar1_logdet_chain()] for the correlation's
#'   chain, [econ_scalars()] for the scale's link derivatives, and
#'   [econ_derivative()], its counterpart for the matrix itself.
#'
#' @keywords internal
econ_logdet_derivative <- function(s, eta, order, terms) {
  sc <- econ_scalars(s, eta, order)
  p <- s@dimension
  sd <- sc$scale
  # p * log(scale), composed through the link
  a_log <- lapply(seq_len(order), function(k) {
    (-1)^(k - 1L) * factorial(k - 1L) / sd[[1L]]^k
  })
  d_scale <- compose_order(a_log, sd[-1L], order)
  d_rho <- terms(eta[2L], order)

  idx <- param_tuple_indices(s, order)
  out <- vapply(idx, function(t) {
    if (all(t == 1L)) return(p * d_scale)
    if (all(t == 2L)) return(d_rho)
    0
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, order))
}


#' The Pattern of a Compound-Symmetric Parameter
#'
#' @description
#' Returns \eqn{P(\rho) = I + \rho(J - I)}, the correlation pattern of a
#' compound-symmetric matrix, together with its derivatives in the second free
#' value to the order that the scalars carry. The pattern is **linear** in the
#' correlation, so every derivative is the matching derivative of \eqn{\rho}
#' times the constant matrix \eqn{J - I}, and every order uses the same
#' formula.
#'
#' @param s A [CompoundSymmetryParam()] object, whose `dimension` supplies
#'   \eqn{I} and \eqn{J}.
#' @param sc The scalars of [econ_scalars()], whose `rho` component supplies the
#'   correlation and its derivatives.
#'
#' @return A list of `s@dimension` by `s@dimension` matrices, one more than the
#'   number of derivatives that `sc` carries: the pattern at index 1 and its derivatives at
#'   the following indices. Each derivative has a zero diagonal, the diagonal
#'   of the pattern being the constant 1.
#'
#' @seealso [econ_derivative()], the caller, and `ar1_pattern()`, the AR(1)
#'   counterpart, which is not linear in the correlation.
#'
#' @keywords internal
cs_pattern <- function(s, sc) {
  p <- s@dimension
  off <- matrix(1, p, p) - diag(p)
  r <- sc$rho
  c(list(diag(p) + r[[1L]] * off),
    lapply(seq_along(r)[-1L], function(k) r[[k]] * off))
}


#' @title Value of a Compound Symmetry Parameter
#' @name param_value.CompoundSymmetryParam
#' @description
#' Returns \eqn{\sigma^2\{(1-\rho)I + \rho J\}}: the common variance on the
#' diagonal and the common covariance \eqn{\sigma^2\rho} everywhere else.
#' Positive definiteness follows from the correlation's link, which is bounded
#' below at \eqn{-1/(p-1)}, so nothing is tested here.
#' @param s A [CompoundSymmetryParam()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` symmetric positive definite numeric
#'   matrix with a constant diagonal and constant off-diagonal entries, with
#'   dimnames `v1`, `v2`, ...
#' @seealso [param_free.CompoundSymmetryParam()] for the inverse map, and
#'   [param_solve.CompoundSymmetryParam()] for the closed matrix inverse.
#' @keywords internal
S7::method(param_value, CompoundSymmetryParam) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta, 0L)
  name_dims(sc$scale[[1L]] * cs_pattern(s, sc)[[1L]], s)
}


#' @title Free Vector of a Compound Symmetry Parameter
#' @name param_free.CompoundSymmetryParam
#' @description
#' Returns the two free values behind a compound-symmetric matrix. The variance
#' is the common diagonal entry and the correlation is the common off-diagonal
#' entry divided by it, both read exactly and then carried onto the free scale by
#' the two links. The round trip closes up to rounding.
#' @details
#' A matrix whose diagonal is not constant, or whose off-diagonal entries are not
#' all equal, is **rejected** with a message naming which, and is never averaged
#' into the nearest compound-symmetric matrix, because a silent projection would
#' hide a mistake in the caller's matrix. A caller who wants the nearest such
#' matrix computes it separately and passes the result in.
#' @param s A [CompoundSymmetryParam()] object.
#' @param m A compound-symmetric numeric matrix of side `s@dimension`, already
#'   checked for shape and symmetry by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length 2, named by `s@free_names`.
#' @seealso [param_value.CompoundSymmetryParam()], the map this inverts.
#' @keywords internal
S7::method(param_free, CompoundSymmetryParam) <- function(s, m, ...) {
  p <- s@dimension
  d <- diag(m)
  scl <- max(1, max(abs(m)))
  if (diff(range(d)) > 1e-8 * scl) {
    stop("'m' does not have a constant diagonal, so it is not compound symmetric.",
      call. = FALSE)
  }
  if (d[1L] <= 0) stop("'m' has a non-positive diagonal.", call. = FALSE)
  off <- m[lower.tri(m)]
  if (diff(range(off)) > 1e-8 * scl) {
    stop(paste0(
      "the off-diagonal entries of 'm' are not all equal, so it is not\n",
      "  compound symmetric. It is not averaged."
    ), call. = FALSE)
  }
  stats::setNames(
    c(linkfunctions7::linkfun(s@param_params$link_scale, d[1L]),
      linkfunctions7::linkfun(s@param_params$link_rho, off[1L] / d[1L])),
    s@free_names
  )
}


#' @title Solve of a Compound Symmetry Parameter
#' @name param_solve.CompoundSymmetryParam
#' @description
#' Exact, by Sherman-Morrison. The inverse of
#' \eqn{\sigma^2\{(1-\rho)I + \rho J\}} is
#'
#' \deqn{\frac{1}{\sigma^2(1-\rho)}
#'   \left[I - \frac{\rho}{1 + (p-1)\rho}J\right],}
#'
#' compound symmetric again, so the inverse is written down and no factorization
#' is performed: the cost is \eqn{O(p^2)} per column of `b`, against the
#' \eqn{O(p^3)} Cholesky factorization of the base method. The factors
#' \eqn{1 + (p-1)\rho} and \eqn{1 - \rho} are evaluated from the second free
#' value, as \eqn{p\omega} and \eqn{p(1-\omega)/(p-1)} with \eqn{\omega} the
#' logistic function of that value, so the inverse stays accurate where the
#' correlation rounds to one of its bounds.
#'
#' An exchangeable covariance has an exchangeable precision, so the family is
#' closed under the choice of side.
#' @param s A [CompoundSymmetryParam()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param b A numeric matrix with `s@dimension` rows, already coerced from a
#'   vector and defaulted to the identity by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric matrix with `s@dimension` rows and as many columns as `b`.
#' @seealso [param_solve()] for the generic, and
#'   [param_solve.matrix_parameter()] for the factorization this avoids.
#' @keywords internal
S7::method(param_solve, CompoundSymmetryParam) <- function(s, eta, b = NULL, ...) {
  p <- s@dimension
  sc <- econ_scalars(s, eta, 0L)
  v <- sc$scale[[1L]]
  r <- sc$rho[[1L]]
  # 1 + (p-1) rho = p sigma and 1 - rho = p (1 - sigma)/(p-1), sigma being the
  # logistic function of the free value: exact at both ends of the interval.
  one_plus <- p * stats::plogis(eta[2L])
  one_minus <- p / (p - 1) * stats::plogis(-eta[2L])
  inv <- (diag(p) - (r / one_plus) * matrix(1, p, p)) / (v * one_minus)
  inv %*% b
}


#' @title Derivatives of a Compound Symmetry Parameter
#' @name param_d1.CompoundSymmetryParam
#' @description
#' Closed form at all four orders. One page covers [param_d1()], [param_d2()],
#' [param_d3()] and [param_d4()] for this family because the four are the same
#' product taken at different orders: the value is
#' \eqn{M = \sigma^2 P(\rho)} with \eqn{P(\rho) = I + \rho(J - I)} linear in the
#' correlation, so a component with \eqn{a} scale indices and \eqn{b} correlation
#' indices is
#'
#' \deqn{\partial^{(a,b)} M
#'   = \frac{\mathrm{d}^a \sigma^2}{\mathrm{d}\eta_1^a}
#'     \cdot \frac{\mathrm{d}^b P}{\mathrm{d}\eta_2^b},}
#'
#' one derivative of each chain. A tuple enters only through the counts \eqn{a}
#' and \eqn{b}, and every order uses the same formula.
#' @details
#' The four methods return lists of `order + 1` matrices, keyed by the tuple names
#' of their own order:
#'
#' - [param_d1()]: 2 entries, keyed by `free_names`.
#' - [param_d2()]: 3 entries, `param_tuple_names(s, 2)`.
#' - [param_d3()]: 4 entries, `param_tuple_names(s, 3)`.
#' - [param_d4()]: 5 entries, `param_tuple_names(s, 4)`.
#'
#' The lists are this short because there are two free values: at \eqn{d = 2} the
#' count \eqn{\binom{d+k-1}{k}} is \eqn{k+1}.
#'
#' Every component with at least one correlation index has a zero diagonal, the
#' diagonal of \eqn{P} being the constant 1; and a component with \eqn{b}
#' correlation indices carries the \eqn{b}-th derivative of the link alone,
#' \eqn{P} being linear in \eqn{\rho}.
#' @param s A [CompoundSymmetryParam()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named list of `order + 1` symmetric `s@dimension` by `s@dimension`
#'   matrices; see **Details** for the keying of each order.
#' @seealso [econ_derivative()], which assembles all four, [cs_pattern()] for the
#'   pattern and its derivatives, and
#'   [param_dlogdet.CompoundSymmetryParam()] for the log-determinant's orders.
#' @keywords internal
S7::method(param_d1, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 1L, cs_pattern)
}

#' @rdname param_d1.CompoundSymmetryParam
#' @name param_d2.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d2, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 2L, cs_pattern)
}

#' @rdname param_d1.CompoundSymmetryParam
#' @name param_d3.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d3, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 3L, cs_pattern)
}

#' @rdname param_d1.CompoundSymmetryParam
#' @name param_d4.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d4, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_derivative(s, eta, 4L, cs_pattern)
}


#' Log-Determinant Chain of a Compound-Symmetric Parameter
#'
#' @description
#' Returns a function of the second free value \eqn{\eta_2} and an order that
#' gives the order-th derivative of \eqn{q(\rho) = \log\{1 + (p-1)\rho\} +
#' (p-1)\log(1-\rho)} in \eqn{\eta_2}, with the correlation carried by
#' `bounded_link(-1/(p-1), 1)`.
#'
#' @details
#' With \eqn{\omega = (1 + e^{-\eta_2})^{-1}} and \eqn{u = 1 - \omega}, the
#' correlation is \eqn{\rho = -1/(p-1) + p\omega/(p-1)}, so
#' \eqn{1 + (p-1)\rho = p\omega} and \eqn{1 - \rho = pu/(p-1)} exactly, and
#'
#' \deqn{q = \log p + \log\omega + (p-1)\{\log(p/(p-1)) + \log u\}.}
#'
#' The two factors are evaluated from \eqn{\eta_2} through
#' `plogis(., log.p = TRUE)`, never from a rounded \eqn{\rho}, which loses all
#' accuracy once \eqn{\rho} is within rounding of either end of its interval.
#' The derivatives follow from \eqn{\omega' = \omega u} and
#' \eqn{u' = -\omega u}: \eqn{u - (p-1)\omega}, then \eqn{-p\omega u},
#' \eqn{-p\omega u(u - \omega)} and
#' \eqn{-p\omega u\{(u - \omega)^2 - 2\omega u\}}.
#'
#' @param s A [CompoundSymmetryParam()] object, whose `dimension` supplies
#'   \eqn{p}.
#'
#' @return A function of `(e, order)`, with `e` a single number and `order` an
#'   integer from 1 to 4, returning a single number.
#'
#' @seealso [econ_logdet_derivative()], which calls it, and
#'   [ar1_logdet_chain()] for the AR(1) counterpart.
#'
#' @keywords internal
cs_logdet_chain <- function(s) {
  p <- s@dimension
  function(e, order) {
    a <- stats::plogis(e)
    b <- stats::plogis(-e)
    switch(order,
      b - (p - 1) * a,
      -p * a * b,
      -p * a * b * (b - a),
      -p * a * b * ((b - a)^2 - 2 * a * b))
  }
}


#' @title Log-Determinant of a Compound Symmetry Parameter
#' @name param_logdet.CompoundSymmetryParam
#' @description
#' Closed form, from the two distinct eigenvalues
#' \eqn{\sigma^2\{1+(p-1)\rho\}} and \eqn{\sigma^2(1-\rho)}:
#'
#' \deqn{\log|M| = p\log\sigma^2 + \log\{1+(p-1)\rho\} + (p-1)\log(1-\rho).}
#'
#' Three logarithms and no factorization, whatever \eqn{p} is, where the base
#' method takes an \eqn{O(p^3)} eigendecomposition.
#'
#' The two correlation logarithms are evaluated from the second free value, as
#' \eqn{\log p + \log\omega} and \eqn{\log(p/(p-1)) + \log(1-\omega)} with
#' \eqn{\omega} the logistic function of that value, through
#' `plogis(., log.p = TRUE)`, and not from the rounded correlation. The value is
#' then accurate at every free value, including those where the correlation
#' rounds to one of its bounds.
#' @param s A [CompoundSymmetryParam()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number, finite at every free vector: both eigenvalues are
#'   strictly positive inside the correlation's link bounds.
#' @seealso [param_dlogdet.CompoundSymmetryParam()] for its four derivative
#'   orders, and [cs_logdet_chain()] for the derivatives of the correlation's
#'   term.
#' @keywords internal
S7::method(param_logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  sc <- econ_scalars(s, eta, 0L)
  p <- s@dimension
  p * log(sc$scale[[1L]]) + log(p) + stats::plogis(eta[2L], log.p = TRUE) +
    (p - 1) * (log(p) - log(p - 1) + stats::plogis(-eta[2L], log.p = TRUE))
}

#' @title Log-Determinant Derivatives of a Compound Symmetry Parameter
#' @name param_dlogdet.CompoundSymmetryParam
#' @description
#' Closed form at all four orders. One page covers [param_dlogdet()],
#' [param_d2logdet()], [param_d3logdet()] and [param_d4logdet()] because the
#' four differ only in the order taken:
#' \eqn{\log|M| = p\log\sigma^2 + q(\rho)} is a **sum** of a function of one free
#' value and a function of the other, so every mixed component is exactly zero
#' and the pure ones are the two chains taken separately.
#' @details
#' The scale's chain is \eqn{p} times the derivatives of \eqn{\log h}, which
#' under the default log link is \eqn{p} at first order and 0 above it. The
#' correlation's is [cs_logdet_chain()], the derivatives of \eqn{q} in the second
#' free value written in closed form in the logistic function of that value, so
#' that they stay accurate where the correlation rounds to one of its bounds.
#'
#' The four methods return vectors of `order + 1` entries, keyed by the tuple
#' names of their own order. At \eqn{p = 4} and \eqn{\eta = (\log 2, 0.8)} the
#' second order is \eqn{(0, -0.856, 0)} over `log_scale:log_scale`,
#' `logit_rho:logit_rho` and `log_scale:logit_rho`: the first zero comes from
#' the log link and the last from the separability.
#'
#' Unlike [log_cholesky()], where the whole log-determinant is linear, the
#' correlation's chain is non-zero at all four orders, so a check of the higher
#' orders compares non-zero numbers for this family, as for [ar1()],
#' [autoregressive()] and [correlation_matrix()].
#' @param s A [CompoundSymmetryParam()] object.
#' @param eta A numeric vector of two free values, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named numeric vector of `order + 1` entries, keyed as
#'   `param_tuple_names(s, order)`, with every mixed entry exactly zero.
#' @seealso [param_logdet.CompoundSymmetryParam()] for the quantity
#'   differentiated, [econ_logdet_derivative()], which assembles all four, and
#'   [param_dlogdet.Ar1Param()] for the other two-value family.
#' @keywords internal
S7::method(param_dlogdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 1L, cs_logdet_chain(s))
}

#' @rdname param_dlogdet.CompoundSymmetryParam
#' @name param_d2logdet.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d2logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 2L, cs_logdet_chain(s))
}

#' @rdname param_dlogdet.CompoundSymmetryParam
#' @name param_d3logdet.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d3logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 3L, cs_logdet_chain(s))
}

#' @rdname param_dlogdet.CompoundSymmetryParam
#' @name param_d4logdet.CompoundSymmetryParam
#' @keywords internal
S7::method(param_d4logdet, CompoundSymmetryParam) <- function(s, eta, ...) {
  econ_logdet_derivative(s, eta, 4L, cs_logdet_chain(s))
}
