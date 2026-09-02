#' @include numerical_fallbacks.R
NULL


#' Scaled Fixed Matrix Parameter
#'
#' @description
#' The S7 class of a fixed symmetric positive semidefinite matrix \eqn{P} carried
#' by a single positive scale, \eqn{M(\eta) = h(\eta) P}. It is the only family
#' in the package that is routinely **rank deficient**: \eqn{P} may be a
#' difference penalty or a basis Gram matrix with a genuine null space, and the
#' class records that rank and null basis at construction.
#'
#' [scaled_matrix()] builds one. With `link = NULL` the object holds \eqn{P}
#' itself and has no free value at all, and `param_name` is then `"fixed"`
#' instead of `"scaled"`.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `ScaledMatrixParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `p`, the fixed matrix; `link`, the link object or `NULL`; and `logdet_p`,
#'   the log pseudo-determinant of \eqn{P} computed once at construction.
#'   `n_free` is 1, or 0 when `link` is `NULL`.
#'
#' @seealso [scaled_matrix()], the constructor, and [matrix_parameter()] for the
#'   properties this inherits.
#'
#' @examples
#' # A ridge is the identity, scaled: full rank, one free value.
#' r <- scaled_matrix(diag(3))
#' c(S7::S7_inherits(r, ScaledMatrixParam), rank = r@rank, n_free = r@n_free)
#'
#' # A second-difference penalty is deficient by two, and says so.
#' q <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
#' c(dimension = q@dimension, rank = q@rank, null = ncol(q@null_basis))
#'
#' # With no link there is nothing to estimate.
#' f <- scaled_matrix(diag(3), link = NULL)
#' c(n_free = f@n_free, name = f@param_name)
#'
#' @export
ScaledMatrixParam <- S7::new_class("ScaledMatrixParam", parent = matrix_parameter)


#' Construct a Scaled Fixed Matrix
#'
#' @description
#' Returns an object holding the map \eqn{M(\eta) = h(\eta)\,P}: a fixed
#' symmetric positive semidefinite matrix \eqn{P}, supplied by the caller, times
#' one positive scale. With `link = NULL` it holds \eqn{P} alone and has no free
#' value.
#'
#' This is the commonest penalty in semiparametric regression, and it is the
#' reason the package admits rank-deficient matrices at all. \eqn{P} may be the
#' Gram matrix of a basis derivative, a difference penalty
#' \eqn{\Delta^\top \Delta}, a neighborhood matrix, or the identity, which makes
#' the object a ridge.
#'
#' @details
#' # Everything is a constant times a function of the scale
#'
#' Nothing here needs deriving. With the default log link, where
#' \eqn{h(\eta) = e^{\eta}} and \eqn{\lambda = h(\eta)},
#'
#' \deqn{\partial_\eta M = M, \qquad \partial^2_\eta M = M,}
#' \deqn{\log|M|_+ = r\,\eta + \log|P|_+, \qquad
#'       \partial_\eta \log|M|_+ = r, \qquad
#'       \partial^2_\eta \log|M|_+ = 0,}
#'
#' with \eqn{r} the rank and \eqn{\log|P|_+} the log pseudo-determinant of
#' \eqn{P}, computed once at construction and stored. Under another link the
#' derivatives carry that link's own, through [diag_dlog()].
#'
#' # Why the derivative of the log pseudo-determinant matters
#'
#' It equals the rank, and that is what leaves the scale estimable. Write a
#' penalty as a negative log prior,
#'
#' \deqn{\tfrac{\lambda}{2}\beta^\top P \beta - \tfrac{r}{2}\log\lambda,}
#'
#' and the stationary point is \eqn{\lambda = r / (\beta^\top P \beta)}. Drop the
#' second term and the derivative keeps one sign, sending the scale to zero. That
#' second term is the normalizing constant of the prior, which is why this
#' package keeps it.
#'
#' # Rank deficiency is admitted, and what it means
#'
#' A deficient \eqn{P} makes the corresponding Gaussian improper, so it is a
#' legitimate **penalty** without being a legitimate density. Both readings of a
#' multivariate Gaussian need full rank: a singular covariance is degenerate on a
#' subspace, and a singular precision does not normalize. The object still
#' answers [param_logdet()], with the pseudo-determinant, and [param_solve()] and
#' [param_factor()] reject it, a consumer of an improper prior needing the
#' quadratic form and the pseudo-determinant instead of an inverse.
#'
#' @section Notation:
#' \eqn{P} is the fixed matrix, \eqn{r} its rank, \eqn{\eta} the single free
#' value, \eqn{h = g^{-1}} the inverse link and \eqn{\lambda = h(\eta)} the
#' scale. \eqn{\log|\cdot|_+} is the log pseudo-determinant, the sum of the
#' logarithms of the \eqn{r} non-zero eigenvalues.
#'
#' @param p A symmetric positive semidefinite numeric matrix. It must be square,
#'   free of `NA`, and symmetric to \eqn{10^{-8}} relative, and it is symmetrized
#'   before use. A matrix whose largest eigenvalue is not positive throws `'p'
#'   must be positive semidefinite and not identically zero.`, and one whose
#'   smallest eigenvalue is below `-tol * max(ev)` throws a message quoting both
#'   eigenvalues.
#' @param link A \pkg{linkfunctions7} link carrying the free value onto the
#'   positive scale, `linkfunctions7::log_link()` by default. It must map onto
#'   the positive half line, so `identity_link()` is rejected, and from the whole
#'   real line, which rules out `sqrt_link()` and its relatives; see
#'   [diagonal_matrix()] for the two conditions. `NULL` means the matrix is fully
#'   known: `n_free` is then 0, `free_names` is empty, and `param_name` is
#'   `"fixed"`.
#' @param tol The relative tolerance below which a singular value of `p` counts
#'   as zero when its rank is determined, `1e-10` by default. It is passed to
#'   [param_null_basis()] and also sets how negative an eigenvalue may be before
#'   `p` is rejected.
#'
#' @return An object of class [ScaledMatrixParam()], with `n_free` 1 or 0,
#'   `free_names` a single link-tagged label or empty, `rank` and `null_basis`
#'   from [param_null_basis()], and `param_params` holding `p`, `link` and
#'   `logdet_p`.
#'
#' @seealso [scalar_matrix()], which is this on the identity from the other
#'   direction, [sum_struct()] for a non-negative combination of several fixed
#'   matrices, [param_null_basis()] for the rank, and [param_logdet()] for the
#'   pseudo-determinant.
#'
#' @examples
#' # A ridge: the identity, scaled. Full rank, and log|M| = p * eta.
#' r <- scaled_matrix(diag(4))
#' c(rank = r@rank, logdet_at_0 = param_logdet(r, 0))
#' all.equal(param_logdet(r, 1.3), 4 * 1.3)
#'
#' # A second-difference penalty on six coefficients: deficient by two, its
#' # null space being the constants and the straight lines.
#' d <- diff(diag(6), differences = 2)
#' s <- scaled_matrix(crossprod(d))
#' c(dimension = s@dimension, rank = s@rank)
#' max(abs(param_value(s, 0.4) %*% s@null_basis))
#'
#' # The derivative of the log pseudo-determinant is the rank, at any scale.
#' c(at_minus_4 = param_dlogdet(s, -4), at_4 = param_dlogdet(s, 4),
#'   rank = s@rank)
#'
#' # Which leaves the scale estimable: the stationary point of
#' # lambda/2 b'Pb - r/2 log(lambda) is r / (b'Pb).
#' set.seed(2)
#' b <- rnorm(6)
#' bPb <- drop(crossprod(b, crossprod(d) %*% b))
#' obj <- function(lam) lam / 2 * bPb - s@rank / 2 * log(lam)
#' c(closed_form = s@rank / bPb, numeric = optimize(obj, c(1e-6, 100))$minimum)
#'
#' # A deficient family has no inverse and says so.
#' try(param_solve(s, 0))
#'
#' # With no link the matrix is fully known and there is nothing to estimate.
#' f <- scaled_matrix(diag(3), link = NULL)
#' c(n_free = f@n_free, logdet = param_logdet(f, numeric(0)))
#'
#' @export
scaled_matrix <- function(p, link = linkfunctions7::log_link(),
                          tol = 1e-10) {

  if (!is.matrix(p) || !is.numeric(p)) {
    stop("'p' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(p) != ncol(p)) stop("'p' must be square.", call. = FALSE)
  if (anyNA(p)) stop("'p' must not contain missing values.", call. = FALSE)
  if (max(abs(p - t(p))) > 1e-8 * max(1, max(abs(p)))) {
    stop("'p' must be symmetric.", call. = FALSE)
  }
  p <- (p + t(p)) / 2
  dim_p <- as.integer(nrow(p))

  ev <- eigen(p, symmetric = TRUE, only.values = TRUE)$values
  if (max(ev) <= 0) {
    stop("'p' must be positive semidefinite and not identically zero.",
      call. = FALSE
    )
  }
  if (min(ev) < -tol * max(ev)) {
    stop(sprintf(paste0(
      "'p' is not positive semidefinite: its smallest eigenvalue is %s\n",
      "  against a largest of %s."
    ), format(min(ev)), format(max(ev))), call. = FALSE)
  }

  ns <- param_null_basis(p, tol = tol)
  # The pseudo-determinant of P: a constant of the parameter, computed once.
  logdet_p <- sum(log(ev[seq_len(ns$rank)]))

  if (!is.null(link)) check_positive_link(link)

  ScaledMatrixParam(
    param_name = if (is.null(link)) "fixed" else "scaled",
    dimension = dim_p,
    n_free = if (is.null(link)) 0L else 1L,
    free_names = if (is.null(link)) character(0) else tagged_name(link, "scale"),
    rank = ns$rank,
    null_basis = ns$null_basis,
    param_params = list(p = p, link = link, logdet_p = logdet_p)
  )
}


#' The Scale Behind a Free Vector, and Its Derivatives
#'
#' @description
#' Returns the inverse link's value and its first two derivatives at the single
#' free value, so the family's [param_value()], [param_d1()] and [param_d2()] all
#' read one call. For a fixed parameter, built with `link = NULL`, it returns
#' \eqn{h = 1} and both derivatives 0, which makes those three methods work
#' unchanged: the matrix is \eqn{P} and its derivatives are empty lists.
#'
#' @param s A [ScaledMatrixParam()] object, whose `param_params$link` is read.
#' @param eta A numeric vector of free values: one value, or `numeric(0)` for a
#'   fixed parameter.
#'
#' @return A list with three single numbers, `h`, `d1` and `d2`.
#'
#' @seealso [diag_dlog()] for the higher orders of \eqn{\log h}, and
#'   [param_value.ScaledMatrixParam()], the first caller.
#'
#' @keywords internal
scaled_scale <- function(s, eta) {
  lk <- s@param_params$link
  if (is.null(lk)) return(list(h = 1, d1 = 0, d2 = 0))
  list(
    h = linkfunctions7::linkinv(lk, eta),
    d1 = linkfunctions7::dlinkinv(lk, eta),
    d2 = linkfunctions7::d2linkinv(lk, eta)
  )
}


#' @title Matrix of a Scaled Parameter
#' @name param_value.ScaledMatrixParam
#' @description
#' Returns \eqn{M = h(\eta)\,P}, the stored fixed matrix times the scale, which
#' is \eqn{P} itself for a fixed parameter. One elementwise multiplication and no
#' decomposition. The result inherits \eqn{P}'s rank exactly, a positive multiple
#' leaving a null space alone, which is why the class can record that rank at
#' construction.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A symmetric positive semidefinite `s@dimension` by `s@dimension`
#'   numeric matrix with dimnames `v1`, `v2`, ..., of rank `s@rank`. Positive
#'   definite only where `s@rank == s@dimension`.
#' @seealso [param_free.ScaledMatrixParam()] for the inverse, and
#'   [scaled_scale()] for the scale.
#' @keywords internal
S7::method(param_value, ScaledMatrixParam) <- function(s, eta, ...) {
  name_dims(scaled_scale(s, eta)$h * s@param_params$p, s)
}


#' @title First Derivative of a Scaled Parameter
#' @name param_d1.ScaledMatrixParam
#' @description
#' Closed form: \eqn{\partial_\eta M = h'(\eta)\,P}, the fixed matrix times the
#' link's own first derivative. Under the default log link \eqn{h' = h}, so the
#' derivative **is** the matrix, which the example on [param_d1()] checks.
#'
#' A fixed parameter, built with `link = NULL`, has no free value, so the list is
#' empty.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list with one symmetric matrix named by `s@free_names`, or an empty
#'   list for a fixed parameter.
#' @seealso [param_d2.ScaledMatrixParam()] for the order above, and
#'   [scaled_scale()] for the link derivatives.
#' @keywords internal
S7::method(param_d1, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(list(), character(0)))
  stats::setNames(
    list(name_dims(scaled_scale(s, eta)$d1 * s@param_params$p, s)),
    s@free_names
  )
}


#' @title Second Derivative of a Scaled Parameter
#' @name param_d2.ScaledMatrixParam
#' @description
#' Closed form: \eqn{\partial^2_\eta M = h''(\eta)\,P}. With one free value there
#' is one component, and under the default log link it is the matrix again,
#' \eqn{h'' = h}. Empty for a fixed parameter.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list with one symmetric matrix keyed as `param_tuple_names(s)`, or
#'   an empty list for a fixed parameter.
#' @seealso [param_d1.ScaledMatrixParam()] and [param_d3.ScaledMatrixParam()] for
#'   the neighboring orders.
#' @keywords internal
S7::method(param_d2, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(list(), character(0)))
  stats::setNames(
    list(name_dims(scaled_scale(s, eta)$d2 * s@param_params$p, s)),
    param_tuple_names(s)
  )
}


#' @title Log-Determinant of a Scaled Parameter
#' @name param_logdet.ScaledMatrixParam
#' @description
#' Closed form,
#'
#' \deqn{\log|M|_+ = r \log h(\eta) + \log|P|_+,}
#'
#' with \eqn{r} the rank and \eqn{\log|P|_+} a constant of the object, computed
#' once at construction. The scale enters once per non-zero eigenvalue, which is
#' where the factor \eqn{r} comes from. One expression covers both cases: it is
#' the log-determinant where the family is of full rank and the log
#' pseudo-determinant otherwise. Under the log link it is \eqn{r\eta + \log|P|_+},
#' so a step of 1 in the free value moves it by \eqn{r}.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic. For a fixed parameter it is `numeric(0)` and the
#'   result is `logdet_p` alone.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.ScaledMatrixParam()] for its gradient, and
#'   [param_null_basis()] for the rank it uses.
#' @keywords internal
S7::method(param_logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  s@rank * log(scaled_scale(s, eta)$h) + s@param_params$logdet_p
}


#' @title Log-Determinant Gradient of a Scaled Parameter
#' @name param_dlogdet.ScaledMatrixParam
#' @description
#' Closed form: \eqn{r\, h'(\eta)/h(\eta)}, which under the default log link is
#' the **rank itself**, at every scale. That is the fact a penalty rests on: it
#' is the term that makes the smoothing parameter estimable, the stationary point
#' of \eqn{\tfrac{\lambda}{2}\beta^\top P\beta - \tfrac{r}{2}\log\lambda} being
#' \eqn{\lambda = r/(\beta^\top P \beta)}.
#'
#' The constant \eqn{\log|P|_+} contributes nothing, so a deficient family
#' answers with its rank, never with its dimension.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`, or
#'   empty for a fixed parameter.
#' @seealso [scaled_matrix()] for the estimability argument in full, and
#'   [param_d2logdet.ScaledMatrixParam()] for the order above.
#' @keywords internal
S7::method(param_dlogdet, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(numeric(0), character(0)))
  sc <- scaled_scale(s, eta)
  stats::setNames(s@rank * sc$d1 / sc$h, s@free_names)
}


#' @title Log-Determinant Hessian of a Scaled Parameter
#' @name param_d2logdet.ScaledMatrixParam
#' @description
#' Closed form: \eqn{r\{h''/h - (h'/h)^2\}}, the rank times the second derivative
#' of \eqn{\log h}. Under the default log link it is exactly zero, the log
#' pseudo-determinant being \eqn{r\eta + \log|P|_+} and so linear in the free
#' value. Under a link with curvature it is not.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector with one entry keyed as `param_tuple_names(s)`, or
#'   empty for a fixed parameter.
#' @seealso [param_dlogdet.ScaledMatrixParam()] for the order below, and
#'   [diag_dlog()] for orders three and four.
#' @keywords internal
S7::method(param_d2logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  if (!s@n_free) return(stats::setNames(numeric(0), character(0)))
  sc <- scaled_scale(s, eta)
  stats::setNames(
    s@rank * (sc$d2 / sc$h - (sc$d1 / sc$h)^2),
    param_tuple_names(s)
  )
}


#' @title Free Vector of a Scaled Parameter
#' @name param_free.ScaledMatrixParam
#' @description
#' Returns \eqn{g(h)} where \eqn{h} is the multiple of \eqn{P} that `m` is, read
#' off the ratio at the entry of \eqn{P} of largest magnitude, which is where it
#' is best determined. Exact, and a true inverse of
#' [param_value.ScaledMatrixParam()].
#' @details
#' Three rejections. The ratio must be positive, or `m` is not a positive
#' multiple. `m` must then agree with `h * p` to \eqn{10^{-8}} relative
#' everywhere, or it is not a multiple of \eqn{P} at all: the ratio at one entry
#' is not enough, since any matrix has *some* ratio there. And for a fixed
#' parameter, built with `link = NULL`, the multiple must be 1 to \eqn{10^{-8}},
#' the object having no free value to absorb anything else; the result is then
#' `numeric(0)`.
#' @param s A [ScaledMatrixParam()] object.
#' @param m A symmetric numeric matrix, a positive multiple of the object's
#'   fixed matrix, already checked for shape and symmetry by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free` named by `s@free_names`, or
#'   `numeric(0)` for a fixed parameter.
#' @seealso [param_value.ScaledMatrixParam()], the map this inverts.
#' @keywords internal
S7::method(param_free, ScaledMatrixParam) <- function(s, m, ...) {
  p <- s@param_params$p
  at <- which.max(abs(p))
  h <- m[at] / p[at]
  if (h <= 0) {
    stop("'m' is not a positive multiple of the parameter's fixed matrix.",
      call. = FALSE
    )
  }
  if (max(abs(m - h * p)) > 1e-8 * max(abs(m))) {
    stop(paste0(
      "'m' is not a multiple of the parameter's fixed matrix, so it is not\n",
      "  in the set this parameter parametrizes."
    ), call. = FALSE)
  }
  if (!s@n_free) {
    if (abs(h - 1) > 1e-8) {
      stop("this parameter is fixed, and 'm' is not its matrix.", call. = FALSE)
    }
    return(stats::setNames(numeric(0), character(0)))
  }
  stats::setNames(
    linkfunctions7::linkfun(s@param_params$link, h), s@free_names
  )
}


#' @title Third Derivatives of a Scaled Parameter
#' @name param_d3.ScaledMatrixParam
#' @description
#' Closed form: \eqn{\partial^3_\eta M = h'''(\eta)\,P}, from
#' `linkfunctions7::d3linkinv()`. Every order is the same fixed matrix times a
#' link derivative, so this family costs nothing at any order and is exact at
#' all of them. Under the default log link it is the matrix again.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector with one free value, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list with one matrix keyed as `param_tuple_names(s, 3)`, or an empty
#'   list for a fixed parameter.
#' @seealso [param_d4.ScaledMatrixParam()] for the order above.
#' @keywords internal
S7::method(param_d3, ScaledMatrixParam) <- function(s, eta, ...) {
  v <- linkfunctions7::d3linkinv(s@param_params$link, eta)
  stats::setNames(
    list(name_dims(v * s@param_params$p, s)),
    param_tuple_names(s, 3L)
  )
}

#' @title Fourth Derivatives of a Scaled Parameter
#' @name param_d4.ScaledMatrixParam
#' @description
#' Closed form: \eqn{\partial^4_\eta M = h''''(\eta)\,P}, from
#' `linkfunctions7::d4linkinv()`. Exact, where a family without a closed form
#' would get a product stencil good to about five digits at this order.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector with one free value, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list with one matrix keyed as `param_tuple_names(s, 4)`, or an empty
#'   list for a fixed parameter.
#' @seealso [param_d3.ScaledMatrixParam()] for the order below, and
#'   [numerical_d4()] for the alternative.
#' @keywords internal
S7::method(param_d4, ScaledMatrixParam) <- function(s, eta, ...) {
  v <- linkfunctions7::d4linkinv(s@param_params$link, eta)
  stats::setNames(
    list(name_dims(v * s@param_params$p, s)),
    param_tuple_names(s, 4L)
  )
}

#' @title Third Log-Determinant Derivatives of a Scaled Parameter
#' @name param_d3logdet.ScaledMatrixParam
#' @description
#' Closed form. The log-(pseudo-)determinant is \eqn{r \log h(\eta) + \log|P|_+}
#' with \eqn{r} the rank, so every derivative is \eqn{r} times the matching
#' derivative of \eqn{\log h}, and the constant contributes nothing. With one
#' free value there is one component.
#'
#' Under the default log link \eqn{\log h(\eta) = \eta}, so the answer is exactly
#' zero and this family cannot exercise the order; a link with curvature can.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector with one free value, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector with one entry keyed as `param_tuple_names(s, 3)`, or
#'   empty for a fixed parameter.
#' @seealso [diag_dlog()], which supplies the derivative of \eqn{\log h}, and
#'   [param_d4logdet.ScaledMatrixParam()] for the order above.
#' @keywords internal
S7::method(param_d3logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  stats::setNames(
    s@rank * diag_dlog(s@param_params$link, eta, 3L),
    param_tuple_names(s, 3L)
  )
}

#' @title Fourth Log-Determinant Derivatives of a Scaled Parameter
#' @name param_d4logdet.ScaledMatrixParam
#' @description
#' Closed form: \eqn{r} times the fourth derivative of \eqn{\log h}, which
#' [diag_dlog()] writes out as
#' \eqn{u_4 - 4u_1u_3 - 3u_2^2 + 12u_1^2u_2 - 6u_1^4} with
#' \eqn{u_m = h^{(m)}/h}. The constant \eqn{\log|P|_+} contributes nothing.
#'
#' Exactly zero under the default log link, the quantity being linear in the free
#' value. This is the order at which a numerical route is least usable, so a
#' family with a curved link gains most from the closed form here.
#' @param s A [ScaledMatrixParam()] object.
#' @param eta A numeric vector with one free value, already checked by the
#'   generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector with one entry keyed as `param_tuple_names(s, 4)`, or
#'   empty for a fixed parameter.
#' @seealso [diag_dlog()] for the expression, and
#'   [param_d4logdet.matrix_parameter()] for the numerical route.
#' @keywords internal
S7::method(param_d4logdet, ScaledMatrixParam) <- function(s, eta, ...) {
  stats::setNames(
    s@rank * diag_dlog(s@param_params$link, eta, 4L),
    param_tuple_names(s, 4L)
  )
}
