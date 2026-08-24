#' @include generics.R
NULL

# What a parameter gets for free when it implements only param_value().
# Every method here is registered on the base class, so a subclass that
# supplies a closed form takes over through dispatch with no registration step.


#' Is This the Package's Own Base Class?
#'
#' @description
#' Asks whether an S7 class is [parameter()] or [matrix_parameter()], the two
#' abstract classes this package registers its fallback methods on. That is how
#' a method a subclass wrote is told apart from one it merely inherited, which is
#' the question [param_is_numerical()] answers and [check_parameter()] acts on.
#'
#' @details
#' Two comparisons, in order. Identity is tried first, since it is the usual
#' case and costs nothing. Then the class name together with the owning package,
#' because identity does not survive a package's code being re-evaluated instead
#' of loaded: a class rebuilt from the same definition is a different object.
#' Coverage tools re-evaluate, so an identity-only test passes every ordinary
#' check and fails under `covr` alone.
#'
#' @param cls An S7 class object, as `S7::S7_class(x)` returns or as
#'   `attr(method, "signature")[[1]]` holds. Anything else returns `FALSE`
#'   rather than throwing, `attr()` on a non-class giving `NULL`.
#'
#' @return `TRUE` when `cls` is [parameter()] or [matrix_parameter()], `FALSE`
#'   for any subclass and for anything that is not an S7 class.
#'
#' @seealso [param_is_numerical()], which calls this on the class a method was
#'   registered against.
#'
#' @keywords internal
is_base_param_class <- function(cls) {
  for (base in list(parameter, matrix_parameter)) {
    if (identical(cls, base)) return(TRUE)
    if (identical(attr(cls, "name"), attr(base, "name")) &&
        identical(attr(cls, "package"), attr(base, "package"))) {
      return(TRUE)
    }
  }
  FALSE
}


#' Which of a Parameter's Quantities Come From the Base Class
#'
#' @description
#' Reports, one derivative quantity at a time, whether the parameter has a
#' method of its own or takes the numerical one registered on [parameter()] and
#' [matrix_parameter()]. Ask it before trusting a fourth-order derivative of a
#' family you did not write: `FALSE` everywhere means every quantity is a closed
#' form, and `TRUE` somewhere marks a component whose accuracy is a stencil's.
#'
#' @details
#' # Why the question is worth asking
#'
#' Whether an independent check exists depends on the answer. A derivative
#' computed by finite differences cannot be checked against a finite difference,
#' and a log-determinant read off an eigendecomposition cannot be checked against
#' an eigendecomposition: the comparison is the same arithmetic twice, and it
#' agrees however wrong the parametrization is. [check_parameter()] calls this
#' and reports such a quantity as **not checked** instead of as passed.
#'
#' # What is not in the list
#'
#' [param_solve()] and [param_factor()] are absent deliberately. Their
#' base-class versions are a Cholesky factorization, which is exact whoever
#' performs it, and [check_parameter()] compares them against `base::solve()`
#' either way. Listing them as numerical would suggest an approximation that is
#' not there. [param_value()] is absent because a family that does not implement
#' it has nothing at all, and [param_free()] because its base method throws.
#'
#' # Every shipped family answers FALSE
#'
#' Measured over all fifteen constructors in this package, at every component:
#' none of them uses a numerical route. The fallbacks exist for a family written
#' elsewhere, and the example below builds one to show what a `TRUE` looks like.
#'
#' @param s An object inheriting from class [parameter()].
#'
#' @return A named logical vector, `TRUE` where the base-class method is in
#'   force. Its length depends on the branch: **nine** entries for a
#'   [matrix_parameter()], over `param_d1`, `param_d2`, `param_d3`, `param_d4`,
#'   `param_logdet`, `param_dlogdet`, `param_d2logdet`, `param_d3logdet` and
#'   `param_d4logdet`; **four** for a family that is not a matrix, over the
#'   derivative orders alone, the log-determinant not existing there.
#'
#' @seealso [check_parameter()], which reports a numerical component as not
#'   checked, and [numerical_d1()] and its higher-order siblings, the routes a
#'   `TRUE` names.
#'
#' @examples
#' # Every family in the package is closed form throughout.
#' param_is_numerical(log_cholesky(3))
#' any(param_is_numerical(diagonal_matrix(3)))
#'
#' # Nine components for a matrix family, four for one that is not.
#' c(matrix = length(param_is_numerical(log_cholesky(3))),
#'   not_matrix = length(param_is_numerical(simplex(3))))
#'
#' # A family written from scratch, with param_value() alone: every derivative
#' # order is then numerical, and so is the whole log-determinant block.
#' Toy <- S7::new_class("Toy", parent = matrix_parameter)
#' S7::method(param_value, Toy) <- function(s, eta, ...) {
#'   m <- diag(rep(exp(eta[1]), 2))
#'   dimnames(m) <- list(c("v1", "v2"), c("v1", "v2"))
#'   m
#' }
#' toy <- Toy(param_name = "toy", n_free = 1L, free_names = "log_s",
#'            param_params = list(), dimension = 2L, rank = 2L,
#'            null_basis = matrix(numeric(0), 2, 0), role = "either")
#' param_is_numerical(toy)
#'
#' # The numerical route is usable, and here it can be checked by hand: the
#' # matrix is exp(eta) times the identity, so every derivative is itself.
#' c(numerical = param_d1(toy, 0.4)[[1]][1, 1], exact = exp(0.4))
#'
#' @export
param_is_numerical <- function(s) {
  cls <- S7::S7_class(s)
  gens <- list(
    param_d1 = param_d1,
    param_d2 = param_d2,
    param_d3 = param_d3,
    param_d4 = param_d4,
    param_logdet = param_logdet,
    param_dlogdet = param_dlogdet,
    param_d2logdet = param_d2logdet,
    param_d3logdet = param_d3logdet,
    param_d4logdet = param_d4logdet
  )
  if (!S7::S7_inherits(s, matrix_parameter)) {
    gens <- gens[c("param_d1", "param_d2", "param_d3", "param_d4")]
  }
  vapply(gens, function(g) {
    m <- tryCatch(S7::method(g, cls), error = function(e) NULL)
    if (is.null(m)) return(TRUE)
    is_base_param_class(attr(m, "signature")[[1L]])
  }, logical(1))
}


#' The Spectral Decomposition a Parameter's Quantities Are Read From
#'
#' @description
#' Returns the eigenvalues and eigenvectors of \eqn{M(\eta)}, together with a
#' flag saying which directions carry the matrix. It is what the base-class
#' log-determinant, the pseudo-inverse and the rank-deficient solve all read, so
#' one decomposition serves them.
#'
#' @details
#' Which eigenvalues count as zero is settled **by position**, never by size.
#' `eigen()` returns them in decreasing order, so the first `s@rank` are kept and
#' the rest dropped, whatever their numerical values are. The rank itself comes
#' from the object and is never re-derived: counting eigenvalues above a relative
#' tolerance is not scale invariant, so a family whose components differ by many
#' orders of magnitude would be assigned a different rank at different \eqn{\eta},
#' and a fitted model with smoothing parameters that far apart is ordinary. The
#' object settled the question once, at construction, from the components. See
#' [param_null_basis()] for the measurement.
#'
#' The matrix is symmetrized as `(m + t(m)) / 2` before the decomposition, so an
#' asymmetry of rounding size does not produce complex eigenvalues.
#'
#' @param s A [matrix_parameter()] object, whose `rank` and `dimension` are read.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list with four components
#'   \describe{
#'     \item{`values`}{numeric of length `s@dimension`, in decreasing order.}
#'     \item{`vectors`}{a `s@dimension` by `s@dimension` orthonormal matrix, one
#'       eigenvector per column, matching `values`.}
#'     \item{`keep`}{logical of length `s@dimension`, `TRUE` in the first
#'       `s@rank` positions.}
#'     \item{`matrix`}{the value itself, as [param_value()] returned it, so a
#'       caller needing both does not evaluate the map twice.}
#'   }
#'
#' @seealso [spectrum_pinv()], which builds the Moore-Penrose inverse from this,
#'   and [param_null_basis()], which fixes the rank at construction.
#'
#' @keywords internal
param_spectrum <- function(s, eta) {
  m <- param_value(s, eta)
  e <- eigen((m + t(m)) / 2, symmetric = TRUE)
  keep <- rep(FALSE, s@dimension)
  if (s@rank > 0L) keep[seq_len(s@rank)] <- TRUE
  list(values = e$values, vectors = e$vectors, keep = keep, matrix = m)
}


#' Moore-Penrose Inverse From a Parameter's Spectrum
#'
#' @description
#' Builds the pseudo-inverse of \eqn{M(\eta)} as
#' \eqn{\sum_{j \in \mathrm{keep}} \lambda_j^{-1} v_j v_j^\top}, over the
#' directions [param_spectrum()]'s `keep` flag marks. It is what the
#' rank-deficient branches of the base-class log-determinant derivatives use in
#' place of \eqn{M^{-1}} in the trace identities.
#'
#' @param sp The result of [param_spectrum()]: a list with `values`, `vectors`
#'   and `keep`.
#'
#' @return A symmetric numeric matrix of the same side as `sp$vectors`. When
#'   `keep` is all `FALSE`, which happens only at rank 0, the zero matrix.
#'
#' @seealso [param_spectrum()] for the input, and [param_solve()], which rejects
#'   a deficient family rather than handing this back.
#'
#' @keywords internal
spectrum_pinv <- function(sp) {
  v <- sp$vectors[, sp$keep, drop = FALSE]
  d <- sp$values[sp$keep]
  if (!length(d)) return(matrix(0, nrow(sp$vectors), nrow(sp$vectors)))
  v %*% (t(v) / d)
}


#' Cholesky Factorization, With the Rank Decided Before It
#'
#' @description
#' Returns the lower triangular Cholesky factor of a symmetric matrix, or `NULL`
#' when the matrix is not positive definite to the given relative tolerance.
#' Positive definiteness is decided **before** the factorization is attempted,
#' from the eigenvalues.
#'
#' @details
#' The verdict comes from the spectrum, because `chol()` is not a rank test. On a matrix with an exactly zero
#' eigenvalue the pivot that ought to be zero comes out positive or negative
#' according to rounding, so `chol()` succeeds on some platforms and fails on
#' others, and a branch that asks it whether a matrix is usable gets a different
#' answer on different machines. `min(ev) <= tol * max(ev)` is a statement about
#' the matrix; a caught error is a statement about the arithmetic.
#'
#' @param m A symmetric numeric matrix. Not checked for symmetry; the callers
#'   have already symmetrized.
#' @param tol The relative tolerance below which the smallest eigenvalue counts
#'   as zero, so `m` is rejected when `min(ev) <= tol * max(ev)`. Defaults to
#'   `1e-12`. At that default `diag(c(1, 1e-12))` returns `NULL`, the test being
#'   an inequality, and `diag(c(1, 1e-11))` factors.
#'
#' @return The lower triangular \eqn{L} with \eqn{M = L L^\top}, or `NULL` when
#'   `m` has no eigenvalues, has an `NA` among them, has a non-positive largest
#'   eigenvalue, fails the relative test, or makes [base::chol()] raise after
#'   all.
#'
#' @seealso [param_factor()], the generic whose base method calls this, and
#'   [param_spectrum()] for the decomposition the deficient branches use
#'   instead.
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
#' Returns the step a central difference should use in one component of
#' \eqn{\eta}, scaled by the size of that component. A one-line forward to
#' [numericals7::fd_step()] at accuracy 2, so the step a fallback here takes and
#' the step the stencil library documents cannot drift apart.
#'
#' @details
#' The rule is \eqn{\varepsilon^{1/(k+2)} \max(1, |\eta_k|)} for a \eqn{k}-th
#' derivative, which balances the truncation error against the rounding the
#' division by \eqn{h^k} amplifies. In doubles that is \eqn{6.1 \times 10^{-6}}
#' at first order and \eqn{1.2 \times 10^{-4}} at second, both scaling up once
#' \eqn{|\eta_k|} passes 1.
#'
#' No clamping happens, and there is nothing to clamp away from. The response
#' and parameter steps in \pkg{distributions7} have to keep a node inside a
#' bounded support; the unconstrained scale has no boundary anywhere, so
#' `bounds` is left at its default.
#'
#' @param eta_k The value of the component, a single number.
#' @param order The derivative order the step is for: 1, 2, 3 or 4.
#'
#' @return A single positive number.
#'
#' @seealso [numericals7::fd_step()] for the rule, and [fd_along()], which
#'   applies a stencil at this step.
#'
#' @keywords internal
fd_step <- function(eta_k, order = 1L) {
  numericals7::fd_step(eta_k, order = order, accuracy = 2L)
}


#' One Stencil Along One Free Value
#'
#' @description
#' Applies one central stencil to a function of the free vector, along a single
#' component of it, and returns the derivative in whatever shape `f` produces.
#' The nodes and weights come from [numericals7::fd_offsets()] and
#' [numericals7::fd_weights()] at accuracy 2, so a difference taken here uses the
#' same stencils as the rest of the toolkit.
#'
#' @details
#' The sum \eqn{h^{-k}\sum_j w_j f(\eta + s_j h e_k)} is accumulated term by
#' term, skipping the nodes whose weight is zero, so `f` is called once per
#' non-zero weight: three times for a first or second derivative, five for a
#' third or fourth. The accumulator starts at `NULL` and takes the shape of the
#' first term, which is why `f` may return a matrix, a vector or a scalar without
#' this function knowing which.
#'
#' [numericals7::fd_derivative()] cannot serve here. Its `f` maps a vector of
#' points to the values at those points, and what is differentiated here maps a
#' whole free vector to a matrix, so the nodes have to be placed along one
#' coordinate by hand.
#'
#' @param f A function of the free vector, returning anything that supports `+`
#'   and scalar multiplication.
#' @param eta The free vector, numeric.
#' @param k The component to differentiate along: a position in `1:length(eta)`.
#' @param order The derivative order: 1, 2, 3 or 4.
#' @param h The step. `NULL`, the default, takes [fd_step()] at `eta[k]` and this
#'   order. Pass a value when several stencils have to share one step, as a mixed
#'   derivative does.
#'
#' @return Whatever `f` returns, differentiated `order` times in component `k`.
#'
#' @seealso [numericals7::fd_weights()] for the weights, [fd_step()] for the
#'   step, and [mixed_stencil()], which composes factors across several
#'   components.
#'
#' @keywords internal
fd_along <- function(f, eta, k, order = 1L, h = NULL) {
  s <- numericals7::fd_offsets(order, accuracy = 2L)$central
  w <- numericals7::fd_weights(s, order)
  if (is.null(h)) h <- fd_step(eta[k], order)
  acc <- NULL
  for (j in seq_along(s)) {
    if (w[j] == 0) next
    e <- eta
    e[k] <- eta[k] + s[j] * h
    v <- w[j] * f(e)
    acc <- if (is.null(acc)) v else acc + v
  }
  acc / h^order
}


#' Numerical First Derivatives of a Parameter's Matrix
#'
#' @description
#' Estimates \eqn{\partial V / \partial \eta_k} for every free value by one
#' three-point central difference of [param_value()] in each component,
#'
#' \deqn{\partial_k V \approx
#'   \frac{V(\eta + h e_k) - V(\eta - h e_k)}{2h},}
#'
#' at the step [fd_step()] gives for a first derivative. This is the method
#' [param_d1()] dispatches to when a family has not written its own, and it is
#' exported so that a family being developed can be compared against it.
#'
#' @details
#' The step is \eqn{\varepsilon^{1/3}\max(1, |\eta_k|)}, about
#' \eqn{6.1 \times 10^{-6}} near the origin, giving a truncation error of order
#' \eqn{h^2}. Measured against the closed form of a \eqn{2 \times 2}
#' log-Cholesky covariance, the agreement is \eqn{7 \times 10^{-11}} absolute on
#' entries of size 3, so eleven digits or so. The example below runs that
#' comparison.
#'
#' The result is symmetrized for a [matrix_parameter()], as `(A + t(A)) / 2`. The
#' derivative of a symmetric matrix is symmetric, so the halving corrects
#' rounding and changes nothing else. `V` costs \eqn{2d} evaluations of the map.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `s@n_free` estimates named by `s@free_names`, each shaped
#'   like [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d1()], the generic this serves, [numerical_d2()],
#'   [numerical_d3()] and [numerical_d4()] for the higher orders, and
#'   [param_is_numerical()] to ask whether a given family reaches this at all.
#'
#' @examples
#' # A scalar matrix is exp(eta) times the identity, so the derivative is the
#' # matrix itself and the error of the difference can be read off.
#' s <- scalar_matrix(2)
#' numerical_d1(s, 0.3)
#' max(abs(numerical_d1(s, 0.3)[[1]] - param_value(s, 0.3)))
#'
#' # Against a family that writes its own: the two agree to about 1e-10.
#' q <- log_cholesky(2)
#' eta <- c(0.2, -0.1, 0.4)
#' max(abs(unlist(numerical_d1(q, eta)) - unlist(param_d1(q, eta))))
#'
#' @export
numerical_d1 <- function(s, eta) {
  out <- vector("list", s@n_free)
  names(out) <- s@free_names
  for (k in seq_len(s@n_free)) {
    out[[k]] <- fd_along(function(e) param_value(s, e), eta, k, 1L)
  }
  out
}


#' Numerical Second Derivatives of a Parameter's Matrix
#'
#' @description
#' Estimates the \eqn{d(d+1)/2} distinct second derivatives
#' \eqn{\partial^2 V / \partial \eta_k \partial \eta_l}, taking exactly one
#' difference for each. Which quantity is differenced depends on what the family
#' already supplies: the analytic first derivative where there is one, and
#' [param_value()] itself where there is not. This is the method [param_d2()]
#' dispatches to when a family has not written its own.
#'
#' @details
#' # Three routes, one layer each
#'
#' Writing \eqn{V(\eta)} for [param_value()] and \eqn{\partial_k V} for an
#' analytic first derivative, the component \eqn{(k, l)} is
#'
#' \deqn{\partial_{kl} V \approx
#'   \frac{\partial_k V(\eta + h e_l) - \partial_k V(\eta - h e_l)}{2h}}
#'
#' whenever [param_d1()] is analytic. Where it is not, an off-diagonal pair takes
#' the four-point mixed stencil
#'
#' \deqn{\partial_{kl} V \approx \frac{V(\eta + h_k e_k + h_l e_l)
#'   - V(\eta + h_k e_k - h_l e_l) - V(\eta - h_k e_k + h_l e_l)
#'   + V(\eta - h_k e_k - h_l e_l)}{4 h_k h_l},}
#'
#' and a diagonal pair the three-point second difference. All three carry a
#' truncation error of order \eqn{h^2}.
#'
#' # Why none of them nests
#'
#' The rule the toolkit follows is that differences are never composed **in the
#' same variable**: a difference of a difference multiplies the error of the
#' first stage into the second, and a fourth derivative built that way is noise.
#' Differencing an analytic first derivative is one layer, since only one stage
#' is numerical. The four-point stencil differences two *different* components,
#' and two differences along different coordinates commute into a single product
#' stencil, so it is one layer as well. The diagonal case uses the
#' second-difference stencil directly and never sees a first difference at all.
#'
#' The steps are the order-2 ones, \eqn{\varepsilon^{1/4}\max(1, |\eta_k|)},
#' about \eqn{1.2 \times 10^{-4}} near the origin, since it is a second
#' derivative being estimated whichever route is taken.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `choose(s@n_free + 1, 2)` estimates keyed as
#'   `param_tuple_names(s)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d2()], the generic this serves, [numerical_d1()] for the
#'   order below, and [mixed_stencil()], which the third and fourth orders use.
#'
#' @examples
#' # A scalar matrix: every order is the matrix again, so the error shows.
#' s <- scalar_matrix(2)
#' numerical_d2(s, 0.3)
#' max(abs(numerical_d2(s, 0.3)[[1]] - param_value(s, 0.3)))
#'
#' # Against a family that writes its own second derivatives.
#' q <- log_cholesky(2)
#' eta <- c(0.2, -0.1, 0.4)
#' max(abs(unlist(numerical_d2(q, eta)) - unlist(param_d2(q, eta))))
#'
#' @export
numerical_d2 <- function(s, eta) {
  idx <- param_tuple_indices(s)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s)
  if (!length(idx)) return(out)

  analytic <- !param_is_numerical(s)[["param_d1"]]

  for (i in seq_along(idx)) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    if (analytic) {
      # One layer on the analytic first derivative in the other component.
      out[[i]] <- fd_along(function(e) param_d1(s, e)[[k]], eta, l, 1L)
    } else if (k == l) {
      out[[i]] <- fd_along(function(e) param_value(s, e), eta, k, 2L)
    } else {
      # Two stencils in two DIFFERENT components, which is one mixed stencil
      # and not a difference of a difference: nesting is forbidden along one
      # variable and is what a mixed derivative is along two. The steps are
      # the order-2 ones, as a second derivative asks for.
      hk <- fd_step(eta[k], 2L)
      hl <- fd_step(eta[l], 2L)
      out[[i]] <- fd_along(function(e) {
        fd_along(function(e2) param_value(s, e2), e, l, 1L, h = hl)
      }, eta, k, 1L, h = hk)
    }
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


#' @title Default First Derivatives
#' @name param_d1.parameter
#' @description Fallback: one central difference of
#'   [param_value()] per component (see
#'   [numerical_d1()]).
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d1, parameter) <- function(s, eta, ...) {
  numerical_d1(s, eta)
}

#' @title Default Second Derivatives
#' @name param_d2.parameter
#' @description Fallback: see [numerical_d2()].
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d2, parameter) <- function(s, eta, ...) {
  numerical_d2(s, eta)
}

#' @title Default Log-Determinant
#' @name param_logdet.parameter
#' @description
#' Fallback: the sum of the logs of the eigenvalues the declared rank keeps,
#' which is the log-determinant for a full-rank family and the log
#' pseudo-determinant otherwise.
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, matrix_parameter) <- function(s, eta, ...) {
  sp <- param_spectrum(s, eta)
  v <- sp$values[sp$keep]
  if (any(v <= 0)) {
    stop(sprintf(paste0(
      "'%s' produced %d non-positive eigenvalue(s) among the %d its rank\n",
      "  keeps, so it is outside the set it claims to parametrize."
    ), s@param_name, sum(v <= 0), s@rank), call. = FALSE)
  }
  sum(log(v))
}

#' @title Default Log-Determinant Gradient
#' @name param_dlogdet.parameter
#' @description
#' Fallback: \eqn{\mathrm{tr}(M^{+} \partial_k M)}, with the pseudo-inverse
#' formed from the directions the declared rank keeps, which is the ordinary
#' inverse when the family is of full rank.
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, matrix_parameter) <- function(s, eta, ...) {
  sp <- param_spectrum(s, eta)
  mi <- spectrum_pinv(sp)
  d <- param_d1(s, eta)
  stats::setNames(
    vapply(d, function(dk) sum(mi * dk), numeric(1)),
    s@free_names
  )
}

#' @title Default Log-Determinant Hessian
#' @name param_d2logdet.parameter
#' @description
#' Fallback: \eqn{\mathrm{tr}(M^{+} \partial_{kl} M) -
#' \mathrm{tr}(M^{+} \partial_k M\, M^{+} \partial_l M)}, the derivative of the
#' identity behind [param_dlogdet()].
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d2logdet, matrix_parameter) <- function(s, eta, ...) {
  sp <- param_spectrum(s, eta)
  mi <- spectrum_pinv(sp)
  d <- param_d1(s, eta)
  d2 <- param_d2(s, eta)
  idx <- param_tuple_indices(s)
  out <- vapply(seq_along(idx), function(i) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    sum(mi * d2[[i]]) - sum(t(mi %*% d[[k]]) * (mi %*% d[[l]]))
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s))
}

#' @title Default Solve
#' @name param_solve.parameter
#' @description Fallback: a Cholesky of [param_value()], with the
#'   definiteness verdict taken spectrally (see [chol_pd()]).
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param b A numeric matrix with `s@dimension` rows.
#' @param ... Unused.
#' @return A numeric matrix.
#' @keywords internal
S7::method(param_solve, matrix_parameter) <- function(s, eta, b = NULL, ...) {
  l <- param_factor(s, eta)
  backsolve(t(l), forwardsolve(l, b))
}

#' @title Default Factor
#' @name param_factor.parameter
#' @description Fallback: the lower Cholesky factor of
#'   [param_value()], rejected when the matrix is not positive
#'   definite spectrally.
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A lower triangular numeric matrix.
#' @keywords internal
S7::method(param_factor, matrix_parameter) <- function(s, eta, ...) {
  l <- chol_pd(param_value(s, eta))
  if (is.null(l)) {
    stop(sprintf(paste0(
      "'%s' is not positive definite at this eta, although it declares full\n",
      "  rank. The verdict is spectral, so this is a statement about the\n",
      "  matrix rather than about a factorization that happened to fail."
    ), s@param_name), call. = FALSE)
  }
  l
}

#' @title Rejection to Invert Without a Closed Form
#' @name param_free.parameter
#' @description
#' The base class rejects rather than inverting the map numerically: an
#' optimization-based inverse would return a plausible \eqn{\eta} for a matrix
#' outside the set the family parametrizes.
#' @param s A [parameter()] object.
#' @param m A symmetric numeric matrix.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @keywords internal
S7::method(param_free, parameter) <- function(s, m, ...) {
  stop(sprintf(paste0(
    "'%s' does not implement param_free(). The inverse map is exact or\n",
    "  rejected, never obtained by optimization, because a numerical inverse\n",
    "  would return a plausible eta for a matrix outside the set."
  ), s@param_name), call. = FALSE)
}


#' Numerical Third Derivatives of a Parameter's Value
#'
#' @description
#' Estimates the distinct third derivatives by applying one product stencil
#' directly to [param_value()], per index tuple. A component the tuple repeats
#' contributes the one-dimensional stencil of the matching order, and a component
#' appearing once contributes a two-point central factor; the product is
#' evaluated in a single pass over the map. This is the method [param_d3()]
#' dispatches to when a family has not written its own.
#'
#' @details
#' # The stencil
#'
#' Let the tuple's distinct components be \eqn{c_1, \dots, c_r} with
#' multiplicities \eqn{m_1, \dots, m_r} summing to three, and let
#' \eqn{(o^{(m)}_i, w^{(m)}_i)} be the offsets and weights of the central stencil
#' for an \eqn{m}-th derivative in one variable. The estimate is their tensor
#' product,
#'
#' \deqn{\partial_{c_1}^{m_1} \cdots \partial_{c_r}^{m_r} V \approx
#'   \frac{1}{\prod_{j} h_j^{m_j}}
#'   \sum_{i_1, \dots, i_r} \left(\prod_{j=1}^{r} w^{(m_j)}_{i_j}\right)
#'   V\Bigl(\eta + \textstyle\sum_{j=1}^{r} o^{(m_j)}_{i_j} h_j e_{c_j}\Bigr),}
#'
#' evaluated at one point per combination of nodes, with the zero-weight nodes
#' skipped. A tuple naming one component three times costs four evaluations of
#' the map, one naming a component twice and another once six, and one naming
#' three distinct components eight.
#'
#' # Why it goes straight to the map
#'
#' A lower-order numerical derivative is never differenced. The rounding of a
#' difference is amplified by \eqn{h^{-k}}, so two stages multiply their errors
#' and a third derivative built from a first and then a second is far worse than
#' one stencil of the order wanted. Going to [param_value()] directly keeps the
#' rounding at a single layer whatever the family implements, at the price of
#' more evaluations of the map.
#'
#' # Accuracy
#'
#' Truncation is of order \eqn{h^2} at a step of
#' \eqn{\varepsilon^{1/5}\max(1, |\eta_k|)}. Measured against the closed form of
#' a \eqn{2 \times 2} log-Cholesky covariance, the agreement is
#' \eqn{7 \times 10^{-6}} absolute on entries of size 12, so about six digits.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `choose(s@n_free + 2, 3)` estimates keyed as
#'   `param_tuple_names(s, 3)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d3()], the generic this serves, [mixed_stencil()], which
#'   builds the product, and [numerical_d4()] for the order above.
#'
#' @examples
#' # A scalar matrix: every order is the matrix again, so the error shows.
#' s <- scalar_matrix(2)
#' numerical_d3(s, 0.3)
#' max(abs(numerical_d3(s, 0.3)[[1]] - param_value(s, 0.3)))
#'
#' # Against a family that writes its own third derivatives.
#' q <- log_cholesky(2)
#' eta <- c(0.2, -0.1, 0.4)
#' ana <- param_d3(q, eta)
#' c(gap = max(abs(unlist(numerical_d3(q, eta)) - unlist(ana))),
#'   scale = max(abs(unlist(ana))))
#'
#' @export
numerical_d3 <- function(s, eta) {
  idx <- param_tuple_indices(s, 3L)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, 3L)
  f <- function(e) param_value(s, e)
  for (i in seq_along(idx)) {
    out[[i]] <- mixed_stencil(f, eta, idx[[i]])
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


#' Numerical Fourth Derivatives of a Parameter's Value
#'
#' @description
#' Estimates the distinct fourth derivatives, one product stencil per index
#' tuple, applied directly to [param_value()]. It is the order-four analogue of
#' [numerical_d3()] and the method [param_d4()] dispatches to when a family has
#' not written its own. Treat it as a starting point: at this order a difference
#' keeps about five digits, which is why every family this package ships carries
#' a closed form instead.
#'
#' @details
#' # The stencil
#'
#' The tensor product is the one written out under [numerical_d3()], with the
#' multiplicities summing to four. A tuple naming one component four times costs
#' five evaluations of the map, one naming two components twice each nine, and
#' one naming four distinct components sixteen.
#'
#' # Accuracy
#'
#' Truncation is of order \eqn{h^2} and rounding of order
#' \eqn{\varepsilon / h^4}, so balancing the two gives a step of
#' \eqn{\varepsilon^{1/6}\max(1, |\eta_k|)}, about \eqn{2.5 \times 10^{-3}} near
#' the origin, and an attainable accuracy of order \eqn{\varepsilon^{1/3}}.
#' Measured against the closed form of a \eqn{2 \times 2} log-Cholesky
#' covariance, the agreement is \eqn{1.2 \times 10^{-4}} absolute on entries of
#' size 24, so five digits. That is enough to catch a transcription error in a
#' closed form. It is not enough to fit with.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `choose(s@n_free + 3, 4)` estimates keyed as
#'   `param_tuple_names(s, 4)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d4()], the generic this serves, [numerical_d3()] for the order
#'   below, and [check_parameter()], which uses a stencil of this kind as the
#'   independent reference for a closed form.
#'
#' @examples
#' # A scalar matrix: every order is the matrix again, so the error shows.
#' s <- scalar_matrix(2)
#' max(abs(numerical_d4(s, 0.3)[[1]] - param_value(s, 0.3)))
#'
#' # Against a family that writes its own: five digits, as stated.
#' q <- log_cholesky(2)
#' eta <- c(0.2, -0.1, 0.4)
#' ana <- param_d4(q, eta)
#' c(gap = max(abs(unlist(numerical_d4(q, eta)) - unlist(ana))),
#'   scale = max(abs(unlist(ana))))
#'
#' # Which is still ample to catch a wrong closed form: a 1 per cent error in
#' # one component is four hundred times the noise.
#' 0.01 * max(abs(unlist(ana)))
#'
#' @export
numerical_d4 <- function(s, eta) {
  idx <- param_tuple_indices(s, 4L)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, 4L)
  f <- function(e) param_value(s, e)
  for (i in seq_along(idx)) {
    out[[i]] <- mixed_stencil(f, eta, idx[[i]])
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


#' One Product Stencil for a Mixed Partial Derivative
#'
#' @description
#' Differentiates `f` once in each component an index tuple names, using a
#' central factor per distinct component of the order that component's
#' multiplicity asks for. The whole tensor product is summed in one pass, so the
#' result is a single stencil, never a composition of lower-order numerical
#' derivatives.
#'
#' @details
#' Each factor's nodes and weights come from [numericals7::fd_offsets()] and
#' [numericals7::fd_weights()] at accuracy 2: two points at order one, three at
#' order two, five at orders three and four. They are read from there instead of
#' being written out here, because a table of stencil coefficients kept in a
#' second place is a table that can come to disagree with the first.
#'
#' The step for each factor is [fd_step()] at that component's value and at the
#' **total** order of the tuple, so all factors of one component share a step and
#' the division at the end is by \eqn{\prod_j h_j^{m_j}}. Nodes whose weight is
#' zero are skipped, so a repeated-index factor costs one evaluation fewer than
#' its node count suggests.
#'
#' @param f A function of the free vector, returning anything that supports `+`
#'   and scalar multiplication.
#' @param eta The point, a numeric vector.
#' @param tuple An integer vector of component indices, possibly with
#'   repetitions, as [param_tuple_indices()] returns. Its length is the
#'   derivative order.
#'
#' @return The stencil's value, shaped like `f(eta)`.
#'
#' @seealso [numericals7::fd_weights()] for the weights, [fd_along()] for the
#'   single-component case, and [numerical_d3()] and [numerical_d4()], the two
#'   callers.
#'
#' @keywords internal
mixed_stencil <- function(f, eta, tuple) {
  counts <- table(tuple)
  ks <- as.integer(names(counts))
  ms <- as.integer(counts)
  hs <- vapply(seq_along(ks), function(j) {
    fd_step(eta[ks[j]], sum(ms))
  }, numeric(1))

  one_dim <- function(m) {
    s <- numericals7::fd_offsets(m, accuracy = 2L)$central
    w <- numericals7::fd_weights(s, m)
    keep <- w != 0
    list(offsets = s[keep], weights = w[keep], power = m)
  }
  facs <- lapply(ms, one_dim)

  combos <- expand.grid(lapply(facs, function(fc) seq_along(fc$offsets)))
  acc <- NULL
  for (r in seq_len(nrow(combos))) {
    e <- eta
    w <- 1
    for (j in seq_along(ks)) {
      pick <- combos[r, j]
      e[ks[j]] <- e[ks[j]] + facs[[j]]$offsets[pick] * hs[j]
      w <- w * facs[[j]]$weights[pick]
    }
    term <- w * f(e)
    acc <- if (is.null(acc)) term else acc + term
  }
  denom <- prod(vapply(seq_along(ks), function(j) {
    hs[j]^facs[[j]]$power
  }, numeric(1)))
  acc / denom
}


#' @title Default Third Derivatives
#' @name param_d3.parameter
#' @description Fallback: see [numerical_d3()].
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list.
#' @keywords internal
S7::method(param_d3, parameter) <- function(s, eta, ...) {
  numerical_d3(s, eta)
}

#' @title Default Fourth Derivatives
#' @name param_d4.parameter
#' @description Fallback: see [numerical_d4()].
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list.
#' @keywords internal
S7::method(param_d4, parameter) <- function(s, eta, ...) {
  numerical_d4(s, eta)
}

#' @title Default Higher Log-Determinant Derivatives
#' @name param_d3logdet.matrix_parameter
#' @description
#' Fallback: one central stencil on [param_d2logdet()], which is
#' the exact trace identity given the matrix derivatives -- a single layer on
#' an analytic quantity, per the toolkit's rule.
#' @param s A [matrix_parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d3logdet, matrix_parameter) <- function(s, eta, ...) {
  idx3 <- param_tuple_indices(s, 3L)
  idx2 <- param_tuple_indices(s, 2L)
  key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
  out <- vapply(seq_along(idx3), function(i) {
    t3 <- idx3[[i]]
    k <- t3[3L]
    pos <- match(paste(sort(t3[1:2]), collapse = ","), key2)
    h <- fd_step(eta[k], 1L)
    up <- eta
    dn <- eta
    up[k] <- eta[k] + h
    dn[k] <- eta[k] - h
    (param_d2logdet(s, up)[[pos]] - param_d2logdet(s, dn)[[pos]]) / (2 * h)
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, 3L))
}

#' @title Default Fourth Log-Determinant Derivatives
#' @name param_d4logdet.matrix_parameter
#' @description
#' Fallback: one second-order stencil on [param_d2logdet()] in the
#' last two components of the tuple -- one layer, mixed across components
#' where they differ.
#' @param s A [matrix_parameter()] object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d4logdet, matrix_parameter) <- function(s, eta, ...) {
  idx4 <- param_tuple_indices(s, 4L)
  idx2 <- param_tuple_indices(s, 2L)
  key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
  out <- vapply(seq_along(idx4), function(i) {
    t4 <- idx4[[i]]
    pos <- match(paste(sort(t4[1:2]), collapse = ","), key2)
    k <- t4[3L]
    l <- t4[4L]
    g <- function(e) param_d2logdet(s, e)[[pos]]
    if (k == l) {
      h <- fd_step(eta[k], 2L)
      up <- eta
      dn <- eta
      up[k] <- eta[k] + h
      dn[k] <- eta[k] - h
      (g(up) - 2 * g(eta) + g(dn)) / h^2
    } else {
      hk <- fd_step(eta[k], 2L)
      hl <- fd_step(eta[l], 2L)
      pp <- pm <- mp <- mm <- eta
      pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
      pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
      mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
      mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
      (g(pp) - g(pm) - g(mp) + g(mm)) / (4 * hk * hl)
    }
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, 4L))
}
