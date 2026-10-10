#' @include generics.R
NULL

# What a parameter gets for free when it implements only param_value().
# Every method here is registered on the base class, so a subclass that
# supplies a closed form takes over through dispatch with no registration step.


#' Test for the Package's Own Base Classes
#'
#' @description
#' Returns whether an S7 class is [parameter()] or [matrix_parameter()], the
#' two base classes on which this package registers its fallback methods. The
#' test distinguishes a method that a subclass wrote from one that it inherited,
#' which is what [param_is_numerical()] reports and what [check_parameter()]
#' acts on.
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
#'   without an error, `attr()` on a non-class giving `NULL`.
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
#' method of its own or takes the one registered on [parameter()] and
#' [matrix_parameter()]. It is the test to run before relying on a fourth-order
#' derivative of a family written elsewhere: `FALSE` everywhere means that every
#' quantity has a method of the family's own, and `TRUE` marks a component
#' computed by a base method, a stencil for the derivatives of the value, an
#' eigendecomposition for the log-determinant and its first two derivatives, and
#' a stencil on the order below for its third and fourth.
#'
#' @details
#' # Why the answer matters
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
#' it has nothing at all, and [param_free()] because its base method signals an
#' error.
#'
#' # Every family in the package returns FALSE
#'
#' Every constructor in this package returns `FALSE` at every component. The
#' fallbacks serve a family written elsewhere, and the example below builds one
#' to show what a `TRUE` looks like.
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
#'   checked, and [numerical_d1()] and its higher-order siblings, the routes that
#'   a `TRUE` names.
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
#'            null_basis = matrix(numeric(0), 2, 0))
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


#' The Spectral Decomposition of a Parameter's Matrix
#'
#' @description
#' Returns the eigenvalues and eigenvectors of \eqn{M(\eta)}, together with a
#' flag marking the directions that carry the matrix. The base-class
#' log-determinant and its first two derivatives read it, so one decomposition
#' serves them.
#'
#' @details
#' Which eigenvalues count as zero is settled **by position**, never by size.
#' `eigen()` returns them in decreasing order, so the first `s@rank` are kept and
#' the rest dropped, whatever their numerical values are. The rank itself comes
#' from the object and is never re-derived: counting eigenvalues above a relative
#' tolerance is not scale invariant, so a family whose components differ by many
#' orders of magnitude would be assigned a different rank at different \eqn{\eta},
#' and a fitted model with smoothing parameters that far apart is ordinary. The
#' rank is fixed once, at construction, from the components; see
#' [param_null_basis()].
#'
#' The matrix is symmetrized as `(m + t(m)) / 2` before the decomposition, so
#' both triangles of an asymmetry of rounding size enter the eigenvalues.
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
#' directions that the `keep` flag of [param_spectrum()] marks. The base-class
#' log-determinant derivatives use it in place of \eqn{M^{-1}} in the trace
#' identities, and for a full-rank family it is the ordinary inverse.
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
#' The verdict comes from the spectrum because `chol()` is not a rank test. On
#' a matrix with an exactly zero eigenvalue the pivot that ought to be zero
#' comes out positive or negative according to rounding, so `chol()` succeeds on
#' some platforms and fails on others. The test `min(ev) <= tol * max(ev)` does
#' not depend on the platform in that way.
#'
#' The default tolerance, \eqn{10^{-14}} (about 45 times the machine epsilon),
#' is above the rounding level of an eigenvalue that is exactly zero and below
#' the conditioning that [param_value()] reaches on ordinary free vectors, so a
#' badly conditioned but positive definite value is accepted.
#'
#' @param m A symmetric numeric matrix. Not checked for symmetry; the callers
#'   pass the value of [param_value()] or a matrix already symmetrized by
#'   [check_matrix()].
#' @param tol The relative tolerance below which the smallest eigenvalue counts
#'   as zero, so `m` is rejected when `min(ev) <= tol * max(ev)`. Defaults to
#'   `1e-14`. At that default `diag(c(1, 1e-14))` returns `NULL`, the test being
#'   an inequality, and `diag(c(1, 1e-13))` factors.
#'
#' @return The lower triangular \eqn{L} with \eqn{M = L L^\top}, or `NULL` when
#'   `m` is empty, has a missing or non-finite entry, has a non-positive largest
#'   eigenvalue, fails the relative test, or makes [base::chol()] signal an
#'   error.
#'
#' @seealso [param_factor()], the generic whose base method calls this, and
#'   [param_spectrum()] for the decomposition that the log-determinant uses.
#'
#' @keywords internal
chol_pd <- function(m, tol = 1e-14) {
  if (!length(m) || !all(is.finite(m))) return(NULL)
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  if (max(ev) <= 0 || min(ev) <= tol * max(ev)) return(NULL)
  r <- tryCatch(chol(m), error = function(e) NULL)
  if (is.null(r)) NULL else t(r)
}


#' Finite-Difference Step for a Free Value
#'
#' @description
#' Returns the step for a central difference in one component of \eqn{\eta},
#' scaled by the size of that component. It forwards to
#' [numericals7::fd_step()] at accuracy 2, so the step of a fallback here and the
#' step that the stencil library documents cannot drift apart.
#'
#' @details
#' The rule is \eqn{\varepsilon^{1/(k+2)} \max(1, |\eta_k|)} for a \eqn{k}-th
#' derivative, which balances the truncation error against the rounding that
#' the division by \eqn{h^k} amplifies. In doubles that is \eqn{6.1 \times 10^{-6}}
#' at first order and \eqn{1.2 \times 10^{-4}} at second, both scaling up once
#' \eqn{|\eta_k|} passes 1.
#'
#' The step is not clamped. The unconstrained scale has no boundary, so the
#' `bounds` argument of [numericals7::fd_step()] is left at its default.
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
#' non-zero weight: two, three, four and five times at orders one to four. The accumulator starts at `NULL` and takes the shape of the
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
#' at the step that [fd_step()] gives for a first derivative. This is the method
#' [param_d1()] dispatches to when a family has not written its own, and it is
#' exported so that a family being developed can be compared against it.
#'
#' @details
#' The step is \eqn{\varepsilon^{1/3}\max(1, |\eta_k|)}, about
#' \eqn{6.1 \times 10^{-6}} near the origin, giving a truncation error of order
#' \eqn{h^2}. The example below compares the result with a closed form.
#'
#' For a matrix family each estimate is symmetrized as `(m + t(m)) / 2`, as
#' the higher orders are, so that an asymmetry of rounding size in
#' [param_value()] does not reach the derivative. It costs \eqn{2d}
#' evaluations of the map.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `s@n_free` estimates named by `s@free_names`, each shaped
#'   like the result of [param_value()] and symmetrized for a matrix family.
#'
#' @seealso [param_d1()], the generic that this serves, [numerical_d2()],
#'   [numerical_d3()] and [numerical_d4()] for the higher orders, and
#'   [param_is_numerical()] to find out whether a given family reaches this.
#'
#' @examples
#' # A scalar matrix is exp(eta) times the identity, so the derivative is the
#' # matrix itself and the error of the difference can be read off.
#' s <- scalar_matrix(2)
#' numerical_d1(s, 0.3)
#' max(abs(numerical_d1(s, 0.3)[[1]] - param_value(s, 0.3)))
#'
#' # Against a family that writes its own.
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
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[k]] <- (out[[k]] + t(out[[k]])) / 2
    }
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
#' # The three routes
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
#' # Single-layer differences
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
#' On the analytic route the step is the order-1 one,
#' \eqn{\varepsilon^{1/3}\max(1, |\eta_l|)}, about \eqn{6.1 \times 10^{-6}} near
#' the origin, because a first derivative of an analytic array is taken. The two
#' routes on [param_value()] use the order-2 step,
#' \eqn{\varepsilon^{1/4}\max(1, |\eta_k|)}, about \eqn{1.2 \times 10^{-4}}.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `choose(s@n_free + 1, 2)` estimates keyed as
#'   `param_tuple_names(s)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d2()], the generic that this serves, [numerical_d1()] for the
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
#' @description
#' The method that every [parameter()] inherits when it registers no [param_d1()] of
#' its own. It estimates \eqn{\partial V/\partial \eta_k} by one three-point
#' central difference of [param_value()] in each component,
#' \eqn{(V(\eta + h e_k) - V(\eta - h e_k)) / 2h}, at the step
#' \eqn{\varepsilon^{1/3}\max(1, |\eta_k|)}, which is about
#' \eqn{6.1 \times 10^{-6}} near the origin. It costs \eqn{2d} evaluations of the
#' map. The families in this package do not reach it.
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` estimates named by `s@free_names`, each shaped
#'   like the result of [param_value()] and symmetrized for a matrix family.
#' @seealso [numerical_d1()], which does the work, and [param_is_numerical()]
#'   to find out whether a family reaches this method.
#' @keywords internal
S7::method(param_d1, parameter) <- function(s, eta, ...) {
  numerical_d1(s, eta)
}

#' @title Default Second Derivatives
#' @name param_d2.parameter
#' @description
#' The method that every [parameter()] inherits when it registers no [param_d2()] of
#' its own. It takes exactly one difference per component, of whichever quantity
#' the family already supplies: the analytic [param_d1()] where there is one, at
#' the order-1 step, and [param_value()] itself where there is not, through a
#' three-point second difference on the diagonal and a four-point mixed stencil
#' off it, at the order-2 step. The truncation error is of order \eqn{h^2} on
#' every route. The families in this package do not reach it.
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` estimates keyed as
#'   `param_tuple_names(s)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#' @seealso [numerical_d2()], which does the work and writes the three routes
#'   out, and [param_d1.parameter()] for the order below.
#' @keywords internal
S7::method(param_d2, parameter) <- function(s, eta, ...) {
  numerical_d2(s, eta)
}

#' @title Default Log-Determinant
#' @name param_logdet.matrix_parameter
#' @description
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_logdet()] of its own. It takes an eigendecomposition of
#' [param_value()], keeps the first `s@rank` eigenvalues, and returns the sum of
#' their logarithms: the log-determinant for a full-rank family and the log
#' pseudo-determinant otherwise. Which eigenvalues are kept is decided by
#' position, from the declared rank, and never from their size.
#'
#' The result is exact up to rounding, so [param_is_numerical()] reporting
#' `TRUE` here means that the answer costs \eqn{O(p^3)} and cannot be checked
#' against an eigendecomposition; it does not mean that the answer is
#' inaccurate.
#' @details
#' A non-positive eigenvalue among those that the rank keeps signals an error
#' naming the family and the counts: the family has declared a rank it does not have at this
#' \eqn{\eta}, so `log()` of a non-positive number would be the wrong thing to
#' return. This is the check that catches a `param_value()` method whose matrix
#' leaves the positive semidefinite cone.
#' @param s A [matrix_parameter()] object, whose `rank` and `dimension` are read.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_logdet()] for the generic, [param_spectrum()] for the
#'   decomposition, and [param_dlogdet.matrix_parameter()] for its derivative.
#' @keywords internal
S7::method(param_logdet, matrix_parameter) <- function(s, eta, ...) {
  sp <- param_spectrum(s, eta)
  v <- sp$values[sp$keep]
  if (any(v <= 0)) {
    stop(sprintf(paste0(
      "'%s' produced %d non-positive eigenvalue(s) among the %d that its\n",
      "  rank keeps, so it is outside the set that it parametrizes."
    ), s@param_name, sum(v <= 0), s@rank), call. = FALSE)
  }
  sum(log(v))
}

#' @title Default Log-Determinant Gradient
#' @name param_dlogdet.matrix_parameter
#' @description
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_dlogdet()] of its own. It evaluates the trace identity
#'
#' \deqn{\partial_k \log|M| = \mathrm{tr}\!\left(M^{+} \partial_k M\right),}
#'
#' with \eqn{M^{+}} the Moore-Penrose inverse formed from the directions that
#' the declared rank keeps, which is the ordinary inverse for a full-rank family. The
#' trace is computed as `sum(mi * dk)`, the elementwise product summed, both
#' matrices being symmetric, so no matrix product is formed.
#'
#' The identity is exact, so the accuracy is entirely the accuracy of
#' [param_d1()]: rounding with an analytic first derivative, and that of a
#' central difference with a numerical one.
#' @param s A [matrix_parameter()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_dlogdet()] for the generic, [param_d1()] for the derivative
#'   arrays this reads, and [spectrum_pinv()] for \eqn{M^{+}}.
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
#' @name param_d2logdet.matrix_parameter
#' @description
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_d2logdet()] of its own. It evaluates the identity that follows from
#' differentiating [param_dlogdet()]'s trace once more, using
#' \eqn{\partial_l M^{-1} = -M^{-1}(\partial_l M)M^{-1}},
#'
#' \deqn{\partial_{kl} \log|M| = \mathrm{tr}\!\left(M^{+} \partial_{kl} M\right)
#'   - \mathrm{tr}\!\left(M^{+} (\partial_k M)\, M^{+} (\partial_l M)\right),}
#'
#' with \eqn{M^{+}} the pseudo-inverse over the directions that the declared
#' rank keeps. The second trace is formed as `sum(t(mi %*% dk) * (mi %*% dl))`, which
#' is the trace of the product without the product being multiplied out.
#'
#' Exact given the derivative arrays, so the accuracy is theirs: rounding with
#' analytic arrays, and that of the stencils with numerical ones.
#' @param s A [matrix_parameter()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
#'   `param_tuple_names(s)` and in that order.
#' @seealso [param_d2logdet()] for the generic, [param_d1()] and [param_d2()] for
#'   the arrays this reads, and [param_d3logdet.matrix_parameter()] for the order
#'   above, which differences this one.
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

#' Reject an Order That Would Difference a Difference
#'
#' @description
#' Signals an error when [param_d3logdet.matrix_parameter()] or
#' [param_d4logdet.matrix_parameter()] is called for a family whose
#' [param_d2logdet()] is itself numerical. Both fallbacks difference
#' [param_d2logdet()]. A family's own [param_d2logdet()] is exact, and so is
#' the base-class one **given** [param_d1()] and [param_d2()]; in either case
#' the differencing is one layer. Where the base-class [param_d2logdet()] reads
#' numerical arrays it is two, which is the nesting that the toolkit avoids
#' everywhere.
#'
#' @details
#' With numerical arrays, the fourth order can come back larger than the
#' quantity that it estimates, and [param_is_numerical()] reports `TRUE` for the
#' order whether the arrays are analytic or not, so a consumer could not
#' distinguish such a number from a usable one.
#'
#' A family that writes its own [param_d2logdet()] passes whether or not its
#' [param_d1()] and [param_d2()] are analytic.
#'
#' @param s A [matrix_parameter()] object.
#' @param order The order requested, 3 or 4, which the message names.
#'
#' @return Invisibly `TRUE`. For a family whose [param_d2logdet()] comes from
#'   the base class and whose [param_d1()] or [param_d2()] comes from the base
#'   class too, an error that names the missing method and the remedy.
#'
#' @seealso [param_is_numerical()], which reports the routes read here, and
#'   [param_d3logdet.matrix_parameter()] and [param_d4logdet.matrix_parameter()],
#'   the two callers.
#'
#' @keywords internal
check_analytic_arrays <- function(s, order) {
  num <- param_is_numerical(s)
  if (!num[["param_d2logdet"]]) return(invisible(TRUE))
  miss <- c("param_d1()", "param_d2()")[c(num[["param_d1"]], num[["param_d2"]])]
  if (!length(miss)) return(invisible(TRUE))
  stop(sprintf(paste0(
    "param_d%dlogdet() would difference a quantity that is itself\n",
    "  numerical. This family has no analytic %s, so\n",
    "  param_d2logdet(), which this order differences, is already one\n",
    "  difference deep, and the two layers would compound.\n",
    "  Write param_d1() and param_d2() out, or register param_d2logdet()\n",
    "  or param_d%dlogdet()."
  ), order, paste(miss, collapse = " and "), order), call. = FALSE)
}

#' @title Default Solve
#' @name param_solve.matrix_parameter
#' @description
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_solve()] of its own. It takes the Cholesky factor \eqn{L} from
#' [param_factor()] and applies it twice,
#' `backsolve(t(l), forwardsolve(l, b))`, which is
#' \eqn{L^{-\top}L^{-1}B = M^{-1}B}. The inverse is not formed, and each
#' triangular system is solved once.
#'
#' The result is exact up to rounding, as `base::solve()` is, which is why
#' [param_is_numerical()] does not list it. The generic has already rejected a rank-deficient family and one whose
#' value is not a symmetric matrix, and filled `b` with the identity when the
#' caller left it out, so by the time this runs `b` is a matrix of the right
#' height.
#' @details
#' Positive definiteness is decided inside [param_factor()], from the
#' eigenvalues, before any factorization is attempted; see [chol_pd()] for why
#' the verdict does not come from whether [base::chol()] signals an error.
#'
#' The cost is one \eqn{O(p^3)} factorization plus \eqn{O(p^2)} per column of
#' `b`. The families that override this method are listed under
#' [param_solve()].
#' @param s A [matrix_parameter()] object, of full rank.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param b A numeric matrix with `s@dimension` rows, already coerced from a
#'   vector and defaulted to the identity by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric matrix with `s@dimension` rows and as many columns as `b`.
#' @seealso [param_solve()] for the generic and the argument handling,
#'   [param_factor.matrix_parameter()] for the factor, and [chol_pd()] for the
#'   definiteness test.
#' @keywords internal
S7::method(param_solve, matrix_parameter) <- function(s, eta, b = NULL, ...) {
  l <- param_factor(s, eta)
  backsolve(t(l), forwardsolve(l, b))
}

#' @title Default Factor
#' @name param_factor.matrix_parameter
#' @description
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_factor()] of its own. It evaluates [param_value()] and returns the
#' lower triangular Cholesky factor through [chol_pd()], which decides positive
#' definiteness from the eigenvalues before attempting the factorization. Exact,
#' at \eqn{O(p^3)}.
#' @details
#' A family that declares full rank and is then not positive definite at this
#' \eqn{\eta} signals an error naming the family. The verdict comes from the
#' eigenvalues and not from a caught `chol()` error, which could differ between
#' platforms on a matrix with an exactly zero eigenvalue. The error means that
#' the family's own [param_value()] has left the cone that the family
#' parametrizes.
#' @param s A [matrix_parameter()] object, of full rank.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` lower triangular numeric matrix with
#'   a positive diagonal, satisfying `L %*% t(L) == param_value(s, eta)`.
#' @seealso [param_factor()] for the generic, which rejects a rank-deficient
#'   family before this runs, and [chol_pd()] for the definiteness test.
#' @keywords internal
S7::method(param_factor, matrix_parameter) <- function(s, eta, ...) {
  l <- chol_pd(param_value(s, eta))
  if (is.null(l)) {
    stop(sprintf(paste0(
      "'%s' is not positive definite at this eta, although it declares full\n",
      "  rank. The test is on the eigenvalues of the matrix."
    ), s@param_name), call. = FALSE)
  }
  l
}

#' @title Default Inverse Map
#' @name param_free.parameter
#' @description
#' The method that every [parameter()] inherits when it registers no [param_free()] of
#' its own. It always signals an error naming the family; it is the one generic
#' whose base method rejects the call instead of computing. An inverse obtained
#' by minimizing \eqn{\lVert V(\eta) - m \rVert} would return an \eqn{\eta} even
#' for a matrix that lies outside the family's set, and the caller could not
#' distinguish that result from a correct one, so the inverse map is either
#' written out exactly or rejected.
#'
#' Every family in this package writes its inverse out, so this method is
#' reached only by a family defined elsewhere.
#' @param s A [parameter()] object, whose `param_name` goes into the message.
#' @param m A value of the family's shape. Never read.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return Never returns; always signals an error.
#' @seealso [param_free()] for the generic and for what each family's own method
#'   rejects, and [param_value()] for the forward map.
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
#' directly to [param_value()], per index tuple. A component that the tuple repeats
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
#' # The stencil is applied to the map
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
#' \eqn{\varepsilon^{1/5}\max(1, |\eta_k|)}, and rounding is of order
#' \eqn{\varepsilon / h^3}.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `choose(s@n_free + 2, 3)` estimates keyed as
#'   `param_tuple_names(s, 3)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d3()], the generic that this serves, [mixed_stencil()], which
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
#' not written its own. At this order a difference is the least accurate of the
#' four, which is why every family in this package carries a closed form
#' instead.
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
#' the origin, and an attainable accuracy of order \eqn{\varepsilon^{1/3}}. That
#' is enough to catch a transcription error in a closed form and not enough for
#' use in a fit.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `choose(s@n_free + 3, 4)` estimates keyed as
#'   `param_tuple_names(s, 4)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d4()], the generic that this serves, [numerical_d3()] for the
#'   order below, and [check_parameter()], which uses this stencil as the
#'   reference for a family that is not a matrix.
#'
#' @examples
#' # A scalar matrix: every order is the matrix again, so the error shows.
#' s <- scalar_matrix(2)
#' max(abs(numerical_d4(s, 0.3)[[1]] - param_value(s, 0.3)))
#'
#' # Against a family that writes its own.
#' q <- log_cholesky(2)
#' eta <- c(0.2, -0.1, 0.4)
#' ana <- param_d4(q, eta)
#' c(gap = max(abs(unlist(numerical_d4(q, eta)) - unlist(ana))),
#'   scale = max(abs(unlist(ana))))
#'
#' # A 1 percent error in the largest entry would be about two thousand times
#' # the gap above, so a wrong closed form is caught.
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
#' Differentiates `f` once in each component that an index tuple names, using a
#' central factor per distinct component, of the order given by that
#' component's multiplicity. The whole tensor product is summed in one pass, so the
#' result is a single stencil, never a composition of lower-order numerical
#' derivatives.
#'
#' @details
#' Each factor's nodes and weights come from [numericals7::fd_offsets()] and
#' [numericals7::fd_weights()] at accuracy 2: three nodes at orders one and two
#' and five at orders three and four, of which two, three, four and five carry a
#' non-zero weight. They are read from there instead of being written out here,
#' so that the two tables cannot disagree.
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
#' @description
#' The method that every [parameter()] inherits when it registers no [param_d3()] of
#' its own. It applies one product stencil per index tuple directly to
#' [param_value()], never to a lower-order numerical derivative, with a
#' one-dimensional factor per distinct component, of the order given by that
#' component's multiplicity. The step is \eqn{\varepsilon^{1/5}\max(1,
#' |\eta_k|)}, about \eqn{7.4 \times 10^{-4}} near the origin, and the
#' truncation error is of order \eqn{h^2}. The families in this package do not
#' reach it.
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` estimates keyed as
#'   `param_tuple_names(s, 3)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#' @seealso [numerical_d3()], which does the work and writes the stencil out, and
#'   [param_d4.parameter()] for the order above.
#' @keywords internal
S7::method(param_d3, parameter) <- function(s, eta, ...) {
  numerical_d3(s, eta)
}

#' @title Default Fourth Derivatives
#' @name param_d4.parameter
#' @description
#' The method that every [parameter()] inherits when it registers no [param_d4()] of
#' its own. It applies one product stencil per index tuple directly to
#' [param_value()], as [param_d3.parameter()] does, with the multiplicities
#' summing to four. The step is \eqn{\varepsilon^{1/6}\max(1, |\eta_k|)}, about
#' \eqn{2.5 \times 10^{-3}} near the origin, and rounding is amplified by
#' \eqn{h^{-4}}, so the result is enough to catch a wrong closed form and not
#' enough for use in a fit. The families in this package do not reach it.
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` estimates keyed as
#'   `param_tuple_names(s, 4)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#' @seealso [numerical_d4()], which does the work, and [param_d3.parameter()] for
#'   the order below.
#' @keywords internal
S7::method(param_d4, parameter) <- function(s, eta, ...) {
  numerical_d4(s, eta)
}

#' @title Default Third Log-Determinant Derivatives
#' @name param_d3logdet.matrix_parameter
#' @description
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_d3logdet()] of its own. For the tuple \eqn{(k, l, m)} it takes the
#' \eqn{(k, l)} component of [param_d2logdet()] and applies one three-point
#' central difference in the remaining component,
#'
#' \deqn{\partial_{klm} \log|M| \approx
#'   \frac{\partial_{kl}\log|M|\,(\eta + h e_m)
#'       - \partial_{kl}\log|M|\,(\eta - h e_m)}{2h},}
#'
#' at the order-1 step \eqn{\varepsilon^{1/3}\max(1, |\eta_m|)}. Since the tuple
#' is sorted, the differenced component is the largest index, and the pair is
#' looked up by matching sorted indices against
#' [param_tuple_indices()]'s order-2 list. It costs two evaluations of the whole
#' second-order block per tuple.
#' @details
#' # Accuracy
#'
#' [param_d2logdet()] is an exact identity **given** the matrix derivative
#' arrays, so differencing it is a single numerical layer whenever
#' [param_d1()] and [param_d2()] are analytic.
#'
#' A family's own [param_d2logdet()] is exact, and differencing it is a single
#' layer too. Where the base-class [param_d2logdet()] reads numerical arrays
#' the layers would compound, and the method signals an error instead of
#' returning a value; [check_analytic_arrays()] is the guard. A family that
#' needs this order writes [param_d1()] and [param_d2()] out, or its own
#' [param_d2logdet()].
#' @param s A [matrix_parameter()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
#'   `param_tuple_names(s, 3)` and in that order.
#' @seealso [param_d3logdet()] for the generic, [param_d2logdet()] for the
#'   quantity differenced, and [param_d4logdet.matrix_parameter()] for the order
#'   above.
#' @keywords internal
S7::method(param_d3logdet, matrix_parameter) <- function(s, eta, ...) {
  check_analytic_arrays(s, 3L)
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
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_d4logdet()] of its own. For the tuple \eqn{(k, l, m, n)} it takes the
#' \eqn{(k, l)} component of [param_d2logdet()] and applies one second-order
#' stencil in the remaining two components: the three-point second difference
#' where \eqn{m = n}, and the four-point mixed stencil
#'
#' \deqn{\frac{g(\eta + h_m e_m + h_n e_n) - g(\eta + h_m e_m - h_n e_n)
#'   - g(\eta - h_m e_m + h_n e_n) + g(\eta - h_m e_m - h_n e_n)}
#'   {4 h_m h_n}}
#'
#' where they differ, both at the order-2 step
#' \eqn{\varepsilon^{1/4}\max(1, |\eta_k|)}. It costs three or four evaluations
#' of the whole second-order block per tuple.
#' @details
#' # Accuracy
#'
#' With analytic [param_d1()] and [param_d2()], or with a family's own
#' [param_d2logdet()], the differencing is a single layer on an exact
#' quantity, at the order-2 step.
#'
#' With only [param_value()] supplied, the numerical arrays would feed a
#' numerical second-order block that is then differenced twice more, and the
#' error can exceed the quantity itself. The method signals an error instead,
#' through [check_analytic_arrays()]. A family that needs a fourth derivative of
#' its log-determinant supplies at least [param_d1()] and [param_d2()] in closed
#' form, or its own [param_d2logdet()]. [param_is_numerical()] reports which components are numerical.
#' @param s A [matrix_parameter()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
#'   `param_tuple_names(s, 4)` and in that order.
#' @seealso [param_d4logdet()] for the generic,
#'   [param_d3logdet.matrix_parameter()] for the order below, and
#'   [param_is_numerical()] to find out which components are numerical.
#' @keywords internal
S7::method(param_d4logdet, matrix_parameter) <- function(s, eta, ...) {
  check_analytic_arrays(s, 4L)
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
