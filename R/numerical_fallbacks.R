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
#' [numericals7::fd_step()] at accuracy 2, so the step of a difference written
#' here and the step that the stencil library documents cannot drift apart.
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
#' @seealso [numericals7::fd_step()] for the rule, and [check_parameter()],
#'   which uses it for its log-determinant reference.
#'
#' @keywords internal
fd_step <- function(eta_k, order = 1L) {
  numericals7::fd_step(eta_k, order = order, accuracy = 2L)
}


#' Which Derivatives of a Parameter Are Analytic
#'
#' @description
#' Returns, for each derivative generic that [param_is_numerical()] reports,
#' `TRUE` where the quantity is computed in closed form. It differs from
#' `!param_is_numerical(s)` in two ways. The base-class [param_dlogdet()] and
#' [param_d2logdet()] are identities in the derivative arrays, so they count
#' as analytic where [param_d1()], and [param_d2()] for the second, are. And a
#' composite family, which registers methods of its own that are built from
#' the parameters held in its `param_params`, counts a quantity as analytic
#' only where every such parameter's quantity is analytic too.
#'
#' @param s A [parameter()] object.
#'
#' @return A named logical vector with the names of [param_is_numerical()].
#'
#' @seealso [param_tensor()], which declares the analytic quantities to the
#'   stencil, and [param_children()].
#'
#' @keywords internal
param_analytic <- function(s) {
  a <- !param_is_numerical(s)
  if ("param_d2logdet" %in% names(a)) {
    if (!a[["param_dlogdet"]]) a[["param_dlogdet"]] <- a[["param_d1"]]
    if (!a[["param_d2logdet"]]) {
      a[["param_d2logdet"]] <- a[["param_d1"]] && a[["param_d2"]]
    }
  }
  for (k in param_children(s)) {
    ak <- param_analytic(k)
    common <- intersect(names(a), names(ak))
    a[common] <- a[common] & ak[common]
  }
  a
}


#' The Parameters Held Inside a Parameter
#'
#' @description
#' Returns the [parameter()] objects found in `s@param_params`, at any depth
#' of nesting in lists: the blocks of a [block_diag()], the inner parameter of
#' a [kron_identity()] or an [inverse_of()], the terms of a [sum_struct()].
#'
#' @param s A [parameter()] object.
#'
#' @return A list of [parameter()] objects, empty for a family that holds
#'   none.
#'
#' @seealso [param_analytic()], which reads them.
#'
#' @keywords internal
param_children <- function(s) {
  walk <- function(x) {
    if (S7::S7_inherits(x, parameter)) return(list(x))
    if (is.list(x)) return(unlist(lapply(x, walk), recursive = FALSE))
    list()
  }
  walk(s@param_params)
}


#' Derivatives of a Parameter by One Tensor Stencil per Component
#'
#' @description
#' Computes the derivatives of order `order` of [param_value()] or, with
#' `what = "logdet"`, of [param_logdet()], through
#' [numericals7::fd_tensor()]. The analytic derivatives of lower order, as
#' [param_analytic()] reads them, are declared to it, so each component is one
#' tensor-product central difference of the highest analytic quantity below
#' it, and of the value itself where there is none.
#'
#' @details
#' The coordinates are the free values. The value is flattened to a vector,
#' and each component is returned with the attributes of [param_value()]'s
#' result and, for a [matrix_parameter()], symmetrized as `(m + t(m)) / 2`.
#' The step of a component whose difference is of order \eqn{r} is
#' \eqn{\varepsilon^{1/(r + a)}\max(1, |\eta_k|)}, with \eqn{a} the accuracy;
#' the unconstrained scale has no bounds.
#'
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order, 1 to 4.
#' @param what `"value"` or `"logdet"`.
#' @param accuracy The order of accuracy of the one-dimensional stencils.
#'
#' @return For `"value"`, a list of estimates keyed as
#'   `param_tuple_names(s, order)` and in that order, each shaped like
#'   [param_value()]'s result. For `"logdet"`, a named numeric vector in the
#'   same order.
#'
#' @seealso [numerical_d1()] to [numerical_d4()],
#'   [param_d3logdet.matrix_parameter()] and
#'   [param_d4logdet.matrix_parameter()], which call this.
#'
#' @keywords internal
param_tensor <- function(s, eta, order, what = c("value", "logdet"),
                         accuracy = 4L) {
  what <- match.arg(what)
  nms <- param_tuple_names(s, order)
  if (!length(nms)) {
    if (what == "logdet") return(stats::setNames(numeric(0), character(0)))
    return(stats::setNames(vector("list", 0L), character(0)))
  }
  a <- param_analytic(s)
  x <- stats::setNames(as.list(eta), paste0("e", seq_len(s@n_free)))
  free <- function(x) unlist(x, use.names = FALSE)
  if (what == "value") {
    shape <- param_value(s, eta)
    f <- function(x) as.vector(param_value(s, free(x)))
    gens <- list(param_d1 = param_d1, param_d2 = param_d2, param_d3 = param_d3)
    wrap <- function(g) {
      force(g)
      function(x) lapply(g(s, free(x)), as.vector)
    }
  } else {
    f <- function(x) param_logdet(s, free(x))
    gens <- list(param_dlogdet = param_dlogdet, param_d2logdet = param_d2logdet,
                 param_d3logdet = param_d3logdet)
    wrap <- function(g) {
      force(g)
      function(x) as.list(unname(g(s, free(x))))
    }
  }
  known <- list()
  for (k in seq_len(min(order - 1L, 3L))) {
    if (a[[names(gens)[k]]]) {
      known[[length(known) + 1L]] <- numericals7::fd_known(wrap(gens[[k]]),
                                                           order = k)
    }
  }
  got <- numericals7::fd_tensor(f, x, order = order, known = known,
                                accuracy = accuracy)
  if (what == "logdet") {
    return(stats::setNames(vapply(got, function(v) v[1L], numeric(1)), nms))
  }
  sym <- S7::S7_inherits(s, matrix_parameter)
  out <- lapply(got, function(v) {
    m <- shape
    m[] <- v
    if (sym) (m + t(m)) / 2 else m
  })
  stats::setNames(out, nms)
}


#' Numerical First Derivatives of a Parameter's Matrix
#'
#' @description
#' Estimates \eqn{\partial V / \partial \eta_k} for every free value by one
#' central difference of [param_value()] in each component, through
#' [param_tensor()]. This is the method [param_d1()] dispatches to when a family
#' has not written its own, and it is exported so that a family being developed
#' can be compared against it.
#'
#' @details
#' At the default accuracy 4 the difference is the five-point rule
#' \eqn{(V(\eta - 2h e_k) - 8V(\eta - h e_k) + 8V(\eta + h e_k)
#' - V(\eta + 2h e_k)) / 12h}, at the step
#' \eqn{\varepsilon^{1/5}\max(1, |\eta_k|)}, about \eqn{7.4 \times 10^{-4}}
#' near the origin, with a truncation error of order \eqn{h^4}. It costs
#' \eqn{4d} evaluations of the map; `accuracy = 2` gives the three-point rule
#' and \eqn{2d}.
#'
#' For a matrix family each estimate is symmetrized as `(m + t(m)) / 2`, as
#' the higher orders are, so that an asymmetry of rounding size in
#' [param_value()] does not reach the derivative.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param accuracy The order of accuracy of the central rule, a positive even
#'   whole number.
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
numerical_d1 <- function(s, eta, accuracy = 4L) {
  param_tensor(s, eta, 1L, accuracy = accuracy)
}


#' Numerical Second Derivatives of a Parameter's Matrix
#'
#' @description
#' Estimates the \eqn{d(d+1)/2} distinct second derivatives
#' \eqn{\partial^2 V / \partial \eta_k \partial \eta_l}, taking exactly one
#' difference for each, through [param_tensor()]. Which quantity is differenced
#' depends on what the family already supplies: the analytic first derivative
#' where there is one, and [param_value()] itself where there is not. This is
#' the method [param_d2()] dispatches to when a family has not written its own.
#'
#' @details
#' # The two routes
#'
#' Where [param_d1()] is analytic, the component \eqn{(k, l)} is one central
#' difference of \eqn{\partial_k V} along \eqn{\eta_l}. Where it is not, a
#' diagonal component is the central second difference of [param_value()] and
#' an off-diagonal one the tensor product of two first differences, in
#' \eqn{\eta_k} and \eqn{\eta_l}. Differences are never composed in the same
#' variable, and two differences along different coordinates form one product
#' stencil, so every route is one layer.
#'
#' The step of a difference of order \eqn{r} is
#' \eqn{\varepsilon^{1/(r + a)}\max(1, |\eta_k|)}, with \eqn{a} the accuracy:
#' \eqn{r = 1} on the analytic route and \eqn{r = 2} on the value. The
#' truncation error is of order \eqn{h^a}.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param accuracy The order of accuracy of the central rules, a positive even
#'   whole number.
#'
#' @return A list of `choose(s@n_free + 1, 2)` estimates keyed as
#'   `param_tuple_names(s)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d2()], the generic that this serves, [numerical_d1()] for the
#'   order below, and [param_tensor()], which builds the stencils.
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
numerical_d2 <- function(s, eta, accuracy = 4L) {
  param_tensor(s, eta, 2L, accuracy = accuracy)
}


#' @title Default First Derivatives
#' @name param_d1.parameter
#' @description
#' The method that every [parameter()] inherits when it registers no [param_d1()] of
#' its own. It estimates \eqn{\partial V/\partial \eta_k} by one central
#' difference of [param_value()] in each component, through [numerical_d1()]
#' at its default accuracy 4: the five-point rule at the step
#' \eqn{\varepsilon^{1/5}\max(1, |\eta_k|)}, at a cost of \eqn{4d}
#' evaluations of the map. The families in this package do not reach it.
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
#' its own. It takes exactly one difference per component, through
#' [numerical_d2()] at its default accuracy 4, of whichever quantity the family
#' already supplies: the analytic [param_d1()] where there is one, and
#' [param_value()] itself where there is not. The truncation error is of order
#' \eqn{h^4} on every route. The families in this package do not reach it.
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` estimates keyed as
#'   `param_tuple_names(s)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#' @seealso [numerical_d2()], which does the work and describes the two
#'   routes, and [param_d1.parameter()] for the order below.
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
#' Estimates the distinct third derivatives, one tensor-product stencil per
#' index tuple, through [param_tensor()]. The stencil is applied to the
#' highest analytic derivative below the third that the family supplies, as
#' [param_analytic()] reads it, and to [param_value()] where there is none.
#' This is the method [param_d3()] dispatches to when a family has not written
#' its own.
#'
#' @details
#' # The stencil
#'
#' Let the tuple be \eqn{\alpha} and the analytic quantity differenced be the
#' derivative \eqn{\beta \le \alpha}, taken from the components that the tuple
#' repeats most. The residual \eqn{\alpha - \beta} has multiplicities
#' \eqn{m_1, \dots, m_r} in the components \eqn{c_1, \dots, c_r}, and with
#' \eqn{(o^{(m)}_i, w^{(m)}_i)} the offsets and weights of the central stencil
#' for an \eqn{m}-th derivative in one variable the estimate is their tensor
#' product,
#'
#' \deqn{\partial^\alpha V \approx
#'   \frac{1}{\prod_{j} h_j^{m_j}}
#'   \sum_{i_1, \dots, i_r} \left(\prod_{j=1}^{r} w^{(m_j)}_{i_j}\right)
#'   \partial^\beta V\Bigl(\eta + \textstyle\sum_{j=1}^{r} o^{(m_j)}_{i_j} h_j
#'   e_{c_j}\Bigr),}
#'
#' with the zero-weight nodes skipped and every point shared by several tuples
#' evaluated once. With an analytic [param_d2()] each tuple is one first
#' difference of a second derivative.
#'
#' # The stencil is never applied to a difference
#'
#' A lower-order numerical derivative is never differenced. The rounding of a
#' difference is amplified by \eqn{h^{-k}}, so two stages multiply their
#' errors. The quantity differenced is either the map or an analytic
#' derivative, so the rounding stays at a single layer whatever the family
#' implements.
#'
#' # Accuracy
#'
#' Truncation is of order \eqn{h^a} at a step of
#' \eqn{\varepsilon^{1/(r + a)}\max(1, |\eta_k|)}, with \eqn{r} the order of
#' the difference and \eqn{a} the accuracy, and rounding is of order
#' \eqn{\varepsilon / h^r}.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param accuracy The order of accuracy of the central rules, a positive even
#'   whole number.
#'
#' @return A list of `choose(s@n_free + 2, 3)` estimates keyed as
#'   `param_tuple_names(s, 3)` and in that order, each shaped like
#'   [param_value()]'s result and symmetrized for a matrix family.
#'
#' @seealso [param_d3()], the generic that this serves, [param_tensor()], which
#'   builds the stencils, and [numerical_d4()] for the order above.
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
numerical_d3 <- function(s, eta, accuracy = 4L) {
  param_tensor(s, eta, 3L, accuracy = accuracy)
}


#' Numerical Fourth Derivatives of a Parameter's Value
#'
#' @description
#' Estimates the distinct fourth derivatives, one tensor-product stencil per
#' index tuple, through [param_tensor()], applied to the highest analytic
#' derivative below the fourth that the family supplies and to [param_value()]
#' where there is none. It is the order-four analogue of [numerical_d3()] and
#' the method [param_d4()] dispatches to when a family has not written its own.
#'
#' @details
#' The tensor product is the one written out under [numerical_d3()]. With an
#' analytic [param_d3()] each tuple is one first difference of a third
#' derivative, and the result is accurate to about \eqn{10^{-12}}. From
#' [param_value()] alone the difference is of order four: rounding is then of
#' order \eqn{\varepsilon / h^4}, and the attainable accuracy at the default
#' accuracy 4 is of order \eqn{\varepsilon^{1/2}}, enough to catch a
#' transcription error in a closed form.
#'
#' @param s A [parameter()] object, of any branch.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param accuracy The order of accuracy of the central rules, a positive even
#'   whole number.
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
#' @export
numerical_d4 <- function(s, eta, accuracy = 4L) {
  param_tensor(s, eta, 4L, accuracy = accuracy)
}


#' @title Default Third Derivatives
#' @name param_d3.parameter
#' @description
#' The method that every [parameter()] inherits when it registers no [param_d3()] of
#' its own. It applies one tensor-product stencil per index tuple, through
#' [numerical_d3()] at its default accuracy 4, to the highest analytic
#' derivative below the third and to [param_value()] where there is none,
#' never to a lower-order numerical derivative. The families in this package
#' do not reach it.
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
#' its own. It applies one tensor-product stencil per index tuple, through
#' [numerical_d4()] at its default accuracy 4, as [param_d3.parameter()] does,
#' with the multiplicities summing to four. From [param_value()] alone the
#' result is enough to catch a wrong closed form and not enough for use in a
#' fit. The families in this package do not reach it.
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
#' [param_d3logdet()] of its own. Each component is one tensor-product central
#' difference, through [param_tensor()] at accuracy 4, of the highest analytic
#' log-determinant derivative below the third, as [param_analytic()] reads
#' it: for the tuple \eqn{(k, l, m)}, one first difference of the
#' \eqn{(k, l)} component of [param_d2logdet()] along \eqn{\eta_m}, at the
#' step \eqn{\varepsilon^{1/5}\max(1, |\eta_m|)}.
#' @details
#' # Accuracy
#'
#' [param_d2logdet()] is an exact identity **given** the matrix derivative
#' arrays, so it is differenced where [param_d1()] and [param_d2()] are
#' analytic or where the family writes it out, and the error is of order
#' \eqn{10^{-11}}. Where the base-class [param_d2logdet()] would read
#' numerical arrays it is not used: the stencil is applied to
#' [param_dlogdet()] where [param_d1()] is analytic, and otherwise to
#' [param_logdet()], which is exact, as a third difference. A numerical
#' quantity is therefore never differenced again, and the third difference of
#' the log-determinant is accurate to about \eqn{10^{-7}}.
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
  param_tensor(s, eta, 3L, what = "logdet")
}

#' @title Default Fourth Log-Determinant Derivatives
#' @name param_d4logdet.matrix_parameter
#' @description
#' The method that every [matrix_parameter()] inherits when it registers no
#' [param_d4logdet()] of its own. Each component is one tensor-product central
#' difference, through [param_tensor()] at accuracy 4, of the highest analytic
#' log-determinant derivative below the fourth: one first difference of
#' [param_d3logdet()] where the family registers it, and otherwise, for the
#' tuple \eqn{(k, l, m, n)}, one second-order stencil of the \eqn{(k, l)}
#' component of [param_d2logdet()] in the remaining two components, the
#' second difference where \eqn{m = n} and the product of two first
#' differences where they differ.
#' @details
#' # Accuracy
#'
#' With analytic [param_d1()] and [param_d2()], or with a family's own
#' [param_d2logdet()], the differencing is a single layer on an exact
#' quantity, accurate to about \eqn{10^{-9}}.
#'
#' With only [param_value()] supplied, the base-class [param_d2logdet()] would
#' read numerical arrays and is not used: the stencil is a fourth difference
#' of [param_logdet()], which is exact, accurate to about \eqn{10^{-6}}. A
#' family that needs this order to more digits supplies [param_d1()] and
#' [param_d2()] in closed form, or its own [param_d2logdet()].
#' [param_is_numerical()] reports which components are numerical.
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
  param_tensor(s, eta, 4L, what = "logdet")
}
