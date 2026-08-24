#' @include parameter_class.R
NULL


#' The Value a Parameter Produces
#'
#' @description
#' Evaluates the map a parametrization stands for: given a free vector
#' \eqn{\eta} on the unconstrained scale, returns the constrained value it
#' produces. The shape depends on the family. A [matrix_parameter()] returns a
#' symmetric \eqn{p \times p} matrix, [simplex()] a probability vector, and
#' [transition_matrix()] a row-stochastic matrix. Whatever \eqn{\eta} is
#' handed in, the value satisfies the family's constraint, so a caller never
#' has to test the result.
#'
#' This is the only method a new family must write. Every derivative order is
#' then available numerically, and for a matrix family so are the
#' log-determinant, the solve and the factor.
#'
#' @details
#' # Validation happens before dispatch
#'
#' The generic checks `eta` in its own body and passes the checked vector on,
#' so every method inherits the check, including one written outside the
#' package. Three things are rejected: a non-numeric `eta`, a length other than
#' `s@n_free` (the message names both counts and lists the family's
#' `free_names`), and any `NA`, `NaN` or infinite entry. The last is a caller
#' error rather than a boundary of the domain, the unconstrained scale having no
#' edge to reach.
#'
#' # What is exact and what is not
#'
#' Every family in this package writes its own [param_value()] out in closed
#' form. None of them falls back to anything, and `param_is_numerical()`
#' returns `FALSE` for all nine components of all fifteen. The numerical
#' methods on [parameter()] exist for a family written elsewhere.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, the point on the unconstrained scale, of
#' length \eqn{d = } `s@n_free`. \eqn{p} is the side of the matrix a matrix
#' family returns.
#'
#' @param s An object inheriting from class [parameter()], from any of the
#'   constructors: [log_cholesky()], [ar1()], [simplex()] and the rest.
#' @param eta A numeric vector of length `s@n_free`, finite in every entry. The
#'   order of its entries is the order of `s@free_names`. Names are ignored and
#'   stripped, so a value that arrives labeled by a link does not leak that
#'   label into the result.
#' @param ... Passed to the method. No method in this package reads it; it is
#'   part of the signature so that a family written elsewhere can take further
#'   arguments of its own.
#'
#' @return The constrained value, shaped as the family declares.
#'   \describe{
#'     \item{a matrix family}{a symmetric `s@dimension` by `s@dimension`
#'       numeric matrix, with `dimnames` `v1`, `v2`, ..., positive definite
#'       unless the family declares a deficient rank.}
#'     \item{[simplex()]}{a numeric vector of length `s@n_free + 1`, named `p1`,
#'       `p2`, ..., with positive entries summing to 1.}
#'     \item{[transition_matrix()]}{a square numeric matrix with positive
#'       entries and rows summing to 1, with `dimnames` `s1`, `s2`, ...}
#'   }
#'
#' @seealso [param_free()] for the inverse map, [param_d1()] through
#'   [param_d4()] for the derivatives of this one, [param_logdet()],
#'   [param_solve()] and [param_factor()] for what a matrix family adds, and
#'   [check_parameter()] to verify a family end to end.
#'
#' @examples
#' # An unstructured covariance. Any six numbers give a positive definite
#' # matrix, which is the property the parametrization exists for.
#' s <- log_cholesky(3)
#' M <- param_value(s, c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2))
#' round(M, 4)
#' eigen(M, only.values = TRUE)$values > 0
#'
#' # Even absurd free values stay inside the cone.
#' eigen(param_value(s, c(-8, 9, -7, 100, -100, 50)),
#'       only.values = TRUE)$values > 0
#'
#' # A probability vector, and a row-stochastic matrix.
#' sum(param_value(simplex(4), c(1.2, -0.7, 0.3)))
#' rowSums(param_value(transition_matrix(3), rep(0.4, 6)))
#'
#' # The three ways a free vector is rejected.
#' try(param_value(s, c(1, 2)))
#' try(param_value(s, c(0, 0, 0, 0, 0, Inf)))
#' try(param_value(s, letters[1:6]))
#'
#' @export
param_value <- S7::new_generic("param_value", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' The Free Vector Behind a Value
#'
#' @description
#' Inverts the map: given a value in the family's set, returns the free vector
#' \eqn{\eta} that [param_value()] would send there. The two are a bijection
#' onto the set, so the round trip closes to machine precision, and every one of
#' the fifteen families in this package inverts exactly, with a worst measured
#' error of \eqn{5 \times 10^{-16}}.
#'
#' Use it to start an optimizer from a matrix rather than from a free vector:
#' fit an unstructured covariance by moments, invert it, and hand the result in
#' as a starting point.
#'
#' @details
#' # Exact or refused, never approximated
#'
#' A family either writes its inverse out or signals an error. The base method
#' on [parameter()] does the second on behalf of a family that has not written
#' the first, naming the family in the message. An inverse found by minimizing
#' \eqn{\lVert V(\eta) - m \rVert} would return a plausible \eqn{\eta} for a
#' matrix that is not in the set at all, and the caller could not tell that
#' answer from a correct one.
#'
#' # A value outside the set is rejected
#'
#' Every method checks `m` against the family's own constraint before
#' inverting, and the message says what failed. `log_cholesky()` rejects a
#' matrix that is not positive definite, with the verdict taken from the
#' eigenvalues; `simplex()` rejects a vector that does not sum to one, and does
#' not renormalize it, because a silent repair would hide the caller's mistake.
#' The shared checks are the shape and the symmetry: a matrix of the wrong side
#' is rejected with a message naming the side required, and an asymmetry above
#' \eqn{10^{-8}} relative is rejected while one below it is averaged away.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = } `s@n_free`, and \eqn{m}
#' the value on the constrained scale.
#'
#' @param s An object inheriting from class [parameter()].
#' @param m A value of the family's shape: a symmetric `s@dimension` by
#'   `s@dimension` numeric matrix on the matrix branch, a probability vector for
#'   [simplex()], a row-stochastic matrix for [transition_matrix()]. It must
#'   satisfy the family's constraint; see **Details** for what each family
#'   rejects.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`, such
#'   that `param_value(s, param_free(s, m))` recovers `m`.
#'
#' @seealso [param_value()], the map this inverts, and [check_parameter()],
#'   which closes the round trip as one of its checks.
#'
#' @examples
#' # The round trip closes exactly, in both directions.
#' s <- log_cholesky(3)
#' eta <- c(0.3, -0.2, 0.5, 0.1, -0.4, 0.2)
#' param_free(s, param_value(s, eta))
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' # Starting an optimizer from a matrix: invert the sample covariance.
#' set.seed(1)
#' y <- matrix(rnorm(300), 100, 3) %*% chol(matrix(c(1, .5, .2, .5, 1, .3,
#'                                                   .2, .3, 1), 3, 3))
#' start <- param_free(s, cov(y))
#' round(start, 3)
#' max(abs(param_value(s, start) - cov(y)))
#'
#' # A matrix outside the set is rejected, and the message says why.
#' try(param_free(s, diag(c(1, 1, -1))))
#'
#' # So is a vector off the simplex; it is not renormalized.
#' try(param_free(simplex(3), c(0.5, 0.6, 0.2)))
#'
#' @export
param_free <- S7::new_generic("param_free", "s", function(s, m, ...) {
  if (S7::S7_inherits(s, matrix_parameter)) m <- check_matrix(s, m)
  S7::S7_dispatch()
})


#' First Derivatives of a Parameter's Value
#'
#' @description
#' Differentiates the map once in each free value, returning
#' \eqn{\partial V / \partial \eta_k} for \eqn{k = 1, \dots, d} as a list of
#' \eqn{d} objects each shaped like the value itself. A likelihood needs these
#' to turn a score in the matrix into a score in the coordinates an optimizer
#' moves, by the chain rule
#' \eqn{\partial \ell / \partial \eta_k = \sum_{ij} (\partial \ell / \partial
#' M_{ij}) (\partial M_{ij} / \partial \eta_k)}.
#'
#' @details
#' # The list is keyed by the free names
#'
#' Entry `k` is the derivative in the free value `s@free_names[k]`, and the list
#' carries those names, so a consumer assembling a gradient can match a
#' coordinate to its label without tracking positions.
#'
#' # Exact, or one stencil
#'
#' Every family here writes this out in closed form. A family that does not gets
#' the method registered on [parameter()], which applies **one** three-point
#' central difference to [param_value()] in each component, at the step
#' [numericals7::fd_step()] gives for a first derivative. It is never a
#' difference of a lower-order difference, so the errors of two stages cannot
#' compound. `param_is_numerical()` reports which route a given object takes, and
#' returns `FALSE` throughout for every family this package ships.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = } `s@n_free`, and \eqn{V}
#' the value [param_value()] returns.
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#'   Checked before dispatch; see [param_value()] for the three rejections.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A list of `s@n_free` derivatives, named by `s@free_names`. Each entry
#'   has the shape of the value: a symmetric `s@dimension` by `s@dimension`
#'   matrix for a matrix family, a numeric vector of length `s@n_free + 1` for
#'   [simplex()]. The derivative of a symmetric matrix is symmetric, so each
#'   entry is too.
#'
#' @seealso [param_d2()], [param_d3()] and [param_d4()] for the higher orders,
#'   [param_dlogdet()] for the derivative of the log-determinant, and
#'   [param_is_numerical()] to ask which route an object takes.
#'
#' @examples
#' s <- log_cholesky(2)
#' eta <- c(0.2, -0.1, 0.4)
#' d1 <- param_d1(s, eta)
#' names(d1)
#' d1[["L2.1"]]
#'
#' # Each entry is a derivative, and a central difference of the map agrees.
#' h <- 1e-5
#' fd <- (param_value(s, eta + c(0, 0, h)) - param_value(s, eta - c(0, 0, h))) /
#'   (2 * h)
#' max(abs(d1[["L2.1"]] - fd))
#'
#' # For a scaled precision the derivative is the matrix itself, the free
#' # value being the log of the scale.
#' P <- crossprod(diff(diag(4)))
#' r <- scaled_matrix(P)
#' all.equal(param_d1(r, 1.3)[[1]], param_value(r, 1.3))
#'
#' @export
param_d1 <- S7::new_generic("param_d1", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Second Derivatives of a Parameter's Value
#'
#' @description
#' Differentiates the map twice, returning the \eqn{d(d+1)/2} distinct
#' components \eqn{\partial^2 V / \partial \eta_k \partial \eta_l}, each shaped
#' like the value. A mixed partial does not depend on the order of
#' differentiation, so only the unordered pairs are returned; the caller reads
#' the \eqn{(l, k)} entry off the \eqn{(k, l)} one. A Newton step in the free
#' vector needs these, and so does the observed information of a model whose
#' covariance is parametrized this way.
#'
#' @details
#' # Keys and their order
#'
#' The list is keyed by [param_tuple_names()], which puts the \eqn{d} diagonal
#' pairs first and then the \eqn{d(d-1)/2} off-diagonal ones in lexicographic
#' order. Diagonal first is what a consumer filling a Hessian wants, and the
#' ordering is part of the interface: [param_tuple_indices()] returns the index
#' pairs in exactly the same order, so a caller can walk the two together.
#'
#' Both come from one enumeration, and neither is produced by taking a key
#' apart. Splitting `"log_L1:log_L2"` on `":"` looks equivalent and is not: a
#' free value whose own label contains the separator splits into the wrong
#' number of pieces, and the failure is silent.
#'
#' # Exact, or one stencil
#'
#' Every family here is closed form. The fallback on [parameter()] takes one
#' difference and no more, and which quantity it differences depends on what the
#' family already supplies. Where [param_d1()] is analytic it differences that
#' once, in the other index. Where [param_d1()] is itself numerical it works on
#' [param_value()] directly: a three-point second difference where the two
#' indices coincide, and one difference in each of the two components where they
#' differ. A mixed derivative in two different variables is one stencil however
#' it is written; the nesting the toolkit forbids is two differences in the same
#' variable.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = } `s@n_free`, and \eqn{V}
#' the value [param_value()] returns.
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A named list of `choose(s@n_free + 1, 2)` entries, keyed as
#'   `param_tuple_names(s)` and in that order. Each entry has the shape of the
#'   value, symmetric for a matrix family.
#'
#' @seealso [param_d1()] for the first order, [param_d3()] and [param_d4()] for
#'   the higher ones, [param_tuple_names()] and [param_tuple_indices()] for the
#'   keys, and [param_d2logdet()] for the second derivative of the
#'   log-determinant, which is not the trace of this.
#'
#' @examples
#' s <- log_cholesky(2)
#' eta <- c(0.2, -0.1, 0.4)
#' d2 <- param_d2(s, eta)
#' names(d2)
#'
#' # The diagonal pairs come first, and the keys match the index tuples.
#' param_tuple_indices(s, 2)[1:3]
#'
#' # A second difference of the map agrees with the L2.1 diagonal component.
#' h <- 1e-4
#' e <- c(0, 0, h)
#' fd <- (param_value(s, eta + e) - 2 * param_value(s, eta) +
#'          param_value(s, eta - e)) / h^2
#' max(abs(d2[["L2.1:L2.1"]] - fd))
#'
#' # log_cholesky is quadratic in a below-diagonal free value, so the second
#' # derivative there is constant and the third vanishes.
#' d2[["L2.1:L2.1"]]
#' max(abs(param_d3(s, eta)[["L2.1:L2.1:L2.1"]]))
#'
#' @export
param_d2 <- S7::new_generic("param_d2", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Log-Determinant of a Parameter's Matrix
#'
#' @description
#' Returns \eqn{\log|M|} when the family is of full rank, and the log
#' pseudo-determinant -- the sum of the logs of the non-zero eigenvalues --
#' when it is not.
#'
#' @details
#' One generic covers both because a consumer asks the same question of either:
#' what normalizing constant does this matrix contribute. The object records
#' which answer is the right one, through its declared rank.
#'
#' The sign in front of the result is the consumer's arithmetic. A gaussian
#' log-density written in the covariance carries \eqn{-\frac{1}{2}\log|\Sigma|}
#' and one written in the precision carries \eqn{+\frac{1}{2}\log|\Omega|};
#' the parameter answers what the log-determinant is, and the likelihood
#' decides where it goes.
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param ... Passed to methods.
#'
#' @return A single number.
#'
#' @seealso [param_dlogdet()], [param_solve()]
#'
#' @examples
#' param_logdet(log_cholesky(2), c(0, 0, 0))
#'
#' # a rank-deficient precision: the pseudo-determinant, over the 4 non-zero
#' # eigenvalues of a second-difference penalty on 6 coefficients
#' p <- crossprod(diff(diag(6), differences = 2))
#' s <- scaled_matrix(p)
#' c(rank = s@rank, logdet = param_logdet(s, 0))
#'
#' @export
param_logdet <- S7::new_generic("param_logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Gradient of the Log-Determinant
#'
#' @description
#' Returns \eqn{\partial \log|M| / \partial \eta_k} for every free value, or
#' the same derivative of the log pseudo-determinant when the family is rank
#' deficient.
#'
#' @details
#' The identity behind it is
#' \eqn{\partial_k \log|M| = \mathrm{tr}(M^{-1} \partial_k M)}, with the
#' Moore-Penrose inverse in place of \eqn{M^{-1}} in the rank-deficient case.
#' A closed form is therefore never an independent claim: it must agree with
#' [param_d1()] through that identity, and
#' [check_parameter()] compares the two routes.
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of length `s@n_free`, named by
#'   `s@free_names`.
#'
#' @seealso [param_logdet()], [param_d2logdet()]
#'
#' @examples
#' # for a scaled precision the derivative is the rank, whatever the scale
#' s <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
#' c(param_dlogdet(s, -3), param_dlogdet(s, 5), rank = s@rank)
#'
#' @export
param_dlogdet <- S7::new_generic("param_dlogdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Hessian of the Log-Determinant
#'
#' @description
#' Returns the distinct second derivatives of the log-determinant, or of the
#' log pseudo-determinant, keyed as [param_tuple_names()].
#'
#' @details
#' Differentiating \eqn{\partial_k \log|M| = \mathrm{tr}(M^{-1}\partial_k M)}
#' once more, and using \eqn{\partial_l M^{-1} = -M^{-1}(\partial_l M)M^{-1}},
#'
#' \deqn{\partial_{kl} \log|M| = \mathrm{tr}\!\left(M^{-1}\partial_{kl}M\right)
#'   - \mathrm{tr}\!\left(M^{-1}(\partial_k M) M^{-1} (\partial_l M)\right),}
#'
#' with the Moore-Penrose inverse in place of \eqn{M^{-1}} when the family is
#' rank deficient. The second term is what makes the Hessian of the
#' log-determinant differ from the trace of the second derivative of the
#' matrix, and it is where a closed form is usually got wrong;
#' [check_parameter()] compares this route against
#' [param_d2()].
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param ... Passed to methods.
#'
#' @return A named numeric vector, keyed as `param_tuple_names(s)`.
#'
#' @seealso [param_dlogdet()]
#'
#' @examples
#' param_d2logdet(log_cholesky(2), c(0.2, -0.1, 0.4))
#'
#' @export
param_d2logdet <- S7::new_generic("param_d2logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Solve Through a Parameter's Matrix
#'
#' @description
#' Returns \eqn{M^{-1} B}, computed through a factor rather than through an
#' explicit inverse.
#'
#' @details
#' A rank-deficient parameter rejects rather than returning a pseudo-inverse.
#' What a consumer of an improper prior needs is the quadratic form and the
#' log pseudo-determinant -- the penalized normal equations invert
#' \eqn{X^\top X + \lambda P}, which is non-singular even when \eqn{P} is not,
#' and is assembled by the consumer -- so a pseudo-inverse would be a plausible
#' matrix answering a question nobody asked.
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param b A numeric matrix or vector with `s@dimension` rows. Defaults to the
#'   identity, which returns the inverse.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with `s@dimension` rows.
#'
#' @seealso [param_factor()], [param_logdet()]
#'
#' @examples
#' s <- log_cholesky(2)
#' round(param_solve(s, c(0, 0, 0.5)), 4)
#'
#' @export
param_solve <- S7::new_generic("param_solve", "s", function(s, eta, b = NULL, ...) {
  eta <- check_eta(s, eta)
  if (!S7::S7_inherits(s, matrix_parameter)) {
    stop(sprintf(
      "param_solve() is not defined for '%s': the value is not a symmetric matrix.",
      s@param_name
    ), call. = FALSE)
  }
  if (s@rank < s@dimension) {
    stop(sprintf(paste0(
      "'%s' is rank deficient (%d of %d), so it has no inverse. A consumer of\n",
      "  an improper prior needs the quadratic form and the log\n",
      "  pseudo-determinant, not a pseudo-inverse: assemble the matrix the\n",
      "  model actually inverts and solve that."
    ), s@param_name, s@rank, s@dimension), call. = FALSE)
  }
  if (is.null(b)) b <- diag(s@dimension)
  if (!is.matrix(b)) b <- as.matrix(b)
  if (nrow(b) != s@dimension) {
    stop(sprintf("'b' must have %d rows.", s@dimension), call. = FALSE)
  }
  S7::S7_dispatch()
})


#' A Factor of a Parameter's Matrix
#'
#' @description
#' Returns a lower triangular \eqn{L} with \eqn{M = L L^\top}.
#'
#' @details
#' Rejected for a rank-deficient family, for the reason given in
#' [param_solve()].
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param ... Passed to methods.
#'
#' @return A lower triangular numeric matrix with `s@dimension` rows and
#'   columns.
#'
#' @seealso [param_solve()]
#'
#' @examples
#' round(param_factor(log_cholesky(2), c(0.1, 0.2, -0.3)), 4)
#'
#' @export
param_factor <- S7::new_generic("param_factor", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  if (!S7::S7_inherits(s, matrix_parameter)) {
    stop(sprintf(
      "param_factor() is not defined for '%s': the value is not a symmetric matrix.",
      s@param_name
    ), call. = FALSE)
  }
  if (s@rank < s@dimension) {
    stop(sprintf(
      "'%s' is rank deficient (%d of %d), so it has no Cholesky factor.",
      s@param_name, s@rank, s@dimension
    ), call. = FALSE)
  }
  S7::S7_dispatch()
})


#' Third Derivatives of a Parameter's Value
#'
#' @description
#' Differentiates the map three times, returning the distinct components
#' \eqn{\partial^3 V / \partial \eta_k \partial \eta_l \partial \eta_m}, each
#' shaped like the value. The derivative is symmetric in its indices, so the
#' \eqn{\binom{d+2}{3}} unordered triples are returned, not the \eqn{d^3}
#' ordered ones. A third-order chain rule through a link needs them, and so does
#' the exact gradient of a marginal criterion, which differentiates a penalized
#' mode with respect to a hyperparameter.
#'
#' @details
#' # Keys
#'
#' The list is keyed by `param_tuple_names(s, 3)`, the lexicographic
#' combinations with repetition, and `param_tuple_indices(s, 3)` gives the index
#' triples in the same order. At \eqn{d = 3} there are ten of them, at
#' \eqn{d = 6} fifty-six.
#'
#' # Exact, or one stencil
#'
#' Every family here is closed form. The fallback on [parameter()] applies one
#' product stencil per triple, of the width each repeated index calls for: an
#' index appearing three times takes the four-point third-difference factor, one
#' appearing twice the three-point second-difference factor, and a distinct
#' index a two-point first-difference factor. The product is evaluated in one
#' pass over [param_value()], so no order is reached by differencing the order
#' below.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = } `s@n_free`, and \eqn{V}
#' the value [param_value()] returns.
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A named list of `choose(s@n_free + 2, 3)` entries, keyed as
#'   `param_tuple_names(s, 3)` and in that order. Each entry has the shape of the
#'   value.
#'
#' @seealso [param_d2()] and [param_d4()] for the neighbouring orders,
#'   [param_tuple_names()] for the keys, and [param_d3logdet()] for the third
#'   derivative of the log-determinant.
#'
#' @examples
#' # A scalar matrix is exp(eta) times the identity, so every order in the one
#' # free value is the matrix again.
#' s <- scalar_matrix(2)
#' param_d3(s, 0.3)
#' all.equal(param_d3(s, 0.3)[[1]], param_value(s, 0.3))
#'
#' # log_cholesky is quadratic in each below-diagonal free value, so a triple
#' # that repeats one of those three times vanishes exactly.
#' q <- log_cholesky(2)
#' d3 <- param_d3(q, c(0.2, -0.1, 0.4))
#' length(d3)
#' max(abs(d3[["L2.1:L2.1:L2.1"]]))
#'
#' # It does not vanish in a diagonal free value, which enters through a log.
#' d3[["log_L1:log_L1:log_L1"]]
#'
#' @export
param_d3 <- S7::new_generic("param_d3", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Fourth Derivatives of a Parameter's Value
#'
#' @description
#' Differentiates the map four times, returning the distinct components
#'
#' \deqn{\frac{\partial^{4} V(\eta)}
#'   {\partial\eta_k\,\partial\eta_l\,\partial\eta_m\,\partial\eta_n},}
#'
#' each shaped like the value. Fourth order is where the contract stops. A
#' fourth-order chain rule through a link needs exactly this much, and nothing
#' in the toolkit asks for a fifth.
#'
#' @details
#' # Keys
#'
#' The derivative is symmetric in its indices, so the \eqn{\binom{d+3}{4}}
#' unordered quadruples are returned, keyed by `param_tuple_names(s, 4)` in
#' lexicographic order, with `param_tuple_indices(s, 4)` giving the index
#' quadruples in the same order. At \eqn{d = 3} there are fifteen, at
#' \eqn{d = 6} a hundred and twenty-six, so the list grows quickly in the size
#' of the matrix.
#'
#' # Exact, or one stencil
#'
#' Every family here is closed form. The fallback on [parameter()] applies one
#' product stencil per quadruple, with a factor per distinct index of the width
#' its multiplicity calls for, evaluated in one pass over [param_value()]. At
#' fourth order the rounding of a difference grows as \eqn{\varepsilon / h^4},
#' so a numerical answer here is the least accurate of the four orders; a family
#' fitted in earnest is better served by writing the closed form out.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = } `s@n_free`, and \eqn{V}
#' the value [param_value()] returns.
#'
#' @param s An object inheriting from class [parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A named list of `choose(s@n_free + 3, 4)` entries, keyed as
#'   `param_tuple_names(s, 4)` and in that order. Each entry has the shape of
#'   the value.
#'
#' @seealso [param_d3()] for the order below, [param_tuple_names()] for the
#'   keys, and [param_d4logdet()] for the fourth derivative of the
#'   log-determinant.
#'
#' @examples
#' # For a scalar matrix, exp(eta) times the identity, every order is the
#' # matrix again.
#' s <- scalar_matrix(2)
#' all.equal(param_d4(s, 0.3)[[1]], param_value(s, 0.3))
#'
#' # The list is over unordered quadruples, so it is shorter than d^4.
#' q <- log_cholesky(3)
#' c(returned = length(param_d4(q, rep(0.1, 6))),
#'   unordered = choose(6 + 3, 4), ordered = 6^4)
#'
#' # log_cholesky is quadratic in a below-diagonal free value, so any quadruple
#' # repeating one of those three or more times is exactly zero.
#' d4 <- param_d4(log_cholesky(2), c(0.2, -0.1, 0.4))
#' max(abs(d4[["L2.1:L2.1:L2.1:L2.1"]]))
#' max(abs(d4[["log_L1:L2.1:L2.1:L2.1"]]))
#'
#' @export
param_d4 <- S7::new_generic("param_d4", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Third and Fourth Derivatives of the Log-Determinant
#'
#' @description
#' The higher derivatives of the log-(pseudo-)determinant, keyed as
#' [param_tuple_names()] of the matching order.
#'
#' @details
#' They follow from the second derivative of [param_d2logdet()]
#' by the same two rules, applied again: the trace is linear, and
#'
#' \deqn{\partial_m M^{-1} = -M^{-1}(\partial_m M)M^{-1}.}
#'
#' Every term of the result is therefore a trace of an alternating product
#' \eqn{M^{-1}(\partial_{I_1}M)M^{-1}(\partial_{I_2}M)\cdots}, one factor
#' per block of a partition of the index set, with the sign and the
#' multiplicity the two rules produce. The expansion is not transcribed:
#' the package differentiates [param_d1()] through
#' [param_d4()] directly, and [check_parameter()] holds
#' the result against a numerical differentiation of the order below, which
#' shares none of its arithmetic.
#'
#' @param s An object inheriting from class [matrix_parameter()].
#' @param eta A numeric vector of length `s@n_free`.
#' @param ... Passed to methods.
#'
#' @return A named numeric vector.
#'
#' @seealso [param_d2logdet()]
#'
#' @examples
#' param_d3logdet(log_cholesky(2), c(0.1, -0.2, 0.4))
#'
#' @export
param_d3logdet <- S7::new_generic("param_d3logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})

#' @rdname param_d3logdet
#' @export
param_d4logdet <- S7::new_generic("param_d4logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Names of the Distinct Derivative Components
#'
#' @description
#' The keys of the derivative lists: one per unordered tuple of free values.
#' At order 2, the keys of [param_d2()] and
#' [param_d2logdet()], diagonal pairs first; at orders 3 and 4,
#' the keys of [param_d3()] and [param_d4()],
#' lexicographic combinations with repetition.
#'
#' @details
#' A mixed partial derivative does not depend on the order in which the
#' free values are differentiated, so the components at order \eqn{k} are
#' the multisets of size \eqn{k} drawn from the \eqn{n} free values, of
#' which there are
#'
#' \deqn{\binom{n + k - 1}{k},}
#'
#' rather than the \eqn{n^k} ordered tuples. That is the length of the
#' vector returned and of every derivative list of that order.
#'
#' This exists so that nothing has to recover a tuple by splitting a key
#' apart. Generating the keys and the index tuples from one enumeration
#' cannot be fooled by a free name that contains the separator, which
#' splitting can. Use [param_tuple_indices()] for the tuples
#' themselves.
#'
#' @param s An object inheriting from class [parameter()].
#' @param order The derivative order: 2 (default), 3 or 4.
#'
#' @return A character vector, one entry per distinct component.
#'
#' @seealso [param_tuple_indices()]
#'
#' @examples
#' param_tuple_names(log_cholesky(2))
#' param_tuple_names(scalar_matrix(2), 3)
#'
#' @export
param_tuple_names <- function(s, order = 2L) {
  idx <- param_tuple_indices(s, order)
  nm <- s@free_names
  vapply(idx, function(t) paste(nm[t], collapse = ":"), character(1))
}


#' Index Tuples Behind the Derivative Component Names
#'
#' @description
#' The unordered index tuples [param_tuple_names()] names, in
#' exactly the same order. At order 2 the diagonal pairs come first and then
#' the off-diagonal ones, because consumers index a Hessian that way; at
#' orders 3 and 4 the tuples are the lexicographic combinations with
#' repetition, matching the enumeration \pkg{distributions7} uses for its
#' higher derivatives.
#'
#' @param s An object inheriting from class [parameter()].
#' @param order The derivative order: 1 to 4.
#'
#' @return A list of integer vectors of length `order`.
#'
#' @seealso [param_tuple_names()]
#'
#' @examples
#' param_tuple_indices(log_cholesky(2), 2)
#' length(param_tuple_indices(log_cholesky(2), 3))
#'
#' @export
param_tuple_indices <- function(s, order = 2L) {
  tuple_indices(s@n_free, order)
}


#' The Index Tuples of a Given Width
#'
#' @description
#' The enumeration behind [param_tuple_indices()], taken over a
#' number of variables rather than over a parameter, so that anything holding
#' derivatives over \eqn{d} variables -- a jet, for instance -- can share it.
#'
#' @param d The number of variables.
#' @param order The derivative order, 1 to 4.
#'
#' @return A list of integer vectors of length `order`.
#'
#' @seealso [param_tuple_indices()]
#'
#' @keywords internal
tuple_indices <- function(d, order = 2L) {
  # One enumeration for the whole toolkit: delegating is what makes it
  # impossible for this copy to disagree with the one a jet is keyed by.
  numericals7::tuple_indices(d, order)
}
