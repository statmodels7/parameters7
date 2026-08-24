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
#' Returns \eqn{\log|M|} for a full-rank family, and the log
#' pseudo-determinant, the sum of the logs of the non-zero eigenvalues, for a
#' rank-deficient one. This is the quantity a Gaussian likelihood needs beside
#' the quadratic form, and every matrix family answers it, in closed form where
#' one exists.
#'
#' One generic covers both cases because a consumer asks the same question of
#' either: what normalizing constant does this matrix contribute. Which answer
#' is the right one follows from the object's declared `rank`, so the caller
#' does not have to branch.
#'
#' @details
#' # The sign belongs to the caller
#'
#' A Gaussian log-density written in the covariance carries
#' \eqn{-\tfrac{1}{2}\log|\Sigma|} and one written in the precision carries
#' \eqn{+\tfrac{1}{2}\log|\Omega|}. The parameter reports what the
#' log-determinant is; the likelihood decides where it goes and with what sign.
#' The object's `role` records which side it was built for, and no method reads
#' it.
#'
#' # Computed from the parametrization, not from the matrix
#'
#' A closed form is usually far cheaper than a decomposition, and often simply
#' linear. In the log-Cholesky parametrization
#' \eqn{\log|M| = 2\sum_{i=1}^{p} \eta_i}, twice the sum of the first \eqn{p}
#' free values, so no factorization happens at all. The [ar1()] family answers
#' \eqn{p\,\eta_1 + (p-1)\log(1 - \rho^2)}, again in constant work whatever
#' \eqn{p} is. The base method on [matrix_parameter()], which a family written
#' elsewhere inherits, takes an eigendecomposition and sums the logs of the
#' eigenvalues above a relative tolerance, at \eqn{O(p^3)}.
#'
#' # Rank deficiency
#'
#' A deficient family returns the log pseudo-determinant, and its `null_basis`
#' says which directions were left out. That is the quantity an improper prior
#' contributes to a marginal likelihood: the penalized normal equations invert
#' \eqn{X^\top X + \lambda P}, which is non-singular even when \eqn{P} is not,
#' so the deficiency never has to be inverted. [param_solve()] and
#' [param_factor()] reject it for that reason.
#'
#' @section Notation:
#' \eqn{M} is the matrix [param_value()] returns, \eqn{\Sigma} read as a
#' covariance and \eqn{\Omega} read as a precision. \eqn{\eta} is the free
#' vector and \eqn{p} the side of the matrix.
#'
#' @param s An object inheriting from class [matrix_parameter()]. A family whose
#'   value is not a symmetric matrix has no method, and the call fails at
#'   dispatch with `Can't find method for param_logdet(<SimplexParam>)`.
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A single number. `-Inf` is not returned: a full-rank family is
#'   positive definite at every finite `eta`, and a deficient one drops its
#'   null directions instead of taking `log(0)`.
#'
#' @seealso [param_dlogdet()], [param_d2logdet()], [param_d3logdet()] and
#'   [param_d4logdet()] for its derivatives, [param_solve()] for the quadratic
#'   form that goes with it, and [param_null_basis()] for the directions a
#'   deficient family omits.
#'
#' @examples
#' # It agrees with the eigenvalues, and is cheaper: no decomposition is taken.
#' s <- log_cholesky(3)
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#' M <- param_value(s, eta)
#' c(param_logdet(s, eta), determinant(M, logarithm = TRUE)$modulus)
#'
#' # In this parametrization it is linear: twice the sum of the first p entries.
#' all.equal(param_logdet(s, eta), 2 * sum(eta[1:3]))
#'
#' # An AR(1) covariance, where it is not linear.
#' a <- ar1(4)
#' a@free_names
#' eta_a <- c(0.3, 0.8)
#' rho <- tanh(eta_a[2])
#' all.equal(param_logdet(a, eta_a), 4 * eta_a[1] + 3 * log(1 - rho^2))
#'
#' # A rank-deficient precision: the pseudo-determinant over the four non-zero
#' # eigenvalues of a second-difference penalty on six coefficients.
#' P <- crossprod(diff(diag(6), differences = 2))
#' r <- scaled_matrix(P)
#' ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
#' c(rank = r@rank, logdet = param_logdet(r, 0), from_ev = sum(log(ev[1:4])))
#'
#' # The scale enters it rank times over, so a step of 1 in the free value
#' # moves the log pseudo-determinant by the rank.
#' param_logdet(r, 1) - param_logdet(r, 0)
#'
#' @export
param_logdet <- S7::new_generic("param_logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Gradient of the Log-Determinant
#'
#' @description
#' Returns \eqn{\partial \log|M| / \partial \eta_k} for every free value, or the
#' same derivative of the log pseudo-determinant for a rank-deficient family. A
#' Gaussian score in the coordinates an optimizer moves is the derivative of the
#' quadratic form plus this, so it is needed at every iteration of a fit.
#'
#' @details
#' # The identity it must satisfy
#'
#' \deqn{\partial_k \log|M| = \mathrm{tr}\!\left(M^{-1} \partial_k M\right),}
#'
#' with the Moore-Penrose inverse in place of \eqn{M^{-1}} when the family is
#' rank deficient. A closed form here is therefore never an independent claim:
#' it has to agree with [param_d1()] through that trace, and
#' [check_parameter()] runs both routes and compares them. The example below
#' does the same in three lines.
#'
#' # What the answer looks like
#'
#' It is often constant. In the log-Cholesky parametrization the
#' log-determinant is \eqn{2\sum_{i\le p}\eta_i}, so the gradient is 2 in the
#' \eqn{p} diagonal directions and 0 in the rest, at every \eqn{\eta}. For a
#' [scaled_matrix()], where the free value is the log of a multiplier on a fixed
#' matrix, the gradient is the rank, again at every \eqn{\eta}: the scale enters
#' the pseudo-determinant once per non-zero eigenvalue.
#'
#' @section Notation:
#' \eqn{M} is the matrix [param_value()] returns and \eqn{\eta} the free vector,
#' of length \eqn{d = } `s@n_free`. \eqn{p} is the side of the matrix.
#'
#' @param s An object inheriting from class [matrix_parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#'
#' @seealso [param_logdet()] for the quantity differentiated here,
#'   [param_d2logdet()] for the next order, and [param_d1()], which this must
#'   agree with through the trace identity.
#'
#' @examples
#' # The trace identity, checked against param_d1() directly.
#' s <- log_cholesky(3)
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#' Minv <- solve(param_value(s, eta))
#' tr <- vapply(param_d1(s, eta), function(A) sum(diag(Minv %*% A)), numeric(1))
#' max(abs(param_dlogdet(s, eta) - tr))
#'
#' # In this parametrization the answer is 2 on the diagonal directions and 0
#' # elsewhere, whatever eta is.
#' param_dlogdet(s, eta)
#'
#' # For a scaled precision it is the rank, whatever the scale.
#' r <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
#' c(at_minus_3 = param_dlogdet(r, -3), at_5 = param_dlogdet(r, 5),
#'   rank = r@rank)
#'
#' # An AR(1) covariance, where the correlation direction is not constant.
#' param_dlogdet(ar1(4), c(0.3, 0.8))
#'
#' @export
param_dlogdet <- S7::new_generic("param_dlogdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Hessian of the Log-Determinant
#'
#' @description
#' Returns the distinct second derivatives of the log-determinant, or of the log
#' pseudo-determinant, keyed as [param_tuple_names()]. A Newton step in the free
#' vector of a Gaussian model needs them, and so does the observed information
#' of a covariance parametrized this way.
#'
#' @details
#' # The identity it must satisfy
#'
#' Differentiating \eqn{\partial_k \log|M| = \mathrm{tr}(M^{-1}\partial_k M)}
#' once more, with \eqn{\partial_l M^{-1} = -M^{-1}(\partial_l M)M^{-1}},
#'
#' \deqn{\partial_{kl} \log|M| = \mathrm{tr}\!\left(M^{-1}\partial_{kl}M\right)
#'   - \mathrm{tr}\!\left(M^{-1}(\partial_k M) M^{-1} (\partial_l M)\right),}
#'
#' and with the Moore-Penrose inverse in place of \eqn{M^{-1}} when the family
#' is rank deficient.
#'
#' The second trace is the term a hand-written closed form usually drops. Without
#' it the answer would be the trace of the second derivative of the matrix, which
#' is a different quantity, and dropping it is invisible on any family whose
#' log-determinant happens to be linear. [check_parameter()] compares this route
#' against [param_d2()] on every family it is given, and the example below runs
#' the same comparison in five lines.
#'
#' @section Notation:
#' \eqn{M} is the matrix [param_value()] returns and \eqn{\eta} the free vector,
#' of length \eqn{d = } `s@n_free`.
#'
#' @param s An object inheriting from class [matrix_parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A named numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
#'   `param_tuple_names(s)` and in that order: the `s@n_free` diagonal pairs
#'   first, then the off-diagonal ones. The \eqn{(l, k)} entry equals the
#'   \eqn{(k, l)} one and is not repeated.
#'
#' @seealso [param_dlogdet()] for the order below, [param_d3logdet()] and
#'   [param_d4logdet()] for the orders above, [param_tuple_names()] for the
#'   keys, and [param_d2()], which this must agree with through the two traces.
#'
#' @examples
#' # The identity, checked against param_d1() and param_d2() directly.
#' s <- ar1(4)
#' eta <- c(0.3, 0.8)
#' Minv <- solve(param_value(s, eta))
#' d1 <- param_d1(s, eta)
#' d2 <- param_d2(s, eta)
#' idx <- param_tuple_indices(s, 2)
#' pred <- vapply(seq_along(idx), function(i) {
#'   k <- idx[[i]][1]; l <- idx[[i]][2]
#'   sum(diag(Minv %*% d2[[i]])) -
#'     sum(diag(Minv %*% d1[[k]] %*% Minv %*% d1[[l]]))
#' }, numeric(1))
#' max(abs(param_d2logdet(s, eta) - pred))
#'
#' # Dropping the second trace would give a different answer here, so the
#' # comparison above has something to catch.
#' vapply(seq_along(idx), function(i) sum(diag(Minv %*% d2[[i]])), numeric(1))
#' param_d2logdet(s, eta)
#'
#' # For log_cholesky the log-determinant is linear in eta, so every second
#' # derivative is exactly zero and the check above could not see a mistake.
#' max(abs(param_d2logdet(log_cholesky(3), rep(0.2, 6))))
#'
#' @export
param_d2logdet <- S7::new_generic("param_d2logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Solve Through a Parameter's Matrix
#'
#' @description
#' Returns \eqn{M^{-1} B} for a right-hand side \eqn{B}, computed through a
#' factorization instead of by forming the inverse. This is what a Gaussian
#' quadratic form needs: `crossprod(r, param_solve(s, eta, r))` is
#' \eqn{r^\top M^{-1} r} at the cost of one triangular solve, where inverting
#' and multiplying would cost more and be less accurate. Called with no `B` it
#' does return the inverse, which is convenient for reading a covariance off a
#' precision.
#'
#' @details
#' # Rank deficiency is rejected
#'
#' A deficient family signals an error and names its rank, instead of returning
#' a pseudo-inverse. What a consumer of an improper prior needs is the quadratic
#' form and the log pseudo-determinant: penalized normal equations invert
#' \eqn{X^\top X + \lambda P}, which is non-singular even where \eqn{P} is not,
#' and the consumer assembles that matrix itself. A pseudo-inverse returned here
#' would be a plausible matrix answering a question nobody asked, and the caller
#' would have no way to tell.
#'
#' # A non-matrix family has no method
#'
#' [simplex()] and [transition_matrix()] inherit [parameter()] alone, not
#' [matrix_parameter()], and the generic checks for that before dispatching, so
#' the message names the family and says its value is not a symmetric matrix.
#'
#' # Cost
#'
#' The base method on [matrix_parameter()], which most families take, works
#' through a Cholesky factor: \eqn{O(p^3)} once plus \eqn{O(p^2)} per column of
#' \eqn{B}. Seven families override it with a closed-form inverse. [ar1()]
#' writes its tridiagonal precision out entry by entry, [compound_symmetry()]
#' uses Sherman-Morrison, and [block_diag()] and [kron_identity()] solve
#' blockwise, so none of them decomposes anything of side \eqn{p}.
#'
#' @section Notation:
#' \eqn{M} is the matrix [param_value()] returns, \eqn{p} its side, and
#' \eqn{\eta} the free vector.
#'
#' @param s An object inheriting from class [matrix_parameter()], of full rank.
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param b A numeric matrix or vector with `s@dimension` rows; a vector is
#'   treated as a one-column matrix. Defaults to `NULL`, which stands for the
#'   identity and returns the inverse. A wrong number of rows throws a message
#'   naming the number required.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A numeric matrix with `s@dimension` rows and as many columns as `b`,
#'   so `s@dimension` by `s@dimension` when `b` is left out. A vector `b`
#'   returns a one-column matrix, not a vector.
#'
#' @seealso [param_factor()] for the factor itself, [param_logdet()] for the
#'   other half of a Gaussian log-density, and [param_null_basis()] for what a
#'   deficient family offers in place of an inverse.
#'
#' @examples
#' s <- log_cholesky(3)
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#'
#' # With no b it is the inverse, and agrees with solve().
#' max(abs(param_solve(s, eta) - solve(param_value(s, eta))))
#'
#' # The quadratic form of a Gaussian log-density, in one solve.
#' r <- c(0.4, -1.1, 0.7)
#' q <- drop(crossprod(r, param_solve(s, eta, r)))
#' c(q, drop(r %*% solve(param_value(s, eta)) %*% r))
#'
#' # A vector b comes back as a one-column matrix.
#' dim(param_solve(s, eta, r))
#'
#' # A rank-deficient family refuses, and says what to do instead.
#' r_def <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
#' try(param_solve(r_def, 0))
#'
#' # So does a family whose value is not a symmetric matrix.
#' try(param_solve(simplex(3), c(0, 0)))
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
#' Returns the lower triangular \eqn{L} with \eqn{M = L L^\top}, the Cholesky
#' factor of the matrix the parametrization produces. Simulation is the usual
#' reason to want it: \eqn{L z} with \eqn{z} standard normal has covariance
#' \eqn{M}, so one factor draws as many vectors as needed. For
#' [log_cholesky()] the factor is what the parametrization holds anyway, so it
#' is returned without any arithmetic.
#'
#' @details
#' # Rank deficiency is rejected
#'
#' A deficient matrix has no Cholesky factor: a triangular \eqn{L} with
#' \eqn{L L^\top = M} would need a zero on the diagonal, and the factor is then
#' not unique. The generic signals an error naming the family and its rank
#' before dispatching. To simulate from a deficient covariance, take an
#' eigendecomposition and use the eigenvectors of the non-zero eigenvalues;
#' [param_null_basis()] gives the directions that carry no variance.
#'
#' # A non-matrix family has no method
#'
#' [simplex()] and [transition_matrix()] are checked for and refused by name,
#' their value not being a symmetric matrix.
#'
#' @section Notation:
#' \eqn{M} is the matrix [param_value()] returns, \eqn{p} its side, and
#' \eqn{\eta} the free vector.
#'
#' @param s An object inheriting from class [matrix_parameter()], of full rank.
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A `s@dimension` by `s@dimension` lower triangular numeric matrix with
#'   a positive diagonal, satisfying `L %*% t(L) == param_value(s, eta)`.
#'
#' @seealso [param_solve()] for a solve through the same factor,
#'   [param_value()] for the matrix it factors, and [log_cholesky()], whose
#'   free values are the logarithms of this factor's diagonal.
#'
#' @examples
#' s <- log_cholesky(3)
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#' L <- param_factor(s, eta)
#' round(L, 4)
#'
#' # It is a factor, exactly.
#' max(abs(L %*% t(L) - param_value(s, eta)))
#'
#' # In this parametrization the diagonal is exp() of the first p free values.
#' all.equal(diag(L), exp(eta[1:3]))
#'
#' # What it is for: simulating with the right covariance.
#' set.seed(1)
#' y <- t(L %*% matrix(rnorm(3 * 20000), 3, 20000))
#' round(cov(y) - param_value(s, eta), 2)
#'
#' # A rank-deficient family has no factor and says so.
#' try(param_factor(scaled_matrix(crossprod(diff(diag(6), differences = 2))), 0))
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
#' @seealso [param_d2()] and [param_d4()] for the neighboring orders,
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


#' Third Derivatives of the Log-Determinant
#'
#' @description
#' Returns the distinct third derivatives of the log-determinant, or of the log
#' pseudo-determinant, keyed as `param_tuple_names(s, 3)`. The exact gradient of
#' a marginal criterion reaches this order: differentiating a Laplace
#' approximation with respect to a hyperparameter moves the penalized mode, and
#' the third derivative is what the movement contracts against.
#'
#' @details
#' # Where it comes from
#'
#' Two rules generate every order. The trace is linear, and
#'
#' \deqn{\partial_m M^{-1} = -M^{-1}(\partial_m M)M^{-1}.}
#'
#' Applying them to \eqn{\partial_{kl}\log|M|} gives a sum of traces of
#' alternating products
#' \eqn{M^{-1}(\partial_{I_1}M)M^{-1}(\partial_{I_2}M)\cdots}, one factor per
#' block of a partition of the index set, with the sign and the multiplicity the
#' two rules produce. The Moore-Penrose inverse replaces \eqn{M^{-1}} for a
#' rank-deficient family.
#'
#' The expansion is not transcribed anywhere. The methods differentiate the
#' derivative arrays [param_d1()] through [param_d4()] instead, and
#' [check_parameter()] holds the result against a numerical differentiation of
#' the order below, which shares none of its arithmetic. A twenty-term expansion
#' written out by hand is exactly the kind of thing that is wrong and looks
#' right.
#'
#' @section Notation:
#' \eqn{M} is the matrix [param_value()] returns and \eqn{\eta} the free vector,
#' of length \eqn{d = } `s@n_free`. \eqn{\partial_I M} is the derivative of
#' \eqn{M} in the free values named by the index set \eqn{I}.
#'
#' @param s An object inheriting from class [matrix_parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A named numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
#'   `param_tuple_names(s, 3)` and in that order.
#'
#' @seealso [param_d2logdet()] for the order below, [param_d4logdet()] for the
#'   order above, [param_tuple_names()] for the keys, and [param_d3()] for the
#'   third derivative of the matrix itself.
#'
#' @examples
#' # An AR(1) covariance. The log-determinant is linear in the log of the
#' # scale, so only the correlation direction survives to third order.
#' s <- ar1(4)
#' param_d3logdet(s, c(0.3, 0.8))
#'
#' # A central difference of the order below agrees, and shares no arithmetic
#' # with the route the method takes.
#' h <- 1e-4
#' e <- c(0, h)
#' fd <- (param_d2logdet(s, c(0.3, 0.8) + e) -
#'          param_d2logdet(s, c(0.3, 0.8) - e)) / (2 * h)
#' fd[["z_rho:z_rho"]]
#' param_d3logdet(s, c(0.3, 0.8))[["z_rho:z_rho:z_rho"]]
#'
#' # For log_cholesky the log-determinant is linear in eta, so every order
#' # above the first is exactly zero.
#' max(abs(param_d3logdet(log_cholesky(3), rep(0.2, 6))))
#'
#' @export
param_d3logdet <- S7::new_generic("param_d3logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})

#' Fourth Derivatives of the Log-Determinant
#'
#' @description
#' Returns the distinct fourth derivatives of the log-determinant, or of the log
#' pseudo-determinant, keyed as `param_tuple_names(s, 4)`. This is the top of
#' the contract: the exact Hessian of a marginal criterion reaches it, and
#' nothing in the toolkit asks for a fifth.
#'
#' @details
#' # Where it comes from
#'
#' The same two rules as at third order. The trace is linear, and
#'
#' \deqn{\partial_m M^{-1} = -M^{-1}(\partial_m M)M^{-1},}
#'
#' so every term is a trace of an alternating product
#' \eqn{M^{-1}(\partial_{I_1}M)M^{-1}(\partial_{I_2}M)\cdots}, one factor per
#' block of a partition of the four indices. The number of terms grows with the
#' number of partitions, which is why the expansion is differentiated rather
#' than written out: the methods differentiate [param_d1()] through
#' [param_d4()], and [check_parameter()] compares the result with a numerical
#' differentiation of [param_d3logdet()].
#'
#' # Accuracy of the check
#'
#' The reference is one stencil on the analytic third order, so the comparison
#' is limited by that stencil rather than by the closed form. A family whose
#' log-determinant is linear in \eqn{\eta}, which includes [log_cholesky()] and
#' [matrix_log()], returns exact zeros here, and a check against them cannot
#' catch a mistake; [ar1()] and [compound_symmetry()] are the families where
#' this order has something to get wrong.
#'
#' @section Notation:
#' \eqn{M} is the matrix [param_value()] returns and \eqn{\eta} the free vector,
#' of length \eqn{d = } `s@n_free`. \eqn{\partial_I M} is the derivative of
#' \eqn{M} in the free values named by the index set \eqn{I}.
#'
#' @param s An object inheriting from class [matrix_parameter()].
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return A named numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
#'   `param_tuple_names(s, 4)` and in that order.
#'
#' @seealso [param_d3logdet()] for the order below, [param_tuple_names()] for
#'   the keys, and [param_d4()] for the fourth derivative of the matrix itself.
#'
#' @examples
#' # An AR(1) covariance: five components at d = 2, of which only the one
#' # purely in the correlation is non-zero, the scale entering linearly.
#' s <- ar1(4)
#' eta <- c(0.3, 0.8)
#' param_d4logdet(s, eta)
#'
#' # A central difference of the third order agrees.
#' h <- 1e-4
#' e <- c(0, h)
#' fd <- (param_d3logdet(s, eta + e) - param_d3logdet(s, eta - e)) / (2 * h)
#' c(analytic = param_d4logdet(s, eta)[["z_rho:z_rho:z_rho:z_rho"]],
#'   difference = fd[["z_rho:z_rho:z_rho"]])
#'
#' # The count is over unordered quadruples.
#' c(returned = length(param_d4logdet(s, eta)), choose(2 + 3, 4))
#'
#' @export
param_d4logdet <- S7::new_generic("param_d4logdet", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  S7::S7_dispatch()
})


#' Names of the Distinct Derivative Components
#'
#' @description
#' Returns the keys of a derivative list: one per unordered tuple of free
#' values, built by pasting the free names together with `":"`. At order 2 these
#' are the keys of [param_d2()] and [param_d2logdet()], with the diagonal pairs
#' first; at orders 3 and 4 the keys of [param_d3()] and [param_d4()], in
#' lexicographic order. Use it to look a component up by name, or to walk a
#' derivative list beside the index tuples [param_tuple_indices()] returns.
#'
#' @details
#' # How many there are
#'
#' A mixed partial does not depend on the order of differentiation, so the
#' components at order \eqn{k} are the multisets of size \eqn{k} drawn from the
#' \eqn{d} free values, of which there are
#'
#' \deqn{\binom{d + k - 1}{k},}
#'
#' against \eqn{d^k} ordered tuples. That is the length of the vector returned
#' and of every derivative list of that order. For \eqn{d = 6}, which is a
#' \eqn{3 \times 3} unstructured covariance, order 4 has 126 components against
#' 1296 ordered ones.
#'
#' # Why the keys and the tuples come from one enumeration
#'
#' Nothing in the toolkit recovers an index by splitting a key apart. Taking
#' `"log_L1:log_L2"` and splitting on `":"` works until a free name contains the
#' separator itself, and then it yields the wrong number of pieces and the
#' failure is silent. Generating the names and the indices from one enumeration
#' cannot be fooled that way. `param_tuple_names()` calls
#' [param_tuple_indices()] and labels what it gets, so the two agree by
#' construction.
#'
#' @section Notation:
#' \eqn{d} is the length of the free vector, `s@n_free`, and \eqn{k} the
#' derivative order.
#'
#' @param s An object inheriting from class [parameter()], whose `free_names`
#'   supply the labels and whose `n_free` supplies \eqn{d}.
#' @param order The derivative order. `2` by default, which is the order
#'   [param_d2()] and [param_d2logdet()] use; `3` and `4` name the components of
#'   [param_d3()] and [param_d4()]. `1` is accepted and returns the free names
#'   themselves. Anything else throws `'order' must be 1, 2, 3 or 4.`, the
#'   contract stopping at fourth order.
#'
#' @return A character vector of `choose(s@n_free + order - 1, order)` names, in
#'   the same order as the list they key. At order 2 the `s@n_free` diagonal
#'   pairs come first.
#'
#' @seealso [param_tuple_indices()] for the index tuples in the same order, and
#'   [param_d2()], [param_d3()], [param_d4()], [param_d2logdet()],
#'   [param_d3logdet()] and [param_d4logdet()] for the lists these key.
#'
#' @examples
#' s <- log_cholesky(2)
#' param_tuple_names(s)
#' param_tuple_names(s, 3)
#'
#' # They really are the keys of the derivative list, in order.
#' identical(param_tuple_names(s, 4), names(param_d4(s, c(0.2, -0.1, 0.4))))
#'
#' # The count is over unordered tuples.
#' q <- log_cholesky(3)
#' vapply(1:4, function(k) length(param_tuple_names(q, k)), integer(1))
#' choose(6 + 1:4 - 1, 1:4)
#'
#' # Order 2 puts the diagonal pairs first, which is how a consumer filling a
#' # Hessian wants them.
#' param_tuple_names(s, 2)
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
#' Returns the unordered index tuples that [param_tuple_names()] labels, in
#' exactly the same order, as integer positions into the free vector. A consumer
#' that has to know *which* free values a component differentiates in reads
#' these; a consumer that only has to display a label reads the names. At order
#' 2 the diagonal pairs come first and the off-diagonal ones after, which is how
#' a Hessian is filled; at orders 3 and 4 the tuples are the lexicographic
#' combinations with repetition.
#'
#' @details
#' The order at orders 3 and 4 matches the enumeration \pkg{distributions7} uses
#' for its own higher derivatives, so a consumer contracting a parameter's
#' derivative array against a distribution's can walk the two lists together
#' without reindexing. Both come from `numericals7::tuple_indices()`, which is
#' the single copy of the enumeration in the toolkit.
#'
#' @section Notation:
#' \eqn{d} is the length of the free vector, `s@n_free`, and \eqn{k} the
#' derivative order.
#'
#' @param s An object inheriting from class [parameter()], whose `n_free`
#'   supplies \eqn{d}.
#' @param order The derivative order: `1`, `2` (the default), `3` or `4`.
#'   Anything else throws `'order' must be 1, 2, 3 or 4.`
#'
#' @return A list of `choose(s@n_free + order - 1, order)` integer vectors, each
#'   of length `order`, holding positions in `1:s@n_free` in non-decreasing
#'   order within a tuple.
#'
#' @seealso [param_tuple_names()] for the same tuples as labels, and
#'   [param_d2()], [param_d3()] and [param_d4()] for the lists they index.
#'
#' @examples
#' s <- log_cholesky(2)
#'
#' # Order 2: the three diagonal pairs first, then the three off-diagonal ones.
#' param_tuple_indices(s, 2)
#'
#' # The tuples and the names are the same enumeration, so they line up.
#' idx <- param_tuple_indices(s, 3)
#' nm <- param_tuple_names(s, 3)
#' identical(nm, vapply(idx, function(t) paste(s@free_names[t], collapse = ":"),
#'                      character(1)))
#'
#' # What they are for: reading a component's indices to contract with.
#' d2 <- param_d2(s, c(0.2, -0.1, 0.4))
#' d1 <- param_d1(s, c(0.2, -0.1, 0.4))
#' i <- which(param_tuple_names(s, 2) == "log_L1:L2.1")
#' param_tuple_indices(s, 2)[[i]]
#' dim(d2[[i]])
#'
#' # Counts, over unordered tuples.
#' vapply(1:4, function(k) length(param_tuple_indices(s, k)), integer(1))
#'
#' @export
param_tuple_indices <- function(s, order = 2L) {
  tuple_indices(s@n_free, order)
}


#' The Index Tuples of a Given Width
#'
#' @description
#' The enumeration behind [param_tuple_indices()], taken over a plain count of
#' variables instead of over a parameter object, so that anything holding
#' derivatives over \eqn{d} variables can use it without constructing a
#' [parameter()]. A one-line forward to `numericals7::tuple_indices()`, which is
#' the toolkit's single copy of this enumeration. Delegating means this copy
#' cannot come to disagree with the one a consumer is keyed by.
#'
#' @param d The number of variables, a single non-negative integer.
#' @param order The derivative order: 1, 2, 3 or 4. Anything else throws
#'   `'order' must be 1, 2, 3 or 4.`
#'
#' @return A list of `choose(d + order - 1, order)` integer vectors, each of
#'   length `order`, holding positions in `1:d` in non-decreasing order within a
#'   tuple. At order 2 the diagonal pairs come first.
#'
#' @seealso [param_tuple_indices()], the parameter-facing wrapper, and
#'   [numericals7::tuple_indices()], the implementation.
#'
#' @keywords internal
tuple_indices <- function(d, order = 2L) {
  # One enumeration for the whole toolkit: delegating is what makes it
  # impossible for this copy to disagree with the one a jet is keyed by.
  numericals7::tuple_indices(d, order)
}
