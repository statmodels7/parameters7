#' @include parameter_class.R generics.R ar1.R compound_symmetry.R
#' @include autoregressive.R scaled_matrix.R simplex.R
NULL

# The free vector is what a fit estimates and not what a reader reads. Nobody
# reads the hyperbolic arctangent of a partial autocorrelation, and the
# quantity behind it is not recoverable from the printed summary of a
# multivariate fit either, whose standard deviations and correlations describe
# the matrix rather than the family that produced it. A family therefore
# declares the quantities it is about, together with the derivatives a
# consumer needs to carry a variance matrix onto them and the scale each
# interval should be built on.


#' Quantities a Family Is About
#'
#' @description
#' Asks a family which interpretable quantities it stands for, and returns them
#' with everything a consumer needs to report them: the values, the Jacobian of
#' the map from the free vector, the scale each interval should be built on, and
#' a label for the block. It exists because **the free vector is what a fit
#' estimates, never what a reader reads**: nobody reads the hyperbolic arc
#' tangent of a partial autocorrelation, and the quantity behind it cannot be
#' recovered from a printed covariance either.
#'
#' Use it to turn an estimate and its variance matrix into a table a reader can
#' use. The base class declares nothing, so a family whose matrix is all there is
#' to say returns `NULL` and a consumer falls back on reporting the matrix.
#'
#' @details
#' # What a consumer does with it
#'
#' Given the free vector's variance matrix \eqn{V}, the standard errors of the
#' declared quantities are the delta method, \eqn{\sqrt{\mathrm{diag}(J V
#' J^\top)}}. Each interval is then built on the scale `transform` names and
#' mapped back, so a variance stays positive and a correlation stays inside
#' \eqn{(-1, 1)}: on the raw scale an interval for a correlation routinely runs
#' past 1.
#'
#' # The Jacobians are closed form
#'
#' For a family whose coordinates are separate scalar links of separate free
#' values, \eqn{J} is the diagonal matrix of the inverse links' first
#' derivatives, which [readable_diagonal()] assembles. For [simplex()] it is the
#' softmax Jacobian. For [autoregressive()] the autoregressive coefficients come
#' out of the Levinson-Durbin recursion the family already propagates derivative
#' arrays through, so their derivatives in every free value are read off the
#' first-order block instead of being computed again.
#'
#' # Which families declare something
#'
#' [ar1()] and [compound_symmetry()] declare a variance and a correlation.
#' [scaled_matrix()] declares its multiplier, or `NULL` when it was built with
#' `link = NULL` and has none. [simplex()] declares the probability vector.
#' [autoregressive()] declares the variance, the partial autocorrelations that
#' parametrize it, **and** the autoregressive coefficients, which appear nowhere
#' in the covariance. Every other family, including [log_cholesky()],
#' [correlation_matrix()] and [transition_matrix()], takes the base method and
#' returns `NULL`.
#'
#' @param s A [parameter()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return `NULL` when the family declares nothing, otherwise a list with
#'   \describe{
#'     \item{`value`}{a named numeric vector of the quantities, on their own
#'       interpretable scale;}
#'     \item{`jacobian`}{a numeric matrix with one row per quantity and one
#'       column per free value, row names matching `value` and no column names;}
#'     \item{`transform`}{a character vector, one entry per quantity and named
#'       like `value`, naming the scale its interval is built on: one of
#'       `"identity"`, `"log"`, `"atanh"` or `"logit"`;}
#'     \item{`label`}{a single string naming the block, for a consumer laying out
#'       a printed summary. The family supplies it because the family is what
#'       holds the reading.}
#'   }
#'
#' @seealso [param_value()] for the matrix itself, [autoregressive()] for the
#'   family with the most to declare, and [readable_diagonal()] for the shared
#'   assembly.
#'
#' @examples
#' # An AR(1) is about a variance and a correlation, not about a log and an
#' # inverse hyperbolic tangent.
#' r <- param_readable(ar1(5), c(log(2), atanh(0.6)))
#' r$value
#' r$transform
#'
#' # What a consumer does with the Jacobian: the delta method. With a variance
#' # of 0.01 on the second free value, the correlation's standard error is
#' # |drho/dz| * 0.1.
#' V <- diag(c(0.04, 0.01))
#' sqrt(diag(r$jacobian %*% V %*% t(r$jacobian)))
#'
#' # An autoregression declares the coefficients as well as the partial
#' # autocorrelations, and at q = 2 they differ.
#' a <- param_readable(autoregressive(8, 2), c(0, 0.9, -0.4))
#' a$value
#'
#' # At q = 1 they coincide: the first partial autocorrelation IS the AR(1)
#' # coefficient.
#' b <- param_readable(autoregressive(5, 1), c(0, 0.9))
#' b$value[c("pacf1", "phi1")]
#'
#' # A family whose matrix is all there is to say declares nothing.
#' param_readable(log_cholesky(3), rep(0, 6))
#'
#' @export
param_readable <- S7::new_generic("param_readable", "s",
  function(s, eta, ...) S7::S7_dispatch()
)


#' No Declared Quantities
#'
#' @description
#' The method every [parameter()] inherits when it declares no interpretable
#' quantities of its own. It returns `NULL`, which is a positive statement and no
#' kind of gap: for [log_cholesky()], [matrix_log()], [correlation_matrix()] and the
#' four compositions the matrix itself is what a reader reads, and a consumer that
#' gets `NULL` reports the matrix as it already would.
#'
#' Returning `NULL` instead of throwing is what a consumer needs in order to loop
#' over every parameter in a model and ask each one, without knowing which kinds
#' declare something.
#'
#' @param s A [parameter()] object. Not read.
#' @param eta A numeric vector of free values. Not read.
#' @param ... Ignored.
#'
#' @return `NULL`.
#'
#' @seealso [param_readable()] for the contract and the list of families that do
#'   declare something.
#'
#' @name param_readable.parameter
#' @keywords internal
S7::method(param_readable, parameter) <- function(s, eta, ...) NULL


#' Quantities That Are Separate Links of Separate Free Values
#'
#' @description
#' Assembles a [param_readable()] declaration for a family whose quantities are
#' one scalar link each of one free value each. The Jacobian is then **diagonal**,
#' its \eqn{k}-th entry the inverse link's first derivative at \eqn{\eta_k}, and
#' there is nothing else to compute.
#'
#' [ar1()], [compound_symmetry()] and [scaled_matrix()] all declare through this,
#' which is why their three methods are one call each.
#'
#' @param links A list of \pkg{linkfunctions7} links, one per quantity, in the
#'   order of the free values they read. `links[[k]]` must be the link of
#'   `eta[k]`; the correspondence is positional and is not checked.
#' @param eta A numeric vector of free values, at least as long as `links`.
#' @param nm A character vector naming the quantities, the same length as
#'   `links`.
#' @param transform A character vector naming the scale each interval is built
#'   on, the same length as `nm`.
#' @param label A single string naming the block.
#'
#' @return A list with `value`, `jacobian`, `transform` and `label`, as
#'   [param_readable()] describes. The Jacobian is `length(nm)` by `length(eta)`
#'   with `nm` as its row names, and is diagonal.
#'
#' @seealso [param_readable()] for the contract, and
#'   [param_readable.AutoregressiveParam()] for the one family whose Jacobian is
#'   not diagonal.
#'
#' @keywords internal
readable_diagonal <- function(links, eta, nm, transform, label) {
  value <- vapply(seq_along(links), function(k) {
    linkfunctions7::linkinv(links[[k]], eta[k])
  }, numeric(1))
  jac <- matrix(0, length(nm), length(eta), dimnames = list(nm, NULL))
  for (k in seq_along(links)) {
    jac[k, k] <- linkfunctions7::dlinkinv(links[[k]], eta[k])
  }
  list(value = stats::setNames(value, nm), jacobian = jac,
       transform = stats::setNames(transform, nm), label = label)
}


#' The Scale and the Correlation of an AR(1)
#'
#' @description
#' Declares two quantities: the marginal variance and the correlation at lag one,
#' the two things an AR(1) covariance is about. Their intervals are built on the
#' log and the inverse hyperbolic tangent, so a variance stays positive and a
#' correlation stays inside \eqn{(-1, 1)}.
#'
#' The Jacobian is diagonal, through [readable_diagonal()]: each quantity is one
#' link of one free value.
#'
#' @param s An [Ar1Param()] object, whose two links are read.
#' @param eta A numeric vector of two free values.
#' @param ... Ignored.
#'
#' @return A list as described in [param_readable()].
#'
#' @name param_readable.Ar1Param
#' @keywords internal
S7::method(param_readable, Ar1Param) <- function(s, eta, ...) {
  readable_diagonal(
    list(s@param_params$link_scale, s@param_params$link_rho), eta,
    c("scale", "rho"), c("log", "atanh"), "Autoregressive structure"
  )
}


#' The Scale and the Common Correlation of a Compound Symmetry
#'
#' @description
#' Declares two quantities: the marginal variance and the correlation every pair
#' shares. The Jacobian is diagonal, through [readable_diagonal()], each quantity
#' being one link of one free value.
#'
#' The correlation's interval is built on the inverse hyperbolic tangent, which is
#' worth knowing here: the family's own link is bounded at \eqn{-1/(p-1)} and
#' never at \eqn{-1}, so an interval built this way can reach below the bound the
#' parametrization enforces. It respects \eqn{(-1, 1)}, which is the interval a
#' reader of a correlation expects.
#'
#' @param s A [CompoundSymmetryParam()] object, whose two links are read.
#' @param eta A numeric vector of two free values.
#' @param ... Ignored.
#'
#' @return A list as described in [param_readable()].
#'
#' @name param_readable.CompoundSymmetryParam
#' @keywords internal
S7::method(param_readable, CompoundSymmetryParam) <- function(s, eta, ...) {
  readable_diagonal(
    list(s@param_params$link_scale, s@param_params$link_rho), eta,
    c("scale", "rho"), c("log", "atanh"), "Compound symmetry structure"
  )
}


#' The Scale of a Fixed Matrix
#'
#' @description
#' Declares the multiplier \eqn{h(\eta)}, the one quantity a
#' [scaled_matrix()] is about, with its interval on the log scale. The fixed
#' matrix itself is the caller's own and needs no reporting.
#'
#' A parameter built with `link = NULL` has no free value and nothing to declare,
#' so this returns `NULL` there, which is the same answer the base class gives and
#' means the same thing.
#'
#' @param s A [ScaledMatrixParam()] object, whose `param_params$link` is read.
#' @param eta A numeric vector of at most one free value; `numeric(0)` for a
#'   fixed parameter.
#' @param ... Ignored.
#'
#' @return A list as described in [param_readable()], or `NULL`
#'   when the matrix is fixed and carries no free value.
#'
#' @name param_readable.ScaledMatrixParam
#' @keywords internal
S7::method(param_readable, ScaledMatrixParam) <- function(s, eta, ...) {
  if (is.null(s@param_params$link)) return(NULL)
  readable_diagonal(list(s@param_params$link), eta, "scale", "log",
                    "Scale of the fixed matrix")
}


#' The Probabilities Behind an Additive Log-Ratio
#'
#' @description
#' Declares the probability vector itself, `p1` ... `pK`, which a reader of a
#' mixture weight or a state distribution needs: the free values are log
#' ratios against the last category and are not interpretable on their own.
#'
#' The Jacobian is the softmax map's,
#' \eqn{\partial p_i/\partial\eta_j = p_i(\delta_{ij} - p_j)}, with the reference
#' category's row contributing \eqn{-p_K p_j}. It is \eqn{K} by \eqn{K-1}: one
#' more quantity than there are free values, since the reference probability is
#' determined by the others. Its **columns** sum to zero, \eqn{\sum_i p_i} being
#' the constant 1, measured at \eqn{7 \times 10^{-18}}.
#'
#' Every interval is built on the logit scale, so it stays inside \eqn{(0, 1)}.
#'
#' @param s A [SimplexParam()] object.
#' @param eta A numeric vector of free values, of length \eqn{K-1}.
#' @param ... Ignored.
#'
#' @return A list as described in [param_readable()].
#'
#' @name param_readable.SimplexParam
#' @keywords internal
S7::method(param_readable, SimplexParam) <- function(s, eta, ...) {
  p <- as.numeric(param_value(s, eta))
  k <- length(p)
  nm <- paste0("p", seq_len(k))
  jac <- matrix(0, k, k - 1L, dimnames = list(nm, NULL))
  for (i in seq_len(k)) {
    for (j in seq_len(k - 1L)) {
      jac[i, j] <- p[i] * ((i == j) - p[j])
    }
  }
  list(value = stats::setNames(p, nm), jacobian = jac,
       transform = stats::setNames(rep("logit", k), nm),
       label = "Probabilities")
}


#' The Scale, the Partial Autocorrelations and the Coefficients
#'
#' @description
#' Declares three groups: the marginal variance, the \eqn{q} partial
#' autocorrelations that parametrize the family, and the \eqn{q} autoregressive
#' coefficients they produce. It is the only family in the package that declares
#' more quantities than it has free values, and the reason is the third group.
#'
#' @details
#' # Why the coefficients are declared
#'
#' \eqn{\phi_1, \dots, \phi_q} are what an autoregression is reported in, and they
#' appear nowhere in the covariance the fit prints, nor among the free values,
#' which are partial autocorrelations. Reporting a coordinate in their place gives
#' a reader the wrong number under the right name: at \eqn{q = 2} with free values
#' \eqn{(0, 0.9, -0.4)} the coefficients are \eqn{(0.988, -0.380)} while the
#' partial autocorrelations are \eqn{(0.716, -0.380)}, so \eqn{\phi_1} and
#' \eqn{\rho_1} differ in the first decimal. At \eqn{q = 1} they coincide exactly,
#' the first partial autocorrelation being the AR(1) coefficient, which is why the
#' distinction is invisible until \eqn{q} passes 1.
#'
#' # The Jacobian is not diagonal
#'
#' A coefficient reads the **whole** chart, so its row of the Jacobian is dense in
#' the partial autocorrelation columns. Those derivatives come out of the
#' Levinson-Durbin recursion the family already propagates derivative arrays
#' through, so they are read off the first-order block instead of recomputed.
#' The variance's row is a single entry, and each partial autocorrelation's is the
#' link's own derivative.
#'
#' # The coefficients' intervals are on the identity scale
#'
#' The stationary region in the coefficients is not a box: it is bounded by the
#' roots of \eqn{1 - \phi_1 z - \cdots - \phi_q z^q} lying outside the unit
#' circle, which no scalar transformation of a single coefficient expresses. An
#' interval that respected the constraint would not be an interval, so the
#' identity scale is the honest choice and a reported coefficient interval may
#' extend outside the stationary region. The partial autocorrelations, whose chart
#' **is** a box, get `"atanh"` and stay inside it.
#'
#' @param s An [AutoregressiveParam()] object.
#' @param eta A numeric vector of free values, of length \eqn{q+1}.
#' @param ... Ignored.
#'
#' @return A list as described in [param_readable()].
#'
#' @name param_readable.AutoregressiveParam
#' @keywords internal
S7::method(param_readable, AutoregressiveParam) <- function(s, eta, ...) {
  q <- s@param_params$order
  n <- q + 1L
  tay <- ar_taylor(s, eta)
  scl <- c(
    linkfunctions7::linkinv(s@param_params$link_scale, eta[1L]),
    linkfunctions7::dlinkinv(s@param_params$link_scale, eta[1L])
  )
  rv <- vapply(seq_len(q), function(k) {
    c(
      linkfunctions7::linkinv(s@param_params$link_pacf, eta[k + 1L]),
      linkfunctions7::dlinkinv(s@param_params$link_pacf, eta[k + 1L])
    )
  }, numeric(2))
  rv <- matrix(rv, nrow = 2L)
  nm <- c("scale", paste0("pacf", seq_len(q)), paste0("phi", seq_len(q)))
  jac <- matrix(0, length(nm), s@n_free, dimnames = list(nm, NULL))
  jac[1L, 1L] <- scl[2L]
  for (k in seq_len(q)) jac[1L + k, 1L + k] <- rv[2L, k]
  for (j in seq_len(q)) jac[1L + q + j, ] <- tay$phi[j, 1L + seq_len(n)]
  list(
    value = stats::setNames(c(scl[1L], rv[1L, ], tay$phi[, 1L]), nm),
    jacobian = jac,
    transform = stats::setNames(
      c("log", rep("atanh", q), rep("identity", q)), nm
    ),
    label = "Autoregressive structure"
  )
}
