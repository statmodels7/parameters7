#' @include numerical_fallbacks.R
NULL


#' Diagonal Parameter
#'
#' @description
#' The S7 class of diagonal matrices with positive entries, each entry carried
#' onto the free scale by a \pkg{linkfunctions7} link. Two constructors return
#' it, and the difference is recorded in `param_params$shared`:
#' [diagonal_matrix()] gives one free value per entry, and [scalar_matrix()] one
#' free value for the whole diagonal. Every method here reads that flag, so the
#' two share one class and one set of formulas.
#'
#' It is the family in which \pkg{parameters7} and \pkg{linkfunctions7} meet: the
#' Jacobian of a diagonal map is diagonal, which is exactly what a scalar link
#' supplies, so the link's own exact derivatives to fourth order are used as they
#' are.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `DiagMatrixParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   two entries: `link`, the link object, and `shared`, `TRUE` for
#'   [scalar_matrix()] and `FALSE` for [diagonal_matrix()]. `rank` is always
#'   \eqn{p}, every entry being positive.
#'
#' @seealso [diagonal_matrix()] and [scalar_matrix()], the two constructors, and
#'   [matrix_parameter()] for the properties this inherits.
#'
#' @examples
#' # One class, two constructors, told apart by `shared`.
#' d <- diagonal_matrix(3)
#' q <- scalar_matrix(3)
#' c(S7::S7_inherits(d, DiagMatrixParam), S7::S7_inherits(q, DiagMatrixParam))
#' c(diagonal = d@param_params$shared, scalar = q@param_params$shared)
#' c(diagonal = d@n_free, scalar = q@n_free)
#'
#' @export
DiagMatrixParam <- S7::new_class("DiagMatrixParam", parent = matrix_parameter)


#' Construct a Diagonal Parameter
#'
#' @description
#' Returns an object holding the map \eqn{M = \mathrm{diag}(h(\eta_1), \dots,
#' h(\eta_p))}, a diagonal matrix whose \eqn{p} entries are positive, each
#' carried from one free value by a scalar link. Use it for a covariance that
#' assumes independence, which is \eqn{p} free values against
#' [log_cholesky()]'s \eqn{p(p+1)/2}, and for a set of independent variance
#' components.
#'
#' The choice of link is yours, and it is a real choice: the derivatives the
#' family reports are the link's own, so the curvature an optimizer sees on the
#' free scale changes with it. [scalar_matrix()] is the same object with one free
#' value shared by every entry.
#'
#' @details
#' # Where the two packages meet
#'
#' The Jacobian of a diagonal map is diagonal, which is exactly the contract a
#' scalar link satisfies, so the \pkg{linkfunctions7} objects are reused as they
#' are and their exact derivatives to fourth order come with them. Nothing here
#' rederives a chain rule.
#'
#' # The free names record the link
#'
#' `free_names` is `log_d1`, `log_d2`, ... under the default, and `sqrt_d1`,
#' `sqrt_d2`, ... under a square-root link. A label names the coordinate, never
#' the quantity it produces, so a consumer that flattens the free vector into
#' scalars with identity links reports a number on the scale it is really on.
#'
#' # The log-determinant, and when it is linear
#'
#' \deqn{\log|M| = \sum_{i=1}^{p} \log h(\eta_i),}
#'
#' a sum of functions of one free value each, so every mixed derivative of it is
#' exactly zero at every order. Under the **log** link \eqn{\log h(\eta) = \eta},
#' so the gradient is a vector of ones and the second, third and fourth
#' derivatives all vanish. Under any other link they do not: a square-root link
#' at \eqn{\eta = (1, 2)} gives a second derivative of \eqn{(-2, -0.5)} on the
#' diagonal.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = p}, and \eqn{p} the side of
#' the matrix. \eqn{h = g^{-1}} is the inverse link, carrying a free value onto a
#' positive diagonal entry, and \eqn{h'}, \eqn{h''} its derivatives in \eqn{\eta}.
#'
#' @param dimension The side \eqn{p} of the matrix. A single positive whole
#'   number, finite and at least 1; anything else throws `'dimension' must be a
#'   single positive integer.`
#' @param link A \pkg{linkfunctions7} link carrying the free scale onto the
#'   positive entries, `linkfunctions7::log_link()` by default. It must map
#'   **onto** the positive half line, which [check_positive_link()] enforces by
#'   reading its lower bound, so `identity_link()` is rejected with its bounds in
#'   the message; and it must map **from** the whole real line, since the free
#'   vector is unconstrained by design. `softplus_link()` satisfies both, and so
#'   does a bounded link such as `logit_link()`, whose range \eqn{(0, 1)} is
#'   positive and which gives a positive definite matrix with entries below 1.
#'   `sqrt_link()`, `inverse_link()` and `power_link()` satisfy only the first:
#'   their predictor scale is \eqn{(0, \infty)}, so the map is even in \eqn{\eta}
#'   and the round trip returns \eqn{\lvert \eta \rvert}.
#' @param role A label recording which side of a model the matrix parametrizes:
#'   `"either"` (the default), `"covariance"` or `"precision"`. No numeric result
#'   depends on it; see [log_cholesky()] for what carries it.
#'
#' @return An object of class [DiagMatrixParam()], with `n_free` equal to
#'   `dimension`, `free_names` tagged by the link, `rank` equal to `dimension`,
#'   an empty `null_basis`, `param_name` `"diag"`, and `param_params` holding
#'   `link` and `shared = FALSE`.
#'
#' @seealso [scalar_matrix()] for one shared value, [log_cholesky()] for an
#'   unstructured matrix, [dr_prod()] to give a correlation its own diagonal
#'   scale, and [param_value()] for the map.
#'
#' @examples
#' # Three independent variances, on the log scale by default.
#' s <- diagonal_matrix(3)
#' s@free_names
#' round(param_value(s, c(0, 0.5, -0.5)), 4)
#'
#' # The entries are exp() of the free values, and the round trip closes.
#' all.equal(unname(diag(param_value(s, c(0, 0.5, -0.5)))),
#'           exp(c(0, 0.5, -0.5)))
#' max(abs(param_free(s, param_value(s, c(0.2, -0.3, 0.7))) -
#'         c(0.2, -0.3, 0.7)))
#'
#' # Any link with a non-negative lower bound serves, and the free names say
#' # which one was used.
#' r <- diagonal_matrix(2, link = linkfunctions7::sqrt_link())
#' r@free_names
#' round(param_value(r, c(1, 2)), 4)
#'
#' # A link onto the whole line is refused, an entry having to be positive.
#' try(diagonal_matrix(2, link = linkfunctions7::identity_link()))
#'
#' # The log-determinant is a sum over the entries, and under the log link it
#' # is linear, so its higher derivatives vanish exactly.
#' eta <- c(0, 0.5, -0.5)
#' all.equal(param_logdet(s, eta), sum(eta))
#' param_dlogdet(s, eta)
#' max(abs(param_d2logdet(s, eta)))
#'
#' # Under another link they do not.
#' param_d2logdet(r, c(1, 2))
#'
#' @export
diagonal_matrix <- function(dimension, link = linkfunctions7::log_link(),
                        role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)
  check_positive_link(link)

  DiagMatrixParam(
    param_name = "diag",
    dimension = p,
    n_free = p,
    free_names = tagged_name(link, paste0("d", seq_len(p))),
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(link = link, shared = FALSE)
  )
}


#' Construct a Scalar Multiple of the Identity
#'
#' @description
#' Returns an object holding the map \eqn{M = \tau I} with \eqn{\tau = h(\eta_1)}
#' positive: a diagonal matrix with **one** free value shared by every entry,
#' whatever \eqn{p} is. It is the simplest parametrization in the package, and it
#' is what a random effect with a single variance component needs.
#'
#' @details
#' # One free value, whatever the dimension
#'
#' `n_free` is 1 at every `dimension`, and `free_names` is `log_scale` under the
#' default link. The multiplicity shows up in the log-determinant, which is
#' \eqn{p \log h(\eta_1)}, so the gradient carries a factor of \eqn{p} and a
#' step of one in the free value moves the log-determinant by \eqn{p} under the
#' log link.
#'
#' # The same matrix from the other direction
#'
#' Under the default log link this is `scaled_matrix(diag(dimension))`, a
#' positive multiple of a fixed matrix that happens to be the identity. The two
#' give the same value at the same free value. Prefer this one when the matrix is
#' the identity and [scaled_matrix()] when it is not, and note that
#' [scaled_matrix()] admits a rank-deficient fixed matrix while this family is
#' always of full rank.
#'
#' @section Notation:
#' \eqn{\eta_1} is the single free value, \eqn{\tau = h(\eta_1)} the shared
#' positive entry, \eqn{p} the side of the matrix and \eqn{h = g^{-1}} the
#' inverse link.
#'
#' @param dimension The side \eqn{p} of the matrix. A single positive whole
#'   number, finite and at least 1; anything else throws.
#' @param link A \pkg{linkfunctions7} link carrying the free value onto the
#'   positive scale, `linkfunctions7::log_link()` by default. It must map onto
#'   the positive half line, so `identity_link()` is rejected, and from the whole
#'   real line, which rules out `sqrt_link()` and its relatives; see
#'   [diagonal_matrix()] for the two conditions.
#' @param role A label recording which side of a model the matrix parametrizes:
#'   `"either"` (the default), `"covariance"` or `"precision"`. No numeric result
#'   depends on it.
#'
#' @return An object of class [DiagMatrixParam()], with `n_free` 1,
#'   `free_names` a single link-tagged label, `rank` equal to `dimension`, an
#'   empty `null_basis`, `param_name` `"scalar"`, and `param_params` holding
#'   `link` and `shared = TRUE`.
#'
#' @seealso [diagonal_matrix()] for one free value per entry, [scaled_matrix()]
#'   for a multiple of an arbitrary fixed matrix, and [kron_identity()] to
#'   replicate any parametrization over independent groups.
#'
#' @examples
#' # One free value at any dimension.
#' s <- scalar_matrix(3)
#' c(n_free = s@n_free, dimension = s@dimension)
#' s@free_names
#' param_value(s, log(2))
#'
#' # It is scaled_matrix() on the identity, reached from the other side.
#' all.equal(param_value(s, log(2)),
#'           param_value(scaled_matrix(diag(3)), log(2)),
#'           check.attributes = FALSE)
#'
#' # The shared value enters the log-determinant p times over, so the gradient
#' # carries a factor of p.
#' c(logdet = param_logdet(s, log(2)), p_log_2 = 3 * log(2))
#' param_dlogdet(s, log(2))
#'
#' @export
scalar_matrix <- function(dimension, link = linkfunctions7::log_link(),
                          role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)
  check_positive_link(link)

  DiagMatrixParam(
    param_name = "scalar",
    dimension = p,
    n_free = 1L,
    free_names = tagged_name(link, "scale"),
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(link = link, shared = TRUE)
  )
}


#' Reject a Link That Does Not Reach the Positive Half Line
#'
#' @description
#' Checks that an object is a \pkg{linkfunctions7} link whose declared range lies
#' in the non-negative half line, and signals an error naming the bounds
#' otherwise. Called by [diagonal_matrix()], [scalar_matrix()] and
#' [scaled_matrix()] at construction.
#'
#' @details
#' A diagonal entry of a positive definite matrix is positive, so a link onto the
#' whole real line would let a caller build a matrix outside the set without
#' anything saying so. Reading the link's own `link_bounds` catches it at
#' construction, which is the only place it can be caught: at evaluation time the
#' free vector is unconstrained by design and no value of it is inadmissible.
#'
#' The test is on the **lower** bound alone, `b[1] >= 0`, so a link with a
#' bounded range is accepted: `logit_link()`, whose range is \eqn{(0, 1)},
#' produces a positive definite matrix with entries below 1, which is a
#' legitimate thing to want.
#'
#' @param link The object to check.
#'
#' @return Invisibly `TRUE`. An object that is not a \pkg{linkfunctions7} link
#'   throws `'link' must be a linkfunctions7 link object.`, and a link whose
#'   lower bound is negative throws a message quoting both bounds.
#'
#' @seealso [diagonal_matrix()], [scalar_matrix()] and [scaled_matrix()], the
#'   three callers. The property read is the link object's own `link_bounds`.
#'
#' @keywords internal
check_positive_link <- function(link) {
  if (!S7::S7_inherits(link, linkfunctions7::link)) {
    stop("'link' must be a linkfunctions7 link object.", call. = FALSE)
  }
  b <- link@link_bounds
  if (!is.numeric(b) || length(b) != 2L || b[1] < 0) {
    stop(sprintf(paste0(
      "'link' maps onto (%s, %s), which is not inside the positive half\n",
      "  line. A diagonal entry of a positive definite matrix is positive."
    ), format(b[1]), format(b[2])), call. = FALSE)
  }
  invisible(TRUE)
}


#' The Diagonal Entries Behind a Free Vector
#'
#' @description
#' Applies the parameter's inverse link to the free vector and returns the
#' \eqn{p} diagonal entries, recycling the single value across the whole diagonal
#' for a [scalar_matrix()]. Every method of the family that needs the entries
#' starts here, so the shared and unshared cases branch in one place.
#'
#' @param s A [DiagMatrixParam()] object, whose `param_params$link`,
#'   `param_params$shared` and `dimension` are read.
#' @param eta A numeric vector of free values, of length `s@n_free`: \eqn{p}
#'   values for [diagonal_matrix()], one for [scalar_matrix()].
#'
#' @return A numeric vector of length `s@dimension`, strictly positive under any
#'   link [check_positive_link()] admits.
#'
#' @seealso [diag_owner()] for which free value each entry belongs to, and
#'   [diag_multiplicity()] for how many entries each free value owns.
#'
#' @keywords internal
diag_entries <- function(s, eta) {
  v <- linkfunctions7::linkinv(s@param_params$link, eta)
  if (s@param_params$shared) rep(v, s@dimension) else v
}


#' The Free Value Each Diagonal Entry Belongs To
#'
#' @description
#' Returns, for each of the \eqn{p} diagonal entries, the index into the free
#' vector of the value that controls it: `1:p` for a [diagonal_matrix()], and
#' `rep(1, p)` for a [scalar_matrix()], whose one value controls them all. The
#' derivative methods use it to place a link derivative in the right entries.
#'
#' @param s A [DiagMatrixParam()] object, whose `param_params$shared` and
#'   `dimension` are read.
#'
#' @return An integer vector of length `s@dimension`, with values in
#'   `1:s@n_free`.
#'
#' @seealso [diag_multiplicity()], the same information counted the other way,
#'   and [diag_entries()].
#'
#' @keywords internal
diag_owner <- function(s) {
  if (s@param_params$shared) rep(1L, s@dimension) else seq_len(s@dimension)
}


#' @title Matrix of a Diagonal Parameter
#' @name param_value.DiagMatrixParam
#' @description
#' Returns \eqn{M = \mathrm{diag}(h(\eta_1), \dots, h(\eta_p))}, the inverse link
#' applied to each free value and placed on the diagonal. For a
#' [scalar_matrix()] the one free value is recycled, giving \eqn{h(\eta_1) I}.
#' Positive definiteness follows from the link, whose range
#' [check_positive_link()] restricted to the non-negative half line at
#' construction, so nothing is tested here.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A diagonal positive definite `s@dimension` by `s@dimension` numeric
#'   matrix with dimnames `v1`, `v2`, ...
#' @seealso [diag_entries()], which applies the link, and
#'   [param_free.DiagMatrixParam()] for the inverse.
#' @keywords internal
S7::method(param_value, DiagMatrixParam) <- function(s, eta, ...) {
  name_dims(diag(diag_entries(s, eta), nrow = s@dimension), s)
}


#' @title Free Vector of a Diagonal Parameter
#' @name param_free.DiagMatrixParam
#' @description
#' Returns \eqn{g} applied to the diagonal of `m`, the link in the forward
#' direction, after checking that `m` really is in the family's set. Exact, and a
#' true inverse of [param_value.DiagMatrixParam()], the link being a bijection.
#' @details
#' Three rejections, each with its own message. `m` must be diagonal, tested as
#' `max(abs(off)) > 1e-10 * max(1, max(abs(d)))` on the off-diagonal part, so an
#' asymmetry of rounding size passes and a real off-diagonal entry does not. Its
#' diagonal entries must all be positive. And for a [scalar_matrix()] they must
#' all be **equal**, tested by `diff(range(d)) > 1e-10 * max(abs(d))`, a diagonal
#' that varies not being a scalar multiple of the identity; only the first entry
#' is then read.
#' @param s A [DiagMatrixParam()] object.
#' @param m A diagonal positive definite `s@dimension` by `s@dimension` numeric
#'   matrix, already checked for shape and symmetry by the generic. It must have
#'   a constant diagonal when `s` came from [scalar_matrix()].
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.DiagMatrixParam()], the map this inverts.
#' @keywords internal
S7::method(param_free, DiagMatrixParam) <- function(s, m, ...) {
  d <- diag(m)
  off <- m - diag(d, nrow = s@dimension)
  if (max(abs(off)) > 1e-10 * max(1, max(abs(d)))) {
    stop("'m' is not diagonal, so it is not in the set this parameter parametrizes.",
      call. = FALSE
    )
  }
  if (any(d <= 0)) {
    stop("'m' has a non-positive diagonal entry.", call. = FALSE)
  }
  if (s@param_params$shared) {
    if (diff(range(d)) > 1e-10 * max(abs(d))) {
      stop(paste0(
        "'m' does not have a constant diagonal, so it is not a scalar\n",
        "  multiple of the identity."
      ), call. = FALSE)
    }
    d <- d[1L]
  }
  stats::setNames(
    linkfunctions7::linkfun(s@param_params$link, d), s@free_names
  )
}


#' @title First Derivatives of a Diagonal Parameter
#' @name param_d1.DiagMatrixParam
#' @description
#' Closed form, and it is the link's own derivative placed on a diagonal:
#'
#' \deqn{\partial_k M = h'(\eta_k)\, \textstyle\sum_{i:\, o_i = k} E_{ii},}
#'
#' where \eqn{o_i} is the free value entry \eqn{i} belongs to. For a
#' [diagonal_matrix()] that is one entry, so \eqn{\partial_k M} is
#' \eqn{h'(\eta_k) E_{kk}}; for a [scalar_matrix()] the one free value owns the
#' whole diagonal and \eqn{\partial_1 M = h'(\eta_1) I}.
#'
#' \eqn{h'} comes from `linkfunctions7::dlinkinv()`, so its accuracy is the
#' link's and nothing is differenced.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` diagonal matrices named by `s@free_names`, each
#'   `s@dimension` by `s@dimension` with dimnames `v1`, `v2`, ...
#' @seealso [diag_owner()] for the ownership map, and
#'   [param_d2.DiagMatrixParam()] for the order above.
#' @keywords internal
S7::method(param_d1, DiagMatrixParam) <- function(s, eta, ...) {
  owner <- diag_owner(s)
  d1 <- linkfunctions7::dlinkinv(s@param_params$link, eta)
  out <- vector("list", s@n_free)
  names(out) <- s@free_names
  for (k in seq_len(s@n_free)) {
    v <- rep(0, s@dimension)
    v[owner == k] <- d1[if (s@param_params$shared) 1L else k]
    out[[k]] <- name_dims(diag(v, nrow = s@dimension), s)
  }
  out
}


#' @title Second Derivatives of a Diagonal Parameter
#' @name param_d2.DiagMatrixParam
#' @description
#' Closed form, from `linkfunctions7::d2linkinv()`. The family is **separable**:
#' each diagonal entry is a function of one free value alone, so a component
#' \eqn{\partial_{kl} M} with \eqn{k \ne l} is exactly the zero matrix, and a
#' pure one carries \eqn{h''(\eta_k)} in the entries that value owns. A
#' [scalar_matrix()] has one free value and so one component.
#'
#' The same separability holds at third and fourth order, which is why those
#' methods keep only the pure components too.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` diagonal matrices keyed as
#'   `param_tuple_names(s)` and in that order, the mixed ones exactly zero.
#' @seealso [param_d1.DiagMatrixParam()] and [param_d3.DiagMatrixParam()] for the
#'   neighboring orders.
#' @keywords internal
S7::method(param_d2, DiagMatrixParam) <- function(s, eta, ...) {
  owner <- diag_owner(s)
  d2 <- linkfunctions7::d2linkinv(s@param_params$link, eta)
  idx <- param_tuple_indices(s)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s)
  for (i in seq_along(idx)) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    v <- rep(0, s@dimension)
    if (k == l) {
      v[owner == k] <- d2[if (s@param_params$shared) 1L else k]
    }
    out[[i]] <- name_dims(diag(v, nrow = s@dimension), s)
  }
  out
}


#' @title Factor of a Diagonal Parameter
#' @name param_factor.DiagMatrixParam
#' @description
#' Returns \eqn{L = \mathrm{diag}(\sqrt{h(\eta_1)}, \dots, \sqrt{h(\eta_p)})}.
#' The Cholesky factor of a diagonal matrix is the diagonal of its square roots,
#' so this is \eqn{p} square roots and no factorization, against the base class's
#' \eqn{O(p^3)}.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` diagonal numeric matrix with a
#'   positive diagonal, satisfying `L %*% t(L) == param_value(s, eta)`, and
#'   carrying no dimnames.
#' @seealso [param_factor()] for the generic and
#'   [param_factor.matrix_parameter()] for what a family without a closed form
#'   pays.
#' @keywords internal
S7::method(param_factor, DiagMatrixParam) <- function(s, eta, ...) {
  diag(sqrt(diag_entries(s, eta)), nrow = s@dimension)
}


#' @title Log-Determinant of a Diagonal Parameter
#' @name param_logdet.DiagMatrixParam
#' @description
#' Returns \eqn{\log|M| = \sum_{i=1}^{p} \log h(\eta_i)}, the sum of the
#' logarithms of the diagonal entries, which for a diagonal matrix is the
#' log-determinant exactly. For a [scalar_matrix()] the one entry is counted
#' \eqn{p} times, giving \eqn{p \log h(\eta_1)}.
#'
#' Under the **log** link \eqn{\log h(\eta) = \eta}, so the result is `sum(eta)`
#' and the quantity is linear in the free vector; under any other link it is not.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number.
#' @seealso [param_dlogdet.DiagMatrixParam()] for its gradient, and [diag_dlog()]
#'   for the derivatives of \eqn{\log h}.
#' @keywords internal
S7::method(param_logdet, DiagMatrixParam) <- function(s, eta, ...) {
  sum(log(diag_entries(s, eta)))
}


#' The Number of Diagonal Entries Each Free Value Owns
#'
#' @description
#' Returns \eqn{m_k}, the number of diagonal entries the \eqn{k}-th free value
#' controls: 1 for every value of a [diagonal_matrix()], and \eqn{p} for the
#' single value of a [scalar_matrix()]. It is the multiplicity that appears in
#' the log-determinant's derivatives, \eqn{\log|M|} counting a shared value once
#' per entry.
#'
#' @param s A [DiagMatrixParam()] object, whose `param_params$shared`,
#'   `dimension` and `n_free` are read.
#'
#' @return An integer vector of length `s@n_free` summing to `s@dimension`.
#'
#' @seealso [diag_owner()], the same information per entry, and
#'   [param_dlogdet.DiagMatrixParam()], which multiplies by it.
#'
#' @keywords internal
diag_multiplicity <- function(s) {
  if (s@param_params$shared) s@dimension else rep(1L, s@n_free)
}


#' @title Log-Determinant Gradient of a Diagonal Parameter
#' @name param_dlogdet.DiagMatrixParam
#' @description
#' Closed form,
#'
#' \deqn{\partial_k \log|M| = m_k \frac{h'(\eta_k)}{h(\eta_k)},}
#'
#' with \eqn{m_k} the number of diagonal entries the \eqn{k}-th free value owns:
#' 1 for a [diagonal_matrix()], and \eqn{p} for the single value of a
#' [scalar_matrix()], which enters the log-determinant once per entry.
#'
#' Under the log link \eqn{h'/h = 1}, so the answer is a vector of ones for a
#' [diagonal_matrix()] and the scalar \eqn{p} for a [scalar_matrix()].
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [diag_multiplicity()] for \eqn{m_k}, and
#'   [param_d2logdet.DiagMatrixParam()] for the order above.
#' @keywords internal
S7::method(param_dlogdet, DiagMatrixParam) <- function(s, eta, ...) {
  lk <- s@param_params$link
  h <- linkfunctions7::linkinv(lk, eta)
  d1 <- linkfunctions7::dlinkinv(lk, eta)
  stats::setNames(diag_multiplicity(s) * d1 / h, s@free_names)
}


#' @title Log-Determinant Hessian of a Diagonal Parameter
#' @name param_d2logdet.DiagMatrixParam
#' @description
#' Closed form. The log-determinant is a sum of \eqn{\log h(\eta_k)} terms, one
#' free value each, so every mixed component is exactly zero and a pure one is
#'
#' \deqn{\partial_{kk} \log|M| = m_k\left\{\frac{h''(\eta_k)}{h(\eta_k)}
#'   - \left(\frac{h'(\eta_k)}{h(\eta_k)}\right)^{2}\right\},}
#'
#' the second derivative of \eqn{\log h} times the number of entries the value
#' owns.
#'
#' Under the **log** link this is identically zero, \eqn{\log h(\eta) = \eta}
#' being linear, so a check of this quantity against a numerical reference on a
#' default `diagonal_matrix()` compares two zeros and would pass whatever was
#' missing. Under a square-root link at \eqn{\eta = (1, 2)} it is
#' \eqn{(-2, -0.5)}, which is where the formula has content.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
#'   `param_tuple_names(s)` and in that order, the mixed ones exactly zero.
#' @seealso [param_dlogdet.DiagMatrixParam()] for the order below, and
#'   [diag_logdet_higher()] for orders three and four.
#' @keywords internal
S7::method(param_d2logdet, DiagMatrixParam) <- function(s, eta, ...) {
  lk <- s@param_params$link
  h <- linkfunctions7::linkinv(lk, eta)
  d1 <- linkfunctions7::dlinkinv(lk, eta)
  d2 <- linkfunctions7::d2linkinv(lk, eta)
  own <- diag_multiplicity(s)
  idx <- param_tuple_indices(s)
  out <- vapply(idx, function(kl) {
    k <- kl[1L]
    if (k != kl[2L]) return(0)
    own[k] * (d2[k] / h[k] - (d1[k] / h[k])^2)
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s))
}


#' Derivatives of the Logarithm of a Link
#'
#' @description
#' Returns \eqn{\mathrm{d}^m \log h(\eta) / \mathrm{d}\eta^m} for \eqn{m} up to
#' four, assembled from the inverse link's own derivatives. Every
#' log-determinant derivative of a diagonal or scaled family is one of these,
#' multiplied by a count of entries.
#'
#' @details
#' Faa di Bruno's formula for the logarithm, written out. With
#' \eqn{u_m = h^{(m)}/h} the four orders are
#'
#' \deqn{u_1, \quad u_2 - u_1^2, \quad u_3 - 3u_1u_2 + 2u_1^3, \quad
#'   u_4 - 4u_1u_3 - 3u_2^2 + 12u_1^2u_2 - 6u_1^4.}
#'
#' The division by \eqn{h} happens once, at the start, so every term is a ratio
#' of order one, never a derivative that can be large on its own. That matters
#' at the ends of a link's range, where \eqn{h} and \eqn{h^{(m)}} can both be
#' extreme while the ratio is ordinary.
#'
#' Under the log link \eqn{h = e^\eta} gives \eqn{u_m = 1} for every \eqn{m}, so
#' the first order is 1 and the second, third and fourth are exactly 0, which is
#' why a default `diagonal_matrix()` reports a linear log-determinant.
#'
#' @param link A \pkg{linkfunctions7} link.
#' @param eta A numeric vector of free values.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A numeric vector the length of `eta`.
#'
#' @seealso [diag_logdet_higher()] and [scaled_dlog()], the two callers, and
#'   `linkfunctions7::dlinkinv()` for the link derivatives read.
#'
#' @keywords internal
diag_dlog <- function(link, eta, order) {
  h <- linkfunctions7::linkinv(link, eta)
  u1 <- linkfunctions7::dlinkinv(link, eta) / h
  if (order == 1L) return(u1)
  u2 <- linkfunctions7::d2linkinv(link, eta) / h
  if (order == 2L) return(u2 - u1^2)
  u3 <- linkfunctions7::d3linkinv(link, eta) / h
  if (order == 3L) return(u3 - 3 * u1 * u2 + 2 * u1^3)
  u4 <- linkfunctions7::d4linkinv(link, eta) / h
  u4 - 4 * u1 * u3 - 3 * u2^2 + 12 * u1^2 * u2 - 6 * u1^4
}

#' Third and Fourth Derivatives of a Diagonal Matrix
#'
#' @description
#' Assembles a whole derivative order of a diagonal family, orders three and
#' four. A diagonal family is separable, so a component is the zero matrix unless
#' every index of the tuple names the **same** free value, and a surviving one
#' carries \eqn{h'''(\eta_k)} or \eqn{h''''(\eta_k)} in the entries that value
#' owns. Most of the list is therefore exactly zero: at \eqn{p = 3} only 3 of the
#' 10 third-order components are non-zero.
#'
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` diagonal matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension` with dimnames `v1`, `v2`, ...
#'
#' @seealso [diag_owner()] for the ownership map, and
#'   [param_d3.DiagMatrixParam()] and [param_d4.DiagMatrixParam()], the two
#'   callers.
#'
#' @keywords internal
diag_higher <- function(s, eta, order) {
  idx <- param_tuple_indices(s, order)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, order)
  link <- s@param_params$link
  owner <- diag_owner(s)
  dfun <- switch(order - 2L, linkfunctions7::d3linkinv, linkfunctions7::d4linkinv)
  for (i in seq_along(idx)) {
    t <- idx[[i]]
    m <- matrix(0, s@dimension, s@dimension)
    if (all(t == t[1L])) {
      k <- t[1L]
      v <- dfun(link, eta[k])
      d <- ifelse(owner == k, v, 0)
      m <- diag(d, nrow = s@dimension)
    }
    out[[i]] <- name_dims(m, s)
  }
  out
}

#' @title Third Derivatives of a Diagonal Parameter
#' @name param_d3.DiagMatrixParam
#' @description
#' Closed form. A diagonal entry depends on one free value through its link, so
#' the only surviving components are the pure ones, and each carries
#' \eqn{h'''(\eta_k)} on the entries that value owns. Every mixed component is
#' the exact zero matrix.
#'
#' \eqn{h'''} comes from `linkfunctions7::d3linkinv()`, so the accuracy is the
#' link's; a family without a closed form here would get a product stencil good
#' to about six digits instead.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` diagonal matrices keyed as
#'   `param_tuple_names(s, 3)` and in that order.
#' @seealso [diag_higher()], which assembles it, and
#'   [param_d4.DiagMatrixParam()] for the order above.
#' @keywords internal
S7::method(param_d3, DiagMatrixParam) <- function(s, eta, ...) {
  diag_higher(s, eta, 3L)
}

#' @title Fourth Derivatives of a Diagonal Parameter
#' @name param_d4.DiagMatrixParam
#' @description
#' Closed form, the same separable structure as at third order with
#' \eqn{h''''(\eta_k)} in place of \eqn{h'''(\eta_k)}. Only the pure components
#' survive; every tuple naming two different free values is the exact zero
#' matrix.
#'
#' \eqn{h''''} comes from `linkfunctions7::d4linkinv()`. This is the order at
#' which a numerical route is least usable, keeping about five digits, so the
#' closed form matters most here.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` diagonal matrices keyed as
#'   `param_tuple_names(s, 4)` and in that order.
#' @seealso [diag_higher()], which assembles it, [param_d3.DiagMatrixParam()] for
#'   the order below, and [numerical_d4()] for the alternative.
#' @keywords internal
S7::method(param_d4, DiagMatrixParam) <- function(s, eta, ...) {
  diag_higher(s, eta, 4L)
}

#' Higher Derivatives of a Diagonal Log-Determinant
#'
#' @description
#' Assembles the derivative components of orders two to four of \eqn{\log|M|} for
#' a diagonal family. The log-determinant is \eqn{\sum_i \log h(\eta_i)}, a sum
#' of functions of one free value each, so every mixed component is exactly zero
#' and a pure one is the matching derivative from [diag_dlog()], counted once per
#' entry the free value owns.
#'
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 2, 3 or 4.
#'
#' @return A numeric vector of `choose(s@n_free + order - 1, order)` entries,
#'   keyed as `param_tuple_names(s, order)` and in that order. Under the log link
#'   it is all zeros at every order above 1, \eqn{\log h(\eta) = \eta} being
#'   linear.
#'
#' @seealso [diag_dlog()] for the derivatives of \eqn{\log h}, and
#'   [param_d3logdet.DiagMatrixParam()] and [param_d4logdet.DiagMatrixParam()],
#'   the two callers.
#'
#' @keywords internal
diag_logdet_higher <- function(s, eta, order) {
  idx <- param_tuple_indices(s, order)
  owner <- diag_owner(s)
  link <- s@param_params$link
  out <- vapply(seq_along(idx), function(i) {
    t <- idx[[i]]
    if (!all(t == t[1L])) return(0)
    k <- t[1L]
    sum(owner == k) * diag_dlog(link, eta[k], order)
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, order))
}

#' @title Third Log-Determinant Derivatives of a Diagonal Parameter
#' @name param_d3logdet.DiagMatrixParam
#' @description
#' Closed form. The log-determinant is a sum of \eqn{\log h(\eta_k)} terms, so
#' each pure component is the third derivative of \eqn{\log h} at that free
#' value, times the number of entries the value owns, and every mixed component
#' is exactly zero.
#'
#' Under the log link the whole vector is zero, \eqn{\log h(\eta) = \eta} being
#' linear. A link with curvature is where this order has content.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
#'   `param_tuple_names(s, 3)`.
#' @seealso [diag_logdet_higher()], which assembles it, and [diag_dlog()] for the
#'   derivatives of \eqn{\log h}.
#' @keywords internal
S7::method(param_d3logdet, DiagMatrixParam) <- function(s, eta, ...) {
  diag_logdet_higher(s, eta, 3L)
}

#' @title Fourth Log-Determinant Derivatives of a Diagonal Parameter
#' @name param_d4logdet.DiagMatrixParam
#' @description
#' Closed form, by the same separability as at third order: each pure component
#' is the fourth derivative of \eqn{\log h} at that free value, times the number
#' of entries it owns, and every mixed component is exactly zero. The fourth
#' derivative of \eqn{\log h} is
#' \eqn{u_4 - 4u_1u_3 - 3u_2^2 + 12u_1^2u_2 - 6u_1^4} with
#' \eqn{u_m = h^{(m)}/h}, written out in [diag_dlog()].
#'
#' Under the log link the whole vector is zero. This is the order where a
#' numerical route is least usable, so the closed form is worth most here.
#' @param s A [DiagMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
#'   `param_tuple_names(s, 4)`.
#' @seealso [diag_logdet_higher()], which assembles it, [diag_dlog()] for the
#'   formula, and [param_d4logdet.matrix_parameter()] for the numerical route.
#' @keywords internal
S7::method(param_d4logdet, DiagMatrixParam) <- function(s, eta, ...) {
  diag_logdet_higher(s, eta, 4L)
}
