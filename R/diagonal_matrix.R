#' @include numerical_fallbacks.R
NULL


#' Diagonal Parameter
#'
#' @description
#' The S7 class of diagonal positive matrices, one free value per entry.
#' Constructed by \code{\link{diagonal_matrix}} or \code{\link{scalar_matrix}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{DiagMatrixParam}.
#'
#' @seealso \code{\link{diagonal_matrix}}, \code{\link{scalar_matrix}}
#'
#' @examples
#' S7::S7_inherits(diagonal_matrix(3), DiagMatrixParam)
#'
#' @export
DiagMatrixParam <- S7::new_class("DiagMatrixParam", parent = matrix_parameter)


#' Construct a Diagonal Parameter
#'
#' @description
#' A diagonal matrix whose entries are positive, each carried through a scalar
#' link from \pkg{linkfunctions7}.
#'
#' @details
#' This is the case in which the two packages meet, and it is the reason
#' \pkg{parameters7} composes with \pkg{linkfunctions7} rather than competing
#' with it. The Jacobian of a diagonal block is diagonal, which is exactly the
#' contract a scalar link satisfies, so the link objects are reused as they
#' are and their exact derivatives come with them.
#'
#' A single free value shared by every entry is \code{\link{scalar_matrix}}.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param link A \pkg{linkfunctions7} link mapping the free scale to the
#'   positive entries. Defaults to \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{DiagMatrixParam}}.
#'
#' @seealso \code{\link{scalar_matrix}}, \code{\link{log_cholesky}}
#'
#' @examples
#' s <- diagonal_matrix(3)
#' round(param_value(s, c(0, 0.5, -0.5)), 4)
#'
#' # any link whose range is positive works
#' round(param_value(diagonal_matrix(2, link = linkfunctions7::sqrt_link()),
#'                     c(1, 2)), 4)
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
    free_names = paste0("d", seq_len(p)),
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(link = link, shared = FALSE)
  )
}


#' Construct a Scalar Multiple of the Identity
#'
#' @description
#' A diagonal matrix with one free value shared by every entry, so
#' \eqn{M = \tau I} with \eqn{\tau} positive.
#'
#' @details
#' The simplest parameter there is, and the one a random effect with a single
#' variance component uses. With the default log link it is the same matrix as
#' \code{scaled_matrix(diag(dimension))}, reached from the other direction.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param link A \pkg{linkfunctions7} link mapping the free scale to the
#'   positive scale. Defaults to \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{DiagMatrixParam}}.
#'
#' @seealso \code{\link{diagonal_matrix}}, \code{\link{scaled_matrix}}
#'
#' @examples
#' s <- scalar_matrix(3)
#' c(n_free = s@n_free)
#' round(param_value(s, log(2)), 4)
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
    free_names = "scale",
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(link = link, shared = TRUE)
  )
}


#' Refuse a Link That Does Not Reach the Positive Half Line
#'
#' @description
#' Checks that an object is a \pkg{linkfunctions7} link whose bounds are
#' contained in \eqn{(0, \infty)}.
#'
#' @details
#' A diagonal entry of a positive definite matrix is positive, so a link onto
#' the whole real line would let a caller build a matrix outside the set
#' silently. Asking the link for its declared bounds catches it at
#' construction, which is the only place it can be caught.
#'
#' @param link The object to check.
#'
#' @return Invisibly \code{TRUE}; raises an error otherwise.
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
#' Applies the parameter's link, recycling a shared value across the diagonal.
#'
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A numeric vector of length \code{s@dimension}.
#'
#' @keywords internal
diag_entries <- function(s, eta) {
  v <- linkfunctions7::linkinv(s@param_params$link, eta)
  if (s@param_params$shared) rep(v, s@dimension) else v
}


#' The Free Value Each Diagonal Entry Belongs To
#'
#' @description
#' The index into the free vector of the value controlling each entry: the
#' identity for an unshared parameter, and all ones for a shared one.
#'
#' @param s A \code{\link{DiagMatrixParam}} object.
#'
#' @return An integer vector of length \code{s@dimension}.
#'
#' @keywords internal
diag_owner <- function(s) {
  if (s@param_params$shared) rep(1L, s@dimension) else seq_len(s@dimension)
}


#' @title Matrix of a Diagonal Parameter
#' @name param_value.DiagMatrixParam
#' @description The link applied entrywise, on the diagonal.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A diagonal positive definite matrix.
#' @keywords internal
S7::method(param_value, DiagMatrixParam) <- function(s, eta, ...) {
  name_dims(diag(diag_entries(s, eta), nrow = s@dimension), s)
}


#' @title Free Vector of a Diagonal Parameter
#' @name param_free.DiagMatrixParam
#' @description
#' The link applied in the forward direction to the diagonal, after refusing a
#' matrix that is not diagonal or not positive.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param m A diagonal positive definite matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(param_free, DiagMatrixParam) <- function(s, m, ...) {
  d <- diag(m)
  off <- m - diag(d, nrow = s@dimension)
  if (max(abs(off)) > 1e-10 * max(1, max(abs(d)))) {
    stop("'m' is not diagonal, so it is not in the set this parameter parametrises.",
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
#' Closed form from the link's own first derivative: the entries a free value
#' owns carry \eqn{h'(\eta_k)} and everything else is zero.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of diagonal matrices.
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
#' Closed form from the link's second derivative. The entries are independent
#' unless the parameter shares one value, so every cross pair vanishes.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of diagonal matrices.
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
#' @description The square roots of the entries, on the diagonal.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A diagonal numeric matrix.
#' @keywords internal
S7::method(param_factor, DiagMatrixParam) <- function(s, eta, ...) {
  diag(sqrt(diag_entries(s, eta)), nrow = s@dimension)
}


#' @title Log-Determinant of a Diagonal Parameter
#' @name param_logdet.DiagMatrixParam
#' @description The sum of the logs of the entries.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, DiagMatrixParam) <- function(s, eta, ...) {
  sum(log(diag_entries(s, eta)))
}


#' The Number of Diagonal Entries Each Free Value Owns
#'
#' @description
#' One for an unshared parameter, and the whole diagonal for a shared one.
#'
#' @param s A \code{\link{DiagMatrixParam}} object.
#'
#' @return An integer vector of length \code{s@n_free}.
#'
#' @keywords internal
diag_multiplicity <- function(s) {
  if (s@param_params$shared) s@dimension else rep(1L, s@n_free)
}


#' @title Log-Determinant Gradient of a Diagonal Parameter
#' @name param_dlogdet.DiagMatrixParam
#' @description
#' Closed form: \eqn{m_k h'(\eta_k)/h(\eta_k)}, with \eqn{m_k} the number of
#' diagonal entries the free value owns -- one, or the whole diagonal when the
#' parameter shares a single value.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
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
#' Closed form: \eqn{m_k\{h''/h - (h'/h)^2\}} on the diagonal and zero
#' elsewhere, the entries being controlled independently unless the parameter
#' shares one value, in which case there is only one free value.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
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


# d^m of log(h(eta)) from the link's own derivatives, orders 1 to 4: the
# Faa di Bruno chain of the logarithm.
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
#' Closed form: a diagonal entry depends on one free value through its link,
#' so the only surviving components are the pure ones, carrying
#' \eqn{h'''(\eta_k)} on the entries the value owns.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of diagonal matrices.
#' @keywords internal
S7::method(param_d3, DiagMatrixParam) <- function(s, eta, ...) {
  diag_higher(s, eta, 3L)
}

#' @title Fourth Derivatives of a Diagonal Parameter
#' @name param_d4.DiagMatrixParam
#' @description Closed form, with \eqn{h''''} in place of \eqn{h'''}.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of diagonal matrices.
#' @keywords internal
S7::method(param_d4, DiagMatrixParam) <- function(s, eta, ...) {
  diag_higher(s, eta, 4L)
}

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

#' @title Higher Log-Determinant Derivatives of a Diagonal Parameter
#' @name param_d3logdet.DiagMatrixParam
#' @description
#' Closed form: the log-determinant is a sum of \eqn{\log h(\eta_k)} terms,
#' so each pure component is the matching derivative of \eqn{\log h} times
#' the number of entries the value owns, and every mixed component is zero.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d3logdet, DiagMatrixParam) <- function(s, eta, ...) {
  diag_logdet_higher(s, eta, 3L)
}

#' @title Fourth Log-Determinant Derivatives of a Diagonal Parameter
#' @name param_d4logdet.DiagMatrixParam
#' @description Closed form; see \code{\link{param_d3logdet.DiagMatrixParam}}.
#' @param s A \code{\link{DiagMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_d4logdet, DiagMatrixParam) <- function(s, eta, ...) {
  diag_logdet_higher(s, eta, 4L)
}
