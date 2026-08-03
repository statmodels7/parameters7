#' @include numerical_fallbacks.R
NULL


#' Diagonal Structure
#'
#' @description
#' The S7 class of diagonal positive matrices, one free value per entry.
#' Constructed by \code{\link{diag_struct}} or \code{\link{scalar_struct}}.
#'
#' @inheritParams covstruct
#'
#' @return An object of class \code{DiagStruct}.
#'
#' @seealso \code{\link{diag_struct}}, \code{\link{scalar_struct}}
#'
#' @examples
#' S7::S7_inherits(diag_struct(3), DiagStruct)
#'
#' @export
DiagStruct <- S7::new_class("DiagStruct", parent = covstruct)


#' Construct a Diagonal Structure
#'
#' @description
#' A diagonal matrix whose entries are positive, each carried through a scalar
#' link from \pkg{linkfunctions7}.
#'
#' @details
#' This is the case in which the two packages meet, and it is the reason
#' \pkg{covstructs7} composes with \pkg{linkfunctions7} rather than competing
#' with it. The Jacobian of a diagonal block is diagonal, which is exactly the
#' contract a scalar link satisfies, so the link objects are reused as they
#' are and their exact derivatives come with them.
#'
#' A single free value shared by every entry is \code{\link{scalar_struct}}.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param link A \pkg{linkfunctions7} link mapping the free scale to the
#'   positive entries. Defaults to \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{DiagStruct}}.
#'
#' @seealso \code{\link{scalar_struct}}, \code{\link{log_cholesky}}
#'
#' @examples
#' s <- diag_struct(3)
#' round(struct_matrix(s, c(0, 0.5, -0.5)), 4)
#'
#' # any link whose range is positive works
#' round(struct_matrix(diag_struct(2, link = linkfunctions7::sqrt_link()),
#'                     c(1, 2)), 4)
#'
#' @export
diag_struct <- function(dimension, link = linkfunctions7::log_link(),
                        role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_struct_args(dimension, role)
  check_positive_link(link)

  DiagStruct(
    struct_name = "diag",
    dimension = p,
    n_free = p,
    free_names = paste0("d", seq_len(p)),
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    struct_params = list(link = link, shared = FALSE)
  )
}


#' Construct a Scalar Multiple of the Identity
#'
#' @description
#' A diagonal matrix with one free value shared by every entry, so
#' \eqn{M = \tau I} with \eqn{\tau} positive.
#'
#' @details
#' The simplest structure there is, and the one a random effect with a single
#' variance component uses. With the default log link it is the same matrix as
#' \code{scaled_struct(diag(dimension))}, reached from the other direction.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param link A \pkg{linkfunctions7} link mapping the free scale to the
#'   positive scale. Defaults to \code{linkfunctions7::log_link()}.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{DiagStruct}}.
#'
#' @seealso \code{\link{diag_struct}}, \code{\link{scaled_struct}}
#'
#' @examples
#' s <- scalar_struct(3)
#' c(n_free = s@n_free)
#' round(struct_matrix(s, log(2)), 4)
#'
#' @export
scalar_struct <- function(dimension, link = linkfunctions7::log_link(),
                          role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_struct_args(dimension, role)
  check_positive_link(link)

  DiagStruct(
    struct_name = "scalar",
    dimension = p,
    n_free = 1L,
    free_names = "scale",
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    struct_params = list(link = link, shared = TRUE)
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
#' Applies the structure's link, recycling a shared value across the diagonal.
#'
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A numeric vector of length \code{s@dimension}.
#'
#' @keywords internal
diag_entries <- function(s, eta) {
  v <- linkfunctions7::linkinv(s@struct_params$link, eta)
  if (s@struct_params$shared) rep(v, s@dimension) else v
}


#' The Free Value Each Diagonal Entry Belongs To
#'
#' @description
#' The index into the free vector of the value controlling each entry: the
#' identity for an unshared structure, and all ones for a shared one.
#'
#' @param s A \code{\link{DiagStruct}} object.
#'
#' @return An integer vector of length \code{s@dimension}.
#'
#' @keywords internal
diag_owner <- function(s) {
  if (s@struct_params$shared) rep(1L, s@dimension) else seq_len(s@dimension)
}


#' @title Matrix of a Diagonal Structure
#' @name struct_matrix.DiagStruct
#' @description The link applied entrywise, on the diagonal.
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A diagonal positive definite matrix.
#' @keywords internal
S7::method(struct_matrix, DiagStruct) <- function(s, eta, ...) {
  name_dims(diag(diag_entries(s, eta), nrow = s@dimension), s)
}


#' @title Free Vector of a Diagonal Structure
#' @name struct_free.DiagStruct
#' @description
#' The link applied in the forward direction to the diagonal, after refusing a
#' matrix that is not diagonal or not positive.
#' @param s A \code{\link{DiagStruct}} object.
#' @param m A diagonal positive definite matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(struct_free, DiagStruct) <- function(s, m, ...) {
  d <- diag(m)
  off <- m - diag(d, nrow = s@dimension)
  if (max(abs(off)) > 1e-10 * max(1, max(abs(d)))) {
    stop("'m' is not diagonal, so it is not in the set this structure parametrises.",
      call. = FALSE
    )
  }
  if (any(d <= 0)) {
    stop("'m' has a non-positive diagonal entry.", call. = FALSE)
  }
  if (s@struct_params$shared) {
    if (diff(range(d)) > 1e-10 * max(abs(d))) {
      stop(paste0(
        "'m' does not have a constant diagonal, so it is not a scalar\n",
        "  multiple of the identity."
      ), call. = FALSE)
    }
    d <- d[1L]
  }
  stats::setNames(
    linkfunctions7::linkfun(s@struct_params$link, d), s@free_names
  )
}


#' @title First Derivatives of a Diagonal Structure
#' @name struct_dmatrix.DiagStruct
#' @description
#' Closed form from the link's own first derivative: the entries a free value
#' owns carry \eqn{h'(\eta_k)} and everything else is zero.
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of diagonal matrices.
#' @keywords internal
S7::method(struct_dmatrix, DiagStruct) <- function(s, eta, ...) {
  owner <- diag_owner(s)
  d1 <- linkfunctions7::dlinkinv(s@struct_params$link, eta)
  out <- vector("list", s@n_free)
  names(out) <- s@free_names
  for (k in seq_len(s@n_free)) {
    v <- rep(0, s@dimension)
    v[owner == k] <- d1[if (s@struct_params$shared) 1L else k]
    out[[k]] <- name_dims(diag(v, nrow = s@dimension), s)
  }
  out
}


#' @title Second Derivatives of a Diagonal Structure
#' @name struct_d2matrix.DiagStruct
#' @description
#' Closed form from the link's second derivative. The entries are independent
#' unless the structure shares one value, so every cross pair vanishes.
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of diagonal matrices.
#' @keywords internal
S7::method(struct_d2matrix, DiagStruct) <- function(s, eta, ...) {
  owner <- diag_owner(s)
  d2 <- linkfunctions7::d2linkinv(s@struct_params$link, eta)
  idx <- struct_pair_indices(s)
  out <- vector("list", length(idx))
  names(out) <- struct_pair_names(s)
  for (i in seq_along(idx)) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    v <- rep(0, s@dimension)
    if (k == l) {
      v[owner == k] <- d2[if (s@struct_params$shared) 1L else k]
    }
    out[[i]] <- name_dims(diag(v, nrow = s@dimension), s)
  }
  out
}


#' @title Factor of a Diagonal Structure
#' @name struct_factor.DiagStruct
#' @description The square roots of the entries, on the diagonal.
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A diagonal numeric matrix.
#' @keywords internal
S7::method(struct_factor, DiagStruct) <- function(s, eta, ...) {
  diag(sqrt(diag_entries(s, eta)), nrow = s@dimension)
}


#' @title Log-Determinant of a Diagonal Structure
#' @name struct_logdet.DiagStruct
#' @description The sum of the logs of the entries.
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(struct_logdet, DiagStruct) <- function(s, eta, ...) {
  sum(log(diag_entries(s, eta)))
}


#' The Number of Diagonal Entries Each Free Value Owns
#'
#' @description
#' One for an unshared structure, and the whole diagonal for a shared one.
#'
#' @param s A \code{\link{DiagStruct}} object.
#'
#' @return An integer vector of length \code{s@n_free}.
#'
#' @keywords internal
diag_multiplicity <- function(s) {
  if (s@struct_params$shared) s@dimension else rep(1L, s@n_free)
}


#' @title Log-Determinant Gradient of a Diagonal Structure
#' @name struct_dlogdet.DiagStruct
#' @description
#' Closed form: \eqn{m_k h'(\eta_k)/h(\eta_k)}, with \eqn{m_k} the number of
#' diagonal entries the free value owns -- one, or the whole diagonal when the
#' structure shares a single value.
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(struct_dlogdet, DiagStruct) <- function(s, eta, ...) {
  lk <- s@struct_params$link
  h <- linkfunctions7::linkinv(lk, eta)
  d1 <- linkfunctions7::dlinkinv(lk, eta)
  stats::setNames(diag_multiplicity(s) * d1 / h, s@free_names)
}


#' @title Log-Determinant Hessian of a Diagonal Structure
#' @name struct_d2logdet.DiagStruct
#' @description
#' Closed form: \eqn{m_k\{h''/h - (h'/h)^2\}} on the diagonal and zero
#' elsewhere, the entries being controlled independently unless the structure
#' shares one value, in which case there is only one free value.
#' @param s A \code{\link{DiagStruct}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(struct_d2logdet, DiagStruct) <- function(s, eta, ...) {
  lk <- s@struct_params$link
  h <- linkfunctions7::linkinv(lk, eta)
  d1 <- linkfunctions7::dlinkinv(lk, eta)
  d2 <- linkfunctions7::d2linkinv(lk, eta)
  own <- diag_multiplicity(s)
  idx <- struct_pair_indices(s)
  out <- vapply(idx, function(kl) {
    k <- kl[1L]
    if (k != kl[2L]) return(0)
    own[k] * (d2[k] / h[k] - (d1[k] / h[k])^2)
  }, numeric(1))
  stats::setNames(out, struct_pair_names(s))
}
