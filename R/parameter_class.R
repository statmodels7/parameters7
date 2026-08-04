#' Constrained Matrix Parameter
#'
#' @description
#' The abstract S7 class of covariance parameters: a map from an unconstrained
#' vector \eqn{\eta \in \mathbb{R}^{d}} to a symmetric matrix in some
#' constrained set, together with its derivatives and the quantities a
#' likelihood asks of the matrix.
#'
#' @details
#' A parameter owns its dimension: \code{log_cholesky(dimension = 3)} and
#' \code{log_cholesky(dimension = 4)} are different objects, with \eqn{d = 6} and
#' \eqn{d = 10} free values respectively. The alternative, a dimensionless
#' recipe applied to whatever arrives, would leave \code{n_free} and
#' \code{free_names} unanswerable before any data exist, and both are needed
#' then.
#'
#' The rank and the null space are properties of the family rather than of a
#' point: scaling a matrix by a positive number and summing positive
#' semidefinite matrices both leave the null space alone. They are computed
#' once, at construction, from the components rather than from an assembled
#' matrix, because the numerical determination of a rank from an assembled
#' matrix is not scale invariant while the null space is. See
#' \code{\link{param_null_basis}}.
#'
#' Only \code{\link{param_value}} is compulsory. Every other generic has a
#' numerical method registered on this class, so a new parameter is a subclass
#' and one method, and a closed form supplied later replaces the corresponding
#' numerical method through dispatch.
#'
#' @param param_name A single character string naming the family.
#' @param n_free The length \eqn{d} of the free vector.
#' @param free_names A character vector of length \code{n_free}, one label per
#'   free value. Fixed at construction: every consumer builds parameter tables
#'   from these. A label names the coordinate rather than the quantity the
#'   coordinate produces, and the families here follow one convention for it.
#'   Where a link carries a constrained quantity onto the free scale, the
#'   label records that link, so a variance appears as \code{"log_scale"} and
#'   a correlation as \code{"z_rho"}; where the coordinate is already
#'   unrestricted the label is the plain name of the quantity, as the
#'   below-diagonal entries \code{"L2.1"} of a Cholesky factor are. The
#'   distinction matters outside the family: a consumer flattens the free
#'   vector into scalar parameters carrying identity links, so a label
#'   promising a bounded quantity reports a number on a scale that number is
#'   not on.
#' @param param_params A list of whatever the family needs to evaluate itself.
#'
#' @return An object of class \code{parameter}. The class is abstract; use one
#'   of the constructors, such as \code{\link{log_cholesky}}.
#'
#' @seealso \code{\link{log_cholesky}}, \code{\link{scaled_matrix}},
#'   \code{\link{param_value}}, \code{\link{check_parameter}}
#'
#' @examples
#' s <- log_cholesky(3)
#' S7::S7_inherits(s, parameter)
#' c(dimension = s@dimension, n_free = s@n_free, rank = s@rank)
#'
#' @export
parameter <- S7::new_class(
  "parameter",
  properties = list(
    param_name = S7::class_character,
    n_free = S7::class_integer,
    free_names = S7::class_character,
    param_params = S7::class_list
  ),
  validator = function(self) {
    errors <- character()
    if (length(self@n_free) != 1L || self@n_free < 0L) {
      errors <- c(errors, "@n_free must be a single non-negative integer")
    }
    if (length(self@free_names) != self@n_free) {
      errors <- c(errors, "@free_names must have one entry per free value")
    }
    if (anyDuplicated(self@free_names)) {
      errors <- c(errors, "@free_names must be unique")
    }
    if (length(errors)) errors else NULL
  }
)


#' Constrained Symmetric Matrix Parameter
#'
#' @description
#' The abstract S7 class of the symmetric positive semidefinite branch: a
#' \code{\link{parameter}} whose value is a symmetric matrix, together with
#' the quantities only a matrix can answer -- the rank, the null space, the
#' log-(pseudo-)determinant, the solve and the factor.
#'
#' @details
#' The rank and the null space are properties of the family rather than of a
#' point: scaling a matrix by a positive number and summing positive
#' semidefinite matrices both leave the null space alone. They are computed
#' once, at construction, from the components rather than from an assembled
#' matrix, because the numerical determination of a rank from an assembled
#' matrix is not scale invariant while the null space is. See
#' \code{\link{param_null_basis}}.
#'
#' A parameter that is not a matrix -- \code{\link{simplex}}, a
#' \code{\link{transition_matrix}} -- inherits from \code{\link{parameter}}
#' directly, so \code{\link{param_logdet}} and \code{\link{param_solve}} do
#' not exist for it by construction rather than by a run-time refusal.
#'
#' @param param_name A single character string naming the family.
#' @param n_free The length \eqn{d} of the free vector.
#' @param free_names A character vector of length \code{n_free}.
#' @param param_params A list of whatever the family needs to evaluate itself.
#' @param dimension The side \eqn{p} of the matrix.
#' @param rank The rank of the matrix the family produces.
#' @param null_basis A \code{dimension} by \code{dimension - rank} matrix
#'   whose columns are an orthonormal basis of the null space.
#' @param role One of \code{"covariance"}, \code{"precision"} or
#'   \code{"either"}. A label: no method reads it and no result depends on it.
#'
#' @return An object of class \code{matrix_parameter}. The class is abstract;
#'   use one of the constructors, such as \code{\link{log_cholesky}}.
#'
#' @seealso \code{\link{parameter}}, \code{\link{log_cholesky}}
#'
#' @examples
#' S7::S7_inherits(log_cholesky(3), matrix_parameter)
#'
#' @export
matrix_parameter <- S7::new_class(
  "matrix_parameter",
  parent = parameter,
  properties = list(
    dimension = S7::class_integer,
    rank = S7::class_integer,
    null_basis = S7::class_numeric,
    role = S7::class_character
  ),
  validator = function(self) {
    errors <- character()
    if (length(self@dimension) != 1L || self@dimension < 1L) {
      errors <- c(errors, "@dimension must be a single positive integer")
    }
    if (length(self@rank) != 1L || self@rank < 0L || self@rank > self@dimension) {
      errors <- c(errors, "@rank must be a single integer in 0:dimension")
    }
    if (!is.matrix(self@null_basis)) {
      errors <- c(errors, "@null_basis must be a matrix")
    } else if (!identical(dim(self@null_basis), c(self@dimension, self@dimension - self@rank))) {
      errors <- c(errors, "@null_basis must be dimension by (dimension - rank)")
    }
    if (!length(self@role) == 1L ||
      !self@role %in% c("covariance", "precision", "either")) {
      errors <- c(errors, "@role must be 'covariance', 'precision' or 'either'")
    }
    if (length(errors)) errors else NULL
  }
)


#' Validate the Arguments Shared by Every Constructor
#'
#' @description
#' Checks the matrix side and the role, and returns the side as an integer.
#'
#' @param dimension The side of the matrix.
#' @param role The role label.
#'
#' @return \code{dimension}, as a single integer.
#'
#' @keywords internal
check_param_args <- function(dimension, role) {
  if (!is.numeric(dimension) || length(dimension) != 1L || !is.finite(dimension) ||
    dimension < 1 || dimension != round(dimension)) {
    stop("'dimension' must be a single positive integer.", call. = FALSE)
  }
  if (!is.character(role) || length(role) != 1L ||
    !role %in% c("covariance", "precision", "either")) {
    stop(
      "'role' must be one of \"covariance\", \"precision\" or \"either\".",
      call. = FALSE
    )
  }
  as.integer(dimension)
}


#' Validate a Free Vector Against a Parameter
#'
#' @description
#' Checks that \code{eta} is a finite numeric vector of the length the
#' parameter declares, and returns it unnamed.
#'
#' @details
#' Called in the body of every generic before dispatch, so that a parameter
#' written outside the package inherits the check without doing anything. The
#' names are stripped for the reason \code{align_theta()} strips them in
#' \pkg{distributions7}: a value that has been through a link comes back
#' carrying its own name, which is meaningless on a number and would leak into
#' the dimnames of the result.
#'
#' @param s A \code{\link{parameter}} object.
#' @param eta The free vector supplied by the caller.
#'
#' @return \code{eta}, as an unnamed numeric vector.
#'
#' @keywords internal
check_eta <- function(s, eta) {
  if (!is.numeric(eta)) {
    stop("'eta' must be numeric.", call. = FALSE)
  }
  if (length(eta) != s@n_free) {
    stop(sprintf(
      "'eta' has %d value(s) but '%s' has %d free value(s)%s.",
      length(eta), s@param_name, s@n_free,
      if (s@n_free) paste0(" (", paste(s@free_names, collapse = ", "), ")") else ""
    ), call. = FALSE)
  }
  if (anyNA(eta) || any(!is.finite(eta))) {
    stop("'eta' must be finite: the free scale has no boundary to reach.",
      call. = FALSE
    )
  }
  unname(eta)
}


#' Validate a Matrix Handed Back to a Parameter
#'
#' @description
#' Checks that \code{m} is a square symmetric numeric matrix of the parameter's
#' dimension.
#'
#' @param s A \code{\link{parameter}} object.
#' @param m The matrix supplied by the caller.
#' @param tol The relative tolerance for the symmetry check.
#'
#' @return \code{m}, symmetrised.
#'
#' @keywords internal
check_matrix <- function(s, m, tol = 1e-8) {
  if (!is.matrix(m) || !is.numeric(m)) {
    stop("'m' must be a numeric matrix.", call. = FALSE)
  }
  if (!identical(dim(m), c(s@dimension, s@dimension))) {
    stop(sprintf(
      "'m' must be %d by %d, matching the parameter.", s@dimension, s@dimension
    ), call. = FALSE)
  }
  if (anyNA(m)) stop("'m' must not contain missing values.", call. = FALSE)
  asym <- max(abs(m - t(m)))
  if (asym > tol * max(1, max(abs(m)))) {
    stop("'m' must be symmetric.", call. = FALSE)
  }
  (m + t(m)) / 2
}


#' Name the Rows and Columns of a Parameter's Matrix
#'
#' @description
#' Applies the dimension labels every matrix a parameter produces carries.
#'
#' @param m A numeric matrix.
#' @param s A \code{\link{parameter}} object.
#'
#' @return \code{m}, with dimnames.
#'
#' @keywords internal
name_dims <- function(m, s) {
  nm <- paste0("v", seq_len(s@dimension))
  dimnames(m) <- list(nm, nm)
  m
}


#' An Orthonormal Basis of a Null Space, and the Rank That Goes With It
#'
#' @description
#' Returns the rank and an orthonormal basis of the common null space of a set
#' of symmetric positive semidefinite matrices.
#'
#' @details
#' The components are normalised individually before being stacked, and the
#' rank is read off the stack rather than off any assembled combination. The
#' reason is that a rank read off a sum is not scale invariant: on a
#' tensor-product penalty of true rank 28 out of 32, counting eigenvalues above
#' a relative tolerance of \eqn{\lambda_1 P_1 + \lambda_2 P_2} gives 28 when
#' the two are equal, 26 at a ratio of \eqn{10^8} and 16 at \eqn{10^{12}},
#' because the smaller contributions sink below the tolerance and are read as
#' zeros. Smoothing parameters that far apart are an ordinary fitted model.
#'
#' The null space of a sum of positive semidefinite matrices is the
#' intersection of their null spaces, which is the null space of the stack, so
#' the stacked route is not an approximation of the assembled one but the exact
#' statement of the same quantity.
#'
#' @param mats A list of symmetric matrices, or a single matrix.
#' @param tol The relative tolerance below which a singular value counts as
#'   zero.
#'
#' @return A list with \code{rank} and \code{null_basis}.
#'
#' @seealso \code{\link{parameter}}
#'
#' @examples
#' # a second-difference penalty: its null space is the constants and the lines
#' d <- diff(diag(6), differences = 2)
#' param_null_basis(crossprod(d))$rank
#'
#' @export
param_null_basis <- function(mats, tol = 1e-10) {
  if (is.matrix(mats)) mats <- list(mats)
  if (!length(mats)) stop("'mats' must not be empty.", call. = FALSE)

  p <- ncol(mats[[1L]])
  stacked <- do.call(rbind, lapply(mats, function(m) {
    sc <- max(abs(m))
    if (sc == 0) matrix(0, nrow(m), ncol(m)) else m / sc
  }))

  sv <- svd(stacked, nu = 0L, nv = p)
  keep <- sv$d > tol * max(sv$d, 0)
  r <- sum(keep)
  list(
    rank = as.integer(r),
    null_basis = sv$v[, !keep, drop = FALSE]
  )
}


#' No Null Space
#'
#' @description
#' The empty basis a full-rank family declares.
#'
#' @param dimension The side of the matrix.
#'
#' @return A \code{dimension} by 0 numeric matrix.
#'
#' @keywords internal
empty_null_basis <- function(dimension) {
  matrix(numeric(0), nrow = dimension, ncol = 0L)
}
