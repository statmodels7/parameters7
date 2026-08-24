#' The Abstract Class of a Constrained Parameter
#'
#' @description
#' A parameter is a map from an unconstrained vector
#' \eqn{\eta \in \mathbb{R}^{d}} onto a quantity that lives in a constrained
#' set: a symmetric positive definite matrix, a probability vector, a
#' stochastic matrix. `parameter` is the abstract S7 class that map belongs
#' to. It records how many free values the map consumes, what each of them is
#' called, and whatever the family needs in order to evaluate itself. Because
#' the domain is the whole of \eqn{\mathbb{R}^{d}}, an optimizer may move
#' \eqn{\eta} anywhere and every point it visits still maps to a valid value,
#' so no constraint has to be policed during a fit.
#'
#' The class is abstract and is not meant to be instantiated. Call a
#' constructor: [log_cholesky()] for an unstructured covariance, [ar1()] or
#' [compound_symmetry()] for a structured one, [simplex()] for a probability
#' vector, [transition_matrix()] for a Markov chain.
#'
#' @details
#' # A parameter owns its dimension
#'
#' `log_cholesky(3)` and `log_cholesky(4)` are different objects, holding
#' \eqn{d = 6} and \eqn{d = 10} free values. Fixing the size at construction
#' means `n_free` and `free_names` can be answered before any data exist, and
#' both are needed then: a model builds its parameter table from the names
#' while it is still assembling the design.
#'
#' # What a subclass must supply
#'
#' [param_value()] is the one compulsory method. Given it, a subclass inherits
#' [param_d1()], [param_d2()], [param_d3()] and [param_d4()], each a product
#' stencil applied to the map itself. A closed form registered later replaces
#' the inherited numerical method through dispatch, and no caller changes.
#'
#' Four generics are not inherited here:
#'
#' - [param_free()], the inverse map, has a base method that throws and names
#'   the family. An inverse found by optimization would hand back a plausible
#'   \eqn{\eta} for a value that is outside the set the family parametrizes, so
#'   the inverse is either written out exactly or refused.
#' - [param_logdet()], [param_solve()] and [param_factor()] are registered on
#'   [matrix_parameter()], not on this class. A probability vector has no
#'   log-determinant, and asking for one fails at dispatch instead of returning
#'   a number.
#'
#' # Free names label the coordinate, not the quantity
#'
#' The families here follow one convention. Where a link carries a constrained
#' quantity onto the free scale, the label records that link: a variance
#' appears as `"log_scale"` and a correlation as `"z_rho"`. Where the
#' coordinate is already unrestricted the label is the plain name of the
#' quantity, as the below-diagonal entries `"L2.1"` of a Cholesky factor are.
#'
#' The distinction matters outside the family. A consumer flattens the free
#' vector into scalar parameters carrying identity links, so a label promising
#' a bounded quantity would report a number that is not on that scale: a free
#' value of 1.4 called `rho` reads as a correlation of 1.4, while called
#' `z_rho` it reads as the correlation \eqn{\tanh(1.4) = 0.885}.
#'
#' @section Notation:
#' \eqn{\eta} is the **free vector**, the point on the unconstrained scale that
#' an optimizer moves, and \eqn{d} its length, the object's `n_free`. Note that
#' \eqn{\eta} means a linear predictor elsewhere in the toolkit; here it is
#' always the free vector of a parametrization. \eqn{p} is the side of the
#' matrix a matrix family produces.
#'
#' @param param_name A single character string naming the family, used in the
#'   error messages the validators raise and in the object's `print` output.
#'   `"log_cholesky"`, `"ar1"` and so on.
#' @param n_free The length \eqn{d} of the free vector: a single non-negative
#'   integer. Zero is legal and describes a family with nothing to estimate.
#'   The validator rejects a vector, a negative value, or a length that
#'   disagrees with `free_names`.
#' @param free_names A character vector of length `n_free`, one label per free
#'   value, in the order the free vector holds them. Fixed at construction and
#'   part of the interface: consumers build their parameter tables from these
#'   labels, so the ordering is not free to change. The validator rejects a
#'   duplicated label and a length other than `n_free`.
#' @param param_params A list of whatever the family needs in order to evaluate
#'   itself, read only by that family's own methods. `log_cholesky()` stores
#'   the row and column index of each free value here; `sum_struct()` stores
#'   its component matrices. Nothing outside the family looks inside it.
#'
#' @return An object of class `parameter`, with properties
#'   \describe{
#'     \item{`param_name`}{character, the family name.}
#'     \item{`n_free`}{integer, the length \eqn{d} of the free vector.}
#'     \item{`free_names`}{character of length `n_free`.}
#'     \item{`param_params`}{list, the family's own data.}
#'   }
#'   The class is abstract, so a useful object comes from a constructor such as
#'   [log_cholesky()] and carries that constructor's subclass. A matrix family
#'   returns a [matrix_parameter()], which adds `dimension`, `rank`,
#'   `null_basis` and `role`.
#'
#' @seealso [matrix_parameter()] for the symmetric matrix branch.
#'   [param_value()] and [param_free()] for the map and its inverse,
#'   [param_d1()] for the derivative arrays, [param_readable()] for the
#'   quantities a family reports, and [check_parameter()] to verify a
#'   parametrization written from scratch.
#'
#'   The constructors: [log_cholesky()], [matrix_log()], [diagonal_matrix()],
#'   [scalar_matrix()], [scaled_matrix()], [correlation_matrix()],
#'   [compound_symmetry()], [ar1()], [autoregressive()], [simplex()],
#'   [transition_matrix()], and the compositions [kron_identity()],
#'   [block_diag()], [dr_prod()] and [sum_struct()].
#'
#' @examples
#' # Every parametrization in the package is a `parameter`.
#' s <- log_cholesky(3)
#' S7::S7_inherits(s, parameter)
#'
#' # The four properties the class itself holds.
#' s@param_name
#' s@n_free
#' s@free_names
#'
#' # The dimension is fixed at construction, so the length of the free vector
#' # is known before any data arrive: p(p+1)/2 for an unstructured covariance.
#' d <- vapply(2:5, function(p) log_cholesky(p)@n_free, integer(1))
#' rbind(p = 2:5, n_free = d, `p(p+1)/2` = (2:5) * (3:6) / 2)
#'
#' # A family whose value is not a symmetric matrix inherits `parameter`
#' # directly, so it has no log-determinant to be asked for.
#' q <- simplex(3)
#' c(parameter = S7::S7_inherits(q, parameter),
#'   matrix_parameter = S7::S7_inherits(q, matrix_parameter))
#' sum(param_value(q, c(0.7, -0.3)))          # a probability vector: sums to 1
#'
#' # The free names record the link, so a free value of 1.4 in a `z_rho`
#' # coordinate is a correlation of tanh(1.4), not a correlation of 1.4.
#' ar1(4)@free_names
#' param_value(ar1(4), c(0, 1.4))[1, 2]
#' tanh(1.4)
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


#' The Abstract Class of a Symmetric Matrix Parameter
#'
#' @description
#' Extends [parameter()] to the branch whose value is a symmetric positive
#' semidefinite \eqn{p \times p} matrix. Beyond the map itself, such a family
#' can answer the four things a Gaussian likelihood asks of a covariance or a
#' precision: the log-determinant, a solve against a right-hand side, a factor,
#' and the rank and null space when the matrix is singular. Those four generics
#' are registered here, so a subclass supplying only [param_value()] inherits
#' all of them.
#'
#' The class is abstract. Every matrix family in the package returns a subclass
#' of it: [log_cholesky()], [matrix_log()], [diagonal_matrix()],
#' [correlation_matrix()], [compound_symmetry()], [ar1()], [autoregressive()],
#' [scaled_matrix()] and the four compositions.
#'
#' @details
#' # Rank and null space are fixed at construction
#'
#' They are properties of the family, the same at every point \eqn{\eta}.
#' Scaling a matrix by a positive number leaves its null space alone, and the
#' null space of a sum of positive semidefinite matrices is the intersection of
#' theirs, so no free value moves it. Recording them once at construction is
#' exact, and it is also the only stable way to obtain them: counting the small
#' eigenvalues of an assembled matrix is not scale invariant, and reads a
#' component whose weight is small as a null direction. [param_null_basis()]
#' computes them from the components and carries the measurement.
#'
#' Most families here are full rank, so `rank` is \eqn{p} and `null_basis` has
#' zero columns. [scaled_matrix()] admits a deficient fixed matrix, and
#' [sum_struct()] a deficient sum, and both then report the deficiency.
#'
#' # What a non-matrix family does instead
#'
#' [simplex()] and [transition_matrix()] inherit [parameter()] directly. They
#' hold no `dimension`, `rank` or `null_basis`, and [param_logdet()] has no
#' method for them, so asking for the log-determinant of a probability vector
#' fails at dispatch. The absence is structural.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, the point on the unconstrained scale, and
#' \eqn{d} its length. \eqn{p} is the side of the matrix. \eqn{M} is the matrix
#' the map produces, \eqn{\Sigma} when it is read as a covariance and
#' \eqn{\Omega} when it is read as a precision.
#'
#' @param param_name A single character string naming the family.
#' @param n_free The length \eqn{d} of the free vector: a single non-negative
#'   integer, agreeing with `length(free_names)`.
#' @param free_names A character vector of length `n_free`, one label per free
#'   value, in the order the free vector holds them. Must be unique.
#' @param param_params A list of whatever the family needs in order to evaluate
#'   itself, read only by that family's own methods.
#' @param dimension The side \eqn{p} of the matrix: a single integer, no
#'   smaller than 1. The validator rejects a vector and a value below 1.
#' @param rank The rank of the matrix the family produces, a single integer in
#'   `0:dimension`. It is a property of the family, so a family whose value is
#'   positive definite at every \eqn{\eta} declares \eqn{p} here.
#' @param null_basis A `dimension` by `dimension - rank` numeric matrix whose
#'   columns are an orthonormal basis of the common null space. Use
#'   [param_null_basis()] to obtain one, or `matrix(numeric(0), dimension, 0)`
#'   for a full-rank family. The validator rejects any other shape, and
#'   reports both the rank and the shape when the two disagree.
#' @param role A single string, one of `"covariance"`, `"precision"` or
#'   `"either"`, recording which side of a model the matrix parametrizes. **No
#'   numeric result depends on it.** It is carried because the family name does
#'   not record it: the same [log_cholesky()] serves either side, and a
#'   consumer that prefixes a free name with the matrix it describes needs to
#'   know which. [block_diag()] reads it to give a composite the common role of
#'   its blocks, or `"either"` when they disagree, and [kron_identity()] copies
#'   it.
#'
#' @return An object of class `matrix_parameter`, which is a [parameter()] with
#'   four further properties
#'   \describe{
#'     \item{`dimension`}{integer, the side \eqn{p}.}
#'     \item{`rank`}{integer in `0:dimension`.}
#'     \item{`null_basis`}{a `dimension` by `dimension - rank` matrix with
#'       orthonormal columns.}
#'     \item{`role`}{character, as supplied.}
#'   }
#'   plus the four it inherits, `param_name`, `n_free`, `free_names` and
#'   `param_params`. The class is abstract, so a useful object comes from a
#'   constructor such as [log_cholesky()].
#'
#' @seealso [parameter()] for the base class and the naming convention for free
#'   values. [param_logdet()], [param_solve()], [param_factor()] and
#'   [param_null_basis()] for the four quantities this branch adds.
#'
#' @examples
#' # A matrix family is both a `parameter` and a `matrix_parameter`.
#' s <- log_cholesky(3)
#' c(parameter = S7::S7_inherits(s, parameter),
#'   matrix_parameter = S7::S7_inherits(s, matrix_parameter))
#'
#' # The four properties this branch adds.
#' c(dimension = s@dimension, rank = s@rank, null_columns = ncol(s@null_basis))
#' s@role
#'
#' # Belonging to this branch is what gives the family a log-determinant, and
#' # it agrees with the eigenvalues of the matrix itself.
#' eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
#' all.equal(param_logdet(s, eta),
#'           sum(log(eigen(param_value(s, eta), only.values = TRUE)$values)))
#'
#' # A rank-deficient family declares the deficiency at construction, and the
#' # null basis really is annihilated by the matrix at any free value.
#' r <- scaled_matrix(diag(c(1, 1, 0, 0)))
#' c(rank = r@rank, null_columns = ncol(r@null_basis))
#' max(abs(param_value(r, 0.6) %*% r@null_basis))
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


#' Validate the Arguments Shared by Every Matrix Constructor
#'
#' @description
#' Checks the two arguments every matrix family's constructor takes, and returns
#' the matrix side coerced to integer so the caller can store it in the class's
#' integer property. Called by [log_cholesky()], [matrix_log()],
#' [diagonal_matrix()], [correlation_matrix()], [compound_symmetry()], [ar1()],
#' [autoregressive()] and [dr_prod()].
#'
#' @param dimension The side of the matrix. Must be a single finite number, at
#'   least 1, equal to its own `round()`. `0`, `2.5`, `c(1, 2)`, `"3"`, `Inf`
#'   and `NA` all throw `'dimension' must be a single positive integer.`
#' @param role The role label. Must be a single string, one of
#'   `"covariance"`, `"precision"` or `"either"`; anything else throws. The
#'   shipped constructors run `match.arg()` first, so their callers see
#'   `match.arg()`'s message and never this one. The check is here for a
#'   constructor written outside the package.
#'
#' @return `dimension`, as a single integer.
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
#' Checks that `eta` is a finite numeric vector of the length the parameter
#' declares, and returns it with its names stripped. Called in the body of
#' every generic before dispatch, so a parameter written outside the package
#' inherits the check without doing anything.
#'
#' @details
#' Three conditions are checked, each with its own message. A non-numeric `eta`
#' throws `'eta' must be numeric.`; a wrong length throws a message naming both
#' counts and listing the family's `free_names`, so a caller who has mismatched
#' two parametrizations can see which is which; and any `NA`, `NaN` or infinite
#' entry throws `'eta' must be finite: the free scale has no boundary to
#' reach.` The last message states the reason a non-finite free value is a
#' caller error: the unconstrained scale has no edge, so an infinity there is a
#' runaway rather than a limit reached.
#'
#' Names are stripped for the reason `align_theta()` strips them in
#' \pkg{distributions7}. A value that has been through a link comes back
#' carrying its own name, which means nothing on a number and would otherwise
#' appear in the dimnames of the matrix built from it.
#'
#' @param s A [parameter()] object, whose `n_free` and `free_names` are read.
#' @param eta The free vector supplied by the caller.
#'
#' @return `eta`, as an unnamed numeric vector of length `s@n_free`.
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
#' Checks that `m` is a square symmetric numeric matrix of the parameter's
#' dimension, and returns it symmetrized. Called by every [param_free()] method
#' on the matrix the caller is asking to invert.
#'
#' @details
#' The symmetry check is relative: `m` is rejected when
#' \eqn{\max_{ij} |m_{ij} - m_{ji}|} exceeds `tol * max(1, max(abs(m)))`. At the
#' default a discrepancy of \eqn{10^{-4}} throws and one of \eqn{10^{-12}}
#' passes. The return value is `(m + t(m)) / 2`, so an asymmetry small enough to
#' be rounding is averaged away instead of being carried into the inverse map.
#'
#' @param s A [parameter()] object, whose `dimension` is read.
#' @param m The matrix supplied by the caller. A non-matrix or non-numeric `m`
#'   throws, as does a wrong shape (with a message naming the dimension
#'   required), an `NA` anywhere, and an asymmetry above `tol`.
#' @param tol The relative tolerance for the symmetry check, a single positive
#'   number. Defaults to `1e-8`, loose enough to accept a matrix assembled from
#'   a Cholesky factor or a Kronecker product and tight enough to reject one
#'   that is genuinely not symmetric.
#'
#' @return `m`, symmetrized as `(m + t(m)) / 2`.
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
#' Applies the dimension labels every matrix a parameter produces carries:
#' `"v1"`, `"v2"`, ..., `"vp"` on both margins, `p` being the parameter's
#' `dimension`. One convention across the families, so a consumer can read a
#' printed covariance without knowing which parametrization built it.
#'
#' @param m A numeric matrix, `p` by `p`. Not checked; the callers are the
#'   package's own [param_value()] methods.
#' @param s A [parameter()] object, whose `dimension` supplies `p`.
#'
#' @return `m`, with `dimnames` set to `list(v1..vp, v1..vp)`. Any dimnames it
#'   arrived with are replaced.
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
#' Returns the rank and an orthonormal basis of the common null space of one or
#' more symmetric positive semidefinite matrices. Each matrix is scaled to unit
#' maximum entry, the scaled matrices are stacked into one tall matrix, and the
#' rank and the null space are read off its singular value decomposition. Every
#' [matrix_parameter()] calls it once at construction to fill its `rank` and
#' `null_basis` properties.
#'
#' @details
#' # Reading the rank off the components
#'
#' The null space of a sum of positive semidefinite matrices is the
#' intersection of their null spaces, and the intersection is the null space of
#' the stack, so reading the answer off the stack is exact. Reading it off an
#' assembled combination \eqn{\lambda_1 P_1 + \lambda_2 P_2} is not, because a
#' count of small eigenvalues is not scale invariant: a component whose weight
#' is small contributes eigenvalues below the tolerance, and they are then
#' counted as null directions that are not there.
#'
#' Measured on the tensor-product penalty of two second-difference penalties
#' over 4 and 8 coefficients, whose true rank is 28 out of 32. Counting the
#' eigenvalues of the assembled sum above a relative tolerance of
#' \eqn{10^{-10}} answers 28 while the two weights are within \eqn{10^{8}} of
#' each other, and 24 once the ratio reaches \eqn{10^{10}}: four directions the
#' penalty does penalize are read as null. Smoothing parameters ten orders of
#' magnitude apart are an ordinary fitted model, not a pathology. The stacked
#' route answers 28 at every ratio, and the basis it returns is annihilated by
#' the assembled matrix to \eqn{3 \times 10^{-16}} relative even at the worst
#' one.
#'
#' # The scaling
#'
#' Each component is divided by its largest absolute entry before stacking, so
#' that a component contributed in different units still gets an equal vote. A
#' component that is identically zero is passed through as zero and contributes
#' nothing, which is correct: its null space is the whole space.
#'
#' # Cost
#'
#' One singular value decomposition of an \eqn{(kp) \times p} matrix for
#' \eqn{k} components of side \eqn{p}, so \eqn{O(k p^3)}. It runs once per
#' constructed parameter and never inside a fit.
#'
#' @section Notation:
#' \eqn{p} is the side of the matrices and \eqn{k} their number. The **null
#' basis** is a matrix whose columns span \eqn{\{v : M v = 0 \text{ for every }
#' M\}}, and the **rank** is \eqn{p} minus the dimension of that space.
#'
#' @param mats A list of symmetric numeric matrices, all of the same side, or a
#'   single matrix, which is wrapped in a list. An empty list throws `'mats'
#'   must not be empty.`; the matrices themselves are not checked for symmetry
#'   or definiteness, the callers being the package's own constructors.
#' @param tol The relative tolerance below which a singular value counts as
#'   zero: a singular value \eqn{s_j} is null when
#'   \eqn{s_j \le \mathrm{tol} \cdot \max_i s_i}. Defaults to `1e-10`. The
#'   default sits well below the singular values a genuine rank carries and
#'   well above the \eqn{10^{-16}} the zero directions of a stacked, normalized
#'   set of matrices reach, so the gap between the two is wide and the exact
#'   value is not delicate.
#'
#' @return A list with two components
#'   \describe{
#'     \item{`rank`}{an integer, the number of singular values above the
#'       tolerance.}
#'     \item{`null_basis`}{a `p` by `p - rank` numeric matrix with orthonormal
#'       columns, spanning the common null space. It has zero columns when the
#'       stack has full rank.}
#'   }
#'
#' @seealso [matrix_parameter()], whose `rank` and `null_basis` properties this
#'   fills, and the two families that can be rank deficient and so are the ones
#'   whose answer is not trivial: [scaled_matrix()] and [sum_struct()].
#'
#' @examples
#' # A second-difference penalty on six coefficients has rank 4, its null space
#' # being the constants and the straight lines.
#' d <- diff(diag(6), differences = 2)
#' nb <- param_null_basis(crossprod(d))
#' nb$rank
#' dim(nb$null_basis)
#'
#' # The claim, checked: a constant and a line are in the span and a quadratic
#' # is not. Projecting onto the basis leaves the first two untouched.
#' proj <- function(v) nb$null_basis %*% (t(nb$null_basis) %*% v)
#' vapply(list(rep(1, 6), 1:6, (1:6)^2),
#'        function(v) sqrt(sum((v - proj(v))^2)) / sqrt(sum(v^2)), numeric(1))
#'
#' # Why the components are taken separately. Two marginal penalties on a
#' # 4 x 8 tensor product: the pair has rank 28 out of 32.
#' P <- function(m) crossprod(diff(diag(m), differences = 2))
#' P1 <- kronecker(P(4), diag(8))
#' P2 <- kronecker(diag(4), P(8))
#' param_null_basis(list(P1, P2))$rank
#'
#' # Counting small eigenvalues of the assembled sum agrees while the two
#' # weights are comparable and loses four directions when they are not.
#' count <- function(M, tol = 1e-10) {
#'   e <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
#'   sum(e > tol * max(e))
#' }
#' c(equal = count(P1 + P2), ratio_1e10 = count(P1 + 1e10 * P2))
#'
#' # The stacked answer stays correct there, and its basis is annihilated.
#' nb2 <- param_null_basis(list(P1, P2))
#' M <- P1 + 1e10 * P2
#' max(abs(M %*% nb2$null_basis)) / max(abs(M))
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
#' The empty basis a full-rank family declares. Every family whose value is
#' positive definite at every free vector passes this to its
#' [matrix_parameter()] constructor. The validator's shape rule, `dimension` by
#' `dimension - rank`, asks for exactly this when `rank` is `dimension`.
#'
#' @param dimension The side of the matrix, a single integer.
#'
#' @return A `dimension` by 0 numeric matrix. `ncol()` is 0, so a consumer that
#'   loops over the null directions does nothing without a special case.
#'
#' @keywords internal
empty_null_basis <- function(dimension) {
  matrix(numeric(0), nrow = dimension, ncol = 0L)
}
