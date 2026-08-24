#' @include simplex.R
NULL


#' Transition Matrix Parameter
#'
#' @description
#' The S7 class of row-stochastic matrices, each row on the open simplex in the
#' additive log-ratio parametrization. Its value is a matrix but **not a
#' symmetric** one, so it inherits [parameter()] directly, never
#' [matrix_parameter()]: there is no `rank`, no `null_basis` and no `role`, and a
#' transition matrix has no log-determinant, solve or factor to be asked for.
#'
#' [transition_matrix()] builds one. The rows are independent in the
#' parametrization, so every derivative array is block diagonal by row and the
#' family reuses [simplex()]'s kernels row by row.
#'
#' @inheritParams parameter
#'
#' @return An object of class `TransitionMatrixParam`, a subclass of
#'   [parameter()] adding no properties of its own: `param_name` is
#'   `"transition_matrix"`, `n_free` is \eqn{K(K-1)}, `free_names` is
#'   `alr{i}.{j}` row by row, and `param_params` holds `n_state`.
#'
#' @seealso [transition_matrix()], the constructor, [simplex()], which is one
#'   row of this, and [parameter()] for the properties this inherits and the
#'   generics it does not get.
#'
#' @examples
#' # A matrix, but not a symmetric one, so not a matrix_parameter.
#' s <- transition_matrix(3)
#' c(parameter = S7::S7_inherits(s, parameter),
#'   matrix_parameter = S7::S7_inherits(s, matrix_parameter))
#'
#' # K(K-1) free values: one simplex per row.
#' rbind(K = 2:5, n_free = vapply(2:5, function(k)
#'   transition_matrix(k)@n_free, integer(1)))
#'
#' @export
TransitionMatrixParam <- S7::new_class("TransitionMatrixParam", parent = parameter)


#' Construct a Transition Matrix Parameter
#'
#' @description
#' Returns an object holding a \eqn{K \times K} row-stochastic matrix, the
#' transition matrix of a Markov chain on \eqn{K} states, with each row an
#' independent [simplex()] in the additive log-ratio chart. \eqn{K(K-1)} free
#' values in all, and every free vector gives positive entries with rows summing
#' to exactly 1, so a chain can be estimated without a constraint.
#'
#' @details
#' # Rows, not columns
#'
#' The rows live on the simplex, so row \eqn{i} holds the probabilities of moving
#' **from** state \eqn{i}. A transition matrix in this convention acts on row
#' vectors of probabilities, \eqn{\pi_{t+1} = \pi_t P}; the column convention is
#' its transpose.
#'
#' # The rows are independent
#'
#' Each row is parametrized by its own \eqn{K-1} free values and no others, so
#' every derivative array is **block diagonal by row**: a component pairing free
#' values of two different rows is exactly the zero matrix, and a component
#' inside row \eqn{i} is zero everywhere outside row \eqn{i}. The implementation
#' evaluates the simplex kernels row by row and never stores those zeros.
#'
#' Measured at \eqn{K = 3}: the second-derivative component
#' `alr1.1:alr2.1`, spanning two rows, is 0 exactly, while `alr1.1:alr1.1` is
#' not.
#'
#' # Reading the free vector
#'
#' It runs row by row, and the names `alr{i}.{j}` say which row and which chart
#' coordinate, the row.column convention [log_cholesky()] uses. The `alr` records
#' the chart, so a free value of 0.5 is not a probability of 0.5: it is a log
#' ratio against the row's last state.
#'
#' # What it has and has not got
#'
#' All four derivative orders are closed form, inherited from [simplex()]'s
#' cumulant recursion. There is no log-determinant, no solve and no factor: the
#' value is not symmetric, and the object is a [parameter()] and never a
#' [matrix_parameter()], so those generics have no method and fail at dispatch.
#' [check_parameter()] runs its seven-check battery instead of the nine-check one,
#' and the identity it tests here is that every derivative component sums to zero
#' along each row.
#'
#' @section Notation:
#' \eqn{K} is the number of states, \eqn{\eta \in \mathbb{R}^{K(K-1)}} the free
#' vector and \eqn{P} the transition matrix, with \eqn{P_{ij}} the probability of
#' moving from state \eqn{i} to state \eqn{j}. The reference category of each row
#' is its last state.
#'
#' @param n_state The number of states \eqn{K}, **at least 2**. A single integer;
#'   `1`, a fraction, `NA` and a vector all throw `'n_state' must be a single
#'   integer of at least 2.` A one-state chain has nothing to estimate.
#'
#' @return An object of class [TransitionMatrixParam()], with `n_free` equal to
#'   \eqn{K(K-1)}, `free_names` `alr1.1`, `alr1.2`, ..., `alr{K}.{K-1}` row by
#'   row, `param_name` `"transition_matrix"`, and `param_params` holding
#'   `n_state`. No `dimension`, `rank`, `null_basis` or `role`.
#'
#' @seealso [simplex()], which is one row of this and carries the derivative
#'   recursion, [param_value()] and [param_free()] for the map and its inverse,
#'   and [check_parameter()] for the battery a non-matrix family gets.
#'
#' @examples
#' # Three states: six free values, two per row.
#' s <- transition_matrix(3)
#' c(n_free = s@n_free)
#' s@free_names
#'
#' set.seed(4)
#' eta <- rnorm(s@n_free)
#' m <- param_value(s, eta)
#' round(m, 4)
#'
#' # Every row is a probability distribution, exactly.
#' rowSums(m)
#'
#' # The round trip closes exactly.
#' max(abs(param_free(s, m) - eta))
#'
#' # The rows are independent, so a derivative in a row-1 free value is zero
#' # everywhere outside row 1.
#' round(param_d1(s, eta)[["alr1.1"]], 4)
#'
#' # And a second derivative spanning two rows is exactly the zero matrix.
#' d2 <- param_d2(s, eta)
#' c(across_rows = max(abs(d2[["alr1.1:alr2.1"]])),
#'   within_row = max(abs(d2[["alr1.1:alr1.1"]])))
#'
#' # There is no log-determinant to ask for: the value is not symmetric.
#' try(param_logdet(s, eta))
#'
#' @export
transition_matrix <- function(n_state) {
  if (!is.numeric(n_state) || length(n_state) != 1L || !is.finite(n_state) ||
    n_state < 2 || n_state != round(n_state)) {
    stop("'n_state' must be a single integer of at least 2.", call. = FALSE)
  }
  k <- as.integer(n_state)
  nm <- as.vector(t(outer(seq_len(k), seq_len(k - 1L), function(i, j) {
    paste0("alr", i, ".", j)
  })))
  TransitionMatrixParam(
    param_name = "transition_matrix",
    n_free = k * (k - 1L),
    free_names = nm,
    param_params = list(n_state = k)
  )
}


#' Row and Chart Coordinate of Each Free Value
#'
#' @description
#' Returns, for each of the \eqn{K(K-1)} free values, the row it belongs to and
#' its position inside that row's simplex chart, in the order the free vector
#' uses. Every method of the family reads it to split the free vector into rows.
#'
#' @details
#' The free vector runs row by row, so at \eqn{K = 3} the rows are
#' `1 1 2 2 3 3` and the coordinates `1 2 1 2 1 2`, which labels the free values
#' `alr1.1 alr1.2 alr2.1 alr2.2 alr3.1 alr3.2`.
#'
#' @param s A [TransitionMatrixParam()] object, whose `param_params$n_state` is
#'   read.
#'
#' @return A list with two integer vectors of length `s@n_free`: `row`, with
#'   values in `1:K`, and `coord`, with values in `1:(K-1)`.
#'
#' @seealso [transition_matrix()], which turns this into `free_names`, and
#'   [tm_derivative()], which uses it to place a row's component.
#'
#' @keywords internal
tm_positions <- function(s) {
  k <- s@param_params$n_state
  list(
    row = rep(seq_len(k), each = k - 1L),
    coord = rep(seq_len(k - 1L), times = k)
  )
}


#' @title Value of a Transition Matrix Parameter
#' @name param_value.TransitionMatrixParam
#' @description
#' Returns the transition matrix: each row is the softmax of that row's own
#' \eqn{K-1} free values, through [simplex_point()], so each row is positive and
#' sums to exactly 1. Nothing is tested and nothing is renormalized; the row sums
#' are a property of the map. Rows are the distributions, so \eqn{P_{ij}} is the
#' probability of moving from state \eqn{i} to state \eqn{j}.
#' @param s A [TransitionMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A \eqn{K \times K} numeric matrix with positive entries and rows
#'   summing to 1, with dimnames `s1`, `s2`, ... on both margins.
#' @seealso [param_free.TransitionMatrixParam()] for the inverse, and
#'   [simplex_point()] for the per-row arithmetic.
#' @keywords internal
S7::method(param_value, TransitionMatrixParam) <- function(s, eta, ...) {
  k <- s@param_params$n_state
  pos <- tm_positions(s)
  m <- matrix(0, k, k)
  for (i in seq_len(k)) {
    m[i, ] <- simplex_point(eta[pos$row == i])
  }
  nm <- paste0("s", seq_len(k))
  dimnames(m) <- list(nm, nm)
  m
}


#' @title Free Vector of a Transition Matrix Parameter
#' @name param_free.TransitionMatrixParam
#' @description
#' Returns the additive log-ratio of each row, \eqn{\log(P_{ij}/P_{iK})}, exact
#' and a true inverse of [param_value.TransitionMatrixParam()]: the round trip
#' closes to \eqn{2 \times 10^{-16}}.
#' @details
#' A row is rejected when it has a non-positive entry, being then outside the
#' **open** simplex where \eqn{\log 0} is not finite, or when it does not sum to
#' 1, in which case it is not a probability distribution and is **not**
#' renormalized. The message names the row.
#'
#' The first rejection is the one a fit meets: [simplex_point()] saturates a large
#' free value to an exact 0 in the reference state, so a matrix produced at the
#' boundary cannot be inverted back.
#' @param s A [TransitionMatrixParam()] object.
#' @param m A \eqn{K \times K} row-stochastic numeric matrix: strictly positive
#'   with every row summing to 1.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.TransitionMatrixParam()], the map this inverts, and
#'   [param_free.SimplexParam()], which does one row's worth.
#' @keywords internal
S7::method(param_free, TransitionMatrixParam) <- function(s, m, ...) {
  k <- s@param_params$n_state
  if (!is.matrix(m) || !identical(dim(m), c(k, k))) {
    stop(sprintf("'m' must be a %d by %d matrix.", k, k), call. = FALSE)
  }
  if (anyNA(m) || any(m <= 0)) {
    stop("every entry of 'm' must be strictly positive.", call. = FALSE)
  }
  if (max(abs(rowSums(m) - 1)) > 1e-8) {
    stop(paste0(
      "the rows of 'm' do not all sum to one, so it is not row stochastic.\n",
      "  It is rejected rather than renormalized."
    ), call. = FALSE)
  }
  out <- numeric(s@n_free)
  pos <- tm_positions(s)
  for (i in seq_len(k)) {
    out[pos$row == i] <- log(m[i, seq_len(k - 1L)] / m[i, k])
  }
  stats::setNames(out, s@free_names)
}


#' Derivative Components of a Transition Matrix Parameter
#'
#' @description
#' Assembles a whole derivative order from the row-wise simplex tensors. A
#' component whose free values span two rows is the zero matrix; one whose free
#' values all belong to row \eqn{i} is that row's [simplex()] component embedded
#' in row \eqn{i} of a matrix of zeros. The four derivative methods of the family
#' are one call each to this function.
#'
#' @details
#' The rows are parametrized independently, so [simplex_tensors()] is evaluated
#' once per row and [simplex_components()] slices it with a `wrap` that does the
#' embedding. Nothing of size \eqn{K^2 (K(K-1))^{\text{order}}} is built: the
#' cross-row components are never computed, only skipped.
#'
#' @param s A [TransitionMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` \eqn{K \times K}
#'   matrices keyed as `param_tuple_names(s, order)` and in that order, most of
#'   them exactly zero. Each non-zero one is supported on a single row.
#'
#' @seealso [simplex_tensors()] and [simplex_components()], which do the per-row
#'   work, and [tm_positions()] for the row map.
#'
#' @keywords internal
tm_derivative <- function(s, eta, order) {
  k <- s@param_params$n_state
  pos <- tm_positions(s)
  tens <- lapply(seq_len(k), function(i) {
    simplex_tensors(simplex_point(eta[pos$row == i]), order)
  })
  arrname <- paste0("d", order)

  idx <- param_tuple_indices(s, order)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, order)
  for (t_i in seq_along(idx)) {
    t <- idx[[t_i]]
    rows <- pos$row[t]
    m <- matrix(0, k, k)
    if (all(rows == rows[1L])) {
      i <- rows[1L]
      co <- pos$coord[t]
      arr <- tens[[i]][[arrname]]
      slice <- switch(order,
        arr[, co[1L]],
        arr[, co[1L], co[2L]],
        arr[, co[1L], co[2L], co[3L]],
        arr[, co[1L], co[2L], co[3L], co[4L]]
      )
      m[i, ] <- slice
    }
    nm <- paste0("s", seq_len(k))
    dimnames(m) <- list(nm, nm)
    out[[t_i]] <- m
  }
  out
}


#' @title First Derivatives of a Transition Matrix Parameter
#' @name param_d1.TransitionMatrixParam
#' @description
#' Closed form. The free value `alr{i}.{j}` belongs to row \eqn{i} alone, so
#' \eqn{\partial P} is zero in every other row and, inside row \eqn{i}, is the
#' [simplex()] first derivative
#' \eqn{\partial_b \pi_a = \pi_a(\delta_{ab} - \pi_b)}: the covariance structure
#' of a categorical indicator.
#'
#' Each component therefore sums to zero along its own row, the row sums of
#' \eqn{P} being the constant 1.
#' @param s A [TransitionMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` \eqn{K \times K} matrices named by
#'   `s@free_names`, each supported on one row.
#' @seealso [tm_derivative()], which assembles it, [param_d1.SimplexParam()] for
#'   the per-row formula, and [param_d2.TransitionMatrixParam()] for the order
#'   above.
#' @keywords internal
S7::method(param_d1, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 1L)
}

#' @title Second Derivatives of a Transition Matrix Parameter
#' @name param_d2.TransitionMatrixParam
#' @description
#' Closed form, and mostly zero. A component pairing free values from two
#' **different** rows is exactly the zero matrix, the rows being parametrized
#' independently; one pairing two free values of the same row is that row's
#' [simplex()] second derivative, embedded in that row.
#'
#' Measured at \eqn{K = 3}: `alr1.1:alr2.1` is 0 exactly and `alr1.1:alr1.1` is
#' not. Of the 21 components at \eqn{K = 3}, only 9 can be non-zero.
#' @param s A [TransitionMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` \eqn{K \times K} matrices keyed as
#'   `param_tuple_names(s)` and in that order, the cross-row ones exactly zero.
#' @seealso [tm_derivative()], which assembles it, and
#'   [param_d1.TransitionMatrixParam()] and [param_d3.TransitionMatrixParam()]
#'   for the neighboring orders.
#' @keywords internal
S7::method(param_d2, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 2L)
}

#' @title Third Derivatives of a Transition Matrix Parameter
#' @name param_d3.TransitionMatrixParam
#' @description
#' Closed form, from [simplex()]'s cumulant recursion at third order, embedded row
#' by row. A component whose three free values do not all belong to one row is
#' exactly the zero matrix, so the great majority of the list is zero and is
#' skipped instead of computed.
#' @param s A [TransitionMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` \eqn{K \times K} matrices keyed as
#'   `param_tuple_names(s, 3)` and in that order.
#' @seealso [tm_derivative()], which assembles it, and
#'   [param_d4.TransitionMatrixParam()] for the order above.
#' @keywords internal
S7::method(param_d3, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 3L)
}

#' @title Fourth Derivatives of a Transition Matrix Parameter
#' @name param_d4.TransitionMatrixParam
#' @description
#' Closed form, from [simplex()]'s cumulant recursion at fourth order, embedded
#' row by row, and the top of the contract. A component whose four free values do
#' not all belong to one row is exactly the zero matrix.
#'
#' Exactness matters most here: a product stencil at fourth order keeps about five
#' digits, and the list is large, so a numerical route would be both slow and
#' poor. The row-wise structure is what keeps it affordable at all, the cross-row
#' components never being evaluated.
#' @param s A [TransitionMatrixParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` \eqn{K \times K} matrices keyed as
#'   `param_tuple_names(s, 4)` and in that order. At \eqn{K = 3} that is 126
#'   components, of which 15 can be non-zero.
#' @seealso [tm_derivative()], which assembles it,
#'   [param_d3.TransitionMatrixParam()] for the order below, and [numerical_d4()]
#'   for the alternative.
#' @keywords internal
S7::method(param_d4, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 4L)
}
