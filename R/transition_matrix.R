#' @include simplex.R
NULL


#' Transition Matrix Parameter
#'
#' @description
#' The S7 class of row-stochastic matrices, each row on the open simplex in
#' the additive log-ratio parametrisation. Constructed by
#' \code{\link{transition_matrix}}.
#'
#' @inheritParams parameter
#'
#' @return An object of class \code{TransitionMatrixParam}.
#'
#' @seealso \code{\link{transition_matrix}}
#'
#' @examples
#' S7::S7_inherits(transition_matrix(3), TransitionMatrixParam)
#'
#' @export
TransitionMatrixParam <- S7::new_class("TransitionMatrixParam", parent = parameter)


#' Construct a Transition Matrix Parameter
#'
#' @description
#' A \eqn{K \times K} row-stochastic matrix -- the transition matrix of a
#' Markov chain on \eqn{K} states -- with each row an independent
#' \code{\link{simplex}} in the additive log-ratio chart, so
#' \eqn{K(K-1)} free values in all.
#'
#' @details
#' The rows are independent in the parametrisation, so every derivative
#' tensor is block diagonal by row: a component pairing free values of two
#' different rows is exactly zero, and the implementation evaluates the
#' simplex kernels row by row rather than storing those zeros.
#'
#' The free vector runs row by row, and the names \code{alr\{i\}.\{j\}} say
#' which row and which chart coordinate, the row.column convention the
#' log-Cholesky names already use. Rows, not columns, live on the simplex: a
#' transition matrix acts on row vectors of probabilities, and the column
#' convention is its transpose.
#'
#' @param n_state The number of states \eqn{K}, at least 2.
#'
#' @return An object of class \code{\link{TransitionMatrixParam}}.
#'
#' @seealso \code{\link{simplex}}
#'
#' @examples
#' s <- transition_matrix(3)
#' eta <- rnorm(s@n_free)
#' rowSums(param_value(s, eta))
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
#' The row each free value belongs to and its position inside that row's
#' chart, in the order the free vector uses.
#'
#' @param s A \code{\link{TransitionMatrixParam}} object.
#'
#' @return A list with integer vectors \code{row} and \code{coord}.
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
#' @description Each row is the softmax of its own free values.
#' @param s A \code{\link{TransitionMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A \eqn{K \times K} row-stochastic matrix.
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
#' The additive log-ratio of each row, exact; refused when any row is outside
#' the open simplex or fails to sum to one.
#' @param s A \code{\link{TransitionMatrixParam}} object.
#' @param m A \eqn{K \times K} row-stochastic matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
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
      "  It is refused rather than renormalised."
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
#' Assembles a derivative order from the row-wise simplex tensors: a
#' component whose free values span two rows is zero, and one inside a row
#' embeds that row's simplex component.
#'
#' @param s A \code{\link{TransitionMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of \eqn{K \times K} matrices.
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
#' @description Closed form, row by row; see \code{\link{tm_derivative}}.
#' @param s A \code{\link{TransitionMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of matrices.
#' @keywords internal
S7::method(param_d1, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 1L)
}

#' @title Second Derivatives of a Transition Matrix Parameter
#' @name param_d2.TransitionMatrixParam
#' @description Closed form, row by row.
#' @param s A \code{\link{TransitionMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of matrices.
#' @keywords internal
S7::method(param_d2, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 2L)
}

#' @title Third Derivatives of a Transition Matrix Parameter
#' @name param_d3.TransitionMatrixParam
#' @description Closed form, row by row.
#' @param s A \code{\link{TransitionMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of matrices.
#' @keywords internal
S7::method(param_d3, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 3L)
}

#' @title Fourth Derivatives of a Transition Matrix Parameter
#' @name param_d4.TransitionMatrixParam
#' @description Closed form, row by row.
#' @param s A \code{\link{TransitionMatrixParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of matrices.
#' @keywords internal
S7::method(param_d4, TransitionMatrixParam) <- function(s, eta, ...) {
  tm_derivative(s, eta, 4L)
}
