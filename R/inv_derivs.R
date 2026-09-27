#' @include generics.R log_cholesky.R
NULL


#' Derivatives of the Inverse of a Matrix Parameter
#'
#' @description
#' `param_inv_d1()` returns the first derivatives of \eqn{M^{-1}} in the free
#' values, and `param_inv_d2()` the second, where \eqn{M} is the matrix
#' [param_value()] returns.
#'
#' @details
#' The default method differentiates the inverse through the derivatives of
#' \eqn{M} itself,
#' \deqn{\partial_k M^{-1} = -M^{-1} A_k M^{-1}, \qquad
#'   \partial_{kl} M^{-1} = M^{-1}(A_l M^{-1} A_k + A_k M^{-1} A_l - A_{kl})
#'   M^{-1},}
#' with \eqn{A_k} and \eqn{A_{kl}} from [param_d1()] and [param_d2()].
#'
#' Where \eqn{M} is nearly singular along a direction that is not a coordinate
#' axis, \eqn{M^{-1}} has large entries of nearly rank one and the products
#' above cancel: at a log-Cholesky free value \eqn{\log L_{22} = -11.5} the
#' default reads \eqn{\partial M^{-1}} with a relative error of
#' \eqn{10^{-7}}, and a trace against it, of products of order \eqn{10^{10}},
#' is out by a number of order one.
#'
#' The method for [log_cholesky()] is exact there. With \eqn{M = L L^\top},
#' \eqn{G = L^{-1}}, \eqn{B_k = G\,\partial_k L} and
#' \eqn{C_k = B_k + B_k^\top},
#' \deqn{\partial_k M^{-1} = -G^\top C_k G, \qquad
#'   \partial_{kl} M^{-1} = G^\top\big(B_l^\top C_k + C_k B_l -
#'   \partial_l C_k\big) G,}
#' where \eqn{\partial_l B_k = -B_l B_k + G\,\partial_{kl} L}. No product
#' carries a cancellation: the entries of \eqn{G} grow as \eqn{1/L_{ii}} and
#' each term is of the size of the result.
#'
#' @param s An object inheriting from class [matrix_parameter()], of full rank.
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. No method in this package reads it.
#'
#' @return `param_inv_d1()`: a list of `s@n_free` matrices, named as
#'   [param_d1()] names its own. `param_inv_d2()`: a list with one matrix per
#'   pair in [param_tuple_indices()] at order 2, named as [param_d2()] names
#'   its own.
#'
#' @seealso [param_solve()] for the inverse itself, [param_d1()] and
#'   [param_d2()] for the derivatives of \eqn{M}.
#'
#' @examples
#' s <- log_cholesky(2)
#' eta <- c(-0.5, -11.5, 0.07)
#' d <- param_inv_d1(s, eta)
#' # against a difference of the inverse, which here is itself accurate
#' h <- 1e-6
#' num <- (param_solve(s, eta + c(0, 0, h)) - param_solve(s, eta - c(0, 0, h))) /
#'   (2 * h)
#' max(abs(d[[3]] - num)) / max(abs(num))
#' names(param_inv_d2(s, eta))
#'
#' @export
param_inv_d1 <- S7::new_generic("param_inv_d1", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  inv_derivs_check(s, "param_inv_d1")
  S7::S7_dispatch()
})


#' @rdname param_inv_d1
#' @export
param_inv_d2 <- S7::new_generic("param_inv_d2", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  inv_derivs_check(s, "param_inv_d2")
  S7::S7_dispatch()
})


#' The Contract Check of the Inverse's Derivatives
#'
#' @description
#' Refuses a parameter that is not a matrix or has no inverse, naming the
#' function that was called.
#'
#' @param s A parameter.
#' @param what The name of the function, for the message.
#'
#' @return `NULL`, invisibly, or an error.
#'
#' @keywords internal
inv_derivs_check <- function(s, what) {
  if (!S7::S7_inherits(s, matrix_parameter)) {
    stop(sprintf(
      "%s() is not defined for '%s': the value is not a symmetric matrix.",
      what, s@param_name
    ), call. = FALSE)
  }
  if (s@rank < s@dimension) {
    stop(sprintf(
      "'%s' is rank deficient (%d of %d), so it has no inverse.",
      s@param_name, s@rank, s@dimension
    ), call. = FALSE)
  }
  invisible(NULL)
}


#' @rdname param_inv_d1
#' @name param_inv_d1.matrix_parameter
#' @keywords internal
S7::method(param_inv_d1, matrix_parameter) <- function(s, eta, ...) {
  si <- unname(param_solve(s, eta))
  a <- param_d1(s, eta)
  out <- lapply(a, function(ak) -(si %*% unname(ak) %*% si))
  names(out) <- names(a)
  out
}


#' @rdname param_inv_d1
#' @name param_inv_d2.matrix_parameter
#' @keywords internal
S7::method(param_inv_d2, matrix_parameter) <- function(s, eta, ...) {
  si <- unname(param_solve(s, eta))
  a <- lapply(param_d1(s, eta), unname)
  a2 <- param_d2(s, eta)
  idx <- param_tuple_indices(s, 2L)
  out <- lapply(seq_along(idx), function(i) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    si %*% (a[[l]] %*% si %*% a[[k]] + a[[k]] %*% si %*% a[[l]] -
              unname(a2[[i]])) %*% si
  })
  names(out) <- names(a2)
  out
}


#' @rdname param_inv_d1
#' @name param_inv_d1.LogCholeskyParam
#' @keywords internal
S7::method(param_inv_d1, LogCholeskyParam) <- function(s, eta, ...) {
  pc <- chol_inv_pieces(s, eta)
  out <- lapply(seq_len(s@n_free), function(k)
    -(t(pc$g) %*% pc$c[[k]] %*% pc$g))
  names(out) <- s@free_names
  out
}


#' @rdname param_inv_d1
#' @name param_inv_d2.LogCholeskyParam
#' @keywords internal
S7::method(param_inv_d2, LogCholeskyParam) <- function(s, eta, ...) {
  pc <- chol_inv_pieces(s, eta)
  idx <- param_tuple_indices(s, 2L)
  out <- lapply(idx, function(kl) {
    k <- kl[1L]
    l <- kl[2L]
    d2l <- chol_dfactor(s, pc$l, c(k, l))
    dbk <- -pc$b[[l]] %*% pc$b[[k]]
    if (!is.null(d2l)) dbk <- dbk + pc$g %*% d2l
    dc <- dbk + t(dbk)
    t(pc$g) %*% (t(pc$b[[l]]) %*% pc$c[[k]] + pc$c[[k]] %*% pc$b[[l]] - dc) %*%
      pc$g
  })
  names(out) <- param_tuple_names(s, 2L)
  out
}


#' The Pieces Behind a Log-Cholesky Inverse's Derivatives
#'
#' @description
#' \eqn{L}, \eqn{G = L^{-1}}, and for each free value \eqn{B_k = G\,\partial_k
#' L} and \eqn{C_k = B_k + B_k^\top}, which [param_inv_d1()] and
#' [param_inv_d2()] read for a [log_cholesky()] parameter.
#'
#' @param s A [LogCholeskyParam()] object.
#' @param eta Its free vector.
#'
#' @return A list with `l`, `g`, `b` and `c`, the last two lists over the
#'   free values.
#'
#' @keywords internal
chol_inv_pieces <- function(s, eta) {
  l <- chol_assemble(s, eta)
  g <- forwardsolve(l, diag(s@dimension))
  b <- lapply(seq_len(s@n_free), function(k) g %*% chol_dfactor(s, l, k))
  list(l = l, g = g, b = b, c = lapply(b, function(bk) bk + t(bk)))
}
