#' @include generics.R log_cholesky.R
NULL


#' Derivatives of the Inverse of a Matrix Parameter
#'
#' @description
#' `param_inv_d1()`, `param_inv_d2()`, `param_inv_d3()` and `param_inv_d4()`
#' return the derivatives of first to fourth order of \eqn{M^{-1}} in the free
#' values, where \eqn{M} is the matrix that [param_value()] returns.
#'
#' @details
#' The default method differentiates the inverse through the derivatives of
#' \eqn{M} itself,
#' \deqn{\partial_k M^{-1} = -M^{-1} A_k M^{-1}, \qquad
#'   \partial_{kl} M^{-1} = M^{-1}(A_l M^{-1} A_k + A_k M^{-1} A_l - A_{kl})
#'   M^{-1},}
#' with \eqn{A_k} and \eqn{A_{kl}} from [param_d1()] and [param_d2()]. At third
#' and fourth order it uses the sum over the ordered set partitions that
#' [inverse_of()] writes out, with the arrays of [param_d3()] and [param_d4()].
#'
#' Where \eqn{M} is nearly singular along a direction that is not a coordinate
#' axis, \eqn{M^{-1}} has large entries of nearly rank one, the products above
#' cancel, and the default method loses accuracy as the condition number of
#' \eqn{M} grows, roughly as its square at the higher orders.
#'
#' The method for [log_cholesky()] is exact there. With \eqn{M = L L^\top},
#' \eqn{G = L^{-1}}, \eqn{B_k = G\,\partial_k L} and
#' \eqn{C_k = B_k + B_k^\top},
#' \deqn{\partial_k M^{-1} = -G^\top C_k G, \qquad
#'   \partial_{kl} M^{-1} = G^\top\big(B_l^\top C_k + C_k B_l -
#'   \partial_l C_k\big) G,}
#' where \eqn{\partial_l B_k = -B_l B_k + G\,\partial_{kl} L}. These products
#' involve no cancellation: the entries of \eqn{G} grow as \eqn{1/L_{ii}} and
#' each term is of the size of the result.
#'
#' At third and fourth order the same method differentiates \eqn{L G = I} and
#' \eqn{M^{-1} = G^\top G} by the Leibniz rule. For a multiset \eqn{S} of free
#' values,
#' \deqn{\partial_S G = -G \sum_{\emptyset \ne U \subseteq S}
#'   (\partial_U L)(\partial_{S \setminus U} G), \qquad
#'   \partial_S M^{-1} = \sum_{U \subseteq S}
#'   (\partial_U G)^\top (\partial_{S \setminus U} G),}
#' the subsets being taken over the positions of \eqn{S}. Every factor is
#' triangular and every derivative of \eqn{L} has a single non-zero entry, so
#' these products do not cancel either.
#'
#' @param s An object inheriting from class [matrix_parameter()], of full rank.
#' @param eta A numeric vector of length `s@n_free`, finite in every entry.
#' @param ... Passed to the method. The methods in this package do not read it.
#'
#' @return `param_inv_d1()`: a list of `s@n_free` matrices, named as
#'   [param_d1()] names its own. `param_inv_d2()`, `param_inv_d3()` and
#'   `param_inv_d4()`: a list with one matrix per tuple in
#'   [param_tuple_indices()] at that order, named as [param_d2()],
#'   [param_d3()] and [param_d4()] name their own.
#'
#' @seealso [param_solve()] for the inverse itself, [param_d1()] to
#'   [param_d4()] for the derivatives of \eqn{M}, and [inverse_of()] for the
#'   family whose value is the inverse.
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
#' # At third order the default route has lost the component that the exact
#' # one keeps: the inverse is quadratic in L2.1, so this one is zero.
#' sandwich <- S7::method(param_inv_d3, matrix_parameter)(s, eta)
#' c(exact = param_inv_d3(s, eta)[["L2.1:L2.1:L2.1"]][2, 2],
#'   sandwich = sandwich[["L2.1:L2.1:L2.1"]][2, 2])
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


#' @rdname param_inv_d1
#' @export
param_inv_d3 <- S7::new_generic("param_inv_d3", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  inv_derivs_check(s, "param_inv_d3")
  S7::S7_dispatch()
})


#' @rdname param_inv_d1
#' @export
param_inv_d4 <- S7::new_generic("param_inv_d4", "s", function(s, eta, ...) {
  eta <- check_eta(s, eta)
  inv_derivs_check(s, "param_inv_d4")
  S7::S7_dispatch()
})


#' Argument Check for the Inverse's Derivatives
#'
#' @description
#' Signals an error when a parameter is not a matrix or has no inverse. The
#' message names the function that was called.
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
#' @name param_inv_d3.matrix_parameter
#' @keywords internal
S7::method(param_inv_d3, matrix_parameter) <- function(s, eta, ...) {
  stats::setNames(inv_partition_sum(s, eta, 3L), param_tuple_names(s, 3L))
}


#' @rdname param_inv_d1
#' @name param_inv_d4.matrix_parameter
#' @keywords internal
S7::method(param_inv_d4, matrix_parameter) <- function(s, eta, ...) {
  stats::setNames(inv_partition_sum(s, eta, 4L), param_tuple_names(s, 4L))
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


#' @rdname param_inv_d1
#' @name param_inv_d3.LogCholeskyParam
#' @keywords internal
S7::method(param_inv_d3, LogCholeskyParam) <- function(s, eta, ...) {
  chol_inv_dk(s, eta, 3L)
}


#' @rdname param_inv_d1
#' @name param_inv_d4.LogCholeskyParam
#' @keywords internal
S7::method(param_inv_d4, LogCholeskyParam) <- function(s, eta, ...) {
  chol_inv_dk(s, eta, 4L)
}


#' Higher Derivatives of a Log-Cholesky Inverse
#'
#' @description
#' Computes the derivatives of one order of \eqn{M^{-1} = G^\top G}, with
#' \eqn{G = L^{-1}}, by the Leibniz rule on \eqn{L G = I} and on
#' \eqn{G^\top G}, as [param_inv_d1()] describes. The derivatives of \eqn{G}
#' are kept by sorted multiset, so that a multiset shared by several
#' components is computed once.
#'
#' @param s A [LogCholeskyParam()] object.
#' @param eta Its free vector, already checked.
#' @param order The derivative order, 3 or 4 in the package's use; any
#'   positive order is computed the same way.
#'
#' @return A list with one matrix per tuple in `param_tuple_indices(s, order)`,
#'   named as `param_tuple_names(s, order)`.
#'
#' @seealso [chol_inv_pieces()] for the first two orders and [leibniz_gram()]
#'   for the Leibniz sum of a Gram product.
#'
#' @keywords internal
chol_inv_dk <- function(s, eta, order) {
  l <- chol_assemble(s, eta)
  p <- s@dimension
  g <- forwardsolve(l, diag(p))
  memo <- new.env(parent = emptyenv())
  dg <- function(ks) {
    if (length(ks) == 0L) return(g)
    key <- paste(sort(ks), collapse = ",")
    hit <- memo[[key]]
    if (!is.null(hit)) return(hit)
    m <- length(ks)
    acc <- matrix(0, p, p)
    # L G = I differentiated over the positions of ks: the term with every
    # position on G is moved to the left and multiplied by G
    for (b in seq_len(2^m - 1L)) {
      take <- bitwAnd(b, bitwShiftL(1L, seq_len(m) - 1L)) > 0L
      dl <- chol_dfactor(s, l, ks[take])
      if (is.null(dl)) next
      acc <- acc + dl %*% dg(ks[!take])
    }
    out <- -(g %*% acc)
    assign(key, out, envir = memo)
    out
  }
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(tp) leibniz_gram(function(ks) t(dg(ks)), tp, p))
  stats::setNames(out, param_tuple_names(s, order))
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
