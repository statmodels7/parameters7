#' @include chain.R generics.R
NULL

# A jet is a number carried together with all its partial derivatives up to
# fourth order in the free values. Adding and multiplying jets propagates
# those derivatives exactly, so a quantity built from a recursion of sums and
# products -- which is what the Levinson-Durbin recursion behind an AR(p)
# parameter is -- comes out analytic at every order with no formula
# transcribed. The alternative, expanding the recursion by hand, is pages of
# algebra whose only check would be the very finite differences the toolkit
# refuses to trust at fourth order.


#' The Bookkeeping a Jet Needs
#'
#' @description
#' The index tuples up to fourth order over \eqn{d} variables, together with
#' the lookup from a tuple to its position, computed once and shared by every
#' jet in a calculation.
#'
#' @details
#' The tuples are the package's own enumeration, so a jet's components are
#' keyed exactly as \code{\link{param_tuple_names}} keys a derivative list and
#' nothing has to be reordered on the way out.
#'
#' @param d The number of variables.
#'
#' @return A list with \code{d}, the tuple lists \code{tuples}, and the
#'   environment \code{pos} mapping a sorted tuple to its order and position.
#'
#' @seealso \code{\link{jet_mul}}
#'
#' @keywords internal
jet_layout <- function(d) {
  tuples <- lapply(1:4, function(o) tuple_indices(d, o))
  # The empty tuple gets no key: R has no zero-length variable name, and the
  # lookups short-circuit on it before reaching the environment anyway.
  pos <- new.env(parent = emptyenv())
  for (o in 1:4) {
    for (i in seq_along(tuples[[o]])) {
      assign(paste(sort(tuples[[o]][[i]]), collapse = ","), c(o, i), envir = pos)
    }
  }
  list(d = d, tuples = tuples, pos = pos)
}


#' A Constant Jet
#'
#' @description
#' A number with every derivative zero.
#'
#' @param v The value.
#' @param lay A layout from \code{\link{jet_layout}}.
#'
#' @return A jet.
#'
#' @keywords internal
jet_const <- function(v, lay) {
  list(v = v, d = lapply(1:4, function(o) numeric(length(lay$tuples[[o]]))))
}


#' A Jet for One Variable
#'
#' @description
#' The jet of a scalar map of a single free value: its value and its four
#' derivatives sit in the pure components of that variable, everything else
#' being zero.
#'
#' @param k The index of the free value.
#' @param dv A list of five numbers: the value and four derivatives.
#' @param lay A layout from \code{\link{jet_layout}}.
#'
#' @return A jet.
#'
#' @keywords internal
jet_var <- function(k, dv, lay) {
  out <- jet_const(dv[[1L]], lay)
  for (o in 1:4) {
    slot <- get(paste(rep(k, o), collapse = ","), envir = lay$pos)
    out$d[[o]][slot[2L]] <- dv[[o + 1L]]
  }
  out
}


#' Sum of Jets
#'
#' @description
#' Adds two jets, or a jet and a number, componentwise.
#'
#' @param a,b Jets, or a jet and a single number.
#'
#' @return A jet.
#'
#' @keywords internal
jet_add <- function(a, b) {
  list(v = a$v + b$v, d = lapply(1:4, function(o) a$d[[o]] + b$d[[o]]))
}


#' Product of Jets
#'
#' @description
#' Multiplies two jets, propagating every derivative exactly.
#'
#' @details
#' The Leibniz rule again, in its scalar form: the component of the product at
#' a tuple \eqn{T} is \eqn{\sum_{S \subseteq T} a_S b_{T \setminus S}}, the sum
#' running over subsets of \emph{positions} so that a repeated variable is
#' counted with the right multiplicity without a bookkeeping of its own. It is
#' the same enumeration \code{\link{leibniz_gram}} uses for a matrix product.
#'
#' @param a,b Jets over the same layout.
#' @param lay A layout from \code{\link{jet_layout}}.
#'
#' @return A jet.
#'
#' @keywords internal
jet_mul <- function(a, b, lay) {
  grab <- function(x, tup) {
    if (!length(tup)) return(x$v)
    slot <- get(paste(sort(tup), collapse = ","), envir = lay$pos)
    x$d[[slot[1L]]][slot[2L]]
  }
  out <- jet_const(a$v * b$v, lay)
  for (o in 1:4) {
    tt <- lay$tuples[[o]]
    vals <- numeric(length(tt))
    positions <- seq_len(o)
    for (i in seq_along(tt)) {
      t <- tt[[i]]
      acc <- 0
      for (m in 0:(2^o - 1L)) {
        take <- as.logical(bitwAnd(m, bitwShiftL(1L, positions - 1L)) > 0L)
        acc <- acc + grab(a, t[take]) * grab(b, t[!take])
      }
      vals[i] <- acc
    }
    out$d[[o]] <- vals
  }
  out
}


#' A Jet Times a Constant
#'
#' @description
#' Scales a jet by a plain number, which is the common case and avoids
#' building a constant jet for it.
#'
#' @param a A jet.
#' @param c A single number.
#'
#' @return A jet.
#'
#' @keywords internal
jet_cmul <- function(a, c) {
  list(v = a$v * c, d = lapply(a$d, function(x) x * c))
}
