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


#' The Set Partitions of the First n Integers
#'
#' @description
#' Every way of splitting \code{1:n} into disjoint non-empty blocks, which is
#' what a chain rule of order \eqn{n} sums over.
#'
#' @details
#' Built by the standard recursion: the partitions of \code{1:n} are obtained
#' from those of \code{1:(n-1)} by placing \code{n} into each existing block in
#' turn and then into a block of its own. There are 1, 2, 5 and 15 of them for
#' \eqn{n = 1, \dots, 4}, the Bell numbers.
#'
#' The blocks index \strong{positions} rather than variables, which is what
#' makes a repeated variable count with the right multiplicity without any
#' bookkeeping of its own -- the same device \code{\link{jet_mul}} uses with
#' subsets of positions.
#'
#' @param n A positive integer, at most four here.
#'
#' @return A list of partitions, each a list of integer vectors.
#'
#' @seealso \code{\link{jet_compose}}
#'
#' @keywords internal
set_partitions <- function(n) {
  if (n == 1L) return(list(list(1L)))
  prev <- set_partitions(n - 1L)
  out <- list()
  for (p in prev) {
    for (k in seq_along(p)) {
      q <- p
      q[[k]] <- c(q[[k]], n)
      out[[length(out) + 1L]] <- q
    }
    out[[length(out) + 1L]] <- c(p, list(n))
  }
  out
}


#' A Smooth Function of a Jet
#'
#' @description
#' Applies a scalar function to a jet, given that function's own derivatives at
#' the jet's value, and propagates every partial derivative exactly.
#'
#' @details
#' This is Faa di Bruno in the form the multivariate case takes when the inner
#' function is scalar valued: for a tuple \eqn{T} of positions,
#' \deqn{(f \circ a)_T = \sum_{\pi} f^{(|\pi|)}(a)\prod_{B \in \pi} a_{T[B]}}
#' the sum running over the set partitions of the positions of \eqn{T}.
#'
#' Every transcendental a jet needs is one call to this with the right five
#' numbers, so the exponential, the logarithm, a power and the gamma function
#' need no derivative machinery of their own. Extending the vocabulary is a
#' matter of writing five derivatives, not of writing a chain rule.
#'
#' @param a A jet.
#' @param fd A numeric vector of five values: \eqn{f} and its first four
#'   derivatives, all evaluated at \code{a$v}.
#' @param lay A layout from \code{\link{jet_layout}}.
#'
#' @return A jet.
#'
#' @seealso \code{\link{jet_mul}}, \code{\link{set_partitions}}
#'
#' @keywords internal
jet_compose <- function(a, fd, lay) {
  grab <- function(tup) {
    slot <- get(paste(sort(tup), collapse = ","), envir = lay$pos)
    a$d[[slot[1L]]][slot[2L]]
  }
  out <- jet_const(fd[[1L]], lay)
  for (o in 1:4) {
    parts <- set_partitions(o)
    tt <- lay$tuples[[o]]
    vals <- numeric(length(tt))
    for (i in seq_along(tt)) {
      t <- tt[[i]]
      acc <- 0
      for (p in parts) {
        term <- fd[[length(p) + 1L]]
        if (term == 0) next
        for (b in p) term <- term * grab(t[b])
        acc <- acc + term
      }
      vals[i] <- acc
    }
    out$d[[o]] <- vals
  }
  out
}


#' The Exponential of a Jet
#'
#' @description \eqn{e^a}, whose four derivatives are all \eqn{e^a}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_exp <- function(a, lay) {
  e <- exp(a$v)
  jet_compose(a, rep(e, 5L), lay)
}


#' The Logarithm of a Jet
#'
#' @description \eqn{\log a}, with derivatives \eqn{(-1)^{k-1}(k-1)!/a^k}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_log <- function(a, lay) {
  v <- a$v
  jet_compose(a, c(log(v), 1 / v, -1 / v^2, 2 / v^3, -6 / v^4), lay)
}


#' A Power of a Jet
#'
#' @description
#' \eqn{a^p} for any real exponent, which covers the square root, the
#' reciprocal and the square without a routine for each.
#' @param a A jet.
#' @param p The exponent.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_pow <- function(a, p, lay) {
  v <- a$v
  fd <- c(v^p,
          p * v^(p - 1),
          p * (p - 1) * v^(p - 2),
          p * (p - 1) * (p - 2) * v^(p - 3),
          p * (p - 1) * (p - 2) * (p - 3) * v^(p - 4))
  jet_compose(a, fd, lay)
}


#' The Reciprocal of a Jet
#'
#' @description \eqn{1/a}, the common case of \code{\link{jet_pow}}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_pow}}
#' @keywords internal
jet_inv <- function(a, lay) jet_pow(a, -1, lay)


#' The Square Root of a Jet
#'
#' @description \eqn{\sqrt{a}}, the common case of \code{\link{jet_pow}}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_pow}}
#' @keywords internal
jet_sqrt <- function(a, lay) jet_pow(a, 0.5, lay)


#' The Quotient of Two Jets
#'
#' @description \eqn{a/b}, as a product with the reciprocal.
#' @param a,b Jets over the same layout.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_mul}}
#' @keywords internal
jet_div <- function(a, b, lay) jet_mul(a, jet_inv(b, lay), lay)


#' The Log-Gamma of a Jet
#'
#' @description
#' \eqn{\log\Gamma(a)}, whose derivatives are the polygamma functions.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_lgamma <- function(a, lay) {
  v <- a$v
  jet_compose(a, c(lgamma(v), digamma(v), trigamma(v),
                   psigamma(v, 2L), psigamma(v, 3L)), lay)
}


#' The Gamma Function of a Jet
#'
#' @description
#' \eqn{\Gamma(a)}, formed as \eqn{\exp\log\Gamma(a)} so that the derivatives
#' come from the polygamma functions rather than from a second table.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_lgamma}}
#' @keywords internal
jet_gamma <- function(a, lay) jet_exp(jet_lgamma(a, lay), lay)


#' The Digamma of a Jet
#'
#' @description \eqn{\psi(a)}, with the polygamma functions above it.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_digamma <- function(a, lay) {
  v <- a$v
  jet_compose(a, c(digamma(v), trigamma(v), psigamma(v, 2L),
                   psigamma(v, 3L), psigamma(v, 4L)), lay)
}
