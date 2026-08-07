#' @include chain.R numerical_fallbacks.R
NULL


#' Correlation Matrix Parameter
#'
#' @description
#' The S7 class of correlation matrices in the spherical parametrization.
#' Constructed by \code{\link{correlation_matrix}}.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class \code{CorrelationParam}.
#'
#' @seealso \code{\link{correlation_matrix}}
#'
#' @examples
#' S7::S7_inherits(correlation_matrix(3), CorrelationParam)
#'
#' @export
CorrelationParam <- S7::new_class("CorrelationParam", parent = matrix_parameter)


#' Construct a Correlation Matrix Parameter
#'
#' @description
#' A correlation matrix -- symmetric positive definite with a unit diagonal --
#' carried by the spherical parametrization: \eqn{R = LL^\top} with each row
#' of \eqn{L} a point on the unit sphere written in angular coordinates, and
#' each angle carried onto the real line by a bounded link.
#'
#' @details
#' Row \eqn{i} of \eqn{L} is built from \eqn{i-1} angles
#' \eqn{\theta_{i1}, \ldots, \theta_{i,i-1}} in \eqn{(0, \pi)}:
#' \deqn{L_{ij} = \cos\theta_{ij} \prod_{k < j} \sin\theta_{ik} \quad (j < i),
#'   \qquad L_{ii} = \prod_{k < i} \sin\theta_{ik}.}
#' The squared entries of the row then telescope to one, so \eqn{R = LL^\top}
#' has a unit diagonal by construction rather than by a correction, and it is
#' positive definite for every value of the angles because \eqn{L} is
#' triangular with a positive diagonal. The construction is that of Rapisarda,
#' Brigo and Mercurio (2007), and the free values are the angles carried onto
#' \eqn{\mathbb{R}} by \code{\link[linkfunctions7]{bounded_link}}, so
#' there is nothing to constrain and no boundary to reach.
#'
#' Derivatives follow from the same Leibniz rule the log-Cholesky family uses,
#' since \eqn{R} is again a Gram product; what changes is the factor, whose
#' entries are products of sines and cosines of angles that each depend on one
#' free value. The rows of \eqn{L} are independent, so a derivative of the
#' \emph{factor} in free values from two different rows vanishes. The same is
#' not true of \eqn{R}: its entry \eqn{(i, j)} is the inner product of rows
#' \eqn{i} and \eqn{j} of \eqn{L}, so a second derivative across those two
#' rows need not vanish. What holds instead is that such a component is
#' supported on exactly the entries \eqn{(i, j)} and \eqn{(j, i)}, and even
#' there it is zero whenever the two angles differentiated sit beyond the
#' columns the two rows share.
#'
#' The free vector runs row by row, and the names \code{z\{i\}.\{j\}} say which
#' row and which angle, the row.column convention \code{\link{log_cholesky}}
#' already uses. The ordering is part of the contract: every consumer builds
#' its parameter tables from these names.
#'
#' A correlation matrix on its own is what a copula or an LKJ prior is written
#' against. A covariance is this matrix conjugated by a diagonal of standard
#' deviations, which is a composition of two parameters rather than a family
#' of its own.
#'
#' @param dimension The side \eqn{p} of the matrix.
#' @param role A label; see \code{\link{log_cholesky}}.
#'
#' @return An object of class \code{\link{CorrelationParam}}.
#'
#' @references
#' Rapisarda, F., Brigo, D. and Mercurio, F. (2007). Parameterizing
#' correlations: a geometric interpretation. \emph{IMA Journal of Management
#' Mathematics} 18, 55-73.
#'
#' @seealso \code{\link{log_cholesky}}, \code{\link{compound_symmetry}},
#'   \code{\link{ar1}}
#'
#' @examples
#' s <- correlation_matrix(3)
#' s@free_names
#'
#' r <- param_value(s, c(0.4, -0.2, 0.6))
#' round(r, 4)
#' diag(r)
#'
#' # the round trip closes exactly
#' eta <- c(0.4, -0.2, 0.6)
#' max(abs(param_free(s, param_value(s, eta)) - eta))
#'
#' @export
correlation_matrix <- function(dimension, role = c("either", "covariance", "precision")) {
  role <- match.arg(role)
  p <- check_param_args(dimension, role)

  rows <- integer(0)
  cols <- integer(0)
  if (p > 1L) {
    for (i in seq.int(2L, p)) {
      rows <- c(rows, rep(i, i - 1L))
      cols <- c(cols, seq_len(i - 1L))
    }
  }
  nm <- paste0("z", rows, ".", cols)

  CorrelationParam(
    param_name = "correlation",
    dimension = p,
    n_free = length(nm),
    free_names = nm,
    rank = p,
    null_basis = empty_null_basis(p),
    role = role,
    param_params = list(
      row = rows, col = cols,
      link = linkfunctions7::bounded_link(lwr = 0, upr = pi)
    )
  )
}


#' Sines and Cosines of a Correlation Parameter's Angles
#'
#' @description
#' For every angle, the value and the first four derivatives, in the free
#' value, of both \eqn{\sin\theta} and \eqn{\cos\theta}.
#'
#' @details
#' Each angle depends on exactly one free value, so these tables are all the
#' derivative machinery the family needs: an entry of \eqn{L} is a product of
#' such factors, and differentiating it replaces each factor by the derivative
#' of the matching order. The composition of the trigonometric function with
#' the link is done by \code{\link{compose4}}.
#'
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with \code{sin} and \code{cos}, each a list indexed by free
#'   value holding five numbers: the value and four derivatives.
#'
#' @keywords internal
corr_tables <- function(s, eta) {
  link <- s@param_params$link
  th <- linkfunctions7::linkinv(link, eta)
  td <- list(
    linkfunctions7::dlinkinv(link, eta),
    linkfunctions7::d2linkinv(link, eta),
    linkfunctions7::d3linkinv(link, eta),
    linkfunctions7::d4linkinv(link, eta)
  )
  sn <- vector("list", s@n_free)
  cs <- vector("list", s@n_free)
  for (k in seq_len(s@n_free)) {
    gd <- lapply(td, `[[`, k)
    t0 <- th[[k]]
    fs <- list(cos(t0), -sin(t0), -cos(t0), sin(t0))
    fc <- list(-sin(t0), -cos(t0), sin(t0), cos(t0))
    sn[[k]] <- c(list(sin(t0)), compose4(fs, gd))
    cs[[k]] <- c(list(cos(t0)), compose4(fc, gd))
  }
  list(sin = sn, cos = cs)
}


#' The Factor of a Correlation Parameter, and Its Derivatives
#'
#' @description
#' Returns \eqn{\partial^S L} for a multiset \eqn{S} of free-value indices, or
#' \code{NULL} when that derivative is identically zero.
#'
#' @details
#' An entry of \eqn{L} depends only on the angles of its own row, so a
#' multiset spanning two rows gives zero; within a row, an entry gives zero
#' unless every differentiated angle appears among its factors.
#'
#' @param s A \code{\link{CorrelationParam}} object.
#' @param tb The tables of \code{\link{corr_tables}}.
#' @param ks A multiset of free-value indices, possibly empty.
#'
#' @return A numeric matrix, or \code{NULL}.
#'
#' @keywords internal
corr_dfactor <- function(s, tb, ks) {
  p <- s@dimension
  pos <- s@param_params
  l <- matrix(0, p, p)
  l[1L, 1L] <- if (length(ks)) 0 else 1

  if (length(ks)) {
    rows <- pos$row[ks]
    if (any(rows != rows[1L])) return(NULL)
  }
  target <- if (length(ks)) pos$row[ks[1L]] else NULL

  # multiplicity of each angle position, within the row being differentiated
  mult <- integer(0)
  if (length(ks)) {
    mult <- integer(target - 1L)
    for (k in ks) mult[pos$col[k]] <- mult[pos$col[k]] + 1L
  }
  # the free-value index of angle (i, k)
  idx_of <- function(i, k) which(pos$row == i & pos$col == k)

  for (i in seq_len(p)) {
    if (i == 1L) next
    if (!is.null(target) && i != target) next
    for (j in seq_len(i)) {
      # the angles this entry is a product of, and the role each plays
      angles <- if (j < i) seq_len(j) else seq_len(i - 1L)
      roles <- if (j < i) c(rep("sin", j - 1L), "cos") else rep("sin", i - 1L)
      if (!length(angles)) {
        l[i, j] <- if (length(ks)) 0 else 1
        next
      }
      if (length(ks) && any(mult[setdiff(seq_len(i - 1L), angles)] > 0L)) {
        l[i, j] <- 0
        next
      }
      v <- 1
      for (a in seq_along(angles)) {
        k <- angles[a]
        m <- if (length(ks)) mult[k] else 0L
        tab <- if (roles[a] == "sin") tb$sin else tb$cos
        v <- v * tab[[idx_of(i, k)]][[m + 1L]]
      }
      l[i, j] <- v
    }
  }
  l
}


#' @title Value of a Correlation Parameter
#' @name param_value.CorrelationParam
#' @description \eqn{R = LL^\top}, with \eqn{L} assembled from the angles.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A correlation matrix.
#' @keywords internal
S7::method(param_value, CorrelationParam) <- function(s, eta, ...) {
  tb <- corr_tables(s, eta)
  l <- corr_dfactor(s, tb, integer(0))
  m <- tcrossprod(l)
  # the diagonal is one by construction; the assignment removes the rounding
  # of a telescoping sum of squares, which would otherwise leave 1 - 2e-16
  diag(m) <- 1
  name_dims(m, s)
}


#' @title Factor of a Correlation Parameter
#' @name param_factor.CorrelationParam
#' @description The factor is the parametrization, so it needs no computing.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A lower triangular numeric matrix.
#' @keywords internal
S7::method(param_factor, CorrelationParam) <- function(s, eta, ...) {
  corr_dfactor(s, corr_tables(s, eta), integer(0))
}


#' @title Free Vector of a Correlation Parameter
#' @name param_free.CorrelationParam
#' @description
#' The angles read off the Cholesky factor, exact: \eqn{\theta_{i1}} is the
#' arc cosine of \eqn{L_{i1}}, and each subsequent angle divides out the sines
#' already recovered. A matrix that is not a correlation matrix, or one whose
#' factor reaches an angle of \eqn{0} or \eqn{\pi}, is rejected.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param m A correlation matrix.
#' @param ... Unused.
#' @return A named numeric vector of free values.
#' @keywords internal
S7::method(param_free, CorrelationParam) <- function(s, m, ...) {
  p <- s@dimension
  if (max(abs(diag(m) - 1)) > 1e-8) {
    stop(paste0(
      "'m' does not have a unit diagonal, so it is not a correlation matrix.\n",
      "  It is rejected rather than rescaled."
    ), call. = FALSE)
  }
  l <- chol_pd(m)
  if (is.null(l)) {
    stop(paste0(
      "'m' is not positive definite, so it is not in the set\n",
      "  correlation_matrix() parametrizes. The verdict is spectral."
    ), call. = FALSE)
  }
  pos <- s@param_params
  th <- numeric(s@n_free)
  for (k in seq_len(s@n_free)) {
    i <- pos$row[k]
    j <- pos$col[k]
    denom <- 1
    if (j > 1L) {
      prev <- which(pos$row == i & pos$col < j)
      denom <- prod(sin(th[prev]))
    }
    if (!is.finite(denom) || abs(denom) < 1e-12) {
      stop(paste0(
        "'m' sits on the boundary of the set: an angle of the factor is 0 or\n",
        "  pi, which no finite free value produces."
      ), call. = FALSE)
    }
    v <- l[i, j] / denom
    th[k] <- acos(min(1, max(-1, v)))
  }
  stats::setNames(
    linkfunctions7::linkfun(s@param_params$link, th), s@free_names
  )
}


#' Derivative Components of a Correlation Parameter
#'
#' @description
#' Assembles one derivative order by the Leibniz rule on \eqn{R = LL^\top}.
#'
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of symmetric matrices.
#'
#' @keywords internal
corr_derivative <- function(s, eta, order) {
  tb <- corr_tables(s, eta)
  dfac <- function(ks) corr_dfactor(s, tb, ks)
  idx <- param_tuple_indices(s, order)
  out <- lapply(idx, function(t) {
    name_dims(leibniz_gram(dfac, t, s@dimension), s)
  })
  stats::setNames(out, param_tuple_names(s, order))
}


#' @title First Derivatives of a Correlation Parameter
#' @name param_d1.CorrelationParam
#' @description
#' Closed form by the Leibniz rule on \eqn{R = LL^\top}, the factor's
#' derivatives coming from the trigonometric tables of the angles.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices with a zero diagonal, the
#'   diagonal of the value being constant.
#' @keywords internal
S7::method(param_d1, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 1L)
}

#' @title Second Derivatives of a Correlation Parameter
#' @name param_d2.CorrelationParam
#' @description Closed form; see \code{\link{param_d1.CorrelationParam}}.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d2, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 2L)
}

#' @title Third Derivatives of a Correlation Parameter
#' @name param_d3.CorrelationParam
#' @description Closed form; see \code{\link{param_d1.CorrelationParam}}.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d3, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 3L)
}

#' @title Fourth Derivatives of a Correlation Parameter
#' @name param_d4.CorrelationParam
#' @description Closed form; see \code{\link{param_d1.CorrelationParam}}.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named list of symmetric matrices.
#' @keywords internal
S7::method(param_d4, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 4L)
}


#' Log-Determinant Chains of a Correlation Parameter
#'
#' @description
#' For each free value, the four derivatives of \eqn{2\log\sin\theta} in that
#' free value: the log-determinant's contribution from one angle.
#'
#' @details
#' The factor is triangular, so \eqn{\lvert R \rvert = \prod_i L_{ii}^2} and
#' \deqn{\log\lvert R \rvert = 2 \sum_{i,k} \log \sin\theta_{ik},}
#' a sum with one term per free value. The log-determinant is therefore
#' separable, every mixed derivative is exactly zero, and each pure one is the
#' logarithm composed with the sine table \code{\link{corr_tables}} already
#' holds.
#'
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#'
#' @return A list with one element per free value, each a list of four
#'   numbers.
#'
#' @keywords internal
corr_logdet_chains <- function(s, eta) {
  tb <- corr_tables(s, eta)
  lapply(seq_len(s@n_free), function(k) {
    sn <- tb$sin[[k]]
    f <- lapply(1:4, function(j) {
      2 * (-1)^(j - 1L) * factorial(j - 1L) / sn[[1L]]^j
    })
    compose4(f, sn[-1L])
  })
}


#' @title Log-Determinant of a Correlation Parameter
#' @name param_logdet.CorrelationParam
#' @description
#' Closed form: twice the sum of the logarithms of the sines of the angles,
#' the factor being triangular with those products on its diagonal. No
#' factorization and no determinant is computed.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(param_logdet, CorrelationParam) <- function(s, eta, ...) {
  if (s@n_free == 0L) return(0)
  tb <- corr_tables(s, eta)
  2 * sum(vapply(tb$sin, function(x) log(x[[1L]]), numeric(1)))
}


#' Log-Determinant Components of a Correlation Parameter
#'
#' @description
#' Assembles one derivative order of the log-determinant from the per-angle
#' chains, every mixed component being exactly zero.
#'
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named numeric vector.
#'
#' @keywords internal
corr_logdet_derivative <- function(s, eta, order) {
  nm <- param_tuple_names(s, order)
  if (!length(nm)) return(stats::setNames(numeric(0), character(0)))
  ch <- corr_logdet_chains(s, eta)
  idx <- param_tuple_indices(s, order)
  out <- vapply(idx, function(t) {
    if (any(t != t[1L])) return(0)
    ch[[t[1L]]][[order]]
  }, numeric(1))
  stats::setNames(out, nm)
}


#' @title Log-Determinant Derivatives of a Correlation Parameter
#' @name param_dlogdet.CorrelationParam
#' @description
#' Closed form at every order. The log-determinant is a sum with one term per
#' angle, so it is separable: every mixed component is exactly zero and each
#' pure one is a logarithm composed with that angle's sine.
#' @param s A \code{\link{CorrelationParam}} object.
#' @param eta A numeric vector of free values.
#' @param ... Unused.
#' @return A named numeric vector.
#' @keywords internal
S7::method(param_dlogdet, CorrelationParam) <- function(s, eta, ...) {
  corr_logdet_derivative(s, eta, 1L)
}

#' @rdname param_dlogdet.CorrelationParam
#' @name param_d2logdet.CorrelationParam
#' @keywords internal
S7::method(param_d2logdet, CorrelationParam) <- function(s, eta, ...) {
  corr_logdet_derivative(s, eta, 2L)
}

#' @rdname param_dlogdet.CorrelationParam
#' @name param_d3logdet.CorrelationParam
#' @keywords internal
S7::method(param_d3logdet, CorrelationParam) <- function(s, eta, ...) {
  corr_logdet_derivative(s, eta, 3L)
}

#' @rdname param_dlogdet.CorrelationParam
#' @name param_d4logdet.CorrelationParam
#' @keywords internal
S7::method(param_d4logdet, CorrelationParam) <- function(s, eta, ...) {
  corr_logdet_derivative(s, eta, 4L)
}
