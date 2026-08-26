#' @include chain.R numerical_fallbacks.R
NULL


#' Correlation Matrix Parameter
#'
#' @description
#' The S7 class of correlation matrices, symmetric positive definite with a unit
#' diagonal, in the spherical parametrization of Rapisarda, Brigo and Mercurio
#' (2007). The unit diagonal holds by construction, never by a correction, so
#' every free vector gives a genuine correlation matrix.
#'
#' [correlation_matrix()] builds one. The free values are angles carried onto the
#' real line, and `param_params` records the row and column each belongs to
#' together with the `bounded_link(0, pi)` that carries it.
#'
#' @inheritParams matrix_parameter
#'
#' @return An object of class `CorrelationParam`, a subclass of
#'   [matrix_parameter()] adding no properties of its own. `param_params` holds
#'   `row` and `col`, the position each angle belongs to, and `link`, a
#'   `linkfunctions7::bounded_link(lwr = 0, upr = pi)`. `n_free` is
#'   \eqn{p(p-1)/2} and `rank` is \eqn{p}.
#'
#' @references
#' Rapisarda, F., Brigo, D. and Mercurio, F. (2007). Parameterizing
#' correlations: a geometric interpretation. *IMA Journal of Management
#' Mathematics* **18**, 55-73.
#'
#' @seealso [correlation_matrix()], the constructor, [dr_prod()] to give this a
#'   diagonal scale and make it a covariance, and [matrix_parameter()] for the
#'   properties this inherits.
#'
#' @examples
#' s <- correlation_matrix(3)
#' S7::S7_inherits(s, CorrelationParam)
#'
#' # The free values are angles: p(p-1)/2 of them, one per below-diagonal entry.
#' c(n_free = s@n_free, p_choose_2 = 3 * 2 / 2)
#' s@free_names
#'
#' # Every free vector gives a unit diagonal exactly, not approximately.
#' diag(param_value(s, c(2.5, -3, 1.7)))
#'
#' @export
CorrelationParam <- S7::new_class("CorrelationParam", parent = matrix_parameter)


#' Construct a Correlation Matrix Parameter
#'
#' @description
#' Returns an object holding the spherical parametrization of a correlation
#' matrix: symmetric, positive definite, and with a unit diagonal that holds
#' exactly, with no correction applied afterwards. Each row of the Cholesky
#' factor is a point on the unit sphere written in angular coordinates, and each
#' angle is carried onto the whole real line by a bounded link, so any free
#' vector in \eqn{\mathbb{R}^{p(p-1)/2}} gives a valid correlation matrix.
#'
#' Use it when the correlations are what the model is about and a scale belongs
#' elsewhere: a copula, an LKJ prior, or the correlation half of a covariance
#' built with [dr_prod()].
#'
#' @details
#' # The construction
#'
#' Row \eqn{i} of \eqn{L} is built from \eqn{i-1} angles
#' \eqn{\theta_{i1}, \dots, \theta_{i,i-1}} in \eqn{(0, \pi)}:
#'
#' \deqn{L_{ij} = \cos\theta_{ij} \prod_{k < j} \sin\theta_{ik} \quad (j < i),
#'   \qquad L_{ii} = \prod_{k < i} \sin\theta_{ik}.}
#'
#' The squared entries of a row telescope to one, so every row of \eqn{L} is a
#' unit vector and \eqn{R = LL^\top} has a unit diagonal by construction. It is
#' positive definite at every value of the angles, \eqn{L} being triangular with
#' a positive diagonal. The angles reach the free scale through
#' `linkfunctions7::bounded_link(lwr = 0, upr = pi)`, so there is nothing to
#' constrain and no boundary to run into.
#'
#' # How the derivatives behave, and a claim that does not hold
#'
#' Derivatives come from the same Leibniz rule the log-Cholesky family uses,
#' \eqn{R} being a Gram product again; what changes is the factor, whose entries
#' are products of sines and cosines of angles that each depend on one free
#' value.
#'
#' The rows of \eqn{L} are independent, so a derivative of the **factor** in free
#' values from two different rows is zero. That is **not** true of \eqn{R}: its
#' entry \eqn{(i, j)} is the inner product of rows \eqn{i} and \eqn{j} of
#' \eqn{L}, so a second derivative across those two rows need not vanish. What
#' does hold is that such a component is supported on exactly the entries
#' \eqn{(i, j)} and \eqn{(j, i)}, and is zero even there when the two angles
#' differentiated sit beyond the columns the two rows share. A structural claim
#' about a product is not a claim about its factors.
#'
#' # The log-determinant is separable
#'
#' The factor is triangular, so \eqn{|R| = \prod_i L_{ii}^2} and
#'
#' \deqn{\log|R| = 2 \sum_{i,k} \log \sin\theta_{ik},}
#'
#' one term per free value. Every mixed derivative of the log-determinant is
#' therefore exactly zero at every order, and each pure one is a logarithm
#' composed with that angle's sine. No factorization and no determinant is
#' computed.
#'
#' # Reading the free vector
#'
#' It runs row by row, and the names `z{i}.{j}` say which row and which angle,
#' the row.column convention [log_cholesky()] uses. The `z` records the link, so
#' a free value of 0.4 is not an angle of 0.4; the angle is
#' `bounded_link(0, pi)`'s inverse of it. The ordering is part of the interface.
#'
#' # Where double precision gives out
#'
#' Measured at \eqn{p = 3} with every free value equal: at 5 the smallest
#' eigenvalue is \eqn{3 \times 10^{-8}} and at 10 it is \eqn{-2 \times 10^{-16}},
#' the matrix having reached its boundary in double precision with a correlation
#' of \eqn{-1} to every printed digit. The parametrization is exact; what gives
#' out is the representation of a correlation one ulp from the edge. A fit that
#' walks a free value past about \eqn{\pm 8} is reporting a boundary, and
#' [check_parameter()] sweeps to \eqn{\pm 2} for that reason.
#'
#' # p = 1 is a constant
#'
#' A \eqn{1 \times 1} correlation matrix has no angles, so `n_free` is 0 and
#' `free_names` is empty. [param_value()] returns the \eqn{1 \times 1} identity
#' at the only free vector there is, `numeric(0)`, [param_d1()] and [param_d2()]
#' return empty lists, and [param_logdet()] is 0. [check_parameter()] passes,
#' reporting its two log-determinant derivative rows as `NOT CHECKED`, there
#' being no free value to differentiate in. The family therefore degenerates to
#' a constant rather than refusing, which is what keeps it composable inside
#' [block_diag()]; it carries no information of its own.
#'
#' @section Notation:
#' \eqn{\eta} is the free vector, of length \eqn{d = p(p-1)/2}, and \eqn{p} the
#' side of the matrix. \eqn{\theta_{ij} \in (0, \pi)} is the \eqn{j}-th angle of
#' row \eqn{i}, \eqn{L} the lower triangular factor and \eqn{R = LL^\top} the
#' correlation matrix. Note that \eqn{\theta} is an angle here and a distribution
#' parameter elsewhere in the toolkit.
#'
#' @param dimension The side \eqn{p} of the matrix. A single positive whole
#'   number, finite and at least 1; \eqn{p = 1} gives a constant with no free
#'   values, as **Details** describes. Anything else throws `'dimension' must be
#'   a single positive integer.`
#' @param role A label recording which side of a model the matrix parametrizes:
#'   `"either"` (the default), `"covariance"` or `"precision"`. No numeric result
#'   depends on it.
#'
#' @return An object of class [CorrelationParam()], with `n_free` equal to
#'   \eqn{p(p-1)/2}, `free_names` `z2.1`, `z3.1`, `z3.2`, ... row by row, `rank`
#'   equal to `dimension`, an empty `null_basis`, `param_name` `"correlation"`,
#'   and `param_params` holding `row`, `col` and `link`.
#'
#' @references
#' Rapisarda, F., Brigo, D. and Mercurio, F. (2007). Parameterizing
#' correlations: a geometric interpretation. *IMA Journal of Management
#' Mathematics* **18**, 55-73.
#'
#' @seealso [dr_prod()] to combine this with a diagonal of standard deviations
#'   into a covariance, [log_cholesky()] for an unstructured matrix,
#'   [compound_symmetry()] and [ar1()] for structured correlations with two free
#'   values, and [param_value()] for the map.
#'
#' @examples
#' # Three angles for a 3 x 3 correlation matrix, named by row and column.
#' s <- correlation_matrix(3)
#' s@free_names
#'
#' eta <- c(0.4, -0.2, 0.6)
#' r <- param_value(s, eta)
#' round(r, 4)
#'
#' # The unit diagonal is exact, and the matrix is positive definite.
#' diag(r)
#' eigen(r, only.values = TRUE)$values > 0
#'
#' # The round trip closes exactly.
#' max(abs(param_free(s, r) - eta))
#'
#' # Absurd free values stay in the set, which is the point.
#' m <- param_value(s, c(-30, 40, -20))
#' c(min_diag = min(diag(m)), max_abs_corr = max(abs(m[upper.tri(m)])))
#'
#' # The log-determinant is twice the sum of the logs of the sines, so it needs
#' # no factorization, and it agrees with the eigenvalues.
#' theta <- linkfunctions7::linkinv(s@param_params$link, eta)
#' c(closed_form = param_logdet(s, eta), from_sines = 2 * sum(log(sin(theta))),
#'   from_eigen = sum(log(eigen(r, only.values = TRUE)$values)))
#'
#' # And it is separable, so every mixed second derivative is exactly zero.
#' round(param_d2logdet(s, eta), 10)
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
  # paste0 recycles a zero-length argument against the length-one literals, so
  # at p = 1, where there are no angles, the unguarded call gives "z." instead
  # of nothing. A one by one correlation matrix is the constant 1.
  nm <- if (length(rows)) paste0("z", rows, ".", cols) else character(0)

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
#' Returns, for every angle, the value and the first four derivatives **in the
#' free value** of both \eqn{\sin\theta} and \eqn{\cos\theta}. These tables are
#' the whole derivative machinery of the family: an entry of \eqn{L} is a product
#' of such factors, so differentiating it replaces each factor by the derivative
#' of the matching order, and nothing else has to be derived.
#'
#' @details
#' Each angle depends on exactly one free value, so a table is indexed by free
#' value and no cross terms appear at this level. The chain from
#' \eqn{\mathrm{d}/\mathrm{d}\theta} to \eqn{\mathrm{d}/\mathrm{d}\eta} is
#' [compose4()]'s, applied to the link's own four derivatives, so the accuracy is
#' the link's and nothing is differenced.
#'
#' @param s A [CorrelationParam()] object, whose `param_params$link` is read.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list with two components, `sin` and `cos`, each a list of
#'   `s@n_free` numeric vectors of length 5: the value at index 1 and the four
#'   derivatives in the free value at indices 2 to 5.
#'
#' @seealso [compose4()] for the chain rule, [corr_dfactor()] which reads these
#'   tables, and [corr_logdet_chains()], which reads the `sin` half.
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
#' `NULL` where that derivative is identically zero, which is most of the time.
#' The empty multiset gives \eqn{L} itself.
#'
#' @details
#' Two rules make almost everything vanish. An entry of \eqn{L} depends only on
#' the angles of **its own row**, so a multiset spanning two rows gives the zero
#' matrix. Within a row, an entry is a product over a prefix of the angles, so it
#' gives zero unless every differentiated angle appears among its factors: a
#' derivative in \eqn{\theta_{ij}} touches only the entries of row \eqn{i} from
#' column \eqn{j} onwards.
#'
#' Note what this does **not** say about \eqn{R}. The factor's derivative across
#' two rows is zero; the derivative of \eqn{R = LL^\top} across those rows is not,
#' \eqn{R_{ij}} being the inner product of two rows. See
#' [correlation_matrix()] for the support that does hold.
#'
#' @param s A [CorrelationParam()] object, whose `param_params$row` and
#'   `param_params$col` are read.
#' @param tb The tables of [corr_tables()], evaluated at the point.
#' @param ks A multiset of free-value indices, possibly empty. Positions in
#'   `1:s@n_free`.
#'
#' @return A `s@dimension` by `s@dimension` lower triangular numeric matrix, or
#'   `NULL` where the derivative is identically zero.
#'
#' @seealso [corr_tables()] for the input, [corr_derivative()] for the Leibniz
#'   sum that consumes it, and [chol_dfactor()], the same idea for the
#'   log-Cholesky family.
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
#' @description
#' Returns \eqn{R = LL^\top}, with \eqn{L} assembled from the angles: row \eqn{i}
#' is a unit vector in spherical coordinates, so the diagonal of \eqn{R} is
#' exactly 1 and its off-diagonal entries are correlations. Nothing is tested and
#' nothing is corrected; the unit diagonal and the positive definiteness are
#' properties of the construction. The cost is one \eqn{O(p^3)} product.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` correlation matrix: symmetric,
#'   positive definite, unit diagonal, dimnames `v1`, `v2`, ...
#' @seealso [param_free.CorrelationParam()] for the inverse, and
#'   [param_factor.CorrelationParam()] for \eqn{L} itself.
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
#' @description
#' Returns \eqn{L} by assembling it from the angles, so no factorization is
#' taken: the factor is what the parametrization holds. Its rows are unit
#' vectors, which is where \eqn{R = LL^\top} gets its unit diagonal.
#'
#' It is the natural way to simulate a correlated standard normal vector:
#' \eqn{Lz} with \eqn{z} standard normal has correlation matrix \eqn{R} and unit
#' variances.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A `s@dimension` by `s@dimension` lower triangular numeric matrix with
#'   a positive diagonal and unit row norms, satisfying
#'   `L %*% t(L) == param_value(s, eta)`, and carrying no dimnames.
#' @seealso [param_value.CorrelationParam()] for the matrix, and
#'   [param_factor.matrix_parameter()] for what a family without a closed form
#'   pays.
#' @keywords internal
S7::method(param_factor, CorrelationParam) <- function(s, eta, ...) {
  corr_dfactor(s, corr_tables(s, eta), integer(0))
}


#' @title Free Vector of a Correlation Parameter
#' @name param_free.CorrelationParam
#' @description
#' Returns the angles behind a correlation matrix, read off its Cholesky factor
#' and carried onto the free scale by the link. Exact:
#' \eqn{\theta_{i1} = \arccos L_{i1}}, and each subsequent angle divides out the
#' sines already recovered before taking an arc cosine. A true inverse of
#' [param_value.CorrelationParam()], the round trip closing to
#' \eqn{2 \times 10^{-16}}.
#' @details
#' Three rejections. `m` must have a unit diagonal, or it is not a correlation
#' matrix. It must be positive definite, tested spectrally through [chol_pd()].
#' And its factor must not reach an angle of exactly 0 or \eqn{\pi}, where the
#' link has no finite value: that happens at a correlation of \eqn{\pm 1}, which
#' is on the boundary of the set without being in it.
#'
#' The last case is worth knowing before it is met. The parametrization reaches a
#' correlation of \eqn{-1} to every printed digit at a free value near 10, and
#' the matrix is then singular in double precision, so a matrix produced by a fit
#' that ran to its boundary cannot be inverted back.
#' @param s A [CorrelationParam()] object.
#' @param m A correlation matrix of side `s@dimension`: symmetric with a unit
#'   diagonal and positive definite, already checked for shape and symmetry by
#'   the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A numeric vector of length `s@n_free`, named by `s@free_names`.
#' @seealso [param_value.CorrelationParam()], the map this inverts, and
#'   [chol_pd()] for the definiteness test.
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
#' Assembles a whole derivative order of \eqn{R = LL^\top} by the Leibniz rule,
#' the factor's derivatives coming from [corr_dfactor()] and the trigonometric
#' tables of [corr_tables()]. The four methods [param_d1()] through
#' [param_d4()] of this family are one call each to this function.
#'
#' @details
#' Every term of the Leibniz sum whose factor derivative is `NULL` is skipped, and
#' most are: a multiset spanning two rows of \eqn{L} contributes nothing. The
#' diagonal of the result is zero at every order above zero, the diagonal of
#' \eqn{R} being the constant 1.
#'
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A list of `choose(s@n_free + order - 1, order)` symmetric matrices
#'   keyed as `param_tuple_names(s, order)` and in that order, each
#'   `s@dimension` by `s@dimension` with a zero diagonal.
#'
#' @seealso [corr_dfactor()] for the factor's derivatives, [leibniz_gram()] for
#'   the sum, and [chol_leibniz()], the same construction for the log-Cholesky
#'   family.
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
#' Closed form, by the Leibniz rule on \eqn{R = LL^\top},
#'
#' \deqn{\partial_k R = L_k L^\top + L L_k^\top,}
#'
#' with the factor's derivative \eqn{L_k} built from the trigonometric tables of
#' the angles: differentiating an entry of \eqn{L} replaces one sine or cosine
#' factor by its own derivative in the free value.
#'
#' The diagonal of every component is **exactly zero**, the diagonal of \eqn{R}
#' being the constant 1. Angle \eqn{\theta_{ij}} belongs to row \eqn{i}, so
#' \eqn{\partial_k L} is supported on that row alone, though
#' \eqn{\partial_k R} is not.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `s@n_free` symmetric matrices named by `s@free_names`, each
#'   `s@dimension` by `s@dimension` with a zero diagonal.
#' @seealso [corr_derivative()], which assembles it, [corr_tables()] for the
#'   trigonometric derivatives, and [param_d2.CorrelationParam()] for the order
#'   above.
#' @keywords internal
S7::method(param_d1, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 1L)
}

#' @title Second Derivatives of a Correlation Parameter
#' @name param_d2.CorrelationParam
#' @description
#' Closed form, by the Leibniz rule on \eqn{R = LL^\top} taken twice,
#'
#' \deqn{\partial_{kl} R = L_{kl} L^\top + L_k L_l^\top + L_l L_k^\top
#'                       + L L_{kl}^\top,}
#'
#' with each \eqn{\partial^S L} from [corr_dfactor()] and zero wherever the
#' multiset \eqn{S} spans two rows of \eqn{L}.
#'
#' What survives is worth knowing before reading a result. A component in two
#' angles from **different** rows \eqn{i} and \eqn{j} does not vanish: the two
#' outer terms drop, but \eqn{L_k L_l^\top + L_l L_k^\top} is supported on
#' exactly the entries \eqn{(i, j)} and \eqn{(j, i)}. It is zero even there when
#' the two angles sit beyond the columns the two rows share. The diagonal is
#' exactly zero at every order.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
#'   `param_tuple_names(s)` and in that order, each with a zero diagonal.
#' @seealso [corr_derivative()], which assembles it, [corr_dfactor()] for the
#'   vanishing rules, and [param_d1.CorrelationParam()] for the order below.
#' @keywords internal
S7::method(param_d2, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 2L)
}

#' @title Third Derivatives of a Correlation Parameter
#' @name param_d3.CorrelationParam
#' @description
#' Closed form, the same Leibniz rule on \eqn{R = LL^\top} with the three
#' differentiations distributed over the two factors. Each \eqn{\partial^S L}
#' comes from [corr_dfactor()], which returns `NULL` whenever \eqn{S} spans two
#' rows of \eqn{L} or reaches past the columns an entry involves, so most terms
#' of the sum are skipped instead of computed and discarded.
#'
#' The angles reach the free scale through a bounded link, so the chain to third
#' order is [compose4()]'s and the accuracy is the link's; nothing is differenced.
#' The diagonal is exactly zero.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
#'   `param_tuple_names(s, 3)` and in that order, each with a zero diagonal.
#' @seealso [corr_derivative()], which assembles it, and
#'   [param_d4.CorrelationParam()] for the order above.
#' @keywords internal
S7::method(param_d3, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 3L)
}

#' @title Fourth Derivatives of a Correlation Parameter
#' @name param_d4.CorrelationParam
#' @description
#' Closed form, the Leibniz rule with four differentiations distributed over the
#' two factors of \eqn{R = LL^\top}. This is the order at which a numerical route
#' is least usable, keeping about five digits, and the spherical construction
#' pays nothing for it: every factor derivative is a product of trigonometric
#' tables that [corr_tables()] has already built.
#'
#' [corr_dfactor()]'s vanishing rules do most of the work. A quadruple spanning
#' three rows of \eqn{L} contributes nothing at all, since a Leibniz term splits
#' the four indices between two factors and each factor must stay within one row.
#' The diagonal is exactly zero.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
#'   `param_tuple_names(s, 4)` and in that order, each with a zero diagonal.
#' @seealso [corr_derivative()], which assembles it,
#'   [param_d3.CorrelationParam()] for the order below, and [numerical_d4()] for
#'   the alternative.
#' @keywords internal
S7::method(param_d4, CorrelationParam) <- function(s, eta, ...) {
  corr_derivative(s, eta, 4L)
}


#' Log-Determinant Chains of a Correlation Parameter
#'
#' @description
#' Returns, for each free value, the four derivatives of \eqn{2\log\sin\theta} in
#' that free value: the log-determinant's whole contribution from one angle. The
#' family's four log-determinant derivative methods read nothing else.
#'
#' @details
#' The factor is triangular, so \eqn{\lvert R \rvert = \prod_i L_{ii}^2} and
#' \deqn{\log\lvert R \rvert = 2 \sum_{i,k} \log \sin\theta_{ik},}
#' a sum with one term per free value. The log-determinant is therefore
#' separable, every mixed derivative is exactly zero, and each pure one is the
#' logarithm composed with the sine table [corr_tables()] already
#' holds.
#'
#' The four coefficients \eqn{2(-1)^{j-1}(j-1)!/\sin^j\theta} are the
#' derivatives of \eqn{2\log u} at \eqn{u = \sin\theta}, and [compose4()] chains
#' them onto the sine's own derivatives in the free value, which
#' [corr_tables()] has already computed.
#'
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#'
#' @return A list of `s@n_free` elements, each a list of four numbers: the first
#'   to fourth derivative of \eqn{2\log\sin\theta_k} in \eqn{\eta_k}.
#'
#' @seealso [corr_tables()] for the sine table, [compose4()] for the chain, and
#'   [corr_logdet_derivative()], the only caller.
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
#' Closed form. The factor is triangular with \eqn{L_{ii} = \prod_{k<i}
#' \sin\theta_{ik}} on its diagonal, so \eqn{|R| = \prod_i L_{ii}^2} and
#'
#' \deqn{\log|R| = 2 \sum_{i,k} \log \sin\theta_{ik},}
#'
#' one term per free value. It is `2 * sum(log(sines))`: no factorization, no
#' determinant, no eigendecomposition. Measured against the eigenvalues of the
#' assembled matrix at \eqn{p = 3}, the two agree to the printed digit.
#'
#' It is always negative or zero, a correlation matrix having determinant at most
#' 1, with 0 reached only at the identity.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A single number, `0` when `s@n_free` is 0.
#' @seealso [param_dlogdet.CorrelationParam()] for its four derivative orders,
#'   and [corr_logdet_chains()] for the per-angle chains they read.
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
#' chains of [corr_logdet_chains()]. Every mixed component is exactly zero, the
#' log-determinant being a sum with one term per free value, and each pure one is
#' that angle's chain at the matching order. The four log-determinant derivative
#' methods of the family are one call each to this function.
#'
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`.
#' @param order The derivative order: 1, 2, 3 or 4.
#'
#' @return A numeric vector of `choose(s@n_free + order - 1, order)` entries,
#'   keyed as `param_tuple_names(s, order)` and in that order, or an empty named
#'   vector when there are no free values.
#'
#' @seealso [corr_logdet_chains()] for the chains, and
#'   [param_dlogdet.CorrelationParam()], which documents all four orders.
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
#' Closed form at all four orders. One page covers
#' [param_dlogdet()], [param_d2logdet()], [param_d3logdet()] and
#' [param_d4logdet()] for this family because the four differ only in the order
#' of one chain: \eqn{\log|R| = 2\sum_{i,k}\log\sin\theta_{ik}} is a sum with one
#' term per free value, so it is **separable**, every mixed component is exactly
#' zero at every order, and each pure one is the matching derivative of
#' \eqn{2\log\sin\theta_k} in \eqn{\eta_k}.
#'
#' [corr_logdet_chains()] computes those four derivatives per angle, and
#' [corr_logdet_derivative()] places them.
#' @details
#' The four methods return vectors of different lengths, keyed by the tuple names
#' of their own order:
#'
#' - [param_dlogdet()]: `s@n_free` entries, keyed by `free_names`.
#' - [param_d2logdet()]: `choose(d + 1, 2)` entries, `param_tuple_names(s, 2)`.
#' - [param_d3logdet()]: `choose(d + 2, 3)` entries, `param_tuple_names(s, 3)`.
#' - [param_d4logdet()]: `choose(d + 3, 4)` entries, `param_tuple_names(s, 4)`.
#'
#' In each, only the components whose indices are all equal can be non-zero. At
#' \eqn{p = 3} and \eqn{\eta = (0.4, -0.2, 0.6)} the second order is
#' \eqn{(-1.161, -1.215, -1.078)} on the three diagonal pairs and exactly 0 on
#' the three off-diagonal ones.
#'
#' Unlike [log_cholesky()], where the log-determinant is linear and every order
#' above the first vanishes, this family gives non-zero answers at all four, so
#' it is one of the families where a check of those orders has content.
#' @param s A [CorrelationParam()] object.
#' @param eta A numeric vector of free values, of length `s@n_free`, already
#'   checked by the generic.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return A named numeric vector, of the length and keying its own order calls
#'   for; see **Details**.
#' @seealso [param_logdet.CorrelationParam()] for the quantity differentiated,
#'   [corr_logdet_chains()] for the per-angle chains, and
#'   [param_d3logdet.LogCholeskyParam()] for a family where these orders vanish.
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
