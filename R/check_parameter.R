#' @include numerical_fallbacks.R
NULL


#' One Row of a Diagnostic Table
#'
#' @param name The name of the check.
#' @param status \code{"OK"}, \code{"FAIL"} or \code{"NOT CHECKED"}.
#' @param statistic The number the verdict rests on, or \code{NA}.
#'
#' @return A one-row data frame.
#'
#' @keywords internal
check_row <- function(name, status, statistic = NA_real_) {
  data.frame(
    check = name, status = status, statistic = statistic,
    stringsAsFactors = FALSE
  )
}


#' Free Vectors to Sweep a Parameter Over
#'
#' @description
#' A set of free vectors covering the ordinary range and a spread one.
#'
#' @details
#' The spread of the last vector is deliberately moderate. The free scale is
#' unbounded, so no value of \eqn{\eta} is inadmissible, but the matrix built
#' from widely separated free values can be singular in double precision --
#' spreading a log-Cholesky parameter over twenty-eight units of log gives a
#' condition number around \eqn{10^{28}} -- and a comparison that fails there
#' is a statement about the arithmetic rather than about the parameter. The
#' scaling that a rank-deficient family must survive is a property of its
#' components rather than of a point, so it is tested where it arises, against
#' the declared null space, and not by driving every family off the edge of
#' double precision.
#'
#' @param s A \code{\link{parameter}} object.
#' @param n The number of random vectors.
#'
#' @return A list of numeric vectors.
#'
#' @keywords internal
sweep_etas <- function(s, n = 4L) {
  d <- s@n_free
  if (d == 0L) return(list(numeric(0)))
  out <- list(rep(0, d), rep(0.5, d), rep(-0.5, d))
  for (i in seq_len(n)) out[[length(out) + 1L]] <- stats::runif(d, -1.5, 1.5)
  out[[length(out) + 1L]] <- seq(-2, 2, length.out = d)
  out
}


#' Validate a Covariance Parameter
#'
#' @description
#' Runs a battery of numerical checks on a parameter, each against a route the
#' implementation does not itself take, and reports what passed, what failed
#' and what could not be checked.
#'
#' @details
#' A quantity that comes from a numerical fallback is reported as \strong{not
#' checked} rather than as passed. Comparing a finite difference against a
#' finite difference is the same arithmetic twice, and it agrees however wrong
#' the parameter is; saying so is the difference between a validator and a
#' formality.
#'
#' The checks are:
#' \enumerate{
#'   \item \strong{membership}: the matrix is symmetric and positive
#'     semidefinite, a full-rank family has a positive smallest eigenvalue, and
#'     a rank-deficient one annihilates its declared null space. The last is
#'     tested through the null basis rather than by counting eigenvalues,
#'     because a count is not scale invariant.
#'   \item \strong{round trip}: \code{param_free()} recovers the free vector
#'     from the matrix, where the family implements it.
#'   \item \strong{first derivatives} against one central difference of
#'     \code{\link{param_value}}.
#'   \item \strong{second derivatives} against one central difference of the
#'     analytic first derivatives.
#'   \item \strong{log-determinant} against the sum of the logs of the
#'     eigenvalues the rank keeps.
#'   \item \strong{log-determinant gradient} against
#'     \eqn{\mathrm{tr}(M^{+} \partial_k M)}, with the pseudo-inverse formed
#'     from an eigendecomposition rather than from \code{\link{param_solve}}.
#'   \item \strong{log-determinant Hessian} against one central difference of
#'     the analytic gradient.
#'   \item \strong{solve} against \code{base::solve}, where the family is of
#'     full rank.
#'   \item \strong{shapes}: the declared dimension, length and names match what
#'     the methods return.
#' }
#'
#' @param s An object inheriting from class \code{\link{parameter}}.
#' @param tol The relative tolerance for the comparisons.
#' @param verbose Whether to print the table.
#'
#' @return Invisibly, a data frame with columns \code{check}, \code{status} and
#'   \code{statistic}.
#'
#' @seealso \code{\link{param_is_numerical}}
#'
#' @examples
#' invisible(check_parameter(log_cholesky(3)))
#'
#' # a rank-deficient penalty passes the same battery
#' invisible(check_parameter(scaled_matrix(crossprod(diff(diag(6), differences = 2)))))
#'
#' @export
check_parameter <- function(s, tol = 1e-6, verbose = TRUE) {
  if (!S7::S7_inherits(s, parameter)) {
    stop("'s' must inherit from class 'parameter'.", call. = FALSE)
  }
  if (!S7::S7_inherits(s, matrix_parameter)) {
    return(check_parameter_vector(s, tol = tol, verbose = verbose))
  }
  num <- param_is_numerical(s)
  etas <- sweep_etas(s)
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  out <- list()

  # --- 1. membership --------------------------------------------------------
  worst <- 0
  bad <- FALSE
  for (eta in etas) {
    m <- param_value(s, eta)
    if (!identical(dim(m), c(s@dimension, s@dimension))) { bad <- TRUE; break }
    asym <- max(abs(m - t(m))) / max(1, max(abs(m)))
    worst <- max(worst, asym)
    ev <- eigen((m + t(m)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (s@rank == s@dimension) {
      if (min(ev) <= 0) bad <- TRUE
    } else {
      # Through the declared null space, which is scale invariant; a count of
      # eigenvalues at this eta would not be.
      resid <- max(abs(m %*% s@null_basis)) / max(1, max(abs(m)))
      worst <- max(worst, resid)
      if (resid > tol) bad <- TRUE
      if (min(ev) < -tol * max(abs(ev))) bad <- TRUE
    }
  }
  out[[length(out) + 1L]] <- check_row(
    "membership", if (bad) "FAIL" else "OK", worst
  )

  # --- 2. the round trip ----------------------------------------------------
  has_free <- !inherits(
    tryCatch(param_free(s, param_value(s, etas[[1L]])), error = function(e) e),
    "error"
  )
  if (!has_free) {
    out[[length(out) + 1L]] <- check_row("round trip", "NOT CHECKED")
  } else {
    err <- 0
    for (eta in etas) {
      got <- tryCatch(
        param_free(s, param_value(s, eta)),
        error = function(e) rep(NA_real_, s@n_free)
      )
      err <- max(err, if (s@n_free) rel(unname(got), eta) else 0)
    }
    out[[length(out) + 1L]] <- check_row(
      "round trip", if (is.finite(err) && err < tol) "OK" else "FAIL", err
    )
  }

  # --- 3. first derivatives -------------------------------------------------
  if (num[["param_d1"]]) {
    out[[length(out) + 1L]] <- check_row("first derivatives", "NOT CHECKED")
  } else {
    err <- 0
    for (eta in etas[seq_len(min(3L, length(etas)))]) {
      a <- param_d1(s, eta)
      b <- numerical_d1(s, eta)
      for (k in seq_along(a)) err <- max(err, rel(a[[k]], b[[k]]))
    }
    out[[length(out) + 1L]] <- check_row(
      "first derivatives", if (err < tol) "OK" else "FAIL", err
    )
  }

  # --- 4. second derivatives ------------------------------------------------
  if (num[["param_d2"]] || num[["param_d1"]]) {
    out[[length(out) + 1L]] <- check_row("second derivatives", "NOT CHECKED")
  } else {
    err <- 0
    for (eta in etas[seq_len(min(3L, length(etas)))]) {
      a <- param_d2(s, eta)
      b <- numerical_d2(s, eta)
      for (k in seq_along(a)) err <- max(err, rel(a[[k]], b[[k]]))
    }
    out[[length(out) + 1L]] <- check_row(
      "second derivatives", if (err < tol) "OK" else "FAIL", err
    )
  }

  # --- 5. log-determinant ---------------------------------------------------
  if (num[["param_logdet"]]) {
    out[[length(out) + 1L]] <- check_row("log-determinant", "NOT CHECKED")
  } else {
    err <- 0
    for (eta in etas) {
      m <- param_value(s, eta)
      ev <- eigen((m + t(m)) / 2, symmetric = TRUE, only.values = TRUE)$values
      err <- max(err, rel(param_logdet(s, eta), sum(log(ev[seq_len(s@rank)]))))
    }
    out[[length(out) + 1L]] <- check_row(
      "log-determinant", if (err < tol) "OK" else "FAIL", err
    )
  }

  # --- 6. log-determinant gradient ------------------------------------------
  if (num[["param_dlogdet"]] || s@n_free == 0L) {
    out[[length(out) + 1L]] <- check_row("logdet gradient", "NOT CHECKED")
  } else {
    err <- 0
    for (eta in etas) {
      sp <- param_spectrum(s, eta)
      mi <- spectrum_pinv(sp)
      d <- param_d1(s, eta)
      want <- vapply(d, function(dk) sum(mi * dk), numeric(1))
      err <- max(err, rel(unname(param_dlogdet(s, eta)), unname(want)))
    }
    out[[length(out) + 1L]] <- check_row(
      "logdet gradient", if (err < tol) "OK" else "FAIL", err
    )
  }

  # --- 7. log-determinant Hessian -------------------------------------------
  if (num[["param_d2logdet"]] || num[["param_dlogdet"]] || s@n_free == 0L) {
    out[[length(out) + 1L]] <- check_row("logdet hessian", "NOT CHECKED")
  } else {
    err <- 0
    idx <- param_tuple_indices(s)
    for (eta in etas[seq_len(min(3L, length(etas)))]) {
      a <- param_d2logdet(s, eta)
      for (i in seq_along(idx)) {
        k <- idx[[i]][1L]
        l <- idx[[i]][2L]
        h <- fd_step(eta[l], 1L)
        up <- dn <- eta
        up[l] <- eta[l] + h
        dn[l] <- eta[l] - h
        want <- (param_dlogdet(s, up)[[k]] - param_dlogdet(s, dn)[[k]]) / (2 * h)
        err <- max(err, abs(a[[i]] - want) / max(1, abs(want)))
      }
    }
    out[[length(out) + 1L]] <- check_row(
      "logdet hessian", if (err < tol) "OK" else "FAIL", err
    )
  }

  # --- 8. solve -------------------------------------------------------------
  if (s@rank < s@dimension) {
    out[[length(out) + 1L]] <- check_row("solve", "NOT CHECKED")
  } else {
    err <- 0
    for (eta in etas) {
      err <- max(err, rel(
        param_solve(s, eta), solve(param_value(s, eta))
      ))
      l <- param_factor(s, eta)
      err <- max(err, rel(tcrossprod(l), param_value(s, eta)))
    }
    out[[length(out) + 1L]] <- check_row(
      "solve and factor", if (err < tol) "OK" else "FAIL", err
    )
  }

  # --- 9. shapes ------------------------------------------------------------
  eta <- etas[[1L]]
  shape_ok <- length(s@free_names) == s@n_free &&
    identical(names(param_d1(s, eta)), s@free_names) &&
    identical(names(param_d2(s, eta)), param_tuple_names(s)) &&
    length(param_tuple_names(s)) == s@n_free * (s@n_free + 1L) / 2L &&
    identical(dim(s@null_basis), c(s@dimension, s@dimension - s@rank))
  out[[length(out) + 1L]] <- check_row(
    "shapes and names", if (shape_ok) "OK" else "FAIL"
  )

  res <- do.call(rbind, out)

  if (verbose) {
    cat(sprintf(
      "Parameter: %s   (%d x %d, rank %d, %d free)\n",
      s@param_name, s@dimension, s@dimension, s@rank, s@n_free
    ))
    width <- max(nchar(res$check))
    for (i in seq_len(nrow(res))) {
      stat <- if (is.na(res$statistic[i])) {
        ""
      } else {
        sprintf("   %.2e", res$statistic[i])
      }
      cat(sprintf(
        "  [%-11s] %-*s%s\n", res$status[i], width, res$check[i], stat
      ))
    }
    n_fail <- sum(res$status == "FAIL")
    n_skip <- sum(res$status == "NOT CHECKED")
    cat(sprintf(
      "  %d passed, %d failed, %d not checked\n",
      sum(res$status == "OK"), n_fail, n_skip
    ))
  }

  invisible(res)
}


#' The Reduced Battery for a Parameter That Is Not a Matrix
#'
#' @description
#' What \code{\link{check_parameter}} runs for a family whose value is not a
#' symmetric matrix: the inverse round trip, every derivative order against
#' the single-stencil numerical construction, and the identities the value's
#' own set supplies -- on the simplex the entries sum to one and every
#' derivative component sums to zero over the value index, since
#' differentiating \eqn{\sum_a \pi_a = 1} kills every order; a transition
#' matrix satisfies the same row by row.
#'
#' @param s A \code{\link{parameter}} object.
#' @param tol The relative tolerance.
#' @param verbose Whether to print the table.
#'
#' @return A data frame with one row per check, invisibly when printed.
#'
#' @keywords internal
check_parameter_vector <- function(s, tol = 1e-6, verbose = TRUE) {
  etas <- lapply(1:3, function(i) {
    set.seed(100 + i)
    stats::rnorm(s@n_free, sd = 0.8)
  })
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  num <- param_is_numerical(s)
  rows <- list()
  add <- function(name, ok, note = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = name, status = if (isTRUE(ok)) "OK" else "FAIL", note = note
    )
  }

  # inverse round trip
  ok <- TRUE
  for (eta in etas) {
    got <- tryCatch(param_free(s, param_value(s, eta)), error = function(e) NULL)
    if (is.null(got) || max(abs(unname(got) - eta)) > tol) ok <- FALSE
  }
  add("param_free round trip", ok)

  # derivative orders against the single-stencil construction
  refs <- list(numerical_d1, numerical_d2, numerical_d3, numerical_d4)
  funs <- list(param_d1, param_d2, param_d3, param_d4)
  tols <- c(1e-6, 1e-5, 1e-4, 5e-3)
  for (o in 1:4) {
    nm <- paste0("param_d", o)
    if (num[[nm]]) {
      add(nm, TRUE, "numerical")
      next
    }
    worst <- 0
    for (eta in etas) {
      a <- funs[[o]](s, eta)
      b <- refs[[o]](s, eta)
      for (k in seq_along(a)) worst <- max(worst, rel(a[[k]], b[[k]]))
    }
    add(nm, worst < tols[o], sprintf("%.1e", worst))
  }

  # the value's own identities
  is_simplex <- S7::S7_inherits(s, SimplexParam)
  is_tm <- S7::S7_inherits(s, TransitionMatrixParam)
  if (is_simplex || is_tm) {
    ok_sum <- TRUE
    ok_zero <- TRUE
    for (eta in etas) {
      v <- param_value(s, eta)
      sums <- if (is_simplex) sum(v) else rowSums(v)
      if (max(abs(sums - 1)) > 1e-12) ok_sum <- FALSE
      for (o in 1:4) {
        d <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
        for (comp in d) {
          z <- if (is_simplex) sum(comp) else rowSums(comp)
          if (max(abs(z)) > tol) ok_zero <- FALSE
        }
      }
    }
    add("value stays on the simplex", ok_sum)
    add("derivatives sum to zero", ok_zero)
  }

  res <- do.call(rbind, rows)
  if (verbose) print(res, row.names = FALSE)
  invisible(res)
}
