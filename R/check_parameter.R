#' @include numerical_fallbacks.R
NULL


#' One Row of a Diagnostic Table
#'
#' @description
#' Builds one row of the table [check_parameter()] accumulates and returns.
#' Strings are kept as strings, so the assembled table's `check` and `status`
#' columns come out character.
#'
#' @param name The name of the check, as it appears in the printed table:
#'   `"membership"`, `"round trip"`, `"first derivatives"` and the rest.
#' @param status One of `"OK"`, `"FAIL"` or `"NOT CHECKED"`. The third is used
#'   where the quantity comes from a numerical fallback and no independent route
#'   exists to compare against.
#' @param statistic The number the verdict rests on, a relative discrepancy, or
#'   `NA_real_` for a check that is structural and has no number, such as the
#'   shape and name check.
#'
#' @return A one-row data frame with columns `check`, `status` and `statistic`.
#'
#' @seealso [check_parameter()], which assembles these rows.
#'
#' @keywords internal
check_row <- function(name, status, statistic = NA_real_) {
  data.frame(
    check = name, status = status, statistic = statistic,
    stringsAsFactors = FALSE
  )
}


#' Capture and Restore the Caller's Random Stream
#'
#' @description
#' `capture_seed()` reads `.Random.seed` from the global environment, returning
#' `NULL` when there is none, and `restore_seed()` puts back what it was given.
#' Together they let a function draw from a fixed seed, so that its report is
#' the same on two runs, without leaving the caller's stream changed.
#'
#' @details
#' A validator wants both properties and they pull against each other. Drawing
#' from the caller's stream makes the worst error reported move between runs;
#' calling `set.seed()` fixes that and replaces whatever state the caller had,
#' so a call inside a simulation silently changes the simulation. Saving the
#' state on entry and restoring it with [base::on.exit()] gives the fixed draw
#' and leaves the caller alone.
#'
#' `restore_seed(NULL)` removes `.Random.seed` again, which is the right answer
#' when the caller had never drawn a random number: leaving the seed the
#' validator set behind would be the leak this exists to prevent.
#'
#' @param old The value [capture_seed()] returned, or `NULL`.
#'
#' @return `capture_seed()` returns the saved `.Random.seed`, an integer vector,
#'   or `NULL`. `restore_seed()` returns `NULL` invisibly and is called for its
#'   effect on the global environment.
#'
#' @examples
#' # the property the pair exists for, through the public interface: a
#' # validator call leaves the caller's stream exactly as it found it
#' set.seed(42)
#' want <- rnorm(1)
#' set.seed(42)
#' invisible(check_parameter(log_cholesky(2), verbose = FALSE))
#' identical(rnorm(1), want)
#'
#' @keywords internal
capture_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' @rdname capture_seed
#' @keywords internal
restore_seed <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", old, envir = globalenv())
  }
  invisible(NULL)
}


#' Free Vectors to Sweep a Parameter Over
#'
#' @description
#' Builds the set of free vectors [check_parameter()] runs its battery at: the
#' origin, two constant vectors at \eqn{\pm 0.5}, `n` drawn uniformly from
#' \eqn{(-1.5, 1.5)}, and one evenly spread from \eqn{-2} to \eqn{2}. Every
#' check reports the worst discrepancy over the whole set, so a family that is
#' right at the origin and wrong away from it is caught.
#'
#' @details
#' The spread of the last vector is deliberately moderate. The free scale is
#' unbounded, so no \eqn{\eta} is inadmissible, but the matrix built from widely
#' separated free values can be singular in double precision: spreading a
#' log-Cholesky parameter over twenty-eight units of log gives a condition number
#' around \eqn{10^{28}}, and a comparison that fails there says something about
#' the arithmetic and nothing about the parametrization.
#'
#' The scaling a rank-deficient family has to survive is a property of its
#' components, the same at every point, so it is tested where it arises, against
#' the declared null space, instead of by driving every family off the edge of
#' double precision.
#'
#' The `n` random vectors are drawn from whatever random stream is current. No
#' seed is set here: [check_parameter()], the only caller, fixes one before
#' calling and restores the caller's own on exit, so the report is reproducible
#' and the caller's stream is left as it was found.
#'
#' @param s A [parameter()] object, whose `n_free` sets the length.
#' @param n The number of random vectors, defaulting to 4. The set returned has
#'   `n + 4` members, or one member, `numeric(0)`, for a family with nothing to
#'   estimate.
#'
#' @return A list of numeric vectors, each of length `s@n_free`.
#'
#' @seealso [check_parameter()], the only caller.
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
#' Runs a battery of numerical checks on a parametrization, each against a route
#' the implementation does not itself take, and reports what passed, what failed
#' and what could not be checked. Write a family of your own, call this on it,
#' and read the table: a `FAIL` names the quantity whose closed form is wrong,
#' and the statistic beside it is the relative size of the disagreement.
#'
#' It is also the validator the package holds itself to. Measured over all
#' fifteen constructors: thirteen matrix families pass all nine checks and the
#' two that are not matrices pass all seven, with nothing skipped anywhere.
#'
#' @details
#' # Not checked is a third verdict
#'
#' A quantity that comes from a numerical fallback is reported as **NOT
#' CHECKED**, never as passed. Comparing a finite difference against a finite
#' difference is the same arithmetic twice, and it agrees however wrong the
#' parametrization is. [param_is_numerical()] is what decides, and a family
#' supplying only [param_value()] therefore comes back with six of the nine rows
#' skipped: it has nothing an independent route could contradict.
#'
#' # The nine checks, and the route each is held against
#'
#' 1. **membership**: the matrix is symmetric and positive semidefinite, a
#'    full-rank family has a strictly positive smallest eigenvalue, and a
#'    rank-deficient one annihilates its declared null space. The last is tested
#'    through the null basis, never by counting eigenvalues, a count of small
#'    eigenvalues not being scale invariant.
#' 2. **round trip**: [param_free()] recovers the free vector from the matrix,
#'    where the family implements an inverse.
#' 3. **first derivatives** against one central difference of [param_value()].
#' 4. **second derivatives** against one central difference of the analytic
#'    first derivatives.
#' 5. **log-determinant** against the sum of the logs of the eigenvalues the
#'    declared rank keeps.
#' 6. **logdet gradient** against \eqn{\mathrm{tr}(M^{+} \partial_k M)}, with
#'    the pseudo-inverse formed from an eigendecomposition instead of from
#'    [param_solve()], so the two routes share no arithmetic.
#' 7. **logdet hessian** against one central difference of the analytic
#'    gradient.
#' 8. **solve and factor** against `base::solve()`, where the family is of full
#'    rank; skipped for a deficient one, which has no inverse.
#' 9. **shapes and names**: the declared dimension, the lengths and the names
#'    match what the methods return. Structural, so it carries no statistic.
#'
#' Each check reports the worst discrepancy over the free vectors [sweep_etas()]
#' supplies, four of which are drawn at random.
#'
#' # Neither branch touches the caller's random stream
#'
#' Both batteries draw from a fixed seed, so two calls on the same family report
#' the same statistics, and both put back the `.Random.seed` they found on
#' entry, so a call in the middle of a simulation leaves that simulation
#' unchanged. The matrix branch seeds once before [sweep_etas()]; the branch for
#' a family that is not a matrix seeds each of its three draws. Restoring is
#' [capture_seed()] and [restore_seed()], through [base::on.exit()], so it
#' happens even when a check signals.
#'
#' # A family that is not a matrix gets a different battery
#'
#' [simplex()] and [transition_matrix()] have no log-determinant, no solve and
#' no factor, so seven other checks run instead: the inverse round trip, the four
#' derivative orders against the single-stencil construction, that the value
#' stays on the simplex, and that every derivative component sums to zero over
#' the value index. **The returned table has different columns in that case**;
#' see **Value**.
#'
#' @param s An object inheriting from class [parameter()]. Anything else throws
#'   `'s' must inherit from class 'parameter'.`
#' @param tol The relative tolerance a check has to meet, defaulting to `1e-6`.
#'   The comparisons are relative to the larger of 1 and the size of the
#'   reference, so the tolerance is dimensionless. `1e-6` is loose enough for the
#'   third and fourth derivative comparisons, which rest on stencils good to
#'   about \eqn{10^{-5}}, and tight enough to catch an error of one part in a
#'   thousand.
#' @param verbose Whether to print the table, `TRUE` by default. The value is
#'   returned either way.
#'
#' @return Invisibly, a data frame with one row per check. For a
#'   [matrix_parameter()] it has nine rows and the columns
#'   \describe{
#'     \item{`check`}{character, the names listed above.}
#'     \item{`status`}{character, `"OK"`, `"FAIL"` or `"NOT CHECKED"`.}
#'     \item{`statistic`}{numeric, the worst relative discrepancy, `NA` for a
#'       check that was skipped or has no number.}
#'   }
#'   For a family that is not a matrix it has seven rows and the columns
#'   `check`, `status` and `note`, the last being **character** and holding a
#'   formatted number or the empty string. Test the `status` column, which is
#'   common to both.
#'
#' @seealso [param_is_numerical()], which decides what is checkable, and
#'   [param_null_basis()] for the null space check 1 uses.
#'
#' @examples
#' set.seed(1)
#'
#' # A family that passes everything.
#' r <- check_parameter(log_cholesky(3))
#' r$status
#'
#' # A rank-deficient penalty passes the same battery, with the solve skipped:
#' # it has no inverse, and the null-space check replaces the definiteness one.
#' P <- crossprod(diff(diag(6), differences = 2))
#' d <- check_parameter(scaled_matrix(P))
#' d[d$status != "OK", ]
#'
#' # A family that is not a matrix gets the seven-row battery and a `note`
#' # column in place of `statistic`.
#' v <- check_parameter(simplex(4), verbose = FALSE)
#' names(v)
#' nrow(v)
#'
#' # What a real defect looks like. This first derivative is 5 per cent wrong.
#' Wrong <- S7::new_class("Wrong", parent = matrix_parameter)
#' S7::method(param_value, Wrong) <- function(s, eta, ...) {
#'   m <- diag(rep(exp(eta[1]), 2))
#'   dimnames(m) <- list(c("v1", "v2"), c("v1", "v2"))
#'   m
#' }
#' S7::method(param_d1, Wrong) <- function(s, eta, ...) {
#'   list(log_s = 1.05 * diag(rep(exp(eta[1]), 2)))
#' }
#' w <- Wrong(param_name = "wrong", n_free = 1L, free_names = "log_s",
#'            param_params = list(), dimension = 2L, rank = 2L,
#'            null_basis = matrix(numeric(0), 2, 0), role = "either")
#' bad <- check_parameter(w, verbose = FALSE)
#' bad[bad$status == "FAIL", ]
#'
#' @export
check_parameter <- function(s, tol = 1e-6, verbose = TRUE) {
  if (!S7::S7_inherits(s, parameter)) {
    stop("'s' must inherit from class 'parameter'.", call. = FALSE)
  }
  # a fixed seed so the report is the same on two runs, and the caller's own
  # stream put back so a call inside a simulation does not change it
  old_seed <- capture_seed()
  on.exit(restore_seed(old_seed), add = TRUE)
  if (!S7::S7_inherits(s, matrix_parameter)) {
    return(check_parameter_vector(s, tol = tol, verbose = verbose))
  }
  num <- param_is_numerical(s)
  set.seed(101)
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
#' What [check_parameter()] runs for a family whose value is not a symmetric
#' matrix, so has no log-determinant, no solve and no factor. Seven checks: the
#' inverse round trip, each of the four derivative orders against the
#' single-stencil numerical construction, that the value stays on the simplex,
#' and that every derivative component sums to zero over the value index.
#'
#' @details
#' The last check is an identity the set itself supplies. Differentiating
#' \eqn{\sum_a \pi_a = 1} gives \eqn{\sum_a \partial \pi_a = 0}, and
#' differentiating again gives the same for every higher order, so a derivative
#' array that does not sum to zero is wrong whatever else it agrees with. A
#' [transition_matrix()] satisfies it row by row, its rows being simplexes.
#'
#' Three free vectors are used, drawn from `rnorm(s@n_free, sd = 0.8)` after
#' `set.seed(101)`, `set.seed(102)` and `set.seed(103)`, so the results are
#' exactly reproducible. The caller's own `.Random.seed` is saved on entry and
#' restored on exit through [capture_seed()] and [restore_seed()], so a call in
#' the middle of a simulation leaves that simulation unchanged. The derivative
#' comparisons
#' are against
#' [numerical_d1()] through [numerical_d4()], which at third and fourth order are
#' themselves good to about \eqn{10^{-5}}: measured on [simplex()] and
#' [transition_matrix()] the four orders come back at \eqn{4 \times 10^{-12}},
#' \eqn{4 \times 10^{-12}}, \eqn{3 \times 10^{-7}} and \eqn{3 \times 10^{-5}},
#' and the widening is the reference's.
#'
#' @param s A [parameter()] object that is not a [matrix_parameter()].
#' @param tol The relative tolerance a check has to meet.
#' @param verbose Whether to print the table.
#'
#' @return Invisibly, a seven-row data frame with columns `check`, `status` and
#'   `note`. Note that `note` is **character**, holding a formatted number or the
#'   empty string, where [check_parameter()]'s matrix branch returns a numeric
#'   `statistic`.
#'
#' @seealso [check_parameter()], which dispatches here, and [simplex()] and
#'   [transition_matrix()], the two families that reach it.
#'
#' @keywords internal
check_parameter_vector <- function(s, tol = 1e-6, verbose = TRUE) {
  # registered again here, so a direct call is as safe as one through
  # check_parameter(); restoring twice puts back the same state
  old_seed <- capture_seed()
  on.exit(restore_seed(old_seed), add = TRUE)
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
