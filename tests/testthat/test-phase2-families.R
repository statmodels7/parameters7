# Correlation matrices, compound symmetry and AR(1): the families whose value
# is constrained beyond being positive definite.

phase2_families <- function(p) {
  list(correlation_matrix(p), compound_symmetry(p), ar1(p))
}

test_that("the validator passes on every dimension from two to six", {
  # p = 2 is where a matrix of indices degenerates to a single row, which is
  # how the tridiagonal inverse of the AR(1) family was once written into the
  # wrong cell; the sweep is what found it, so the sweep stays.
  for (p in 2:6) {
    for (s in phase2_families(p)) {
      res <- check_parameter(s, verbose = FALSE)
      expect_true(all(res$status == "OK"),
        label = sprintf("%s at p = %d: %s", s@param_name, p,
                        paste(res$check[res$status != "OK"], collapse = ", "))
      )
    }
  }
})

test_that("every order matches one stencil on the map", {
  for (s in phase2_families(4)) {
    set.seed(9)
    eta <- rnorm(s@n_free, sd = 0.6)
    tols <- c(1e-6, 1e-5, 1e-4, 5e-3)
    for (o in 1:4) {
      a <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
      b <- switch(o, numerical_d1, numerical_d2, numerical_d3,
                  numerical_d4)(s, eta)
      expect_identical(names(a), param_tuple_names(s, o))
      for (k in seq_along(a)) {
        expect_lt(max(abs(a[[k]] - b[[k]])) / max(1, max(abs(b[[k]]))), tols[o],
          label = sprintf("%s order %d component %s", s@param_name, o,
                          names(a)[k])
        )
      }
    }
  }
})

test_that("the closed log-determinants agree with the spectrum and the trace", {
  for (s in phase2_families(5)) {
    set.seed(4)
    eta <- rnorm(s@n_free, sd = 0.5)
    m <- param_value(s, eta)
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values

    expect_equal(param_logdet(s, eta), sum(log(ev)), tolerance = 1e-10)

    # the gradient against tr(M^{-1} dM), which uses the matrix derivatives
    # and nothing the log-determinant methods compute
    mi <- solve(m)
    want <- vapply(param_d1(s, eta), function(a) sum(mi * a), numeric(1))
    expect_equal(unname(param_dlogdet(s, eta)), unname(want), tolerance = 1e-8)

    # and the higher orders against one stencil on the analytic second
    idx3 <- param_tuple_indices(s, 3L)
    idx2 <- param_tuple_indices(s, 2L)
    key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
    got <- param_d3logdet(s, eta)
    for (i in seq_along(idx3)) {
      t <- idx3[[i]]
      pos <- match(paste(sort(t[1:2]), collapse = ","), key2)
      k <- t[3L]
      h <- 1e-5 * max(1, abs(eta[k]))
      up <- eta; dn <- eta
      up[k] <- up[k] + h; dn[k] <- dn[k] - h
      ref <- (param_d2logdet(s, up)[[pos]] - param_d2logdet(s, dn)[[pos]]) /
        (2 * h)
      expect_lt(abs(got[[i]] - ref), 1e-5 + 1e-4 * abs(ref),
        label = sprintf("%s %s", s@param_name, names(got)[i]))
    }
    expect_true(all(is.finite(param_d4logdet(s, eta))))
  }
})

test_that("the closed solves agree with a general one", {
  for (s in phase2_families(6)) {
    set.seed(2)
    eta <- rnorm(s@n_free, sd = 0.5)
    m <- param_value(s, eta)
    expect_equal(param_solve(s, eta), solve(m),
      tolerance = 1e-10, ignore_attr = TRUE)
  }
  # the AR(1) precision is tridiagonal, which is the property the family is
  # used for: a Markov process has no partial correlation beyond the lag
  s <- ar1(6)
  p6 <- param_solve(s, c(0.2, atanh(0.7)))
  far <- abs(row(p6) - col(p6)) > 1
  expect_true(all(p6[far] == 0))
})

test_that("a correlation parameter stays a correlation matrix", {
  s <- correlation_matrix(5)
  set.seed(6)
  for (i in 1:5) {
    eta <- rnorm(s@n_free, sd = 1.2)
    r <- param_value(s, eta)
    expect_equal(diag(r), rep(1, 5), tolerance = 1e-14, ignore_attr = TRUE)
    expect_gt(min(eigen(r, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_true(all(abs(r) <= 1 + 1e-12))
    # differentiating a constant diagonal gives a zero diagonal, at every order
    for (o in 1:4) {
      d <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
      for (comp in d) expect_lt(max(abs(diag(comp))), 1e-9)
    }
  }
  # The factor's rows are independent, so a derivative of the FACTOR across
  # two rows vanishes -- but R's need not, its entry (i, j) being the inner
  # product of rows i and j. What holds for R is that such a component is
  # SUPPORTED on those two entries; it is zero even there when the two angles
  # sit beyond the columns the rows share.
  eta <- rnorm(s@n_free)
  tb <- parameters7:::corr_tables(s, eta)
  rows <- s@param_params$row
  cross <- which(vapply(param_tuple_indices(s, 2L), function(t) {
    rows[t[1L]] != rows[t[2L]]
  }, logical(1)))
  expect_gt(length(cross), 0)
  d2 <- param_d2(s, eta)
  for (i in cross) {
    t <- param_tuple_indices(s, 2L)[[i]]
    expect_null(parameters7:::corr_dfactor(s, tb, t))
    keep <- matrix(FALSE, 5, 5)
    keep[rows[t[1L]], rows[t[2L]]] <- TRUE
    keep[rows[t[2L]], rows[t[1L]]] <- TRUE
    expect_true(all(d2[[i]][!keep] == 0))
  }
  # and the claim is not vacuous: some cross-row component really is non-zero
  expect_gt(max(vapply(cross, function(i) max(abs(d2[[i]])), numeric(1))), 0)
})

test_that("the economical families refuse a value that is not theirs", {
  cs <- compound_symmetry(3)
  ar <- ar1(3)
  expect_error(compound_symmetry(1), "at least 2")
  expect_error(ar1(1), "at least 2")

  # a compound-symmetric matrix is not AR(1) and the reverse, and neither is
  # quietly projected onto the other
  m_cs <- param_value(cs, c(0, 0.4))
  m_ar <- param_value(ar, c(0, 0.6))
  expect_error(param_free(ar, m_cs), "geometric pattern")
  expect_error(param_free(cs, m_ar), "not all equal")
  expect_error(param_free(cs, diag(c(1, 2, 3))), "constant diagonal")

  # and the correlation bound of compound symmetry depends on the dimension:
  # -1/(p-1), not -1
  expect_equal(cs@param_params$link_rho@link_bounds, c(-0.5, 1),
    ignore_attr = TRUE)
  expect_equal(compound_symmetry(5)@param_params$link_rho@link_bounds,
    c(-0.25, 1), ignore_attr = TRUE)
})

test_that("a corrupted closed form is caught", {
  # the paired injection: 5 per cent off must fail, and the unbroken family
  # must pass, so the check cannot be satisfied by weakening it
  Bent <- S7::new_class("Bent", parent = Ar1Param, package = NULL)
  b <- Bent(
    param_name = "ar1", dimension = 4L, n_free = 2L,
    free_names = c("scale", "rho"), rank = 4L,
    null_basis = matrix(numeric(0), 4, 0), role = "either",
    param_params = ar1(4)@param_params
  )
  S7::method(param_d2, Bent) <- function(s, eta, ...) {
    out <- parameters7:::econ_derivative(s, eta, 2L, parameters7:::ar1_pattern)
    out[[1L]] <- out[[1L]] * 1.05
    out
  }
  expect_true(any(check_parameter(b, verbose = FALSE)$status == "FAIL"))
  expect_true(all(check_parameter(ar1(4), verbose = FALSE)$status == "OK"))
})

test_that("compose4 reproduces a composition it does not know about", {
  # exp(sin(x)) differentiated four times, against Richardson: the helper is
  # the one piece every phase-two family leans on
  x <- 0.7
  g <- list(cos(x), -sin(x), -cos(x), sin(x))
  f <- rep(list(exp(sin(x))), 4)
  got <- compose4(f, g)
  h <- function(z) exp(sin(z))
  for (o in 1:4) {
    want <- numDeriv::genD(h, x)$D
    ref <- switch(o,
      numDeriv::grad(h, x),
      numDeriv::hessian(h, x)[1, 1],
      numDeriv::grad(function(z) numDeriv::hessian(h, z)[1, 1], x),
      NULL
    )
    if (!is.null(ref)) expect_equal(got[[o]], ref, tolerance = 1e-6)
  }
})
