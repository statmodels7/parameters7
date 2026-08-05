# Autoregressions of general order: the partial-autocorrelation chart, the
# Levinson-Durbin recursion carried in jets, and the banded precision.

test_that("order one reproduces the closed-form ar1 family exactly", {
  # Two independent implementations of the same object -- ar1() writes its
  # formulas out, autoregressive() runs the recursion -- so agreement to
  # machine precision checks both at once, which no finite difference can.
  s <- autoregressive(6, order = 1)
  a <- ar1(6)
  for (eta in list(c(log(2), atanh(0.6)), c(-0.4, atanh(-0.8)), c(1.1, 0))) {
    expect_equal(param_value(s, eta), param_value(a, eta), tolerance = 1e-12)
    expect_equal(param_logdet(s, eta), param_logdet(a, eta), tolerance = 1e-12)
    expect_equal(param_solve(s, eta), param_solve(a, eta), tolerance = 1e-12)
    for (o in 1:4) {
      f <- switch(o, param_d1, param_d2, param_d3, param_d4)
      g <- switch(o, param_dlogdet, param_d2logdet, param_d3logdet,
                  param_d4logdet)
      expect_equal(unname(f(s, eta)), unname(f(a, eta)), tolerance = 1e-10)
      expect_equal(unname(g(s, eta)), unname(g(a, eta)), tolerance = 1e-10)
    }
  }
})

test_that("every order matches one stencil on the map", {
  for (q in 2:3) {
    s <- autoregressive(7, order = q)
    set.seed(3 + q)
    eta <- rnorm(s@n_free, sd = 0.5)
    tols <- c(1e-6, 1e-5, 1e-4, 5e-3)
    for (o in 1:4) {
      a <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
      b <- switch(o, numerical_d1, numerical_d2, numerical_d3,
                  numerical_d4)(s, eta)
      expect_identical(names(a), param_tuple_names(s, o))
      for (k in seq_along(a)) {
        expect_lt(max(abs(a[[k]] - b[[k]])) / max(1, max(abs(b[[k]]))), tols[o],
          label = sprintf("q = %d order %d component %s", q, o, names(a)[k]))
      }
    }
  }
})

test_that("the chart lands inside the stationary region", {
  # The reason the family is parametrized by partial autocorrelations at all:
  # the region is not a box, so this is a claim about the chart rather than
  # about a bound. The roots of the autoregressive polynomial must lie
  # outside the unit circle for every free vector.
  set.seed(21)
  for (q in 2:4) {
    s <- autoregressive(q + 4L, order = q)
    for (i in 1:20) {
      eta <- c(rnorm(1), rnorm(q, sd = 2.5))
      phi <- vapply(parameters7:::ar_jets(s, eta)$phi, function(x) x$v,
                    numeric(1))
      expect_gt(min(Mod(polyroot(c(1, -phi)))), 1)
      m <- param_value(s, eta)
      expect_gt(min(eigen(m, symmetric = TRUE, only.values = TRUE)$values), 0)
    }
  }
})

test_that("the log-determinant is closed and separable", {
  for (q in 1:3) {
    s <- autoregressive(q + 4L, order = q)
    set.seed(q)
    eta <- rnorm(s@n_free, sd = 0.6)
    m <- param_value(s, eta)
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    expect_equal(param_logdet(s, eta), sum(log(ev)), tolerance = 1e-10)

    mi <- solve(m)
    want <- vapply(param_d1(s, eta), function(a) sum(mi * a), numeric(1))
    expect_equal(unname(param_dlogdet(s, eta)), unname(want), tolerance = 1e-8)

    # one term per free value, so every mixed component vanishes
    for (o in 2:4) {
      d <- switch(o - 1L, param_d2logdet, param_d3logdet,
                  param_d4logdet)(s, eta)
      idx <- param_tuple_indices(s, o)
      for (i in seq_along(idx)) {
        if (any(idx[[i]] != idx[[i]][1L])) expect_identical(d[[i]], 0)
      }
    }
  }
})

test_that("the precision is exact and banded of the order's width", {
  for (q in 1:3) {
    s <- autoregressive(8, order = q)
    set.seed(10 + q)
    eta <- rnorm(s@n_free, sd = 0.5)
    pr <- param_solve(s, eta)
    expect_equal(pr, solve(param_value(s, eta)),
      tolerance = 1e-10, ignore_attr = TRUE)
    # an order-q Markov process has no partial correlation beyond the lag,
    # and the entries are exactly zero rather than small
    expect_true(all(pr[abs(row(pr) - col(pr)) > q] == 0))
  }
})

test_that("the inverse map is exact and refuses what is not in the set", {
  s <- autoregressive(7, order = 2)
  set.seed(5)
  eta <- rnorm(3, sd = 0.5)
  expect_equal(unname(param_free(s, param_value(s, eta))), eta,
    tolerance = 1e-10)

  # a matrix that is not Toeplitz is not the covariance of a stationary
  # process, and is refused rather than averaged along its diagonals
  m <- param_value(s, eta)
  m[1L, 2L] <- m[2L, 1L] <- m[1L, 2L] * 1.2
  expect_error(param_free(s, m), "not Toeplitz")

  # an AR(3) covariance is Toeplitz but does not follow an order-two
  # recursion beyond lag two
  m3 <- param_value(autoregressive(7, order = 3), c(0, 0.6, -0.4, 0.5))
  expect_error(param_free(s, m3), "Yule-Walker")

  expect_error(param_free(s, diag(7)), NA)  # white noise IS an AR(2) with r = 0
})

test_that("the dimension must exceed the order", {
  expect_error(autoregressive(3, order = 3), "must exceed")
  expect_error(autoregressive(5, order = 0), "positive integer")
  expect_error(autoregressive(5, order = 2.5), "positive integer")
  expect_equal(autoregressive(4, order = 3)@n_free, 4L)
})

test_that("the validator passes over orders and dimensions", {
  for (q in 1:3) {
    for (p in (q + 1L):(q + 3L)) {
      s <- autoregressive(p, order = q)
      res <- check_parameter(s, verbose = FALSE)
      expect_true(all(res$status == "OK"),
        label = sprintf("q = %d, p = %d: %s", q, p,
                        paste(res$check[res$status != "OK"], collapse = ", ")))
    }
  }
})

test_that("jet arithmetic reproduces a polynomial's derivatives", {
  # The machinery the whole family rests on now lives in numericals7 and is
  # consumed through its exported seam and the operators; this checks the
  # integration against a polynomial whose derivatives can be written down:
  # f(x, y) = x^3 y^2.
  lay <- numericals7::jet_layout(2L)
  x <- 1.3
  y <- -0.7
  jx <- numericals7::jet_var(1L, list(x, 1, 0, 0, 0), lay)
  jy <- numericals7::jet_var(2L, list(y, 1, 0, 0, 0), lay)
  f <- jx * jx * jx * (jy * jy)

  at <- function(tup) {
    o <- length(tup)
    i <- which(vapply(numericals7::tuple_indices(2L, o), function(t) {
      identical(sort(t), sort(as.integer(tup)))
    }, logical(1)))
    f$d[[o]][i]
  }
  expect_equal(f$v, x^3 * y^2)
  expect_equal(at(1L), 3 * x^2 * y^2)
  expect_equal(at(2L), 2 * x^3 * y)
  expect_equal(at(c(1L, 1L)), 6 * x * y^2)
  expect_equal(at(c(1L, 2L)), 6 * x^2 * y)
  expect_equal(at(c(2L, 2L)), 2 * x^3)
  expect_equal(at(c(1L, 1L, 1L)), 6 * y^2)
  expect_equal(at(c(1L, 1L, 2L)), 12 * x * y)
  expect_equal(at(c(1L, 2L, 2L)), 6 * x^2)
  expect_equal(at(c(2L, 2L, 2L)), 0)
  expect_equal(at(c(1L, 1L, 1L, 1L)), 0)
  expect_equal(at(c(1L, 1L, 1L, 2L)), 12 * y)
  expect_equal(at(c(1L, 1L, 2L, 2L)), 12 * x)
  expect_equal(at(c(1L, 2L, 2L, 2L)), 0)
})
