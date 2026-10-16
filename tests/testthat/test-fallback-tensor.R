# The numerical fallbacks on numericals7::fd_tensor(): against the R
# implementations of 0.24.0 (helper-fallback-twins.R) and against the closed
# forms.

normgap <- function(a, b) {
  a <- unlist(a)
  b <- unlist(b)
  sqrt(sum((a - b)^2)) / max(sqrt(sum(b^2)), 1e-300)
}

test_that("at accuracy 2 the value routes reproduce the R implementations", {
  eta <- c(0.3, -0.2, 0.25)
  s <- toy_value()
  expect_identical(numerical_d1(s, eta, accuracy = 2L), twin_numerical_d1(s, eta))
  expect_identical(numerical_d2(s, eta, accuracy = 2L), twin_numerical_d2(s, eta))
  expect_identical(numerical_d3(s, eta, accuracy = 2L), twin_numerical_d3(s, eta))
  # the five-point fourth-derivative weights, which fd_weights() solves to
  # about 1e-15 and Fornberg's recursion gives exactly; the gap is the
  # rounding of the stencil
  expect_lt(normgap(numerical_d4(s, eta, accuracy = 2L),
                    twin_numerical_d4(s, eta)), 1e-4)
  # the analytic first-derivative route of the second order
  d <- toy_d12()
  expect_identical(numerical_d2(d, eta, accuracy = 2L), twin_numerical_d2(d, eta))
})

test_that("at accuracy 2 the first two orders of every structure are unchanged", {
  set.seed(7)
  st <- shipped_structures()
  for (nm in names(st)) {
    s <- st[[nm]]
    e <- stats::rnorm(s@n_free, sd = 0.6)
    expect_identical(numerical_d1(s, e, accuracy = 2L), twin_numerical_d1(s, e),
                     label = nm)
    # the values are identical; the 0.24.0 route on an analytic first
    # derivative took that derivative's names, which simplex() leaves out,
    # and the result now carries the names of param_value()
    strip <- function(l) lapply(l, function(m) { names(m) <- NULL; m })
    expect_identical(strip(numerical_d2(s, e, accuracy = 2L)),
                     strip(twin_numerical_d2(s, e)), label = nm)
  }
})

test_that("every order meets the closed forms of every structure", {
  # measured 2026-10-16 (numericals7 0.23.0): at most 1.0e-12, 1.9e-12,
  # 5.0e-12 and 9.1e-12 at orders one to four
  set.seed(7)
  st <- shipped_structures()
  for (nm in names(st)) {
    s <- st[[nm]]
    e <- stats::rnorm(s@n_free, sd = 0.6)
    expect_lt(normgap(numerical_d1(s, e), param_d1(s, e)), 1e-10, label = nm)
    expect_lt(normgap(numerical_d2(s, e), param_d2(s, e)), 1e-10, label = nm)
    expect_lt(normgap(numerical_d3(s, e), param_d3(s, e)), 1e-10, label = nm)
    expect_lt(normgap(numerical_d4(s, e), param_d4(s, e)), 1e-10, label = nm)
  }
})

test_that("the third and fourth orders difference the analytic order below", {
  eta <- c(0.3, -0.2, 0.25)
  d <- toy_d12()
  v <- toy_value()
  # d3 of the toys: only a_a_a and b_b_b are nonzero
  exact <- numerical_d3(v, eta)
  exact[] <- list(matrix(0, 2L, 2L))
  exact[["a:a:a"]] <- matrix(c(exp(eta[1]), 0, 0, 0), 2L, 2L)
  exact[["b:b:b"]] <- matrix(c(0, 0, 0, exp(eta[2])), 2L, 2L)
  e_d <- normgap(numerical_d3(d, eta), exact)
  e_v <- normgap(numerical_d3(v, eta), exact)
  # from the analytic second order, one first difference; from the value
  # alone, a third difference
  expect_lt(e_d, 1e-11)
  expect_lt(e_d, e_v)
  expect_lt(e_v, 1e-7)
})

test_that("a composite over a value-only parameter differences its own value", {
  # block_diag() registers its own methods, so param_is_numerical() reads
  # them as analytic; they are built from the block's fallbacks
  bd <- block_diag(toy_value(), ar1(3))
  expect_false(any(param_is_numerical(bd)))
  a <- param_analytic(bd)
  expect_false(any(a[c("param_d1", "param_d2", "param_d3", "param_d4")]))
  expect_true(all(param_analytic(block_diag(log_cholesky(2), ar1(3)))))
  set.seed(9)
  e <- stats::rnorm(bd@n_free, sd = 0.5)
  expect_lt(normgap(numerical_d3(bd, e), param_d3(bd, e)), 1e-6)
})

test_that("the log-determinant orders three and four", {
  eta <- c(0.3, -0.2, 0.25)
  d <- toy_d12()
  # at accuracy 2 the third order is the 0.24.0 route; the fourth takes a
  # different pair of the second-order block for tuples with three distinct
  # components, and agrees to the stencil's accuracy
  expect_identical(param_tensor(d, eta, 3L, "logdet", 2L), twin_d3logdet(d, eta))
  expect_lt(normgap(param_tensor(d, eta, 4L, "logdet", 2L), twin_d4logdet(d, eta)),
            1e-6)
  # against the exact derivatives of log(exp(a + b) - c^2): measured 9.3e-12
  # and 1.7e-9 at the default accuracy, 9.8e-11 and 3.8e-8 on the 0.24.0 route
  expect_lt(normgap(param_d3logdet(d, eta), toy_logdet_exact(eta, 3)), 1e-10)
  expect_lt(normgap(param_d4logdet(d, eta), toy_logdet_exact(eta, 4)), 1e-8)
  expect_identical(names(param_d3logdet(d, eta)), param_tuple_names(d, 3L))
  expect_identical(names(param_d4logdet(d, eta)), param_tuple_names(d, 4L))
  # without analytic arrays the stencil is applied to param_logdet(), which
  # is exact: measured 9e-9 to 3.5e-7 at order three and 3.5e-7 to 7.9e-6 at
  # order four over four points
  v <- toy_value()
  expect_lt(normgap(param_d3logdet(v, eta), toy_logdet_exact(eta, 3)), 1e-5)
  expect_lt(normgap(param_d4logdet(v, eta), toy_logdet_exact(eta, 4)), 1e-4)
})

test_that("the shapes, names and symmetry of the results", {
  eta <- c(0.3, -0.2, 0.25)
  s <- toy_value()
  for (o in 1:4) {
    got <- list(numerical_d1, numerical_d2, numerical_d3, numerical_d4)[[o]](s, eta)
    expect_identical(names(got), param_tuple_names(s, o))
    for (m in got) {
      expect_identical(dim(m), c(2L, 2L))
      expect_identical(m, t(m))
    }
  }
  sx <- simplex(4)
  g <- numerical_d3(sx, c(0.1, -0.3, 0.2))
  expect_identical(unname(lengths(g)),
                   rep(length(param_value(sx, c(0.1, -0.3, 0.2))), length(g)))
  expect_identical(names(g[[1L]]), names(param_value(sx, c(0.1, -0.3, 0.2))))
})
