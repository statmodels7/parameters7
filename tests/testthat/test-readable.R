# The quantities a family declares, and the Jacobian a consumer carries a
# variance matrix onto them with.

test_that("the declared Jacobian is the derivative of the declared values", {
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(s = ar1(5), eta = c(log(4), atanh(0.6))),
    list(s = compound_symmetry(4), eta = c(log(2), 0.3)),
    list(s = autoregressive(8, 2), eta = c(log(4), atanh(0.75), atanh(-0.35))),
    list(s = autoregressive(9, 3), eta = c(0.2, 0.5, -0.3, 0.4)),
    list(s = scaled_matrix(diag(3)), eta = 0.7),
    list(s = simplex(4), eta = c(0.3, -0.5, 1.1))
  )
  for (cs in cases) {
    r <- param_readable(cs$s, cs$eta)
    # Richardson extrapolation, which shares no code with either the links or
    # the recursion the analytic Jacobian comes from.
    num <- numDeriv::jacobian(function(e) param_readable(cs$s, e)$value, cs$eta)
    expect_equal(unname(r$jacobian), num, tolerance = 1e-8,
      label = cs$s@param_name)

    expect_identical(nrow(r$jacobian), length(r$value))
    expect_identical(ncol(r$jacobian), cs$s@n_free)
    expect_identical(names(r$transform), names(r$value))
    expect_true(all(r$transform %in%
                      c("identity", "log", "atanh", "logit")))
    expect_true(is.character(r$label) && length(r$label) == 1L)
  }
})

test_that("a family with nothing to declare says so", {
  expect_null(param_readable(log_cholesky(3), rep(0.1, 6)))
  expect_null(param_readable(matrix_log(3), rep(0.1, 6)))
  # a fixed matrix carries no free value and therefore no quantity
  expect_null(param_readable(scaled_matrix(diag(3), link = NULL), numeric(0)))
})

test_that("the autoregressive coefficients are the ones the matrix implies", {
  # Yule-Walker read off the assembled covariance is an independent route to
  # the same numbers: the recursion produces them, the matrix confirms them.
  for (q in 1:3) {
    s <- autoregressive(9, order = q)
    set.seed(40 + q)
    eta <- c(0.3, rnorm(q, sd = 0.7))
    r <- param_readable(s, eta)
    m <- param_value(s, eta)
    rho <- m[1, ] / m[1, 1]
    yw <- solve(toeplitz(rho[seq_len(q)]), rho[1L + seq_len(q)])
    expect_equal(unname(r$value[paste0("phi", seq_len(q))]), unname(yw),
      tolerance = 1e-10, label = sprintf("q = %d", q))

    # the marginal variance and the partial autocorrelations are what the
    # coordinates are the transforms of
    expect_equal(unname(r$value[["scale"]]), m[1, 1], tolerance = 1e-12)
    expect_equal(unname(r$value[paste0("pacf", seq_len(q))]), tanh(eta[-1L]),
      tolerance = 1e-12)
    # the last coefficient IS the last partial autocorrelation
    expect_equal(unname(r$value[[paste0("phi", q)]]),
                 unname(r$value[[paste0("pacf", q)]]), tolerance = 1e-12)
  }
})

test_that("the declared quantities are the family's own reading", {
  r <- param_readable(ar1(5), c(log(4), atanh(0.6)))
  expect_equal(unname(r$value), c(4, 0.6), tolerance = 1e-12)
  expect_identical(unname(r$transform), c("log", "atanh"))

  # compound symmetry is definite only above -1/(p-1), and its link says so
  cs <- param_readable(compound_symmetry(5), c(0, -20))
  expect_gt(cs$value[["rho"]], -1 / 4)

  # the simplex reports probabilities, which sum to one
  sp <- param_readable(simplex(4), c(0.3, -0.5, 1.1))
  expect_equal(sum(sp$value), 1, tolerance = 1e-12)
  expect_true(all(sp$value > 0))
})
