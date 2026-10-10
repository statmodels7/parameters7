rel_err <- function(a, b) {
  max(mapply(function(x, y) max(abs(unname(x) - unname(y))) / max(abs(unname(y))),
             a, b))
}

test_that("inverse_of(log_cholesky()) derivatives match the exact inverse derivatives when ill-conditioned", {
  s <- log_cholesky(2)
  io <- inverse_of(s)
  eta <- c(0, -11.5, 0.7)
  expect_lt(rel_err(param_d1(io, eta), param_inv_d1(s, eta)), 1e-10)
  expect_lt(rel_err(param_d2(io, eta), param_inv_d2(s, eta)), 1e-10)

  s3 <- log_cholesky(3)
  io3 <- inverse_of(s3)
  eta3 <- c(0, -11.5, -0.5, 0.7, 0.3, -0.4)
  expect_lt(rel_err(param_d1(io3, eta3), param_inv_d1(s3, eta3)), 1e-10)
  expect_lt(rel_err(param_d2(io3, eta3), param_inv_d2(s3, eta3)), 1e-10)
})

test_that("the log_cholesky route has the structure of the partition sum", {
  for (p in 2:4) {
    s <- log_cholesky(p)
    io <- inverse_of(s)
    eta <- seq(-0.3, 0.4, length.out = s@n_free)
    for (order in 1:2) {
      fast <- if (order == 1L) param_d1(io, eta) else param_d2(io, eta)
      slow <- inverse_derivs(io, eta, order)
      expect_identical(names(fast), names(slow))
      expect_identical(lapply(fast, dimnames), lapply(slow, dimnames))
      expect_identical(lapply(fast, dim), lapply(slow, dim))
      expect_equal(fast, slow, tolerance = 1e-10)
    }
  }
})

test_that("the dimension-one log_cholesky inverse keeps its derivative structure", {
  s <- log_cholesky(1)
  io <- inverse_of(s)
  eta <- 0.3
  expect_identical(names(param_d1(io, eta)), names(inverse_derivs(io, eta, 1L)))
  expect_equal(param_d1(io, eta), inverse_derivs(io, eta, 1L), tolerance = 1e-10)
  expect_equal(param_d2(io, eta), inverse_derivs(io, eta, 2L), tolerance = 1e-10)
})

test_that("orders three and four of inverse_of(log_cholesky()) are unchanged", {
  s <- log_cholesky(2)
  io <- inverse_of(s)
  eta <- c(0.1, -0.4, 0.5)
  expect_identical(param_d3(io, eta), inverse_derivs(io, eta, 3L))
  expect_identical(param_d4(io, eta), inverse_derivs(io, eta, 4L))
})
