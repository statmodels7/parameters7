# The derivatives of the inverse: the default sandwich and the exact
# log-Cholesky method, against a Richardson difference of param_solve(),
# which shares no arithmetic with either.

rich_inv <- function(s, eta, k, l = NULL) {
  f <- function(e) as.numeric(param_solve(s, e))
  if (is.null(l)) {
    g <- function(t) f(replace(eta, k, eta[k] + t))
    return(matrix(numDeriv::jacobian(g, 0), s@dimension))
  }
  g <- function(t) {
    e <- replace(eta, l, eta[l] + t)
    numDeriv::jacobian(function(u) f(replace(e, k, e[k] + u)), 0)
  }
  matrix(numDeriv::jacobian(g, 0), s@dimension)
}

rel <- function(a, b) max(abs(a - b)) / max(abs(b))

test_that("the two methods agree at ordinary free values", {
  for (p in 2:4) {
    s <- log_cholesky(p)
    set.seed(p)
    eta <- rnorm(s@n_free, 0, 0.4)
    e1 <- param_inv_d1(s, eta)
    e2 <- param_inv_d2(s, eta)
    d1 <- S7::method(param_inv_d1, matrix_parameter)(s, eta)
    d2 <- S7::method(param_inv_d2, matrix_parameter)(s, eta)
    expect_identical(names(e1), names(param_d1(s, eta)))
    expect_identical(names(e2), names(param_d2(s, eta)))
    for (k in seq_along(e1)) expect_lt(rel(e1[[k]], d1[[k]]), 1e-12)
    for (k in seq_along(e2)) expect_lt(rel(e2[[k]], d2[[k]]), 1e-12)
  }
})

test_that("the exact method holds where the covariance is nearly singular", {
  skip_if_not_installed("numDeriv")
  s <- log_cholesky(2)
  eta <- c(-0.5, -11.5, 0.07)
  e1 <- param_inv_d1(s, eta)
  e2 <- param_inv_d2(s, eta)
  d2 <- S7::method(param_inv_d2, matrix_parameter)(s, eta)
  for (k in 1:3) expect_lt(rel(e1[[k]], rich_inv(s, eta, k)), 1e-9)
  idx <- param_tuple_indices(s, 2L)
  for (i in seq_along(idx)) {
    ref <- rich_inv(s, eta, idx[[i]][1], idx[[i]][2])
    expect_lt(rel(e2[[i]], ref), 1e-6)
  }
  # the negative control: the sandwich's second derivative is out there
  worst <- max(vapply(seq_along(idx), function(i)
    rel(d2[[i]], rich_inv(s, eta, idx[[i]][1], idx[[i]][2])), numeric(1)))
  expect_gt(worst, 1e-2)
})

test_that("a parameter with no inverse is rejected by name", {
  expect_error(param_inv_d1(simplex(3), c(0, 0)), "not a symmetric matrix")
  s <- scaled_matrix(diag(c(1, 0)))
  expect_error(param_inv_d1(s, 0), "rank deficient")
})
