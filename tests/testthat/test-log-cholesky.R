# The unstructured positive definite family. Every closed form is checked
# against a route the implementation does not take: the matrix against an
# independently assembled factor, the derivatives against finite differences,
# and the log-determinant against an eigendecomposition.

test_that("the constructor validates its arguments", {
  expect_error(log_cholesky(0), "positive integer")
  expect_error(log_cholesky(2.5), "positive integer")
  expect_error(log_cholesky(c(2, 3)), "positive integer")
  expect_error(log_cholesky(3, role = "nonsense"))
})

test_that("the shape and the names follow the frozen ordering", {
  s <- log_cholesky(3)

  expect_identical(s@dimension, 3L)
  expect_identical(s@n_free, 6L)
  expect_identical(s@rank, 3L)
  expect_identical(
    s@free_names,
    c("log_L1", "log_L2", "log_L3", "L2.1", "L3.1", "L3.2")
  )
  expect_identical(dim(s@null_basis), c(3L, 0L))
  expect_identical(s@role, "either")

  # p = 1 is the degenerate case and must still work
  expect_identical(log_cholesky(1)@free_names, "log_L1")
})

test_that("the matrix is the outer product of the factor it declares", {
  s <- log_cholesky(3)
  eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)

  # Assembled here, from the ordering the documentation states, rather than
  # taken from the package.
  l <- matrix(0, 3, 3)
  diag(l) <- exp(eta[1:3])
  l[2, 1] <- eta[4]
  l[3, 1] <- eta[5]
  l[3, 2] <- eta[6]

  expect_equal(param_value(s, eta), tcrossprod(l), ignore_attr = TRUE)
  expect_equal(param_factor(s, eta), l, ignore_attr = TRUE)
})

test_that("the matrix is positive definite over a sweep", {
  set.seed(3)
  for (p in 1:4) {
    s <- log_cholesky(p)
    for (i in 1:5) {
      eta <- stats::runif(s@n_free, -2, 2)
      m <- param_value(s, eta)
      expect_equal(m, t(m))
      expect_gt(min(eigen(m, symmetric = TRUE, only.values = TRUE)$values), 0)
    }
  }
})

test_that("the round trip closes exactly and refuses what is not in the set", {
  set.seed(4)
  s <- log_cholesky(4)
  for (i in 1:5) {
    eta <- stats::runif(s@n_free, -2, 2)
    expect_equal(unname(param_free(s, param_value(s, eta))), eta)
  }
  expect_named(param_free(s, param_value(s, rep(0, s@n_free))), s@free_names)

  # not positive definite
  m <- diag(4)
  m[1, 1] <- -1
  expect_error(param_free(s, m), "not positive definite")

  # not symmetric, refused before dispatch
  ns <- diag(4)
  ns[1, 2] <- 1
  expect_error(param_free(s, ns), "symmetric")
})

test_that("the derivatives agree with finite differences", {
  set.seed(5)
  s <- log_cholesky(3)
  for (i in 1:3) {
    eta <- stats::runif(s@n_free, -1.5, 1.5)

    a1 <- param_d1(s, eta)
    n1 <- parameters7:::numerical_d1(s, eta)
    expect_named(a1, s@free_names)
    for (k in seq_along(a1)) {
      expect_equal(a1[[k]], n1[[k]], tolerance = 1e-6, ignore_attr = TRUE)
    }

    a2 <- param_d2(s, eta)
    n2 <- parameters7:::numerical_d2(s, eta)
    expect_named(a2, param_tuple_names(s))
    for (k in seq_along(a2)) {
      expect_equal(a2[[k]], n2[[k]], tolerance = 1e-5, ignore_attr = TRUE)
    }
  }
})

test_that("every derivative matrix is symmetric", {
  s <- log_cholesky(3)
  eta <- c(0.4, -0.3, 0.2, 1, -0.5, 0.7)
  for (m in c(param_d1(s, eta), param_d2(s, eta))) {
    expect_equal(m, t(m))
  }
})

test_that("the log-determinant is linear in the free vector", {
  set.seed(6)
  s <- log_cholesky(4)
  for (i in 1:4) {
    eta <- stats::runif(s@n_free, -2, 2)
    m <- param_value(s, eta)
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    expect_equal(param_logdet(s, eta), sum(log(ev)))
    expect_equal(param_logdet(s, eta), 2 * sum(eta[1:4]))

    # the gradient is 2 on the diagonal directions and 0 elsewhere
    expect_equal(
      unname(param_dlogdet(s, eta)), c(rep(2, 4), rep(0, 6))
    )
    expect_true(all(param_d2logdet(s, eta) == 0))
  }
})

test_that("the log-determinant gradient satisfies the trace identity", {
  # tr(M^{-1} dM/deta_k), computed with base::solve, which shares nothing with
  # the closed form.
  set.seed(7)
  s <- log_cholesky(3)
  eta <- stats::runif(s@n_free, -1, 1)
  mi <- solve(param_value(s, eta))
  d <- param_d1(s, eta)
  want <- vapply(d, function(dk) sum(diag(mi %*% dk)), numeric(1))
  expect_equal(unname(param_dlogdet(s, eta)), unname(want))
})

test_that("solve and factor agree with base R", {
  set.seed(8)
  s <- log_cholesky(4)
  eta <- stats::runif(s@n_free, -1, 1)
  m <- param_value(s, eta)

  expect_equal(param_solve(s, eta), solve(m), ignore_attr = TRUE)
  b <- matrix(stats::rnorm(8), 4, 2)
  expect_equal(param_solve(s, eta, b), solve(m, b), ignore_attr = TRUE)
  expect_equal(tcrossprod(param_factor(s, eta)), m, ignore_attr = TRUE)
})

test_that("everything is closed form", {
  expect_false(any(param_is_numerical(log_cholesky(3))))
})
