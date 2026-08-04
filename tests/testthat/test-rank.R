# The rank contract, which is the part of the design that a naive
# implementation gets wrong. A rank read off an assembled matrix is not scale
# invariant; a null space is. These tests pin the difference, because the whole
# reason param_null_basis() takes components rather than a matrix is that the
# obvious alternative fails on an ordinary fitted model.

test_that("param_null_basis reports the rank of a single matrix", {
  expect_identical(param_null_basis(diag(5))$rank, 5L)
  expect_identical(dim(param_null_basis(diag(5))$null_basis), c(5L, 0L))

  p <- crossprod(diff(diag(6), differences = 2))
  ns <- param_null_basis(p)
  expect_identical(ns$rank, 4L)
  expect_identical(dim(ns$null_basis), c(6L, 2L))
  expect_lt(max(abs(p %*% ns$null_basis)), 1e-10)

  # a first-difference penalty leaves only the constants
  expect_identical(param_null_basis(crossprod(diff(diag(6))))$rank, 5L)

  expect_error(param_null_basis(list()), "must not be empty")
})

test_that("the rank of a sum is that of the stack, not the sum of the ranks", {
  p1 <- kronecker(crossprod(diff(diag(8), differences = 2)), diag(4))
  p2 <- kronecker(diag(8), crossprod(diff(diag(4), differences = 2)))

  r1 <- param_null_basis(p1)$rank
  r2 <- param_null_basis(p2)$rank
  rs <- param_null_basis(list(p1, p2))$rank

  expect_identical(r1, 24L)
  expect_identical(r2, 16L)
  expect_identical(rs, 28L)
  # two singular penalties summing to something of larger rank than either,
  # and smaller than their total: neither component nor sum of components
  expect_gt(rs, max(r1, r2))
  expect_lt(rs, r1 + r2)
})

test_that("counting eigenvalues of the assembled sum is not scale invariant", {
  # This is the measurement the design rests on. The mathematical rank of
  # l1 * P1 + l2 * P2 does not depend on the positive scalars; the count of
  # eigenvalues above a relative tolerance does, and badly. The null space
  # does not move at all.
  p1 <- kronecker(crossprod(diff(diag(8), differences = 2)), diag(4))
  p2 <- kronecker(diag(8), crossprod(diff(diag(4), differences = 2)))
  ns <- param_null_basis(list(p1, p2))
  expect_identical(ns$rank, 28L)

  count_ev <- function(m, tol = 1e-8) {
    e <- svd(m, nu = 0, nv = 0)$d
    sum(e > tol * max(e))
  }

  balanced <- p1 + p2
  skewed <- 1e-6 * p1 + 1e6 * p2

  expect_identical(count_ev(balanced), 28L)
  # the count collapses once the two are far apart, although nothing about the
  # matrix's rank has changed
  expect_lt(count_ev(skewed), 28L)

  # the null-space residual, which is what the package uses, does not move
  for (m in list(balanced, skewed, 1e-12 * p1 + 1e12 * p2)) {
    expect_lt(max(abs(m %*% ns$null_basis)) / max(abs(m)), 1e-12)
  }
})

test_that("the components are normalised before being stacked", {
  # Without the individual normalisation the stack itself would lose the
  # smaller component, which is the same defect one level down.
  p1 <- crossprod(diff(diag(6), differences = 2))
  p2 <- crossprod(diff(diag(6)))
  expect_identical(
    param_null_basis(list(p1, p2))$rank,
    param_null_basis(list(1e10 * p1, 1e-10 * p2))$rank
  )
})

test_that("a parameter declares a rank consistent with what it produces", {
  for (s in list(
    log_cholesky(3), diagonal_matrix(4), scalar_matrix(3),
    scaled_matrix(diag(5)), scaled_matrix(crossprod(diff(diag(7), differences = 2)))
  )) {
    eta <- rep(0.3, s@n_free)
    m <- param_value(s, eta)
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    expect_identical(sum(ev > 1e-8 * max(ev)), s@rank, label = s@param_name)
    if (s@rank < s@dimension) {
      expect_lt(max(abs(m %*% s@null_basis)) / max(abs(m)), 1e-10)
    }
  }
})
