# kron_identity(): m identical blocks of one structure, every quantity a
# linear lift of the inner parameter's.

test_that("the lift agrees with the assembled Kronecker product", {
  s <- kron_identity(log_cholesky(2), 3)
  inner <- log_cholesky(2)
  eta <- c(0.2, -0.4, 0.1)

  expect_identical(s@dimension, 6L)
  expect_identical(s@n_free, inner@n_free)
  expect_identical(s@free_names, inner@free_names)

  V <- param_value(s, eta)
  expect_equal(V, kronecker(diag(3), param_value(inner, eta)))

  # logdet, its derivatives, and the inverse against the assembled matrix
  expect_equal(param_logdet(s, eta), determinant(V)$modulus[[1]],
               tolerance = 1e-12)
  expect_equal(param_solve(s, eta), solve(V), tolerance = 1e-10)
  F <- param_factor(s, eta)
  expect_equal(F %*% t(F), V, tolerance = 1e-12)

  d1 <- param_d1(s, eta)
  d1_inner <- param_d1(inner, eta)
  for (k in names(d1_inner)) {
    expect_equal(d1[[k]], kronecker(diag(3), d1_inner[[k]]))
  }
})

test_that("the whole contract passes check_parameter", {
  for (case in list(kron_identity(log_cholesky(2), 3),
                    kron_identity(diagonal_matrix(2), 2))) {
    res <- check_parameter(case, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(case@param_name,
                             paste(res$check[res$status != "OK"],
                                   collapse = ", ")))
  }
})

test_that("param_free inverts exactly and refuses what the blocks cannot be", {
  s <- kron_identity(log_cholesky(2), 2)
  eta <- c(0.3, 0.5, -0.2)
  M <- param_value(s, eta)
  expect_equal(param_free(s, M), stats::setNames(eta, s@free_names),
               tolerance = 1e-10)

  # blocks that differ are not representable
  bad <- M
  bad[1, 1] <- bad[1, 1] * 2
  expect_error(param_free(s, bad), "identical diagonal copies")
})

test_that("rank and null basis replicate blockwise", {
  p <- crossprod(diff(diag(4), differences = 2))
  inner <- scaled_matrix(p)
  s <- kron_identity(inner, 3)
  expect_identical(s@rank, 3L * inner@rank)
  M <- param_value(s, 0)
  # the stored null basis annihilates the assembled matrix
  expect_lt(max(abs(M %*% s@null_basis)), 1e-12)
})

test_that("degenerate inputs are rejected", {
  expect_error(kron_identity(diag(2), 3), "matrix_parameter")
  expect_error(kron_identity(log_cholesky(2), 0), "at least 1")
  expect_error(kron_identity(log_cholesky(2), 2.5), "at least 1")
})
