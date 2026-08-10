# A block diagonal of distinct matrix parameters.
#
# The composition rederives nothing: every quantity is the blocks' own, placed
# where the block sits, and everything spanning two blocks is exactly zero
# because the free values of one block do not enter another. The tests
# therefore check the assembly against an independent route at orders one to
# four, and then check the structural claim -- separability -- in a way that
# cannot be satisfied vacuously.

bd_cases <- function() {
  list(
    block_diag(log_cholesky(2), ar1(3)),
    block_diag(a = diagonal_matrix(2), b = compound_symmetry(3),
               c = log_cholesky(2)),
    block_diag(correlation_matrix(3), scalar_matrix(2))
  )
}

test_that("the validator passes on every combination", {
  for (s in bd_cases()) {
    res <- check_parameter(s, verbose = FALSE)
    expect_true(all(res$status == "OK"),
      label = sprintf("%s: %s", s@param_name,
                      paste(res$check[res$status != "OK"], collapse = ", "))
    )
  }
})

test_that("the smallest case is covered", {
  # two blocks of side one is where every range degenerates to a single index,
  # the case a loop written for several is least likely to be tried on
  s <- block_diag(scalar_matrix(1), scalar_matrix(1))
  expect_identical(s@dimension, 2L)
  expect_identical(s@n_free, 2L)
  res <- check_parameter(s, verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("every order matches one stencil on the map", {
  for (s in bd_cases()) {
    set.seed(11)
    eta <- rnorm(s@n_free, sd = 0.5)
    tols <- c(1e-6, 1e-5, 1e-4, 5e-3)
    for (o in 1:4) {
      a <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
      b <- switch(o, numerical_d1, numerical_d2, numerical_d3,
                  numerical_d4)(s, eta)
      expect_identical(names(a), unname(param_tuple_names(s, o)))
      for (k in seq_along(a)) {
        expect_lt(max(abs(a[[k]] - b[[k]])) / max(1, max(abs(b[[k]]))), tols[o],
          label = sprintf("%s order %d component %s", s@param_name, o,
                          names(a)[k])
        )
      }
    }
  }
})

test_that("a component spanning two blocks is exactly zero, and one inside is not", {
  # the second half is what stops the first from being satisfied by returning
  # zero everywhere
  s <- block_diag(log_cholesky(2), ar1(3))
  owner <- c(1L, 1L, 1L, 2L, 2L)
  set.seed(3)
  eta <- rnorm(s@n_free, sd = 0.5)
  for (o in 2:4) {
    d <- switch(o, NULL, param_d2, param_d3, param_d4)(s, eta)
    idx <- param_tuple_indices(s, o)
    crossed <- vapply(idx, function(t) length(unique(owner[t])) > 1L, logical(1))
    expect_true(any(crossed))
    expect_true(all(vapply(d[crossed], function(m) all(m == 0), logical(1))))
    expect_true(any(vapply(d[!crossed], function(m) max(abs(m)) > 1e-8,
                           logical(1))))
  }
})

test_that("the log-determinant is separable across blocks", {
  s <- block_diag(log_cholesky(2), ar1(3))
  owner <- c(1L, 1L, 1L, 2L, 2L)
  set.seed(7)
  eta <- rnorm(s@n_free, sd = 0.4)
  # it is the sum of the blocks', which is what a block-diagonal determinant is
  blocks <- list(log_cholesky(2), ar1(3))
  want <- param_logdet(blocks[[1]], eta[1:3]) + param_logdet(blocks[[2]], eta[4:5])
  expect_equal(param_logdet(s, eta), want, tolerance = 1e-13)
  for (o in 2:4) {
    v <- switch(o, NULL, param_d2logdet, param_d3logdet, param_d4logdet)(s, eta)
    idx <- param_tuple_indices(s, o)
    crossed <- vapply(idx, function(t) length(unique(owner[t])) > 1L, logical(1))
    expect_true(all(v[crossed] == 0))
  }
  # the AR(1) block's own logdet curvature is not zero, so the check above is
  # not passing because everything is
  expect_gt(max(abs(param_d2logdet(s, eta))), 1e-8)
})

test_that("one block reproduces the block itself", {
  inner <- correlation_matrix(3)
  s <- block_diag(inner)
  set.seed(5)
  eta <- rnorm(inner@n_free, sd = 0.5)
  # ignore_attr, because a block's own dimnames cannot be carried onto a
  # composite that in general holds several blocks
  expect_equal(param_value(s, eta), param_value(inner, eta), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(param_logdet(s, eta), param_logdet(inner, eta), tolerance = 1e-14)
  expect_equal(unname(param_d3(s, eta)), unname(param_d3(inner, eta)),
               tolerance = 1e-14, ignore_attr = TRUE)
  # the names differ by the block label, which is the point of the prefix
  expect_identical(s@free_names, paste0("b1_", inner@free_names))
})

test_that("identical blocks agree with kron_identity on the value", {
  # the two compositions differ in their free vectors and not in what they
  # build, so feeding the same values to every block must reproduce the
  # replication -- a check of each against the other
  inner <- log_cholesky(2)
  eta <- c(0.3, -0.1, 0.2)
  k <- kron_identity(inner, 3)
  b <- block_diag(inner, inner, inner)
  expect_equal(param_value(b, rep(eta, 3)), param_value(k, eta),
               tolerance = 1e-14)
  expect_equal(param_logdet(b, rep(eta, 3)), param_logdet(k, eta),
               tolerance = 1e-13)
  expect_identical(b@n_free, 3L * k@n_free)
})

test_that("the rank and the null basis come from the components", {
  # a rank-deficient block makes the composite deficient by the same amount,
  # and the null basis is the blocks' placed in their own rows
  p <- crossprod(diff(diag(5), differences = 2))
  s <- block_diag(full = log_cholesky(2), deficient = scaled_matrix(p))
  expect_identical(s@rank, 2L + 3L)
  expect_identical(dim(s@null_basis), c(7L, 2L))
  expect_true(all(abs(s@null_basis[1:2, ]) < 1e-14))
  eta <- c(0.1, -0.2, 0.3, 0.5)
  m <- param_value(s, eta)
  expect_lt(max(abs(m %*% s@null_basis)) / max(abs(m)), 1e-12)
  # and the solve is rejected, as it is for a deficient family of its own
  expect_error(param_solve(s, eta), "rank deficient")
})

test_that("a matrix that is not block diagonal is rejected", {
  s <- block_diag(log_cholesky(2), ar1(3))
  eta <- c(0.2, -0.1, 0.3, 0.1, 0.4)
  m <- param_value(s, eta)
  expect_equal(param_free(s, m), eta, tolerance = 1e-10, ignore_attr = TRUE)
  m[1, 4] <- m[4, 1] <- 0.5
  expect_error(param_free(s, m), "not block diagonal")
})

test_that("the constructor states what it rejects", {
  expect_error(block_diag(), "at least one block")
  expect_error(block_diag(log_cholesky(2), simplex(3)), "matrix_parameter")
  expect_error(block_diag(a = log_cholesky(2), a = ar1(2)), "unique")
  # a list is accepted in place of the dots, since a caller assembling blocks
  # in a loop has one
  s <- block_diag(list(log_cholesky(2), ar1(3)))
  expect_identical(s@dimension, 5L)
})
