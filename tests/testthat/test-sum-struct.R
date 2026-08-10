# A non-negative combination of fixed matrices: the variance-components
# covariance, and the matrix a penalty with one smoothing parameter per
# component assembles.
#
# The value is linear in the weights, so its derivatives are trivial and the
# work is in the log-determinant, whose cyclic trace expansion is checked
# against one stencil on param_logdet at every order.

fd_scalar <- function(f, eta, tuple, h) {
  used <- sort(unique(tuple))
  mult <- tabulate(match(tuple, used), length(used))
  fac <- list(list(o = c(-1, 1), w = c(-0.5, 0.5)),
              list(o = c(-1, 0, 1), w = c(1, -2, 1)),
              list(o = c(-2, -1, 1, 2), w = c(-0.5, 1, -1, 0.5)),
              list(o = c(-2, -1, 0, 1, 2), w = c(1, -4, 6, -4, 1)))
  fs <- fac[mult]
  grid <- expand.grid(lapply(fs, function(x) seq_along(x$o)))
  acc <- 0
  for (r in seq_len(nrow(grid))) {
    e <- eta
    w <- 1
    for (j in seq_along(used)) {
      pick <- grid[r, j]
      e[used[j]] <- e[used[j]] + fs[[j]]$o[pick] * h
      w <- w * fs[[j]]$w[pick]
    }
    acc <- acc + w * f(e)
  }
  acc / h^length(tuple)
}

ss_cases <- function() {
  set.seed(9)
  Z <- matrix(stats::rnorm(8), 4, 2)
  list(
    sum_struct(list(between = tcrossprod(Z), within = diag(4))),
    sum_struct(list(diag(3), matrix(1, 3, 3), diag(c(1, 0, 2))))
  )
}

test_that("the validator passes", {
  for (s in ss_cases()) {
    res <- check_parameter(s, verbose = FALSE)
    expect_true(all(res$status == "OK"),
      label = paste(res$check[res$status != "OK"], collapse = ", "))
  }
})

test_that("a single component is the scaled matrix", {
  # the smallest legal case, and it has an independent twin already in the
  # package
  p <- crossprod(matrix(stats::rnorm(9), 3))
  a <- sum_struct(list(p))
  b <- scaled_matrix(p)
  eta <- 0.4
  expect_equal(param_value(a, eta), param_value(b, eta), tolerance = 1e-13,
               ignore_attr = TRUE)
  expect_equal(param_logdet(a, eta), param_logdet(b, eta), tolerance = 1e-12)
  expect_equal(unname(param_dlogdet(a, eta)), unname(param_dlogdet(b, eta)),
               tolerance = 1e-12)
})

test_that("every order matches one stencil on the map", {
  for (s in ss_cases()) {
    set.seed(4)
    eta <- stats::rnorm(s@n_free, sd = 0.4)
    tols <- c(1e-6, 1e-5, 1e-4, 5e-3)
    for (o in 1:4) {
      a <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
      b <- switch(o, numerical_d1, numerical_d2, numerical_d3,
                  numerical_d4)(s, eta)
      expect_identical(names(a), unname(param_tuple_names(s, o)))
      for (k in seq_along(a)) {
        expect_lt(max(abs(a[[k]] - b[[k]])) / max(1, max(abs(b[[k]]))), tols[o])
      }
    }
  }
})

test_that("the log-determinant derivatives match one stencil on it", {
  for (s in ss_cases()) {
    set.seed(6)
    eta <- stats::rnorm(s@n_free, sd = 0.35)
    f <- function(e) param_logdet(s, e)
    hs <- c(1e-4, 1e-3, 3e-3, 1e-2)
    tols <- c(1e-7, 1e-6, 1e-4, 5e-3)
    for (o in 1:4) {
      got <- switch(o, param_dlogdet, param_d2logdet, param_d3logdet,
                    param_d4logdet)(s, eta)
      idx <- param_tuple_indices(s, o)
      for (k in seq_along(idx)) {
        ref <- fd_scalar(f, eta, idx[[k]], hs[o])
        expect_lt(abs(got[[k]] - ref), tols[o] * (1 + abs(ref)),
          label = sprintf("logdet order %d component %s", o, names(got)[k])
        )
      }
    }
  }
})

test_that("the cyclic orderings are counted with multiplicity", {
  # deduplicating them makes a component in one weight repeated three times too
  # small by a factor of two and four times by six, which is what the first
  # version of the expansion did
  expect_length(multiset_orderings(c(1L, 1L)), 2L)
  expect_length(multiset_orderings(c(1L, 1L, 1L)), 6L)
  expect_length(multiset_orderings(c(1L, 2L, 3L)), 6L)
  expect_identical(multiset_orderings(integer(0)), list(integer(0)))
  # and the consequence, measured against the stencil
  s <- ss_cases()[[1]]
  eta <- c(0.2, -0.1)
  f <- function(e) param_logdet(s, e)
  idx <- param_tuple_indices(s, 3L)
  pos <- which(vapply(idx, function(t) all(t == 1L), logical(1)))
  got3 <- param_d3logdet(s, eta)[[pos]]
  ref3 <- fd_scalar(f, eta, c(1L, 1L, 1L), 3e-3)
  expect_lt(abs(got3 - ref3), 1e-4 * (1 + abs(ref3)))
  expect_gt(abs(got3 - ref3 / 2), 1e-3)
})

test_that("the rank is the components' and does not move with the weights", {
  # second and first differences on five coefficients: both annihilate the
  # constants and nothing else in common, so the intersection is one
  # dimensional
  s <- sum_struct(list(crossprod(diff(diag(5), differences = 2)),
                       crossprod(diff(diag(5), differences = 1))))
  expect_identical(s@rank, 4L)
  expect_identical(ncol(s@null_basis), 1L)
  # the null residual stays put however far apart the weights are, which a
  # count of small eigenvalues of the assembled matrix would not
  for (r in c(0, 5, 10, 14)) {
    m <- param_value(s, c(0, r * log(10)))
    expect_lt(max(abs(m %*% s@null_basis)) / max(abs(m)), 1e-12)
  }
  expect_error(param_solve(s, c(0, 0)), "rank deficient")
})

test_that("the constructor states what it rejects", {
  expect_error(sum_struct(list()), "non-empty list")
  expect_error(sum_struct(list(diag(2), diag(3))), "same side")
  expect_error(sum_struct(list(matrix(c(1, 2, 3, 4), 2))), "symmetric")
  expect_error(sum_struct(list(diag(c(1, -1)))), "positive semidefinite")
  expect_error(sum_struct(list(a = diag(2), a = diag(2))), "unique")
  expect_error(sum_struct(list(diag(2)),
                          link = linkfunctions7::identity_link()),
               "not inside the positive half")
})

test_that("a matrix outside the family is rejected by param_free", {
  s <- ss_cases()[[1]]
  eta <- c(0.3, -0.2)
  expect_equal(param_free(s, param_value(s, eta)), eta, tolerance = 1e-10,
               ignore_attr = TRUE)
  bad <- param_value(s, eta)
  bad[1, 2] <- bad[2, 1] <- bad[1, 2] + 5
  expect_error(param_free(s, bad), "not a combination")
  # and a combination needing a non-positive weight, which the link cannot
  # carry
  comp <- s@param_params$components
  expect_error(param_free(s, comp[[1]] - 0.5 * comp[[2]]), "non-positive")
})
