# A covariance written as scales times a correlation.
#
# The point of the parametrization is that the coordinates are the quantities
# a reader reads, so the first thing to check is exactly that: the standard
# deviations come back off the diagonal and the correlation block sees a
# matrix with a unit diagonal. Everything else is the factorization
# d(D R D)_ij = d(d_i d_j) * d(R_ij), checked against one stencil on the map.

# one product stencil of the requested order on a scalar function of eta, the
# independent reference for the log-determinant derivatives
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

test_that("the coordinates are the standard deviations and the correlations", {
  s <- dr_prod(3)
  eta <- c(0, log(2), log(0.5), 0.6, 0.9, 1.2)
  v <- param_value(s, eta)
  expect_equal(sqrt(diag(v)), c(1, 2, 0.5), tolerance = 1e-14)
  r <- v / outer(sqrt(diag(v)), sqrt(diag(v)))
  expect_equal(diag(r), rep(1, 3), tolerance = 1e-14)
  expect_equal(unname(r), unname(param_value(correlation_matrix(3), eta[4:6])),
               tolerance = 1e-14)
  expect_equal(param_free(s, v), eta, tolerance = 1e-12, ignore_attr = TRUE)
  expect_identical(s@free_names[1:3], c("log_sd1", "log_sd2", "log_sd3"))
})

test_that("the validator passes on every dimension from two to six", {
  # p = 2 is the smallest legal case, where the correlation carries a single
  # angle and the scale factor's two-index branch has one off-diagonal pair
  for (p in 2:6) {
    s <- dr_prod(p)
    res <- check_parameter(s, verbose = FALSE)
    expect_true(all(res$status == "OK"),
      label = sprintf("p = %d: %s", p,
                      paste(res$check[res$status != "OK"], collapse = ", "))
    )
  }
})

test_that("every order matches one stencil on the map", {
  for (p in c(2L, 4L)) {
    s <- dr_prod(p)
    set.seed(13)
    eta <- rnorm(s@n_free, sd = 0.4)
    tols <- c(1e-6, 1e-5, 1e-4, 5e-3)
    for (o in 1:4) {
      a <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
      b <- switch(o, numerical_d1, numerical_d2, numerical_d3,
                  numerical_d4)(s, eta)
      expect_identical(names(a), unname(param_tuple_names(s, o)))
      for (k in seq_along(a)) {
        expect_lt(max(abs(a[[k]] - b[[k]])) / max(1, max(abs(b[[k]]))), tols[o],
          label = sprintf("p = %d order %d component %s", p, o, names(a)[k])
        )
      }
    }
  }
})

test_that("the log-determinant and its derivatives are right", {
  s <- dr_prod(4)
  set.seed(6)
  eta <- rnorm(s@n_free, sd = 0.4)
  v <- param_value(s, eta)
  expect_equal(param_logdet(s, eta),
               as.numeric(determinant(v, logarithm = TRUE)$modulus),
               tolerance = 1e-12)
  # 2 sum log d + log|R|, the scales counting twice because they multiply on
  # both sides
  expect_equal(param_logdet(s, eta),
               2 * sum(eta[1:4]) + param_logdet(correlation_matrix(4), eta[-(1:4)]),
               tolerance = 1e-12)
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
})

test_that("a log-determinant component mixing two groups is exactly zero", {
  # separability is the structural claim; the second half stops it from being
  # satisfied by everything being zero
  s <- dr_prod(3)
  set.seed(8)
  eta <- rnorm(s@n_free, sd = 0.4)
  p <- 3L
  for (o in 2:4) {
    v <- switch(o, NULL, param_d2logdet, param_d3logdet, param_d4logdet)(s, eta)
    idx <- param_tuple_indices(s, o)
    mixed <- vapply(idx, function(t) {
      any(t <= p) && any(t > p)
    }, logical(1))
    twoscales <- vapply(idx, function(t) {
      all(t <= p) && length(unique(t)) > 1L
    }, logical(1))
    expect_true(any(mixed))
    expect_true(any(twoscales))
    expect_true(all(v[mixed] == 0))
    expect_true(all(v[twoscales] == 0))
    expect_true(any(abs(v[!mixed & !twoscales]) > 1e-8))
  }
})

test_that("a value derivative is supported where the scale indices say", {
  # d(d_i d_j) vanishes unless every scale index of the tuple names i or j, so
  # a component carrying three distinct scale indices is identically zero
  s <- dr_prod(4)
  set.seed(2)
  eta <- rnorm(s@n_free, sd = 0.4)
  d3 <- param_d3(s, eta)
  idx <- param_tuple_indices(s, 3L)
  three <- vapply(idx, function(t) length(unique(t[t <= 4L])) == 3L, logical(1))
  expect_true(any(three))
  expect_true(all(vapply(d3[three], function(m) all(m == 0), logical(1))))
  # and a component in one scale index alone is not zero
  one <- vapply(idx, function(t) all(t <= 4L) && length(unique(t)) == 1L,
                logical(1))
  expect_true(all(vapply(d3[one], function(m) max(abs(m)) > 1e-10, logical(1))))
})

test_that("another positive link on the scales is served, and a signed one is not", {
  # the derivatives of the inverse link enter every component, so a link other
  # than the default exercises the factorization rather than the log's own
  # algebra; and the identity is rejected, a standard deviation being positive
  expect_error(dr_prod(3, link = linkfunctions7::identity_link()),
               "not inside the positive half")
  s <- dr_prod(3, link = linkfunctions7::softplus_link())
  expect_identical(s@free_names[1:3], c("softplus_sd1", "softplus_sd2",
                                        "softplus_sd3"))
  set.seed(21)
  eta <- c(rnorm(3, 0.5, 0.4), rnorm(3, 0, 0.4))
  expect_equal(sqrt(diag(param_value(s, eta))),
               linkfunctions7::linkinv(linkfunctions7::softplus_link(), eta[1:3]),
               tolerance = 1e-13, ignore_attr = TRUE)
  for (o in 1:4) {
    a <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
    b <- switch(o, numerical_d1, numerical_d2, numerical_d3,
                numerical_d4)(s, eta)
    for (k in seq_along(a)) {
      expect_lt(max(abs(a[[k]] - b[[k]])) / max(1, max(abs(b[[k]]))),
                c(1e-6, 1e-5, 1e-4, 5e-3)[o])
    }
  }
})

test_that("the solve and the factor scale the correlation's own", {
  s <- dr_prod(4)
  set.seed(1)
  eta <- rnorm(s@n_free, sd = 0.3)
  v <- param_value(s, eta)
  expect_equal(param_solve(s, eta), solve(v), tolerance = 1e-11,
               ignore_attr = TRUE)
  b <- matrix(rnorm(8), 4, 2)
  expect_equal(param_solve(s, eta, b), solve(v, b), tolerance = 1e-11,
               ignore_attr = TRUE)
  f <- param_factor(s, eta)
  expect_equal(f %*% t(f), v, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("the constructor states what it rejects", {
  expect_error(dr_prod(1), "at least 2")
  expect_error(dr_prod(3, correlation = correlation_matrix(4)), "side 4")
  expect_error(dr_prod(3, correlation = simplex(3)), "matrix_parameter")
  # a rank-deficient correlation would give the product a null space that
  # moves with the free vector
  p <- crossprod(diff(diag(3), differences = 2))
  expect_error(dr_prod(3, correlation = scaled_matrix(p)), "full rank")
})
