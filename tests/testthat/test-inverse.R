test_that("inverse_of is the inverse, at every free vector swept", {
  for (s in list(ar1(4), correlation_matrix(4), log_cholesky(3),
                 compound_symmetry(4), autoregressive(5, order = 2))) {
    g <- inverse_of(s)
    set.seed(7)
    for (i in 1:3) {
      eta <- stats::rnorm(s@n_free, 0, 0.4)
      expect_lt(max(abs(unclass(param_value(g, eta)) %*%
                          unclass(param_value(s, eta)) -
                          diag(s@dimension))), 1e-10)
    }
  }
})

test_that("inverse_of passes check_parameter on every family it wraps", {
  for (s in list(ar1(4), correlation_matrix(4), log_cholesky(3),
                 compound_symmetry(4), autoregressive(5, order = 2),
                 dr_prod(3), matrix_log(3), diagonal_matrix(4))) {
    r <- check_parameter(inverse_of(s), verbose = FALSE)
    expect_identical(sum(r$status == "FAIL"), 0L,
                     info = paste("inverse_of(", s@param_name, ")"))
  }
})

test_that("the log-determinant is the inner one negated", {
  s <- ar1(5)
  g <- inverse_of(s)
  eta <- c(0.3, atanh(0.55))
  expect_equal(param_logdet(g, eta), -param_logdet(s, eta))
  expect_equal(unlist(param_dlogdet(g, eta)), -unlist(param_dlogdet(s, eta)))
  expect_equal(unlist(param_d4logdet(g, eta)), -unlist(param_d4logdet(s, eta)))
})

test_that("the log-determinant derivatives keep the inner container", {
  # a list where the inner returns a numeric vector makes check_parameter
  # fail several frames away, on a subtraction it cannot do
  g <- inverse_of(ar1(4))
  eta <- c(0.2, atanh(0.4))
  expect_false(is.list(param_dlogdet(g, eta)))
  expect_false(is.list(param_d2logdet(g, eta)))
})

test_that("inverting twice returns the original matrix", {
  s <- ar1(4)
  eta <- c(0.3, atanh(0.55))
  expect_equal(unclass(param_value(inverse_of(inverse_of(s)), eta)),
               unclass(param_value(s, eta)))
})

test_that("a rank-deficient family has no inverse and is rejected", {
  expect_error(inverse_of(scaled_matrix(diag(c(1, 1, 0)))), "rank")
  expect_error(inverse_of(simplex(3)), "matrix_parameter")
})

test_that("the ordered set partitions are the Fubini numbers", {
  expect_identical(lengths(lapply(1:4, ordered_set_partitions)),
                   c(1L, 3L, 13L, 75L))
})


test_that("ar1_inv is tridiagonal and is the AR(1) inverse", {
  s <- ar1_inv(6)
  eta <- c(log(2), atanh(0.6))
  om <- unclass(param_value(s, eta))
  lag <- abs(outer(1:6, 1:6, "-"))
  expect_identical(max(abs(om[lag > 1])), 0)
  expect_lt(max(abs(om %*% unclass(param_value(ar1(6), eta)) - diag(6))), 1e-12)
})

test_that("ar1_inv's written-out derivatives agree with the partition sum", {
  # two independent implementations of one object: the twin writes the three
  # distinct entries out, the wrapper sums over ordered block partitions, and
  # they share no arithmetic
  for (p in c(2, 3, 5, 8)) {
    s <- ar1_inv(p)
    g <- inverse_of(ar1(p))
    for (rho in c(-0.9, -0.3, 0, 0.3, 0.9)) {
      eta <- c(0.37, atanh(rho))
      for (k in 1:4) {
        f <- switch(k, param_d1, param_d2, param_d3, param_d4)
        expect_lt(max(abs(unlist(f(s, eta)) - unlist(f(g, eta)))), 1e-9)
      }
    }
  }
})

test_that("autoregressive_inv is banded of the process's own order", {
  for (q in 1:3) {
    s <- autoregressive_inv(7, order = q)
    eta <- c(0.2, atanh(seq(0.5, by = -0.15, length.out = q)))
    om <- unclass(param_value(s, eta))
    lag <- abs(outer(1:7, 1:7, "-"))
    expect_lt(max(abs(om[lag > q])), 1e-10)
  }
})

test_that("both twins pass check_parameter", {
  for (s in list(ar1_inv(5), ar1_inv(2), autoregressive_inv(6, order = 2),
                 autoregressive_inv(5, order = 1))) {
    r <- check_parameter(s, verbose = FALSE)
    expect_identical(sum(r$status == "FAIL"), 0L, info = s@param_name)
  }
})

test_that("the twins carry the free vector of the family they invert", {
  expect_identical(ar1_inv(5)@free_names, ar1(5)@free_names)
  expect_identical(autoregressive_inv(6, 2)@free_names,
                   autoregressive(6, 2)@free_names)
})

test_that("recip_derivs is the reciprocal's chain, not a truncated power", {
  # power_derivs() truncates beyond the exponent, which for -1 would make
  # every order zero; this is the negative control for using it here
  h <- function(e) exp(e)
  e <- 0.37
  v <- list(h(e), h(e), h(e), h(e), h(e))
  got <- recip_derivs(v)
  want <- lapply(0:4, function(k) (-1)^k * exp(-e))
  expect_equal(unlist(got), unlist(want))
  expect_true(all(unlist(power_derivs(h(e), -1L)) == 0))
})

test_that("autoregressive_inv's written-out derivatives agree with the sum", {
  # two independent implementations of one object: the twin factorizes
  # Omega = U' diag(tau) U and differentiates the factors, the wrapper sums
  # over ordered block partitions of the whole matrix
  for (cfg in list(c(5, 1), c(6, 2), c(9, 2), c(8, 3), c(10, 4))) {
    p <- cfg[1]; q <- cfg[2]
    s <- autoregressive_inv(p, order = q)
    g <- inverse_of(autoregressive(p, order = q))
    for (seed in 1:2) {
      set.seed(seed)
      eta <- c(stats::rnorm(1, 0, 0.5), stats::rnorm(q, 0, 0.6))
      expect_identical(max(abs(unclass(param_value(s, eta)) -
                                 unclass(param_value(g, eta)))), 0)
      for (k in 1:4) {
        f <- switch(k, param_d1, param_d2, param_d3, param_d4)
        expect_lt(max(abs(unlist(f(s, eta)) - unlist(f(g, eta)))), 1e-9)
      }
    }
  }
})

test_that("a coefficient does not depend on a correlation its order misses", {
  # U's row k + 1 carries the order-k coefficients, which read the first k
  # partial autocorrelations only, so differentiating in a later one is
  # EXACTLY zero and not merely small.
  #
  # The claim is about the FACTOR and not about Omega: Omega[2, 1] sums over
  # every row of U, including the later ones that do carry the higher orders,
  # so it moves with all of them. Asserting it there would assert something
  # false.
  s <- autoregressive_inv(7, order = 3)
  eta <- c(0.3, atanh(c(0.6, -0.3, 0.4)))
  n <- s@n_free
  cd <- parameters7:::ar_inv_codes(param_tuple_indices(s, 1L), n)
  f <- parameters7:::ar_inv_factors(s, eta, cd)
  code <- function(v) {
    as.integer(sum(tabulate(v, nbins = n) * 5^(seq_len(n) - 1L))) + 1L
  }
  # row 2 holds the order-one coefficient: it moves with pacf1 (free value 2)
  # and with nothing later
  expect_true(abs(f$u[[code(2L)]][2L, 1L]) > 1e-8)
  expect_identical(f$u[[code(3L)]][2L, 1L], 0)
  expect_identical(f$u[[code(4L)]][2L, 1L], 0)
  # row 3 holds the order-two coefficients: pacf2 reaches it, pacf3 does not
  expect_true(abs(f$u[[code(3L)]][3L, 1L]) > 1e-8)
  expect_identical(f$u[[code(4L)]][3L, 1L], 0)
  # and the scale reaches none of them
  expect_identical(max(abs(f$u[[code(1L)]])), 0)
})

test_that("the lower-order rows of U are the sub-family's own coefficients", {
  # the identity the written-out route rests on, asserted rather than assumed
  p <- 8; q <- 3
  s <- autoregressive(p, order = q)
  eta <- c(0.3, atanh(c(0.6, -0.3, 0.4)))
  pr <- parameters7:::ar_prediction(s, eta)
  for (k in seq_len(q)) {
    sk <- autoregressive(p, order = k)
    phik <- parameters7:::ar_taylor(sk, eta[seq_len(k + 1L)])$phi[, 1L]
    expect_identical(max(abs(-pr$u[k + 1L, seq.int(k, 1L)] - phik)), 0)
  }
  # and the factorization itself
  om <- t(pr$u) %*% diag(1 / pr$v) %*% pr$u
  expect_lt(max(abs(om - solve(unclass(param_value(s, eta))))), 1e-12)
})

test_that("the empty sub-multiset has a usable key", {
  # a list cannot be subscripted by an empty name: x[[""]] is an error, not a
  # lookup, so the undifferentiated factor keys as code 1 and never as ""
  cd <- parameters7:::ar_inv_codes(list(c(2L, 2L)), 3L)
  expect_true(1L %in% cd$all)
  expect_identical(cd$cnt[[match(1L, cd$all)]], integer(3))
})
