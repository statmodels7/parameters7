# The compiled Leibniz assembly against its R twin. The kernel exploits that
# every derivative of the factor is a single-entry matrix; the twin builds
# the same components through dense products and shares no arithmetic with
# it beyond the factor itself.

test_that("the compiled assembly matches the R twin at every order", {
  for (p in c(2L, 4L, 5L)) {
    s <- log_cholesky(p)
    set.seed(p)
    eta <- rnorm(s@n_free, sd = 0.4)
    for (ord in 1:4) {
      a <- parameters7:::chol_leibniz(s, eta, ord)
      b <- parameters7:::.chol_leibniz_r(s, eta, ord)
      expect_identical(names(a), names(b))
      worst <- max(vapply(seq_along(a), function(i) {
        max(abs(a[[i]] - b[[i]]))
      }, numeric(1)))
      expect_lt(worst, 1e-14)
      expect_identical(dimnames(a[[1L]]), dimnames(b[[1L]]))
    }
  }
})

test_that("orders one and two keep their contract through the kernel", {
  s <- log_cholesky(3)
  eta <- c(0.2, -0.1, 0.4, 0.3, -0.2, 0.1)
  d1 <- param_d1(s, eta)
  expect_identical(names(d1), s@free_names)
  # d/dk M against a central difference of param_value
  h <- 1e-6
  for (k in seq_len(s@n_free)) {
    up <- eta; up[k] <- up[k] + h
    dn <- eta; dn[k] <- dn[k] - h
    fd <- (param_value(s, up) - param_value(s, dn)) / (2 * h)
    expect_lt(max(abs(d1[[k]] - fd)), 1e-7)
  }
})
