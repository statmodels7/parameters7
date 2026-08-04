# The renamed API, the higher orders, and the three new families: simplex,
# transition matrix, matrix logarithm.

test_that("every family's third and fourth derivatives match one stencil", {
  fams <- list(
    log_cholesky(3), matrix_log(3), diagonal_matrix(3), scalar_matrix(3),
    scaled_matrix(crossprod(diff(diag(5), differences = 2))),
    simplex(4), transition_matrix(3)
  )
  for (s in fams) {
    set.seed(7)
    eta <- rnorm(s@n_free, sd = 0.5)
    a3 <- param_d3(s, eta)
    b3 <- numerical_d3(s, eta)
    expect_identical(names(a3), param_tuple_names(s, 3L))
    for (k in seq_along(a3)) {
      expect_lt(max(abs(a3[[k]] - b3[[k]])) / max(1, max(abs(b3[[k]]))),
        1e-4, label = paste(s@param_name, "d3", k)
      )
    }
    a4 <- param_d4(s, eta)
    b4 <- numerical_d4(s, eta)
    for (k in seq_along(a4)) {
      expect_lt(max(abs(a4[[k]] - b4[[k]])) / max(1, max(abs(b4[[k]]))),
        5e-3, label = paste(s@param_name, "d4", k)
      )
    }
  }
})

test_that("higher log-determinant derivatives agree with one stencil on d2", {
  for (s in list(log_cholesky(3), matrix_log(3), diagonal_matrix(3),
                 scalar_matrix(3))) {
    set.seed(3)
    eta <- rnorm(s@n_free, sd = 0.5)
    got <- param_d3logdet(s, eta)
    idx3 <- param_tuple_indices(s, 3L)
    idx2 <- param_tuple_indices(s, 2L)
    key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
    for (i in seq_along(idx3)) {
      t3 <- idx3[[i]]
      k <- t3[3L]
      pos <- match(paste(sort(t3[1:2]), collapse = ","), key2)
      h <- 1e-5 * max(1, abs(eta[k]))
      up <- eta; dn <- eta
      up[k] <- up[k] + h; dn[k] <- dn[k] - h
      want <- (param_d2logdet(s, up)[[pos]] - param_d2logdet(s, dn)[[pos]]) / (2 * h)
      expect_lt(abs(got[[i]] - want), 1e-5 + 1e-4 * abs(want),
        label = paste(s@param_name, names(got)[i])
      )
    }
    expect_true(all(is.finite(param_d4logdet(s, eta))))
  }
})

test_that("the simplex stays a simplex through every order", {
  s <- simplex(5)
  set.seed(11)
  eta <- rnorm(4)
  v <- param_value(s, eta)
  expect_equal(sum(v), 1, tolerance = 1e-14)
  expect_true(all(v > 0))
  # differentiating sum(pi) = 1 kills every order
  for (o in 1:4) {
    d <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
    for (comp in d) expect_lt(abs(sum(comp)), 1e-12)
  }
  # the inverse is exact, and refuses off-simplex input
  expect_equal(unname(param_free(s, v)), eta, tolerance = 1e-12)
  expect_error(param_free(s, v * 1.1), "sum to one")
  expect_error(param_free(s, c(-0.1, 0.4, 0.3, 0.2, 0.2)), "positive")
})

test_that("a transition matrix is row-stochastic with row-local derivatives", {
  s <- transition_matrix(3)
  set.seed(13)
  eta <- rnorm(s@n_free, sd = 0.7)
  m <- param_value(s, eta)
  expect_equal(rowSums(m), rep(1, 3), tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(unname(param_free(s, m)), eta, tolerance = 1e-12)

  # a component pairing two rows is exactly zero
  d2 <- param_d2(s, eta)
  pos <- parameters7:::tm_positions(s)
  idx <- param_tuple_indices(s, 2L)
  for (i in seq_along(idx)) {
    rows <- pos$row[idx[[i]]]
    if (rows[1L] != rows[2L]) expect_true(all(d2[[i]] == 0))
  }

  # and solve/factor are not defined off the matrix branch
  expect_error(param_solve(s, eta), "not defined")
  expect_error(param_factor(s, eta), "not defined")
  expect_error(param_logdet(s, eta))
})

test_that("matrix_log is exact where its chart makes it exact", {
  s <- matrix_log(3)
  set.seed(17)
  eta <- rnorm(6, sd = 0.6)
  m <- param_value(s, eta)
  expect_equal(param_logdet(s, eta), as.numeric(determinant(m)$modulus),
    tolerance = 1e-12
  )
  expect_equal(param_solve(s, eta), solve(m),
    tolerance = 1e-12, ignore_attr = TRUE
  )
  expect_equal(unname(param_free(s, m)), eta, tolerance = 1e-10)
  # near-repeated eigenvalues: Opitz does not cancel
  eta0 <- c(0.5, 0.5 + 1e-9, 0.5, 1e-10, 1e-10, 1e-10)
  d1 <- param_d1(s, eta0)
  n1 <- numerical_d1(s, eta0)
  for (k in seq_along(d1)) {
    expect_lt(max(abs(d1[[k]] - n1[[k]])), 1e-6)
  }
})

test_that("check_parameter runs the reduced battery and catches an injection", {
  res <- check_parameter(simplex(4), verbose = FALSE)
  expect_true(all(res$status == "OK"))
  res_tm <- check_parameter(transition_matrix(3), verbose = FALSE)
  expect_true(all(res_tm$status == "OK"))

  # a 5% corruption of a closed form must be caught, and removing it must
  # restore the pass -- the injection discipline of the whole toolkit
  Broken <- S7::new_class("Broken", parent = SimplexParam,
    package = NULL
  )
  b <- Broken(
    param_name = "simplex", n_free = 3L,
    free_names = paste0("alr", 1:3), param_params = list(n_cat = 4L)
  )
  S7::method(param_d3, Broken) <- function(s, eta, ...) {
    out <- parameters7:::simplex_components(
      s, parameters7:::simplex_tensors(parameters7:::simplex_point(eta), 3L), 3L
    )
    out[[1L]] <- out[[1L]] * 1.05
    out
  }
  res_b <- check_parameter(b, verbose = FALSE)
  expect_true(any(res_b$status == "FAIL"))
})

test_that("tuple enumeration matches its own names at every order", {
  s <- log_cholesky(2)
  for (o in 1:4) {
    idx <- param_tuple_indices(s, o)
    nm <- param_tuple_names(s, o)
    expect_length(idx, choose(s@n_free + o - 1, o))
    expect_identical(
      nm,
      vapply(idx, function(t) paste(s@free_names[t], collapse = ":"), character(1))
    )
  }
})
