# Interface defects fixed together: a block with no free values inside
# block_diag(), the orders three and four of a fully known scaled_matrix(),
# the names of param_free(), the membership check of a bounded diagonal link,
# and the messages of the argument and matrix checks.

fx <- function() scaled_matrix(diag(2), link = NULL)

test_that("block_diag() accepts a block with no free values", {
  for (s in list(block_diag(a = fx(), b = ar1(3)),
                 block_diag(a = ar1(3), b = fx()),
                 block_diag(a = log_cholesky(2), f = fx(), b = ar1(3)),
                 block_diag(fx()))) {
    expect_identical(s@n_free, length(s@free_names))
    expect_true(anyDuplicated(s@free_names) == 0L)
  }
  s <- block_diag(a = fx(), b = ar1(3))
  expect_identical(s@free_names, paste0("b_", ar1(3)@free_names))
  expect_identical(s@dimension, 5L)
  expect_identical(s@rank, 5L)
})

test_that("block_diag() with a fixed block agrees with the assembled matrix", {
  set.seed(5)
  fixed_cases <- list(
    block_diag(a = fx(), b = ar1(3)),
    block_diag(a = log_cholesky(2), f = fx(), b = ar1(3)),
    block_diag(a = ar1(3), b = fx())
  )
  for (s in fixed_cases) {
    eta <- rnorm(s@n_free, sd = 0.4)
    M <- param_value(s, eta)
    expect_equal(M, name_dims(M, s))
    expect_lt(max(abs(M - t(M))), 1e-14)
    # the fixed block is the identity of side 2, wherever it sits
    j <- which(vapply(s@param_params$blocks, function(b) b@n_free == 0L,
                      logical(1)))
    r <- s@param_params$rows[[j]]
    expect_equal(unname(M[r, r]), diag(2))
    # logdet, solve, factor
    expect_equal(param_logdet(s, eta), as.numeric(determinant(M)$modulus))
    expect_lt(max(abs(param_solve(s, eta) - solve(M))), 1e-12)
    L <- param_factor(s, eta)
    expect_lt(max(abs(L %*% t(L) - M)), 1e-12)
    # the null basis has no columns and the rank is the dimension
    expect_identical(dim(s@null_basis), c(s@dimension, 0L))
    # round trip
    expect_lt(max(abs(param_free(s, M) - eta)), 1e-10)
    # derivatives against one stencil on the map, and the log-determinant's
    # against differences of the log-determinant
    for (o in 1:4) {
      a <- switch(o, param_d1, param_d2, param_d3, param_d4)(s, eta)
      b <- switch(o, numerical_d1, numerical_d2, numerical_d3,
                  numerical_d4)(s, eta)
      expect_identical(length(a), length(param_tuple_names(s, o)))
      expect_identical(names(a), unname(param_tuple_names(s, o)))
      tol <- c(1e-6, 1e-5, 1e-4, 5e-3)[o]
      for (k in seq_along(a)) {
        expect_lt(max(abs(a[[k]] - b[[k]])) / max(1, max(abs(b[[k]]))), tol)
      }
      dl <- switch(o, param_dlogdet, param_d2logdet, param_d3logdet,
                   param_d4logdet)(s, eta)
      expect_identical(names(dl), unname(param_tuple_names(s, o)))
    }
    # the first order of the log-determinant against a central difference
    g <- vapply(seq_len(s@n_free), function(k) {
      h <- 1e-5
      e1 <- eta; e1[k] <- e1[k] + h
      e2 <- eta; e2[k] <- e2[k] - h
      (param_logdet(s, e1) - param_logdet(s, e2)) / (2 * h)
    }, numeric(1))
    expect_lt(max(abs(unname(param_dlogdet(s, eta)) - g)), 1e-6)
  }
})

test_that("a composition made only of fixed blocks has no free values", {
  s <- block_diag(fx())
  expect_identical(s@n_free, 0L)
  expect_identical(s@free_names, character(0))
  expect_equal(unname(param_value(s, numeric(0))), diag(2))
  expect_length(param_free(s, diag(2)), 0L)
  res <- check_parameter(block_diag(a = fx(), b = ar1(3)), verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("a fully known scaled_matrix() has empty derivatives at every order", {
  f <- scaled_matrix(diag(3), link = NULL)
  expect_identical(param_d3(f, numeric(0)), stats::setNames(list(), character(0)))
  expect_identical(param_d4(f, numeric(0)), stats::setNames(list(), character(0)))
  expect_identical(param_d3logdet(f, numeric(0)),
                   stats::setNames(numeric(0), character(0)))
  expect_identical(param_d4logdet(f, numeric(0)),
                   stats::setNames(numeric(0), character(0)))
  # the shapes agree with the lower orders
  expect_identical(param_d2(f, numeric(0)), param_d3(f, numeric(0)))
  expect_identical(param_d2logdet(f, numeric(0)),
                   param_d3logdet(f, numeric(0)))
})

test_that("param_free() names its result by free_names", {
  s <- block_diag(subject = log_cholesky(2), time = ar1(3))
  eta <- c(0.1, -0.2, 0.3, 0.4, -0.5)
  expect_identical(names(param_free(s, param_value(s, eta))), s@free_names)
  s2 <- block_diag(a = fx(), b = ar1(3))
  expect_identical(names(param_free(s2, param_value(s2, c(0.4, -0.5)))),
                   s2@free_names)
  d <- dr_prod(3)
  m <- param_value(d, c(0.1, -0.2, 0.3, 0.2, -0.1, 0.3))
  expect_identical(names(param_free(d, m)), d@free_names)
  expect_identical(names(param_free(d, m))[1:3], c("log_sd1", "log_sd2", "log_sd3"))
})

test_that("param_free() of a bounded diagonal rejects entries outside the range", {
  d <- diagonal_matrix(2, link = linkfunctions7::logit_link())
  expect_no_error(param_free(d, diag(c(0.2, 0.7))))
  expect_error(param_free(d, diag(c(2, 3))), "range")
  expect_error(param_free(d, diag(c(0.5, 1))), "range")
  expect_no_warning(try(param_free(d, diag(c(2, 3))), silent = TRUE))
  q <- scalar_matrix(2, link = linkfunctions7::logit_link())
  expect_no_error(param_free(q, diag(c(0.3, 0.3))))
  expect_error(param_free(q, diag(c(1.5, 1.5))), "range")
  # the default link has an unbounded range and is unaffected
  expect_equal(unname(param_free(diagonal_matrix(2), diag(c(2, 3)))),
               log(c(2, 3)))
})

test_that("kron_identity() rejects a block count that is not a number", {
  x <- log_cholesky(2)
  msg <- "'m' must be a single integer of at least 1"
  for (bad in list(Inf, -Inf, 1e10, TRUE, "3", NA, 0, 2.5, c(1, 2), NULL,
                   NaN, .Machine$integer.max + 1)) {
    expect_error(kron_identity(x, bad), msg, fixed = TRUE)
  }
  expect_no_warning(try(kron_identity(x, Inf), silent = TRUE))
  expect_identical(kron_identity(x, 3L)@param_params$m, 3L)
  expect_identical(kron_identity(x, 3)@param_params$m, 3L)
})

test_that("param_free() rejects a non-finite matrix with a plain message", {
  s <- log_cholesky(2)
  m <- diag(2)
  m[1, 1] <- Inf
  expect_error(param_free(s, m), "finite", fixed = TRUE)
  m[1, 1] <- -Inf
  expect_error(param_free(s, m), "finite", fixed = TRUE)
  m <- diag(2)
  m[1, 2] <- m[2, 1] <- Inf
  expect_error(param_free(s, m), "finite", fixed = TRUE)
  m <- diag(2)
  m[1, 1] <- NA
  expect_error(param_free(s, m), "missing", fixed = TRUE)
})

test_that("the validators return a message for a non-scalar slot", {
  expect_error(
    parameter(param_name = "x", n_free = 1:2, free_names = c("a", "b"),
              param_params = list()),
    "n_free must be a single non-negative integer", fixed = TRUE
  )
  expect_error(
    parameter(param_name = "x", n_free = NA_integer_, free_names = "a",
              param_params = list()),
    "n_free must be a single non-negative integer", fixed = TRUE
  )
  mk <- function(dimension = 2L, rank = 2L, null_basis = matrix(0, 2, 0)) {
    matrix_parameter(param_name = "x", n_free = 0L, free_names = character(0),
                     param_params = list(), dimension = dimension,
                     rank = rank, null_basis = null_basis)
  }
  expect_error(mk(dimension = c(2L, 3L)),
               "dimension must be a single positive integer", fixed = TRUE)
  expect_error(mk(rank = c(1L, 2L)),
               "rank must be a single integer in 0:dimension", fixed = TRUE)
  expect_error(mk(dimension = integer(0)),
               "dimension must be a single positive integer", fixed = TRUE)
  expect_error(mk(dimension = NA_integer_),
               "dimension must be a single positive integer", fixed = TRUE)
  expect_error(mk(rank = NA_integer_),
               "rank must be a single integer in 0:dimension", fixed = TRUE)
  expect_s3_class(mk(), "S7_object")
})
