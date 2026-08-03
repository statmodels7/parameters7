# The contract every structure inherits: the checks in the bodies of the
# generics, the shared enumeration behind the second-derivative keys, the
# diagonal families, and what a structure defined by its map alone gets.

test_that("the generics validate the free vector before dispatching", {
  s <- log_cholesky(2)

  expect_error(struct_matrix(s, c(1, 2)), "3 free value")
  expect_error(struct_matrix(s, rep(0, 4)), "3 free value")
  expect_error(struct_matrix(s, "a"), "must be numeric")
  expect_error(struct_matrix(s, c(1, NA, 3)), "must be finite")
  expect_error(struct_matrix(s, c(1, Inf, 3)), "must be finite")

  # the message names the free values, so a caller knows what to supply
  expect_error(struct_matrix(s, 1), "log_L1, log_L2, L2.1")

  # names on the values are stripped rather than leaking into the result
  eta <- c(log_L1 = 0.1, log_L2 = 0.2, L2.1 = 0.3)
  expect_identical(
    dimnames(struct_matrix(s, eta)),
    list(c("v1", "v2"), c("v1", "v2"))
  )
})

test_that("the matrix handed back to struct_free is validated", {
  s <- log_cholesky(3)
  expect_error(struct_free(s, 1:9), "numeric matrix")
  expect_error(struct_free(s, diag(2)), "3 by 3")
  m <- diag(3)
  m[1, 2] <- 1
  expect_error(struct_free(s, m), "symmetric")
  m2 <- diag(3)
  m2[1, 1] <- NA
  expect_error(struct_free(s, m2), "missing")
})

test_that("the second-derivative keys and their index pairs agree", {
  s <- log_cholesky(2)
  nm <- struct_pair_names(s)
  idx <- struct_pair_indices(s)

  expect_length(nm, 6L) # d(d+1)/2 with d = 3
  expect_length(idx, 6L)
  # diagonal first, then the off-diagonal pairs
  expect_identical(idx[[1L]], c(1L, 1L))
  expect_identical(idx[[4L]], c(1L, 2L))
  # and the keys are built from the same enumeration
  for (i in seq_along(idx)) {
    expect_identical(
      nm[i],
      paste0(s@free_names[idx[[i]][1L]], ":", s@free_names[idx[[i]][2L]])
    )
  }

  # a structure with no free values has no pairs
  f <- scaled_struct(diag(2), link = NULL)
  expect_length(struct_pair_names(f), 0L)
  expect_length(struct_pair_indices(f), 0L)
})

test_that("the diagonal families behave", {
  d <- diag_struct(3)
  expect_identical(d@free_names, c("d1", "d2", "d3"))
  expect_equal(
    struct_matrix(d, log(c(1, 2, 3))), diag(c(1, 2, 3)),
    ignore_attr = TRUE
  )
  expect_equal(unname(struct_free(d, diag(c(2, 3, 4)))), log(c(2, 3, 4)))
  expect_equal(struct_logdet(d, log(c(2, 3, 4))), sum(log(c(2, 3, 4))))

  sc <- scalar_struct(4)
  expect_identical(sc@n_free, 1L)
  expect_identical(sc@free_names, "scale")
  expect_equal(struct_matrix(sc, log(2.5)), diag(2.5, 4), ignore_attr = TRUE)
  expect_equal(struct_logdet(sc, log(2.5)), 4 * log(2.5))
  expect_equal(unname(struct_dlogdet(sc, log(2.5))), 4)

  # a matrix that is not of the family is refused rather than projected onto it
  expect_error(struct_free(d, matrix(1, 3, 3)), "not diagonal")
  expect_error(struct_free(d, diag(c(1, -1, 1))), "non-positive")
  expect_error(struct_free(sc, diag(c(1, 2, 3, 4))), "constant diagonal")

  # Any link whose range lies inside the positive half line works, and the
  # rule is about the range rather than about the family: a logit reaches only
  # (0, 1), which is a restriction on the entries and not a violation, so it is
  # accepted. What is refused is a range that includes a non-positive value,
  # because a diagonal entry of a positive definite matrix is positive.
  sq <- diag_struct(2, link = linkfunctions7::sqrt_link())
  expect_equal(struct_matrix(sq, c(2, 3)), diag(c(4, 9)), ignore_attr = TRUE)

  lg <- diag_struct(2, link = linkfunctions7::logit_link())
  expect_true(all(diag(struct_matrix(lg, c(0, 1))) > 0))

  expect_error(diag_struct(2, link = linkfunctions7::identity_link()), "positive half")
  expect_error(diag_struct(2, link = "log"), "link object")
})

test_that("the diagonal families are closed form throughout", {
  expect_false(any(struct_is_numerical(diag_struct(3))))
  expect_false(any(struct_is_numerical(scalar_struct(3))))
})

test_that("a structure defined by its map alone is complete", {
  # The bargain: one method, and everything else arrives.
  Bumps <- S7::new_class("Bumps", parent = covstruct, package = NULL)
  gen <- struct_matrix
  S7::method(gen, Bumps) <- function(s, eta, ...) {
    p <- s@dimension
    a <- exp(eta[1L])
    rho <- tanh(eta[2L])
    a * rho^abs(outer(seq_len(p), seq_len(p), "-"))
  }
  b <- Bumps(
    struct_name = "ar1_by_hand", dimension = 4L, n_free = 2L,
    free_names = c("scale", "rho"), rank = 4L,
    null_basis = matrix(numeric(0), 4, 0), role = "covariance",
    struct_params = list()
  )

  eta <- c(0.2, 0.6)
  expect_length(struct_dmatrix(b, eta), 2L)
  expect_length(struct_d2matrix(b, eta), 3L)
  expect_true(is.finite(struct_logdet(b, eta)))
  expect_length(struct_dlogdet(b, eta), 2L)
  expect_equal(
    struct_solve(b, eta), solve(struct_matrix(b, eta)),
    tolerance = 1e-8, ignore_attr = TRUE
  )

  # the numerical first derivative agrees with the closed form nothing in the
  # object knows about
  a <- exp(eta[1L])
  rho <- tanh(eta[2L])
  k <- abs(outer(1:4, 1:4, "-"))
  d_scale <- a * rho^k
  d_rho <- a * k * rho^pmax(k - 1, 0) * (1 - rho^2)
  d_rho[k == 0] <- 0
  expect_equal(struct_dmatrix(b, eta)[[1L]], d_scale, tolerance = 1e-6)
  expect_equal(struct_dmatrix(b, eta)[[2L]], d_rho, tolerance = 1e-6)
})

test_that("the print method states what the object is", {
  expect_output(print(log_cholesky(3)), "Structure: log_cholesky")
  expect_output(print(log_cholesky(3)), "Rank:      3 of 3")
  expect_output(print(log_cholesky(3)), "log_L1")
  expect_output(print(log_cholesky(3)), "From the base class: none")

  p <- crossprod(diff(diag(6), differences = 2))
  expect_output(print(scaled_struct(p)), "null space of dimension 2")
  expect_output(print(scaled_struct(p)), "Role:      precision")

  # a long free vector is truncated rather than filling the screen
  expect_output(print(log_cholesky(8)), "more")
})

test_that("the class validator refuses an inconsistent object", {
  expect_error(
    covstruct(
      struct_name = "x", dimension = 3L, n_free = 2L, free_names = "a",
      rank = 3L, null_basis = matrix(numeric(0), 3, 0), role = "either",
      struct_params = list()
    ),
    "one entry per free value"
  )
  expect_error(
    covstruct(
      struct_name = "x", dimension = 3L, n_free = 2L, free_names = c("a", "a"),
      rank = 3L, null_basis = matrix(numeric(0), 3, 0), role = "either",
      struct_params = list()
    ),
    "unique"
  )
  expect_error(
    covstruct(
      struct_name = "x", dimension = 3L, n_free = 1L, free_names = "a",
      rank = 2L, null_basis = matrix(numeric(0), 3, 0), role = "either",
      struct_params = list()
    ),
    "dimension by"
  )
  expect_error(
    covstruct(
      struct_name = "x", dimension = 3L, n_free = 1L, free_names = "a",
      rank = 3L, null_basis = matrix(numeric(0), 3, 0), role = "sideways",
      struct_params = list()
    ),
    "role"
  )
})
