# The validator, and the injections that give it meaning. Every probe is
# paired: one asserts that a deliberate defect is caught, one asserts that the
# unbroken structure still passes, so the first cannot be satisfied by simply
# weakening the check.

pen2 <- function(k) crossprod(diff(diag(k), differences = 2))

status_of <- function(res, name) res$status[res$check == name]

test_that("the shipped families pass every check they can be given", {
  set.seed(11)
  for (s in list(
    log_cholesky(3), diag_struct(3), scalar_struct(4),
    scaled_struct(diag(4)), scaled_struct(pen2(6))
  )) {
    res <- check_covstruct(s, verbose = FALSE)
    expect_false(any(res$status == "FAIL"), label = s@struct_name)
  }
})

test_that("a rank-deficient family is checked, not excused", {
  set.seed(12)
  res <- check_covstruct(scaled_struct(pen2(6)), verbose = FALSE)

  # membership goes through the null space and really runs
  expect_identical(status_of(res, "membership"), "OK")
  # the solve is the one thing it cannot be asked for
  expect_identical(status_of(res, "solve"), "NOT CHECKED")
  # everything else is live
  expect_identical(status_of(res, "log-determinant"), "OK")
  expect_identical(status_of(res, "logdet gradient"), "OK")
})

test_that("a wrong first derivative is caught", {
  Wrong <- S7::new_class("Wrong", parent = LogCholeskyStruct, package = NULL)
  gen <- struct_dmatrix
  S7::method(gen, Wrong) <- function(s, eta, ...) {
    lapply(S7::method(struct_dmatrix, LogCholeskyStruct)(s, eta), function(m) 1.05 * m)
  }
  good <- log_cholesky(3)
  bad <- Wrong(
    struct_name = "wrong", dimension = good@dimension, n_free = good@n_free,
    free_names = good@free_names, rank = good@rank,
    null_basis = good@null_basis, role = good@role,
    struct_params = good@struct_params
  )

  set.seed(13)
  expect_identical(
    status_of(check_covstruct(bad, verbose = FALSE), "first derivatives"),
    "FAIL"
  )
  set.seed(13)
  expect_identical(
    status_of(check_covstruct(good, verbose = FALSE), "first derivatives"),
    "OK"
  )
})

test_that("a wrong log-determinant gradient is caught", {
  Wrong <- S7::new_class("Wrong2", parent = LogCholeskyStruct, package = NULL)
  gen <- struct_dlogdet
  S7::method(gen, Wrong) <- function(s, eta, ...) {
    1.05 * S7::method(struct_dlogdet, LogCholeskyStruct)(s, eta)
  }
  good <- log_cholesky(3)
  bad <- Wrong(
    struct_name = "wrong", dimension = good@dimension, n_free = good@n_free,
    free_names = good@free_names, rank = good@rank,
    null_basis = good@null_basis, role = good@role,
    struct_params = good@struct_params
  )

  set.seed(14)
  expect_identical(
    status_of(check_covstruct(bad, verbose = FALSE), "logdet gradient"), "FAIL"
  )
  set.seed(14)
  expect_identical(
    status_of(check_covstruct(good, verbose = FALSE), "logdet gradient"), "OK"
  )
})

test_that("a declared rank that does not match the matrix is caught", {
  # The membership check exists for this: a structure that claims a null space
  # its matrix does not have.
  p <- pen2(6)
  s <- scaled_struct(p)
  liar <- ScaledStruct(
    struct_name = "liar", dimension = s@dimension, n_free = s@n_free,
    free_names = s@free_names,
    rank = 3L,
    null_basis = cbind(s@null_basis, diag(6)[, 1]),
    role = s@role, struct_params = s@struct_params
  )

  set.seed(15)
  expect_identical(
    status_of(check_covstruct(liar, verbose = FALSE), "membership"), "FAIL"
  )
  set.seed(15)
  expect_identical(
    status_of(check_covstruct(s, verbose = FALSE), "membership"), "OK"
  )
})

test_that("a quantity with no independent reference is reported unchecked", {
  # A structure that supplies only its matrix: the derivatives are finite
  # differences, and comparing them with finite differences would be the same
  # arithmetic twice.
  Bare <- S7::new_class("Bare", parent = covstruct, package = NULL)
  gen <- struct_matrix
  S7::method(gen, Bare) <- function(s, eta, ...) {
    d <- exp(eta)
    m <- diag(d, nrow = s@dimension)
    m[1, 2] <- m[2, 1] <- 0.3 * sqrt(d[1] * d[2])
    m
  }
  b <- Bare(
    struct_name = "bare", dimension = 3L, n_free = 3L,
    free_names = c("a", "b", "c"), rank = 3L,
    null_basis = matrix(numeric(0), 3, 0), role = "either",
    struct_params = list()
  )

  expect_true(all(struct_is_numerical(b)[c("struct_dmatrix", "struct_d2matrix")]))

  set.seed(16)
  res <- check_covstruct(b, verbose = FALSE)
  expect_identical(status_of(res, "first derivatives"), "NOT CHECKED")
  expect_identical(status_of(res, "second derivatives"), "NOT CHECKED")
  expect_identical(status_of(res, "log-determinant"), "NOT CHECKED")
  # but the membership and the solve, which have independent references, run
  expect_identical(status_of(res, "membership"), "OK")
  expect_identical(status_of(res, "solve and factor"), "OK")
  # and the inverse map is refused rather than approximated
  expect_identical(status_of(res, "round trip"), "NOT CHECKED")
  expect_error(struct_free(b, diag(3)), "exact or")
})

test_that("the validator refuses something that is not a structure", {
  expect_error(check_covstruct(diag(3)), "must inherit")
})

test_that("verbose printing states the shape and the verdicts", {
  set.seed(17)
  expect_output(check_covstruct(log_cholesky(2)), "Structure: log_cholesky")
  set.seed(17)
  expect_output(check_covstruct(log_cholesky(2)), "rank 2")
  set.seed(17)
  expect_output(check_covstruct(scaled_struct(pen2(5))), "NOT CHECKED")
})
