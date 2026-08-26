# check_parameter() draws random free vectors. Until 0.13.0 the matrix branch
# took them from the caller's stream and the branch for a family that is not a
# matrix called set.seed() and never put back what it found, so a call inside a
# simulation silently changed the simulation.

test_that("the matrix branch leaves the caller's stream alone", {
  set.seed(42)
  want <- rnorm(1)
  set.seed(42)
  invisible(check_parameter(log_cholesky(2), verbose = FALSE))
  expect_identical(rnorm(1), want)
})

test_that("the branch for a family that is not a matrix leaves it alone too", {
  set.seed(42)
  want <- rnorm(1)
  set.seed(42)
  invisible(check_parameter(simplex(3), verbose = FALSE))
  expect_identical(rnorm(1), want)

  set.seed(42)
  invisible(check_parameter(transition_matrix(3), verbose = FALSE))
  expect_identical(rnorm(1), want)
})

test_that("the report does not depend on the state it was called from", {
  a <- check_parameter(log_cholesky(3), verbose = FALSE)
  set.seed(999)
  invisible(runif(17))
  b <- check_parameter(log_cholesky(3), verbose = FALSE)
  expect_identical(a, b)
})

test_that("no seed is left behind when the caller had none", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(check_parameter(log_cholesky(2), verbose = FALSE))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the stream is restored even when the call signals", {
  set.seed(42)
  want <- rnorm(1)
  set.seed(42)
  expect_error(check_parameter("not a parameter"))
  expect_identical(rnorm(1), want)
})

test_that("capture_seed and restore_seed are a round trip", {
  set.seed(42)
  before <- capture_seed()
  invisible(runif(5))
  expect_false(identical(capture_seed(), before))
  restore_seed(before)
  expect_identical(capture_seed(), before)

  # the negative control: without the restore the stream really has moved,
  # so the round trip above is testing something
  set.seed(42)
  invisible(runif(5))
  expect_false(identical(capture_seed(), before))

  # and NULL means remove, not "leave whatever is there"
  restore_seed(NULL)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})
