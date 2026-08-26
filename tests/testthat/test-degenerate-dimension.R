# Every family swept down to its minimum dimension, which is where a matrix of
# indices degenerates and where paste0 recycles a zero-length argument against a
# string literal. correlation_matrix(1) reported n_free = 1 with the free name
# "z." until 0.12.0, and the object could not be used.

test_that("correlation_matrix(1) is a constant with no free values", {
  s <- correlation_matrix(1)
  expect_identical(s@n_free, 0L)
  expect_identical(s@free_names, character(0))

  # the whole contract answers at the only free vector there is
  e0 <- numeric(0)
  expect_equal(unname(param_value(s, e0)), matrix(1, 1, 1))
  expect_identical(length(param_d1(s, e0)), 0L)
  expect_identical(length(param_d2(s, e0)), 0L)
  expect_equal(param_logdet(s, e0), 0)
  expect_identical(length(param_free(s, matrix(1, 1, 1))), 0L)
  expect_equal(param_factor(s, e0), matrix(1, 1, 1))
  expect_equal(param_solve(s, e0, matrix(2, 1, 1)), matrix(2, 1, 1))
})

test_that("the validator passes at p = 1, with the two logdet rows not checked", {
  res <- check_parameter(correlation_matrix(1), verbose = FALSE)
  expect_false(any(res$status == "FAIL"))
  # there is no free value to differentiate in, so those two rows report
  # NOT CHECKED rather than a number
  nc <- res$check[res$status == "NOT CHECKED"]
  expect_setequal(nc, c("logdet gradient", "logdet hessian"))
})

test_that("the guard is what fixes it: the unguarded expression still recycles", {
  # the negative control. Without the length test, paste0 recycles the
  # zero-length index against the length-one literals and gives one name.
  rows <- integer(0)
  cols <- integer(0)
  expect_identical(paste0("z", rows, ".", cols), "z.")
  expect_identical(if (length(rows)) paste0("z", rows, ".", cols) else character(0),
                   character(0))
})

test_that("p = 1 does not disturb any larger dimension", {
  for (p in 2:6) {
    s <- correlation_matrix(p)
    expect_identical(s@n_free, as.integer(p * (p - 1L) / 2L))
    expect_identical(length(s@free_names), s@n_free)
  }
  expect_identical(correlation_matrix(2)@free_names, "z2.1")
  expect_identical(correlation_matrix(3)@free_names, c("z2.1", "z3.1", "z3.2"))
})

test_that("every family either works or refuses at its minimum dimension", {
  # the sweep that found this one. A family added later that neither works nor
  # refuses at its smallest legal input fails here.
  builds <- list(
    "log_cholesky(1)"      = function() log_cholesky(1),
    "matrix_log(1)"        = function() matrix_log(1),
    "diagonal_matrix(1)"   = function() diagonal_matrix(1),
    "scalar_matrix(1)"     = function() scalar_matrix(1),
    "correlation_matrix(1)" = function() correlation_matrix(1),
    "compound_symmetry(1)" = function() compound_symmetry(1),
    "ar1(1)"               = function() ar1(1),
    "autoregressive(1, 1)" = function() autoregressive(1, order = 1),
    "scaled_matrix(1x1)"   = function() scaled_matrix(matrix(1, 1, 1)),
    "transition_matrix(1)" = function() transition_matrix(1),
    "simplex(1)"           = function() simplex(1)
  )
  for (nm in names(builds)) {
    s <- tryCatch(builds[[nm]](), error = function(e) e)
    if (inherits(s, "error")) {
      # a refusal is a legitimate answer, provided the message says why
      expect_true(nzchar(conditionMessage(s)), label = nm)
      next
    }
    res <- tryCatch(check_parameter(s, verbose = FALSE), error = function(e) e)
    expect_false(inherits(res, "error"),
                 label = sprintf("%s: check_parameter signalled", nm))
    if (!inherits(res, "error")) {
      expect_false(any(res$status == "FAIL"), label = nm)
    }
  }
})
