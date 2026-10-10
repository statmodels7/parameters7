# Three defects of the validator: the derivative rows read only the three
# constant free vectors, the membership row never failed on an asymmetric
# matrix, and the battery for a family that is not a matrix ignored `tol` and
# reported a numerical fallback as passed.

status_of <- function(res, name) res$status[res$check == name]

two_by_two <- function(cls_name, value, d1 = NULL) {
  Cls <- S7::new_class(cls_name, parent = matrix_parameter, package = NULL)
  S7::method(param_value, Cls) <- value
  if (!is.null(d1)) S7::method(param_d1, Cls) <- d1
  Cls(
    param_name = "toy", n_free = 2L, free_names = c("a", "b"),
    param_params = list(), dimension = 2L, rank = 2L,
    null_basis = matrix(numeric(0), 2, 0)
  )
}

toy_value <- function(s, eta, ...) {
  m <- diag(exp(eta))
  dimnames(m) <- list(c("v1", "v2"), c("v1", "v2"))
  m
}

test_that("the derivative rows run over every free vector, not three", {
  # right at the three constant vectors, wrong at any other
  d1 <- function(s, eta, ...) {
    f <- if (eta[1] != eta[2]) 1.05 else 1
    list(
      a = f * diag(c(exp(eta[1]), 0)),
      b = f * diag(c(0, exp(eta[2])))
    )
  }
  bad <- two_by_two("Fix5Bad", toy_value, d1)
  good <- two_by_two("Fix5Good", toy_value, function(s, eta, ...) {
    list(a = diag(c(exp(eta[1]), 0)), b = diag(c(0, exp(eta[2]))))
  })

  res <- check_parameter(bad, verbose = FALSE)
  expect_identical(status_of(res, "first derivatives"), "FAIL")
  expect_identical(
    status_of(check_parameter(good, verbose = FALSE), "first derivatives"),
    "OK"
  )
})

test_that("an asymmetric matrix fails the membership check", {
  asym <- function(s, eta, ...) {
    m <- matrix(c(2, 0.7, 0.3, 2), 2) * exp(eta[1])
    dimnames(m) <- list(c("v1", "v2"), c("v1", "v2"))
    m
  }
  bad <- two_by_two("Fix6Bad", asym)
  good <- two_by_two("Fix6Good", function(s, eta, ...) {
    m <- matrix(c(2, 0.5, 0.5, 2), 2) * exp(eta[1])
    dimnames(m) <- list(c("v1", "v2"), c("v1", "v2"))
    m
  })

  expect_identical(
    status_of(check_parameter(bad, verbose = FALSE), "membership"), "FAIL"
  )
  expect_identical(
    status_of(check_parameter(good, verbose = FALSE), "membership"), "OK"
  )
})

test_that("a numerical fallback is NOT CHECKED in the non-matrix battery", {
  Bare <- S7::new_class("Fix7Bare", parent = parameter, package = NULL)
  S7::method(param_value, Bare) <- function(s, eta, ...) {
    e <- exp(c(eta, 0))
    e / sum(e)
  }
  S7::method(param_free, Bare) <- function(s, m, ...) {
    unname(log(m[-length(m)] / m[length(m)]))
  }
  b <- Bare(
    param_name = "bare", n_free = 2L, free_names = c("a", "b"),
    param_params = list()
  )
  expect_true(all(param_is_numerical(b)))

  res <- check_parameter(b, verbose = FALSE)
  for (o in 1:4) {
    expect_identical(status_of(res, paste0("param_d", o)), "NOT CHECKED")
  }
  # a shipped family is unchanged
  ok <- check_parameter(simplex(4), verbose = FALSE)
  expect_true(all(ok$status == "OK"))
})

test_that("tol moves the derivative thresholds of the non-matrix battery", {
  s <- simplex(4)
  expect_true(all(check_parameter(s, verbose = FALSE)$status == "OK"))
  res <- check_parameter(s, tol = 1e-30, verbose = FALSE)
  for (o in 1:4) {
    expect_identical(status_of(res, paste0("param_d", o)), "FAIL")
  }
  # a looser tolerance loosens them, and the default is as before
  loose <- check_parameter(s, tol = 1e-3, verbose = FALSE)
  expect_false(any(loose$status == "FAIL"))
})
