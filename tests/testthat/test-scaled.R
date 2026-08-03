# The scaled fixed matrix, and with it the rank-deficient path the package
# exists to admit. The identities the family rests on are checked against
# eigendecompositions and against numerical optimisation, neither of which the
# implementation uses.

pen2 <- function(k) crossprod(diff(diag(k), differences = 2))

test_that("the constructor validates its matrix", {
  expect_error(scaled_struct(1:4), "numeric matrix")
  expect_error(scaled_struct(matrix(1:6, 2, 3)), "square")
  expect_error(scaled_struct(matrix(c(1, 2, 3, 4), 2, 2)), "symmetric")
  expect_error(scaled_struct(matrix(NA_real_, 2, 2)), "missing")
  expect_error(scaled_struct(diag(c(1, -1))), "positive semidefinite")
  expect_error(scaled_struct(matrix(0, 3, 3)), "positive semidefinite")

  # a link whose range is not the positive half line
  expect_error(
    scaled_struct(diag(3), link = linkfunctions7::identity_link()),
    "positive half"
  )
})

test_that("the rank and the null space are those of the fixed matrix", {
  p <- pen2(6)
  s <- scaled_struct(p)

  expect_identical(s@dimension, 6L)
  expect_identical(s@rank, 4L)
  expect_identical(dim(s@null_basis), c(6L, 2L))
  expect_identical(s@role, "precision")

  # the null space of a second-difference penalty is the constants and the
  # straight lines, which is what it means for the penalty to leave them alone
  lin <- cbind(rep(1, 6), 1:6)
  expect_lt(max(abs(p %*% lin)), 1e-10)
  # the declared basis spans exactly that space
  expect_lt(max(abs(qr.resid(qr(lin), s@null_basis))), 1e-10)
  expect_equal(crossprod(s@null_basis), diag(2), ignore_attr = TRUE)

  # ...and the matrix really does annihilate it, at any scale
  for (eta in c(-8, 0, 8)) {
    m <- struct_matrix(s, eta)
    expect_lt(max(abs(m %*% s@null_basis)) / max(abs(m)), 1e-12)
  }

  # a full-rank fixed matrix declares no null space
  expect_identical(scaled_struct(diag(4))@rank, 4L)
  expect_identical(dim(scaled_struct(diag(4))@null_basis), c(4L, 0L))
})

test_that("the matrix is the fixed one times the scale", {
  p <- pen2(5)
  s <- scaled_struct(p)
  for (eta in c(-2, 0, 1.5)) {
    expect_equal(struct_matrix(s, eta), exp(eta) * p, ignore_attr = TRUE)
  }
})

test_that("the log pseudo-determinant is the rank times the log scale", {
  p <- pen2(7)
  s <- scaled_struct(p)
  r <- s@rank

  logdet_plus <- function(m) {
    e <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    sum(log(e[e > 1e-8 * max(e)]))
  }

  for (eta in c(-3, 0, 2.5)) {
    expect_equal(struct_logdet(s, eta), logdet_plus(struct_matrix(s, eta)))
    expect_equal(
      struct_logdet(s, eta),
      r * eta + logdet_plus(p)
    )
  }
})

test_that("the derivative of the log pseudo-determinant is the rank", {
  # Under the log link, at every scale, to the last digit; and against a
  # central difference of the log-determinant itself.
  p <- pen2(8)
  s <- scaled_struct(p)
  for (eta in c(-6, -1, 0, 4, 9)) {
    expect_equal(unname(struct_dlogdet(s, eta)), as.numeric(s@rank))
    h <- 1e-5
    fd <- (struct_logdet(s, eta + h) - struct_logdet(s, eta - h)) / (2 * h)
    expect_equal(unname(struct_dlogdet(s, eta)), fd, tolerance = 1e-8)
    expect_equal(unname(struct_d2logdet(s, eta)), 0)
  }
})

test_that("the constant is what makes the scale estimable", {
  # The identity the package keeps the normalising constant for: minimising
  # (lambda/2) b'Pb - (r/2) log lambda gives lambda = r / (b'Pb), while
  # dropping the constant sends it to zero. Checked against optimize(), which
  # knows nothing about the closed form.
  set.seed(9)
  p <- pen2(9)
  s <- scaled_struct(p)
  b <- stats::rnorm(9)
  q <- drop(t(b) %*% p %*% b)

  with_const <- function(lambda) lambda / 2 * q - s@rank / 2 * log(lambda)
  num <- stats::optimize(with_const, c(1e-8, 1e8), tol = 1e-12)$minimum
  expect_equal(num, s@rank / q, tolerance = 1e-6)

  without <- stats::optimize(function(lambda) lambda / 2 * q, c(1e-8, 1e8),
    tol = 1e-12
  )$minimum
  expect_lt(without, 1e-6)

  # and the constant in question is exactly the structure's log-determinant
  eta <- log(num)
  expect_equal(struct_logdet(s, eta) - struct_logdet(s, 0), s@rank * eta)
})

test_that("the derivatives are the fixed matrix times the link's own", {
  p <- pen2(5)
  s <- scaled_struct(p)
  for (eta in c(-1, 0.7)) {
    expect_equal(struct_dmatrix(s, eta)[[1L]], exp(eta) * p, ignore_attr = TRUE)
    expect_equal(struct_d2matrix(s, eta)[[1L]], exp(eta) * p, ignore_attr = TRUE)
  }

  # against finite differences, which share no code with the closed form
  n1 <- covstructs7:::numerical_dmatrix(s, 0.3)
  expect_equal(struct_dmatrix(s, 0.3)[[1L]], n1[[1L]], tolerance = 1e-6)
})

test_that("a rank-deficient structure refuses to be inverted", {
  s <- scaled_struct(pen2(6))
  expect_error(struct_solve(s, 0), "rank deficient")
  expect_error(struct_factor(s, 0), "rank deficient")

  # a full-rank one does not
  expect_silent(struct_solve(scaled_struct(diag(3)), 0))
})

test_that("a fixed matrix has no free value at all", {
  s <- scaled_struct(diag(3), link = NULL)

  expect_identical(s@n_free, 0L)
  expect_identical(s@free_names, character(0))
  expect_identical(s@struct_name, "fixed")
  expect_equal(struct_matrix(s, numeric(0)), diag(3), ignore_attr = TRUE)
  expect_length(struct_dmatrix(s, numeric(0)), 0L)
  expect_length(struct_d2matrix(s, numeric(0)), 0L)
  expect_length(struct_dlogdet(s, numeric(0)), 0L)
  expect_equal(struct_logdet(s, numeric(0)), 0)
  expect_length(struct_pair_names(s), 0L)

  # and a free vector is refused, since there is none to supply
  expect_error(struct_matrix(s, 0), "0 free value")
})

test_that("the round trip recovers the scale and refuses another matrix", {
  p <- pen2(6)
  s <- scaled_struct(p)
  expect_equal(unname(struct_free(s, 3.7 * p)), log(3.7))
  expect_error(struct_free(s, pen2(6) + diag(6)), "not a multiple")
  expect_error(struct_free(s, -p), "positive multiple")
})

test_that("the ridge is the identity, scaled", {
  s <- scaled_struct(diag(4))
  expect_identical(s@rank, 4L)
  expect_equal(struct_logdet(s, 1.2), 4 * 1.2)
  expect_equal(struct_matrix(s, log(3)), diag(3, 4), ignore_attr = TRUE)
})
