test_that("a family's own param_d2logdet() reaches the higher log-determinant orders", {
  # The guard read only param_d1() and param_d2(), so a family that wrote
  # param_d2logdet() in closed form, the quantity the default orders three
  # and four difference, was rejected as if that quantity were numerical.
  OwnD2 <- S7::new_class("OwnD2", parent = matrix_parameter, package = NULL)
  gv <- param_value
  S7::method(gv, OwnD2) <- function(s, eta, ...) diag(exp(eta), 2L)
  gd <- param_d2logdet
  # log|M| = eta1 + eta2, so every second derivative is zero
  S7::method(gd, OwnD2) <- function(s, eta, ...) {
    stats::setNames(numeric(3L), param_tuple_names(s, 2L))
  }
  s <- OwnD2(param_name = "own_d2", dimension = 2L, n_free = 2L,
             free_names = c("a", "b"), rank = 2L,
             null_basis = matrix(numeric(0), 2, 0), param_params = list())
  eta <- c(0.3, -0.2)
  expect_equal(unname(param_d3logdet(s, eta)), numeric(4L))
  expect_equal(unname(param_d4logdet(s, eta)), numeric(5L))
})

test_that("numerical_d1() is symmetric for a matrix family", {
  # A value that is symmetric only to rounding gave a first derivative with
  # an asymmetry of the same size, where the higher orders were symmetrized.
  Lopsided <- S7::new_class("Lopsided", parent = matrix_parameter,
                            package = NULL)
  gv <- param_value
  S7::method(gv, Lopsided) <- function(s, eta, ...) {
    m <- matrix(c(exp(eta[1L]), eta[1L], eta[1L], 2), 2L, 2L)
    m[1L, 2L] <- m[1L, 2L] * (1 + 1e-9 * eta[1L])
    m
  }
  s <- Lopsided(param_name = "lopsided", dimension = 2L, n_free = 1L,
                free_names = "a", rank = 2L,
                null_basis = matrix(numeric(0), 2, 0), param_params = list())
  d <- numerical_d1(s, 0.4)[[1L]]
  expect_identical(d, t(d))
  expect_identical(param_d1(s, 0.4)[[1L]], d)
})

test_that("the derivatives of correlation_matrix() have an exactly zero diagonal", {
  s <- correlation_matrix(4)
  eta <- c(0.3, -1.1, 0.8, 2.4, -0.6, 1.7)
  gens <- list(param_d1, param_d2, param_d3, param_d4)
  for (n in 1:4) {
    diags <- unlist(lapply(gens[[n]](s, eta), diag))
    expect_identical(unname(diags), numeric(length(diags)),
                     label = sprintf("order %d diagonals", n))
  }
})
