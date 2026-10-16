# Two toy matrix families for the tests of the numerical fallbacks:
# V(eta) = [exp(e1), e3; e3, exp(e2)], with only the value (toy_value()) or
# with the value and its first and second derivatives (toy_d12()).

toy_classes <- local({
  ToyValue <- S7::new_class("ToyValue", parent = matrix_parameter,
                            package = NULL)
  ToyD12 <- S7::new_class("ToyD12", parent = matrix_parameter, package = NULL)
  v <- function(eta) {
    matrix(c(exp(eta[1]), eta[3], eta[3], exp(eta[2])), 2L, 2L)
  }
  gv <- param_value
  S7::method(gv, ToyValue) <- function(s, eta, ...) v(eta)
  S7::method(gv, ToyD12) <- function(s, eta, ...) v(eta)
  g1 <- param_d1
  S7::method(g1, ToyD12) <- function(s, eta, ...) {
    stats::setNames(list(
      matrix(c(exp(eta[1]), 0, 0, 0), 2L, 2L),
      matrix(c(0, 0, 0, exp(eta[2])), 2L, 2L),
      matrix(c(0, 1, 1, 0), 2L, 2L)), s@free_names)
  }
  g2 <- param_d2
  S7::method(g2, ToyD12) <- function(s, eta, ...) {
    z <- matrix(0, 2L, 2L)
    stats::setNames(list(
      matrix(c(exp(eta[1]), 0, 0, 0), 2L, 2L),
      matrix(c(0, 0, 0, exp(eta[2])), 2L, 2L), z, z, z, z),
      param_tuple_names(s, 2L))
  }
  list(value = ToyValue, d12 = ToyD12)
})

toy_make <- function(cls, name) {
  cls(param_name = name, dimension = 2L, n_free = 3L,
      free_names = c("a", "b", "c"), rank = 2L,
      null_basis = matrix(numeric(0), 2, 0), param_params = list())
}
toy_value <- function() toy_make(toy_classes$value, "toy_value")
toy_d12 <- function() toy_make(toy_classes$d12, "toy_d12")

# the exact log-determinant derivatives of the toys, log(exp(a + b) - c^2)
toy_logdet_exact <- function(eta, order) {
  e <- quote(log(exp(a + b) - c^2))
  vars <- c("a", "b", "c")
  idx <- numericals7::tuple_indices(3L, order)
  vapply(idx, function(t) {
    ex <- e
    for (v in vars[t]) ex <- stats::D(ex, v)
    eval(ex, list(a = eta[1], b = eta[2], c = eta[3]))
  }, numeric(1))
}

# the 16 shipped structures of the measurement of 2026-10-16
shipped_structures <- function() {
  list(
    log_cholesky3 = log_cholesky(3), correlation4 = correlation_matrix(4),
    diagonal3 = diagonal_matrix(3), compound4 = compound_symmetry(4),
    ar1_4 = ar1(4), ar1_inv5 = ar1_inv(5),
    autoregressive5_2 = autoregressive(5, order = 2),
    scalar3 = scalar_matrix(3), scaled3 = scaled_matrix(diag(3)),
    matrix_log3 = matrix_log(3), dr_prod3 = dr_prod(3), simplex4 = simplex(4),
    transition3 = transition_matrix(3),
    block_diag = block_diag(log_cholesky(2), ar1(3)),
    kron = kron_identity(log_cholesky(2), 2),
    inverse = inverse_of(log_cholesky(2)))
}
