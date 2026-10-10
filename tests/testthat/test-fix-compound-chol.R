# Two decided repairs. The compound-symmetric log-determinant and inverse are
# evaluated from the free value of the correlation, because the correlation
# itself rounds to one of its bounds long before the free value is large. And
# chol_pd() rejects only below a relative eigenvalue ratio of 1e-14, so a badly
# conditioned but positive definite matrix is not rejected.

# mpmath references (80 digits) at eta = (0, e2): the log-determinant's
# correlation part q(e2) and its first four derivatives in e2, and the diagonal
# and off-diagonal entries of the inverse at scale 1.
cs_ref <- function() {
  data.frame(
    p = c(2, 2, 4, 4, 4, 4, 10, 10),
    e2 = c(-40, 40, -40, -20, 20, 40, -40, 40),
    ld0 = c(-38.613705638880109, -38.613705638880109, -37.750659421524767, -17.750659429769381, -57.750659429769381, -117.75065942152477, -36.749170266085518, -356.74917026608552),
    ld1 = c(0.99999999999999999, -0.99999999999999999, 0.99999999999999998, 0.99999999175538553, -2.9999999917553855, -3.0, 0.99999999999999996, -9.0),
    ld2 = c(-8.4967085105831779e-18, -8.4967085105831779e-18, -1.6993417021166356e-17, -8.2446144557673974e-9, -8.2446144557673974e-9, -1.6993417021166356e-17, -4.248354255291589e-17, -4.248354255291589e-17),
    ld3 = c(-8.4967085105831778e-18, 8.4967085105831778e-18, -1.6993417021166356e-17, -8.2446144217805635e-9, 8.2446144217805635e-9, 1.6993417021166356e-17, -4.2483542552915889e-17, 4.2483542552915889e-17),
    ld4 = c(-8.4967085105831777e-18, -8.4967085105831777e-18, -1.6993417021166355e-17, -8.2446143538068961e-9, -8.2446143538068961e-9, -1.6993417021166355e-17, -4.2483542552915889e-17, -4.2483542552915889e-17),
    invdiag = c(58846316709254997.0, 58846316709254997.0, 14711579177313750.0, 30322825.338111894, 272905423.04300703, 1.3240421259582374e+17, 2353852668370200.7, 1.9066206613798619e+17),
    invoff = c(58846316709254996.0, -58846316709254996.0, 14711579177313749.0, 30322824.588111892, -90968474.264335677, -44134737531941247.0, 2353852668370199.8, -21184674015331799.0)
  )
}

test_that("the compound-symmetric log-determinant is exact at large free values", {
  # the case recorded when the defect was found: p = 4, eta = (0, e2)
  s <- compound_symmetry(4)
  for (e2 in c(20, 30, 36, 40)) {
    truth <- log(4) + stats::plogis(e2, log.p = TRUE) +
      3 * (log(4 / 3) + stats::plogis(-e2, log.p = TRUE))
    expect_equal(param_logdet(s, c(0, e2)), truth, tolerance = 1e-13)
  }

  ref <- cs_ref()
  for (i in seq_len(nrow(ref))) {
    s <- compound_symmetry(ref$p[i])
    eta <- c(0, ref$e2[i])
    got <- c(
      param_logdet(s, eta),
      param_dlogdet(s, eta)[["logit_rho"]],
      param_d2logdet(s, eta)[["logit_rho:logit_rho"]],
      param_d3logdet(s, eta)[["logit_rho:logit_rho:logit_rho"]],
      param_d4logdet(s, eta)[["logit_rho:logit_rho:logit_rho:logit_rho"]]
    )
    want <- unlist(ref[i, paste0("ld", 0:4)])
    expect_equal(unname(got), unname(want), tolerance = 1e-12,
      info = sprintf("p = %d, e2 = %g", ref$p[i], ref$e2[i]))
  }
})

test_that("the compound-symmetric inverse is exact at large free values", {
  ref <- cs_ref()
  for (i in seq_len(nrow(ref))) {
    s <- compound_symmetry(ref$p[i])
    inv <- param_solve(s, c(0, ref$e2[i]))
    expect_equal(inv[1L, 1L], ref$invdiag[i], tolerance = 1e-12,
      info = sprintf("p = %d, e2 = %g", ref$p[i], ref$e2[i]))
    expect_equal(inv[1L, 2L], ref$invoff[i], tolerance = 1e-12,
      info = sprintf("p = %d, e2 = %g", ref$p[i], ref$e2[i]))
  }
})

test_that("chol_pd() accepts a positive definite matrix above a ratio of 1e-14", {
  # a positive definite matrix whose extreme eigenvalues are 8e-13 apart in
  # ratio, produced by correlation_matrix(4) at the recorded free vector
  s <- correlation_matrix(4)
  eta <- c(3.061797, -6.284970, -3.104448, 0.445654, -2.267904, 6.961788)
  m <- param_value(s, eta)
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(ev) / max(ev), 1e-12)
  l <- chol_pd(m)
  expect_false(is.null(l))
  expect_equal(l %*% t(l), m, tolerance = 1e-14)
  expect_equal(unname(param_free(s, m)), eta, tolerance = 1e-4)

  # the matrices that are not positive definite are still rejected
  expect_null(chol_pd(diag(c(1, 0))))
  expect_null(chol_pd(diag(c(-1, 1))))
  expect_null(chol_pd(matrix(1, 3, 3)))
  expect_null(chol_pd(matrix(c(1, 2, 2, 1), 2)))
  expect_null(chol_pd(matrix(c(1, NA, NA, 1), 2)))
  expect_null(chol_pd(matrix(c(Inf, 0, 0, 1), 2)))
  expect_null(chol_pd(matrix(0, 0, 0)))

  # the threshold itself, and an exactly singular matrix whose pivot rounds
  # positive, which chol() alone would accept
  expect_null(chol_pd(diag(c(1, 1e-14))))
  expect_false(is.null(chol_pd(diag(c(1, 1e-13)))))
  expect_null(chol_pd(crossprod(matrix(c(1, 2, 3, 2, 4, 6.0000000001), 2))))

  # and the two param_free methods keep their message
  expect_error(param_free(s, matrix(1, 4, 4)), "not positive definite")
  expect_error(param_free(log_cholesky(3), matrix(1, 3, 3)),
    "not positive definite")
  expect_error(param_free(log_cholesky(3), diag(c(1, 1, 0))),
    "not positive definite")
})
