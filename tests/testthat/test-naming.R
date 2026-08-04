# Free names name coordinates, not the quantities the coordinates produce.

test_that("link_tag names the transformation and not the link's parameters", {
  lk <- linkfunctions7::bounded_link(lwr = -0.25, upr = 1)
  expect_identical(parameters7:::link_tag(linkfunctions7::identity_link()), "")
  expect_identical(parameters7:::link_tag(linkfunctions7::log_link()), "log")
  expect_identical(parameters7:::link_tag(linkfunctions7::rhobit_link()), "z")
  expect_identical(parameters7:::link_tag(lk), "logit")
  expect_identical(parameters7:::link_tag(linkfunctions7::sqrt_link()), "sqrt")

  # a parametric link names itself with its parameters, and none of that may
  # reach an identifier
  expect_identical(parameters7:::link_tag(linkfunctions7::power_link(0.3)),
                   "power")
  expect_false(grepl("[^a-z0-9_]",
                     parameters7:::link_tag(linkfunctions7::softplus_link(2))))

  # the identity link leaves the quantity's own name alone, which is what a
  # coordinate that is already free deserves
  expect_identical(
    parameters7:::tagged_name(linkfunctions7::identity_link(), "rho"), "rho")
  expect_identical(
    parameters7:::tagged_name(linkfunctions7::log_link(), c("d1", "d2")),
    c("log_d1", "log_d2"))
})

test_that("every family names its coordinates and not its quantities", {
  expect_identical(diagonal_matrix(2)@free_names, c("log_d1", "log_d2"))
  expect_identical(scalar_matrix(3)@free_names, "log_scale")
  expect_identical(scaled_matrix(diag(3))@free_names, "log_scale")
  expect_identical(compound_symmetry(4)@free_names, c("log_scale", "logit_rho"))
  expect_identical(ar1(5)@free_names, c("log_scale", "z_rho"))
  expect_identical(autoregressive(8, 2)@free_names,
                   c("log_scale", "z_pacf1", "z_pacf2"))

  # the families that already named coordinates are untouched
  expect_identical(log_cholesky(2)@free_names, c("log_L1", "log_L2", "L2.1"))
  expect_identical(simplex(3)@free_names, c("alr1", "alr2"))
  expect_identical(correlation_matrix(3)@free_names, c("z2.1", "z3.1", "z3.2"))
})

test_that("a swapped link is a different coordinate and says so", {
  expect_identical(ar1(5, link_scale = linkfunctions7::sqrt_link())@free_names,
                   c("sqrt_scale", "z_rho"))
  expect_identical(
    diagonal_matrix(2, link = linkfunctions7::softplus_link(2))@free_names,
    c("softplus_d1", "softplus_d2"))
})

test_that("what a free name promises is what the coordinate does", {
  # The claim behind the convention: a free value ranges over the whole line.
  # A name carrying no tag would say the coordinate IS the quantity, and a
  # quantity that is free of any constraint is what the tag records the
  # absence of -- so the property is worth asserting rather than the spelling.
  fams <- list(
    log_cholesky(3), matrix_log(3), correlation_matrix(3), diagonal_matrix(3),
    scalar_matrix(3), compound_symmetry(4), ar1(5), autoregressive(6, 2),
    simplex(4), transition_matrix(3)
  )
  for (s in fams) {
    for (k in seq_len(s@n_free)) {
      for (v in c(-20, -12, -6, 6, 12, 20)) {
        lbl <- sprintf("%s: coordinate %s at %g", s@param_name,
                       s@free_names[k], v)
        eta <- numeric(s@n_free)
        eta[k] <- v
        m <- param_value(s, eta)
        expect_true(all(is.finite(m)), label = lbl)
        if (!S7::S7_inherits(s, matrix_parameter) ||
            s@rank != s@dimension) next

        ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
        # A chart onto an open set approaches the boundary without reaching
        # it, and in double precision approaching it far enough is reaching
        # it: at a free value of twenty the smallest eigenvalue of a
        # correlation matrix is 1e-33 of the largest. What survives out
        # there is that none of them is negative beyond rounding, and
        # definiteness proper is asserted where the arithmetic can still
        # carry it.
        expect_gt(min(ev), -1e-13 * max(ev), label = lbl)
        if (abs(v) <= 6) expect_gt(min(ev) / max(ev), 1e-12, label = lbl)
      }
    }
  }
})
