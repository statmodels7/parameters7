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
  expect_identical(ar1(5, link_scale = linkfunctions7::logit_link())@free_names,
                   c("logit_scale", "z_rho"))
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
test_that("every family labels its matrices the same way", {
  # The convention name_dims() states, tested on the margins the families
  # return rather than on the names they declare. It covers the value and the
  # four derivative orders. Before it was pinned, the four composition wrappers
  # returned matrices carrying no labels at all while check_parameter()'s
  # shapes-and-names row passed, that row reading only the declared names.
  fam <- list(
    log_cholesky       = log_cholesky(3),
    matrix_log         = matrix_log(3),
    diagonal_matrix    = diagonal_matrix(3),
    scalar_matrix      = scalar_matrix(3),
    correlation_matrix = correlation_matrix(3),
    compound_symmetry  = compound_symmetry(3),
    ar1                = ar1(3),
    autoregressive     = autoregressive(4, order = 2),
    scaled_matrix      = scaled_matrix(diag(3)),
    kron_identity      = kron_identity(log_cholesky(2), 2),
    block_diag         = block_diag(list(log_cholesky(2), ar1(2))),
    dr_prod            = dr_prod(3, correlation = correlation_matrix(3)),
    sum_struct         = sum_struct(list(diag(3), crossprod(diff(diag(3)))))
  )

  for (nm in names(fam)) {
    s <- fam[[nm]]
    eta <- seq(0.15, by = 0.05, length.out = s@n_free)
    want <- rep(list(paste0("v", seq_len(s@dimension))), 2L)

    expect_identical(dimnames(param_value(s, eta)), want, info = nm)
    for (k in 1:4) {
      d <- switch(k, param_d1(s, eta), param_d2(s, eta),
                  param_d3(s, eta), param_d4(s, eta))
      expect_identical(dimnames(d[[1L]]), want,
                       info = paste(nm, "order", k))
    }
  }

  # A family whose rows are states rather than variables keeps its own prefix.
  t3 <- transition_matrix(3)
  expect_identical(dimnames(param_value(t3, rep(0.1, t3@n_free)))[[1L]],
                   c("s1", "s2", "s3"))

  # A factor and a solve are left as each family builds them. sum_struct is the
  # one whose two answers come from its own value, so it unnames them to match
  # the families that never label them.
  ss <- fam$sum_struct
  eta <- c(0.15, 0.2)
  expect_null(dimnames(param_factor(ss, eta)))
  expect_null(dimnames(param_solve(ss, eta, diag(3))))
})
