# Repairs to sum_struct(): the log-determinant derivatives of a rank deficient
# family, the names of param_free(), and components that are linearly dependent.

fd_scalar_fix <- function(f, eta, tuple, h) {
  used <- sort(unique(tuple))
  mult <- tabulate(match(tuple, used), length(used))
  fac <- list(list(o = c(-1, 1), w = c(-0.5, 0.5)),
              list(o = c(-1, 0, 1), w = c(1, -2, 1)),
              list(o = c(-2, -1, 1, 2), w = c(-0.5, 1, -1, 0.5)),
              list(o = c(-2, -1, 0, 1, 2), w = c(1, -4, 6, -4, 1)))
  fs <- fac[mult]
  grid <- expand.grid(lapply(fs, function(x) seq_along(x$o)))
  acc <- 0
  for (r in seq_len(nrow(grid))) {
    e <- eta
    w <- 1
    for (j in seq_along(used)) {
      pick <- grid[r, j]
      e[used[j]] <- e[used[j]] + fs[[j]]$o[pick] * h
      w <- w * fs[[j]]$w[pick]
    }
    acc <- acc + w * f(e)
  }
  acc / h^length(tuple)
}

ss_deficient <- function() {
  list(
    sum_struct(list(matrix(1, 3, 3))),
    sum_struct(list(crossprod(diff(diag(5), differences = 2)),
                    crossprod(diff(diag(5)))))
  )
}

test_that("a rank deficient family has log-determinant derivatives", {
  for (s in ss_deficient()) {
    eta <- seq(-0.3, 0.4, length.out = s@n_free)
    f <- function(e) param_logdet(s, e)
    hs <- c(1e-4, 1e-3, 3e-3, 1e-2)
    tols <- c(1e-6, 1e-5, 1e-3, 5e-2)
    for (o in 1:4) {
      got <- switch(o, param_dlogdet, param_d2logdet, param_d3logdet,
                    param_d4logdet)(s, eta)
      expect_length(got, choose(s@n_free + o - 1, o))
      idx <- param_tuple_indices(s, o)
      for (k in seq_along(idx)) {
        ref <- fd_scalar_fix(f, eta, idx[[k]], hs[o])
        expect_lt(abs(got[[k]] - ref), tols[o] * (1 + abs(ref)),
                  label = sprintf("order %d component %d", o, k))
      }
    }
  }
})

test_that("check_parameter completes on a rank deficient family", {
  for (s in ss_deficient()) {
    res <- check_parameter(s, verbose = FALSE)
    # solve is not checked on a rank deficient family, by design
    bad <- !res$status %in% c("OK", "NOT CHECKED")
    expect_false(any(bad), label = paste(res$check[bad], collapse = ", "))
    expect_true(all(res$status[grepl("logdet", res$check)] == "OK"))
  }
})

test_that("param_free returns a vector named by the free names", {
  s <- sum_struct(list(between = matrix(1, 3, 3), within = diag(3)))
  eta <- c(0.3, -0.2)
  got <- param_free(s, param_value(s, eta))
  expect_identical(names(got), s@free_names)
  expect_equal(unname(got), eta, tolerance = 1e-12)
})

test_that("linearly dependent components are rejected by the constructor", {
  expect_error(sum_struct(list(diag(2), 2 * diag(2))), "linearly dependent")
  expect_error(sum_struct(list(a = diag(3), b = matrix(1, 3, 3),
                               c = diag(3) + matrix(1, 3, 3))),
               "linearly dependent")
  # the test is scale free
  expect_error(sum_struct(list(1e-9 * diag(2), 1e9 * diag(2))),
               "linearly dependent")
  # independent components are accepted at any scale
  expect_no_error(sum_struct(list(1e-9 * diag(2), 1e9 * matrix(1, 2, 2))))
})
