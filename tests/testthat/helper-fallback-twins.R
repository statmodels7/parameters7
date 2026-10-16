# The R implementations of the numerical fallbacks as they stood in
# parameters7 0.24.0, before they became adapters on numericals7::fd_tensor().
# They are kept as twins (test-fallback-tensor.R): at accuracy 2, on a route
# that differences the same quantity, the adapters evaluate the same nodes
# with the same steps.

twin_fd_along <- function(f, eta, k, order = 1L, h = NULL) {
  s <- numericals7::fd_offsets(order, accuracy = 2L)$central
  w <- numericals7::fd_weights(s, order)
  if (is.null(h)) h <- fd_step(eta[k], order)
  acc <- NULL
  for (j in seq_along(s)) {
    if (w[j] == 0) next
    e <- eta
    e[k] <- eta[k] + s[j] * h
    v <- w[j] * f(e)
    acc <- if (is.null(acc)) v else acc + v
  }
  acc / h^order
}


twin_numerical_d1 <- function(s, eta) {
  out <- vector("list", s@n_free)
  names(out) <- s@free_names
  for (k in seq_len(s@n_free)) {
    out[[k]] <- twin_fd_along(function(e) param_value(s, e), eta, k, 1L)
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[k]] <- (out[[k]] + t(out[[k]])) / 2
    }
  }
  out
}


twin_numerical_d2 <- function(s, eta) {
  idx <- param_tuple_indices(s)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s)
  if (!length(idx)) return(out)

  analytic <- !param_is_numerical(s)[["param_d1"]]

  for (i in seq_along(idx)) {
    k <- idx[[i]][1L]
    l <- idx[[i]][2L]
    if (analytic) {
      # One layer on the analytic first derivative in the other component.
      out[[i]] <- twin_fd_along(function(e) param_d1(s, e)[[k]], eta, l, 1L)
    } else if (k == l) {
      out[[i]] <- twin_fd_along(function(e) param_value(s, e), eta, k, 2L)
    } else {
      # Two stencils in two DIFFERENT components, which is one mixed stencil
      # and not a difference of a difference: nesting is forbidden along one
      # variable and is what a mixed derivative is along two. The steps are
      # the order-2 ones, as a second derivative asks for.
      hk <- fd_step(eta[k], 2L)
      hl <- fd_step(eta[l], 2L)
      out[[i]] <- twin_fd_along(function(e) {
        twin_fd_along(function(e2) param_value(s, e2), e, l, 1L, h = hl)
      }, eta, k, 1L, h = hk)
    }
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


twin_numerical_d3 <- function(s, eta) {
  idx <- param_tuple_indices(s, 3L)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, 3L)
  f <- function(e) param_value(s, e)
  for (i in seq_along(idx)) {
    out[[i]] <- twin_mixed_stencil(f, eta, idx[[i]])
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


twin_numerical_d4 <- function(s, eta) {
  idx <- param_tuple_indices(s, 4L)
  out <- vector("list", length(idx))
  names(out) <- param_tuple_names(s, 4L)
  f <- function(e) param_value(s, e)
  for (i in seq_along(idx)) {
    out[[i]] <- twin_mixed_stencil(f, eta, idx[[i]])
    if (S7::S7_inherits(s, matrix_parameter)) {
      out[[i]] <- (out[[i]] + t(out[[i]])) / 2
    }
  }
  out
}


twin_mixed_stencil <- function(f, eta, tuple) {
  counts <- table(tuple)
  ks <- as.integer(names(counts))
  ms <- as.integer(counts)
  hs <- vapply(seq_along(ks), function(j) {
    fd_step(eta[ks[j]], sum(ms))
  }, numeric(1))

  one_dim <- function(m) {
    s <- numericals7::fd_offsets(m, accuracy = 2L)$central
    w <- numericals7::fd_weights(s, m)
    keep <- w != 0
    list(offsets = s[keep], weights = w[keep], power = m)
  }
  facs <- lapply(ms, one_dim)

  combos <- expand.grid(lapply(facs, function(fc) seq_along(fc$offsets)))
  acc <- NULL
  for (r in seq_len(nrow(combos))) {
    e <- eta
    w <- 1
    for (j in seq_along(ks)) {
      pick <- combos[r, j]
      e[ks[j]] <- e[ks[j]] + facs[[j]]$offsets[pick] * hs[j]
      w <- w * facs[[j]]$weights[pick]
    }
    term <- w * f(e)
    acc <- if (is.null(acc)) term else acc + term
  }
  denom <- prod(vapply(seq_along(ks), function(j) {
    hs[j]^facs[[j]]$power
  }, numeric(1)))
  acc / denom
}


twin_d3logdet <- function(s, eta) {
  # (0.24.0 called check_analytic_arrays() here)
  idx3 <- param_tuple_indices(s, 3L)
  idx2 <- param_tuple_indices(s, 2L)
  key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
  out <- vapply(seq_along(idx3), function(i) {
    t3 <- idx3[[i]]
    k <- t3[3L]
    pos <- match(paste(sort(t3[1:2]), collapse = ","), key2)
    h <- fd_step(eta[k], 1L)
    up <- eta
    dn <- eta
    up[k] <- eta[k] + h
    dn[k] <- eta[k] - h
    (param_d2logdet(s, up)[[pos]] - param_d2logdet(s, dn)[[pos]]) / (2 * h)
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, 3L))
}


twin_d4logdet <- function(s, eta) {
  # (0.24.0 called check_analytic_arrays() here)
  idx4 <- param_tuple_indices(s, 4L)
  idx2 <- param_tuple_indices(s, 2L)
  key2 <- vapply(idx2, function(t) paste(sort(t), collapse = ","), character(1))
  out <- vapply(seq_along(idx4), function(i) {
    t4 <- idx4[[i]]
    pos <- match(paste(sort(t4[1:2]), collapse = ","), key2)
    k <- t4[3L]
    l <- t4[4L]
    g <- function(e) param_d2logdet(s, e)[[pos]]
    if (k == l) {
      h <- fd_step(eta[k], 2L)
      up <- eta
      dn <- eta
      up[k] <- eta[k] + h
      dn[k] <- eta[k] - h
      (g(up) - 2 * g(eta) + g(dn)) / h^2
    } else {
      hk <- fd_step(eta[k], 2L)
      hl <- fd_step(eta[l], 2L)
      pp <- pm <- mp <- mm <- eta
      pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
      pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
      mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
      mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
      (g(pp) - g(pm) - g(mp) + g(mm)) / (4 * hk * hl)
    }
  }, numeric(1))
  stats::setNames(out, param_tuple_names(s, 4L))
}
