# Hessian of the Log-Determinant

Returns the distinct second derivatives of the log-determinant, or of
the log pseudo-determinant, keyed as
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md).
A Newton step in the free vector of a Gaussian model needs them, and so
does the observed information of a covariance parametrized this way.

## Usage

``` r
param_d2logdet(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- ...:

  Passed to the method. The methods in this package do not read it.

## Value

A named numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
`param_tuple_names(s)` and in that order: the `s@n_free` diagonal pairs
first, then the off-diagonal ones. The \\(l, k)\\ entry equals the \\(k,
l)\\ one and is not repeated.

## The identity it must satisfy

Differentiating \\\partial_k \log\|M\| = \mathrm{tr}(M^{-1}\partial_k
M)\\ once more, with \\\partial_l M^{-1} = -M^{-1}(\partial_l
M)M^{-1}\\,

\$\$\partial\_{kl} \log\|M\| =
\mathrm{tr}\\\left(M^{-1}\partial\_{kl}M\right) -
\mathrm{tr}\\\left(M^{-1}(\partial_k M) M^{-1} (\partial_l
M)\right),\$\$

and with the Moore-Penrose inverse in place of \\M^{-1}\\ when the
family is rank deficient.

Without the second trace the result would be the trace of the second
derivative of the matrix, which is a different quantity and is not zero
even for
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
whose log-determinant is linear.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
compares this method with one central difference of
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
and the example below compares it with the two traces.

## Notation

\\M\\ is the matrix
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns and \\\eta\\ the free vector, of length \\d = \\ `s@n_free`.

## See also

[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
for the order below,
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for the orders above,
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
for the keys, and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
which this must agree with through the two traces.

## Examples

``` r
# The identity, checked against param_d1() and param_d2() directly.
s <- ar1(4)
eta <- c(0.3, 0.8)
Minv <- solve(param_value(s, eta))
d1 <- param_d1(s, eta)
d2 <- param_d2(s, eta)
idx <- param_tuple_indices(s, 2)
pred <- vapply(seq_along(idx), function(i) {
  k <- idx[[i]][1]; l <- idx[[i]][2]
  sum(diag(Minv %*% d2[[i]])) -
    sum(diag(Minv %*% d1[[k]] %*% Minv %*% d1[[l]]))
}, numeric(1))
max(abs(param_d2logdet(s, eta) - pred))
#> [1] 1.776357e-15

# Dropping the second trace would give a different answer here, so the
# comparison above has something to catch.
vapply(seq_along(idx), function(i) sum(diag(Minv %*% d2[[i]])), numeric(1))
#> [1]  4.000000  5.291338 -3.984221
param_d2logdet(s, eta)
#> log_scale:log_scale         z_rho:z_rho     log_scale:z_rho 
#>            0.000000           -3.354331            0.000000 

# For log_cholesky the log-determinant is linear in eta, so every second
# derivative is exactly zero.
max(abs(param_d2logdet(log_cholesky(3), rep(0.2, 6))))
#> [1] 0
```
