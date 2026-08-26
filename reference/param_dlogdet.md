# Gradient of the Log-Determinant

Returns \\\partial \log\|M\| / \partial \eta_k\\ for every free value,
or the same derivative of the log pseudo-determinant for a
rank-deficient family. A Gaussian score in the coordinates an optimizer
moves is the derivative of the quadratic form plus this, so it is needed
at every iteration of a fit.

## Usage

``` r
param_dlogdet(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## The identity it must satisfy

\$\$\partial_k \log\|M\| = \mathrm{tr}\\\left(M^{-1} \partial_k
M\right),\$\$

with the Moore-Penrose inverse in place of \\M^{-1}\\ when the family is
rank deficient. A closed form here is therefore never an independent
claim: it has to agree with
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
through that trace, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs both routes and compares them. The example below does the same in
three lines.

## What the answer looks like

It is often constant. In the log-Cholesky parametrization the
log-determinant is \\2\sum\_{i\le p}\eta_i\\, so the gradient is 2 in
the \\p\\ diagonal directions and 0 in the rest, at every \\\eta\\. For
a
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
where the free value is the log of a multiplier on a fixed matrix, the
gradient is the rank, again at every \\\eta\\: the scale enters the
pseudo-determinant once per non-zero eigenvalue.

## Notation

\\M\\ is the matrix
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns and \\\eta\\ the free vector, of length \\d = \\ `s@n_free`.
\\p\\ is the side of the matrix.

## See also

[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
for the quantity differentiated here,
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
for the next order, and
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
which this must agree with through the trace identity.

## Examples

``` r
# The trace identity, checked against param_d1() directly.
s <- log_cholesky(3)
eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
Minv <- solve(param_value(s, eta))
tr <- vapply(param_d1(s, eta), function(A) sum(diag(Minv %*% A)), numeric(1))
max(abs(param_dlogdet(s, eta) - tr))
#> [1] 4.440892e-16

# In this parametrization the answer is 2 on the diagonal directions and 0
# elsewhere, whatever eta is.
param_dlogdet(s, eta)
#> log_L1 log_L2 log_L3   L2.1   L3.1   L3.2 
#>      2      2      2      0      0      0 

# For a scaled precision it is the rank, whatever the scale.
r <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
c(at_minus_3 = param_dlogdet(r, -3), at_5 = param_dlogdet(r, 5),
  rank = r@rank)
#> at_minus_3.log_scale       at_5.log_scale                 rank 
#>                    4                    4                    4 

# An AR(1) covariance, where the correlation direction is not constant.
param_dlogdet(ar1(4), c(0.3, 0.8))
#> log_scale     z_rho 
#>  4.000000 -3.984221 
```
