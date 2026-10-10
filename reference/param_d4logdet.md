# Fourth Derivatives of the Log-Determinant

Returns the distinct fourth derivatives of the log-determinant, or of
the log pseudo-determinant, keyed as `param_tuple_names(s, 4)`. This is
the highest order that the package provides, and the exact Hessian of a
marginal criterion reaches it.

## Usage

``` r
param_d4logdet(s, eta, ...)
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

A named numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
`param_tuple_names(s, 4)` and in that order.

## Where it comes from

The same two rules as at third order. The trace is linear, and

\$\$\partial_m M^{-1} = -M^{-1}(\partial_m M)M^{-1},\$\$

so every term is a trace of an alternating product
\\M^{-1}(\partial\_{I_1}M)M^{-1}(\partial\_{I_2}M)\cdots\\, one factor
per block of a partition of the four indices. The number of terms grows
with the number of partitions, so the expansion is not written out: each
family in this package computes this order without differencing, and the
base method on
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
takes one central difference of
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md).

## Checking this order

One central difference of
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md),
as in the example below, is a reference for this order, and the
comparison is limited by that stencil rather than by the closed form.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
does not run this comparison. A family whose log-determinant is linear
in \\\eta\\, which includes
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
and
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
returns exact zeros here;
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
and
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
are among the families where this order is not zero.

## Notation

\\M\\ is the matrix
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns and \\\eta\\ the free vector, of length \\d = \\ `s@n_free`.
\\\partial_I M\\ is the derivative of \\M\\ in the free values named by
the index set \\I\\.

## See also

[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
for the order below,
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
for the keys, and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for the fourth derivative of the matrix itself.

## Examples

``` r
# An AR(1) covariance: five components at d = 2, of which only the one
# purely in the correlation is non-zero, the scale entering linearly.
s <- ar1(4)
eta <- c(0.3, 0.8)
param_d4logdet(s, eta)
#> log_scale:log_scale:log_scale:log_scale     log_scale:log_scale:log_scale:z_rho 
#>                            7.105427e-15                            0.000000e+00 
#>         log_scale:log_scale:z_rho:z_rho             log_scale:z_rho:z_rho:z_rho 
#>                            0.000000e+00                            0.000000e+00 
#>                 z_rho:z_rho:z_rho:z_rho 
#>                           -2.165788e+00 

# A central difference of the third order agrees.
h <- 1e-4
e <- c(0, h)
fd <- (param_d3logdet(s, eta + e) - param_d3logdet(s, eta - e)) / (2 * h)
c(analytic = param_d4logdet(s, eta)[["z_rho:z_rho:z_rho:z_rho"]],
  difference = fd[["z_rho:z_rho:z_rho"]])
#>   analytic difference 
#>  -2.165788  -2.165787 

# The count is over unordered quadruples.
c(returned = length(param_d4logdet(s, eta)), choose(2 + 3, 4))
#> returned          
#>        5        5 
```
