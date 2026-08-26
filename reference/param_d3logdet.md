# Third Derivatives of the Log-Determinant

Returns the distinct third derivatives of the log-determinant, or of the
log pseudo-determinant, keyed as `param_tuple_names(s, 3)`. The exact
gradient of a marginal criterion reaches this order: differentiating a
Laplace approximation with respect to a hyperparameter moves the
penalized mode, and the third derivative is what the movement contracts
against.

## Usage

``` r
param_d3logdet(s, eta, ...)
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

A named numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
`param_tuple_names(s, 3)` and in that order.

## Where it comes from

Two rules generate every order. The trace is linear, and

\$\$\partial_m M^{-1} = -M^{-1}(\partial_m M)M^{-1}.\$\$

Applying them to \\\partial\_{kl}\log\|M\|\\ gives a sum of traces of
alternating products
\\M^{-1}(\partial\_{I_1}M)M^{-1}(\partial\_{I_2}M)\cdots\\, one factor
per block of a partition of the index set, with the sign and the
multiplicity the two rules produce. The Moore-Penrose inverse replaces
\\M^{-1}\\ for a rank-deficient family.

The expansion is not transcribed anywhere. The methods differentiate the
derivative arrays
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
through
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
instead, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
holds the result against a numerical differentiation of the order below,
which shares none of its arithmetic. A twenty-term expansion written out
by hand is exactly the kind of thing that is wrong and looks right.

## Notation

\\M\\ is the matrix
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns and \\\eta\\ the free vector, of length \\d = \\ `s@n_free`.
\\\partial_I M\\ is the derivative of \\M\\ in the free values named by
the index set \\I\\.

## See also

[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
for the order below,
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for the order above,
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
for the keys, and
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
for the third derivative of the matrix itself.

## Examples

``` r
# An AR(1) covariance. The log-determinant is linear in the log of the
# scale, so only the correlation direction survives to third order.
s <- ar1(4)
param_d3logdet(s, c(0.3, 0.8))
#> log_scale:log_scale:log_scale     log_scale:log_scale:z_rho 
#>                  8.881784e-16                  0.000000e+00 
#>         log_scale:z_rho:z_rho             z_rho:z_rho:z_rho 
#>                  0.000000e+00                  4.454798e+00 

# A central difference of the order below agrees, and shares no arithmetic
# with the route the method takes.
h <- 1e-4
e <- c(0, h)
fd <- (param_d2logdet(s, c(0.3, 0.8) + e) -
         param_d2logdet(s, c(0.3, 0.8) - e)) / (2 * h)
fd[["z_rho:z_rho"]]
#> [1] 4.454798
param_d3logdet(s, c(0.3, 0.8))[["z_rho:z_rho:z_rho"]]
#> [1] 4.454798

# For log_cholesky the log-determinant is linear in eta, so every order
# above the first is exactly zero.
max(abs(param_d3logdet(log_cholesky(3), rep(0.2, 6))))
#> [1] 0
```
