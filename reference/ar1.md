# Construct an AR(1) Parameter

Returns an object holding the first-order autoregressive covariance

\$\$M(\eta)\_{ij} = \sigma^2 \rho^{\|i - j\|},\$\$

with \\\sigma^2\\ positive and \\\rho \in (-1, 1)\\: two free values at
every dimension. Use it when the measurements are ordered and the
correlation should fall with the distance between them, which is a time
series, a spatial transect, or repeated measures where an exchangeable
correlation is too strong an assumption.

## Usage

``` r
ar1(dimension, link_scale = linkfunctions7::log_link())
```

## Arguments

- dimension:

  The side \\p\\ of the matrix, **at least 2**. A one by one matrix has
  no correlation, so `ar1(1)` throws a message saying the family would
  carry a free value with no effect.

- link_scale:

  A linkfunctions7 link carrying the first free value onto the positive
  variance,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. It must map onto the positive half line, so
  `identity_link()` is rejected, and from the whole real line, which
  rules out `sqrt_link()` and its relatives; see
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  for the two conditions.

## Value

An object of class
[`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md),
with `n_free` 2, `free_names` `log_scale` and `z_rho` under the
defaults, `rank` equal to `dimension`, an empty `null_basis`,
`param_name` `"ar1"`, and `param_params` holding `link_scale` and
`link_rho`.

## The correlation is unrestricted

Unlike
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
the bound is \\\|\rho\| \< 1\\ at every dimension, so the correlation
rides
[`linkfunctions7::rhobit_link()`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.html),
the inverse hyperbolic tangent, and every free value gives a positive
definite matrix. A free value of 0 is a correlation of 0 here, where in
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
it is the midpoint of a dimension-dependent interval.

## The inverse is tridiagonal

This is the property the family is used for. An AR(1) process is Markov,
so its precision has no entries beyond the first off-diagonal:

\$\$M^{-1} = \frac{1}{\sigma^2(1-\rho^2)} \begin{pmatrix} 1 & -\rho & &
\\ -\rho & 1+\rho^2 & -\rho & \\ & \ddots & \ddots & \ddots \\ & & -\rho
& 1 \end{pmatrix}.\$\$

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
returns it from that form, so no factorization is performed and the cost
is \\O(p)\\ entries against the base class's \\O(p^3)\\.

Note what this does **not** say: the precision is tridiagonal without
being AR(1). Its own correlation at lag 1 is not constant along the
diagonal (measured at \\p = 4\\, \\\rho = 0.6\\: \\-0.514\\ at the ends
and \\-0.441\\ in the middle), and its lag-2 correlation is 0 where an
AR(1) would have \\\rho^2\\. So an AR(1) covariance and an AR(1)
precision are different models. This family is the AR(1) pattern itself,
whichever side a consumer puts it on;
[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
is the family whose value is the matrix above, so that the AR(1) process
is the one written on the other side.

## The log-determinant

The determinant of the correlation pattern is \\(1-\rho^2)^{p-1}\\, so

\$\$\log\|M\| = p\log\sigma^2 + (p-1)\log(1-\rho^2),\$\$

a sum of a function of one free value and a function of the other. Every
mixed derivative of it is exactly zero, and no determinant is computed.

## Why the derivatives are harder than compound symmetry's

The pattern is **not** linear in the correlation. An entry is
\\\rho^{m}\\ for the lag \\m\\, so its derivatives in the free value are
a power composed with the link, taken to fourth order by
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md).
Each distinct lag is composed once and written into every entry that
carries it, the matrix having only \\p\\ distinct values.

## Notation

\\\eta = (\eta_1, \eta_2)\\ is the free vector, \\\sigma^2\\ the common
variance and \\\rho\\ the lag-one correlation. \\p\\ is the side of the
matrix and \\m = \|i-j\|\\ the lag of an entry.

## See also

[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
for an exchangeable correlation,
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for orders above one,
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
for an unrestricted correlation, and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for the tridiagonal inverse.

## Examples

``` r
# A 4 x 4 AR(1) covariance at a variance of 2 and a correlation of 0.6.
s <- ar1(4)
eta <- c(log(2), atanh(0.6))
m <- param_value(s, eta)
round(m, 4)
#>       v1   v2   v3    v4
#> v1 2.000 1.20 0.72 0.432
#> v2 1.200 2.00 1.20 0.720
#> v3 0.720 1.20 2.00 1.200
#> v4 0.432 0.72 1.20 2.000

# The entries really are sigma^2 rho^|i-j|.
all.equal(m[1, 4], 2 * 0.6^3, check.attributes = FALSE)
#> [1] TRUE

# A free value of 0 is a correlation of 0, the link being rhobit.
param_value(s, c(0, 0))[1, 2]
#> [1] 0

# The inverse is tridiagonal, and exact: an AR(1) process is Markov.
om <- param_solve(s, eta)
round(om, 4)
#>         [,1]    [,2]    [,3]    [,4]
#> [1,]  0.7812 -0.4688  0.0000  0.0000
#> [2,] -0.4688  1.0625 -0.4688  0.0000
#> [3,]  0.0000 -0.4688  1.0625 -0.4688
#> [4,]  0.0000  0.0000 -0.4688  0.7812
max(abs(om[abs(row(om) - col(om)) > 1]))
#> [1] 0
max(abs(om - solve(m)))
#> [1] 2.220446e-16

# But it is not itself AR(1): its lag-1 correlation varies along the
# diagonal, and its lag-2 correlation is 0 where rho^2 would be 0.36.
d <- sqrt(diag(om))
round(om / outer(d, d), 4)
#>         [,1]    [,2]    [,3]    [,4]
#> [1,]  1.0000 -0.5145  0.0000  0.0000
#> [2,] -0.5145  1.0000 -0.4412  0.0000
#> [3,]  0.0000 -0.4412  1.0000 -0.5145
#> [4,]  0.0000  0.0000 -0.5145  1.0000

# The log-determinant is closed form and agrees with the eigenvalues.
c(closed = param_logdet(s, eta),
  from_eigen = sum(log(eigen(m, only.values = TRUE)$values)))
#>     closed from_eigen 
#>   1.433727   1.433727 

# The round trip closes exactly, and every mixed second derivative of the
# log-determinant is zero.
max(abs(param_free(s, m) - eta))
#> [1] 0
param_d2logdet(s, eta)
#> log_scale:log_scale         z_rho:z_rho     log_scale:z_rho 
#>                0.00               -3.84                0.00 
```
