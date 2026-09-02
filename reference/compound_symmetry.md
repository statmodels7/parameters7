# Construct a Compound Symmetry Parameter

Returns an object holding the exchangeable covariance

\$\$M(\eta) = \sigma^2\\(1-\rho)I + \rho J\\,\$\$

with \\\sigma^2\\ positive and \\\rho\\ the correlation every pair
shares. Two free values at every dimension, against
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)'s
\\p(p+1)/2\\, so it is the parametrization to use when the measurements
are exchangeable: repeated measures on a subject, items in a block, the
covariance a random intercept induces.

## Usage

``` r
compound_symmetry(dimension, link_scale = linkfunctions7::log_link())
```

## Arguments

- dimension:

  The side \\p\\ of the matrix, **at least 2**. A one by one matrix has
  no correlation, so `compound_symmetry(1)` throws a message saying the
  family would carry a free value with no effect.

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
[`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md),
with `n_free` 2, `free_names` `log_scale` and `logit_rho` under the
defaults, `rank` equal to `dimension`, an empty `null_basis`,
`param_name` `"compound_symmetry"`, and `param_params` holding
`link_scale` and `link_rho`.

## The correlation is bounded below as well as above

The eigenvalues are \\\sigma^2\\1 + (p-1)\rho\\\\ once and
\\\sigma^2(1-\rho)\\ with multiplicity \\p-1\\, so the matrix is
positive definite exactly when

\$\$-\frac{1}{p-1} \< \rho \< 1.\$\$

The correlation therefore rides
`linkfunctions7::bounded_link(-1/(p-1), 1)`. A `rhobit_link()` onto
\\(-1, 1)\\ would let a caller build an indefinite matrix at a perfectly
ordinary free value: at \\p = 4\\ a correlation of \\-0.5\\ is inside
\\(-1, 1)\\ and outside the cone. The bound depends on the dimension, so
two objects of different sizes carry different links.

One consequence to expect: a free value of 0 is the **midpoint** of that
interval, so at \\p = 4\\ it is a correlation of \\1/3\\, not of 0.

## Two quantities that cost nothing

From the eigenvalues,

\$\$\log\|M\| = p\log\sigma^2 + \log\\1 + (p-1)\rho\\ +
(p-1)\log(1-\rho),\$\$

a sum of a function of one free value and a function of the other, so
every mixed derivative of the log-determinant is exactly zero and no
determinant is computed.

The inverse is compound symmetric again, by Sherman-Morrison, so
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
returns it in closed form and factorizes nothing. An exchangeable
covariance has an exchangeable precision, which leaves the family closed
under the choice of side.

## Why the derivatives are easy

The value is the scale times a pattern **linear** in the correlation,
\\I + \rho(J - I)\\, so a component with \\a\\ scale indices and \\b\\
correlation indices is the \\a\\-th derivative of the scale times the
\\b\\-th derivative of the pattern. All four orders follow from the two
links' own derivatives with no further algebra, and any component with
three or more correlation indices reduces to the third derivative of the
link alone.

## Notation

\\\eta = (\eta_1, \eta_2)\\ is the free vector, \\\sigma^2\\ the common
variance and \\\rho\\ the common correlation. \\p\\ is the side of the
matrix, \\I\\ the identity and \\J\\ the matrix of ones.

## See also

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
the other two-value family, where the correlation decays with distance
instead of being shared;
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
for an unrestricted correlation;
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
to replicate this over independent groups; and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for the closed inverse.

## Examples

``` r
# Two free values at p = 4: a variance and one correlation.
s <- compound_symmetry(4)
c(n_free = s@n_free)
#> n_free 
#>      2 
s@free_names
#> [1] "log_scale" "logit_rho"

eta <- c(log(2), 0.8)
m <- param_value(s, eta)
round(m, 4)
#>        v1     v2     v3     v4
#> v1 2.0000 1.1733 1.1733 1.1733
#> v2 1.1733 2.0000 1.1733 1.1733
#> v3 1.1733 1.1733 2.0000 1.1733
#> v4 1.1733 1.1733 1.1733 2.0000

# The eigenvalues are the two the bound comes from.
sig2 <- m[1, 1]
rho <- m[1, 2] / sig2
round(eigen(m, only.values = TRUE)$values, 6)
#> [1] 5.519796 0.826735 0.826735 0.826735
round(c(sig2 * (1 + 3 * rho), rep(sig2 * (1 - rho), 3)), 6)
#> [1] 5.519796 0.826735 0.826735 0.826735

# The log-determinant is closed form and agrees with them.
c(closed = param_logdet(s, eta),
  from_eigen = sum(log(eigen(m, only.values = TRUE)$values)))
#>     closed from_eigen 
#>   1.137527   1.137527 

# The inverse is compound symmetric too, and is not factorized.
round(param_solve(s, eta), 4)
#>         [,1]    [,2]    [,3]    [,4]
#> [1,]  0.9525 -0.2571 -0.2571 -0.2571
#> [2,] -0.2571  0.9525 -0.2571 -0.2571
#> [3,] -0.2571 -0.2571  0.9525 -0.2571
#> [4,] -0.2571 -0.2571 -0.2571  0.9525
max(abs(param_solve(s, eta) - solve(m)))
#> [1] 1.110223e-16

# The round trip closes exactly.
max(abs(param_free(s, m) - eta))
#> [1] 3.330669e-16

# A free value of 0 is the midpoint of (-1/(p-1), 1), not a correlation of 0.
param_value(s, c(0, 0))[1, 2]
#> [1] 0.3333333

# The bound tightens with the dimension.
t(vapply(c(2, 4, 10),
         function(p) compound_symmetry(p)@param_params$link_rho@link_bounds,
         numeric(2)))
#>            [,1] [,2]
#> [1,] -1.0000000    1
#> [2,] -0.3333333    1
#> [3,] -0.1111111    1
```
