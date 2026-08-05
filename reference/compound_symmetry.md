# Construct a Compound Symmetry Parameter

The exchangeable covariance \$\$M(\eta) = \sigma^2\\(1-\rho)I + \rho
J\\,\$\$ with \\\sigma^2\\ positive and \\\rho\\ the common correlation,
so two free values whatever the dimension.

## Usage

``` r
compound_symmetry(
  dimension,
  link_scale = linkfunctions7::log_link(),
  role = c("either", "covariance", "precision")
)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix, at least 2.

- link_scale:

  A linkfunctions7 link onto the positive scale. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Value

An object of class
[`CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md).

## Details

Positive definiteness bounds the correlation below as well as above: the
eigenvalues are \\\sigma^2\\1 + (p-1)\rho\\\\ once and
\\\sigma^2(1-\rho)\\ with multiplicity \\p-1\\, so the matrix is
definite exactly when \\-1/(p-1) \< \rho \< 1\\. The correlation is
therefore carried by
[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.html)
onto that interval and not by
[`rhobit_link`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.html),
which would let a caller build an indefinite matrix at a perfectly
ordinary free value.

Two quantities are closed form and cost nothing. The eigenvalues above
give \$\$\log\lvert M \rvert = p\log\sigma^2 + \log\\1 + (p-1)\rho\\ +
(p-1)\log(1-\rho),\$\$ a sum of a function of one free value and a
function of the other, so every mixed derivative of the log-determinant
is exactly zero. And the inverse is compound symmetric again, by the
Sherman-Morrison identity, so it is returned in closed form rather than
factorized: the inverse of an exchangeable covariance is an exchangeable
precision, which is what makes this family closed under the choice of
side.

The value is the scale times a pattern that is *linear* in the
correlation, \\I + \rho(J - I)\\, so all four derivative orders follow
from the two links' own derivatives with no further algebra.

## See also

[`ar1`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`correlation_matrix`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)

## Examples

``` r
s <- compound_symmetry(4)
c(n_free = s@n_free)
#> n_free 
#>      2 

eta <- c(log(2), 0.8)
round(param_value(s, eta), 4)
#>        v1     v2     v3     v4
#> v1 2.0000 1.1733 1.1733 1.1733
#> v2 1.1733 2.0000 1.1733 1.1733
#> v3 1.1733 1.1733 2.0000 1.1733
#> v4 1.1733 1.1733 1.1733 2.0000

# the inverse is compound symmetric too, and is not factorized
round(param_solve(s, eta), 4)
#>         [,1]    [,2]    [,3]    [,4]
#> [1,]  0.9525 -0.2571 -0.2571 -0.2571
#> [2,] -0.2571  0.9525 -0.2571 -0.2571
#> [3,] -0.2571 -0.2571  0.9525 -0.2571
#> [4,] -0.2571 -0.2571 -0.2571  0.9525

# the round trip closes exactly
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 3.330669e-16
```
