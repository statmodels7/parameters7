# Construct an AR(1) Parameter

The first-order autoregressive covariance \$\$M(\eta)\_{ij} = \sigma^2
\rho^{\lvert i - j \rvert},\$\$ with \\\sigma^2\\ positive and \\\rho
\in (-1, 1)\\, so two free values whatever the dimension.

## Usage

``` r
ar1(
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
[`Ar1Param`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md).

## Details

Unlike compound symmetry, the correlation is bounded only by
\\\lvert\rho\rvert \< 1\\ at every dimension, so it is carried by
[`rhobit_link`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.html)
– the inverse hyperbolic tangent – and every free value gives a positive
definite matrix.

Two quantities are closed form. The determinant of the correlation
pattern is \\(1-\rho^2)^{p-1}\\, so \$\$\log\lvert M \rvert =
p\log\sigma^2 + (p-1)\log(1-\rho^2),\$\$ again a sum of a function of
one free value and a function of the other, with every mixed derivative
exactly zero. And the inverse is **tridiagonal**, which is the property
the family is used for: an AR(1) process is Markov, so its precision has
no entries beyond the first off-diagonal, and the solve is returned from
that form rather than from a factorization.

The pattern is not linear in the correlation as compound symmetry's is:
an entry is \\\rho^{m}\\ for the lag \\m\\, so its derivatives in the
free value are the composition of a power with the link, taken to fourth
order by
[`compose4`](https://statmodels7.github.io/parameters7/reference/compose4.md).
This is also why an AR(1) covariance and an AR(1) precision are
different models – the inverse is tridiagonal but not AR(1) – which is
what the `role` label exists to record.

## See also

[`compound_symmetry`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`correlation_matrix`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)

## Examples

``` r
s <- ar1(4)
eta <- c(log(2), atanh(0.6))
round(param_value(s, eta), 4)
#>       v1   v2   v3    v4
#> v1 2.000 1.20 0.72 0.432
#> v2 1.200 2.00 1.20 0.720
#> v3 0.720 1.20 2.00 1.200
#> v4 0.432 0.72 1.20 2.000

# the inverse is tridiagonal: an AR(1) process is Markov
round(param_solve(s, eta), 4)
#>         [,1]    [,2]    [,3]    [,4]
#> [1,]  0.7812 -0.4688  0.0000  0.0000
#> [2,] -0.4688  1.0625 -0.4688  0.0000
#> [3,]  0.0000 -0.4688  1.0625 -0.4688
#> [4,]  0.0000  0.0000 -0.4688  0.7812

# two free values whatever the dimension
c(p4 = ar1(4)@n_free, p20 = ar1(20)@n_free)
#>  p4 p20 
#>   2   2 
```
