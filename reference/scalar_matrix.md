# Construct a Scalar Multiple of the Identity

A diagonal matrix with one free value shared by every entry, so \\M =
\tau I\\ with \\\tau\\ positive.

## Usage

``` r
scalar_matrix(
  dimension,
  link = linkfunctions7::log_link(),
  role = c("either", "covariance", "precision")
)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix.

- link:

  A linkfunctions7 link mapping the free scale to the positive scale.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Value

An object of class
[`DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md).

## Details

The simplest parameter there is, and the one a random effect with a
single variance component uses. With the default log link it is the same
matrix as `scaled_matrix(diag(dimension))`, reached from the other
direction.

## See also

[`diagonal_matrix`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scaled_matrix`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)

## Examples

``` r
s <- scalar_matrix(3)
c(n_free = s@n_free)
#> n_free 
#>      1 
round(param_value(s, log(2)), 4)
#>    v1 v2 v3
#> v1  2  0  0
#> v2  0  2  0
#> v3  0  0  2
```
