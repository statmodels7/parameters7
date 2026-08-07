# A Factor of a Parameter's Matrix

Returns a lower triangular \\L\\ with \\M = L L^\top\\.

## Usage

``` r
param_factor(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A lower triangular numeric matrix with `s@dimension` rows and columns.

## Details

Rejected for a rank-deficient family, for the reason given in
[`param_solve`](https://statmodels7.github.io/parameters7/reference/param_solve.md).

## See also

[`param_solve`](https://statmodels7.github.io/parameters7/reference/param_solve.md)

## Examples

``` r
round(param_factor(log_cholesky(2), c(0.1, 0.2, -0.3)), 4)
#>         [,1]   [,2]
#> [1,]  1.1052 0.0000
#> [2,] -0.3000 1.2214
```
