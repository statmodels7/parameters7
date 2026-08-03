# A Factor of a Structure's Matrix

Returns a lower triangular \\L\\ with \\M = L L^\top\\.

## Usage

``` r
struct_factor(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A lower triangular numeric matrix with `s@dimension` rows and columns.

## Details

Refused for a rank-deficient family, for the reason given in
[`struct_solve`](https://statmodels7.github.io/covstructs7/reference/struct_solve.md).

## See also

[`struct_solve`](https://statmodels7.github.io/covstructs7/reference/struct_solve.md)

## Examples

``` r
round(struct_factor(log_cholesky(2), c(0.1, 0.2, -0.3)), 4)
#>         [,1]   [,2]
#> [1,]  1.1052 0.0000
#> [2,] -0.3000 1.2214
```
