# Higher Derivatives of a Diagonal Log-Determinant

The derivative components of orders two to four of \\\log\lvert
M\rvert\\ for a diagonal family.

## Usage

``` r
diag_logdet_higher(s, eta, order)
```

## Arguments

- s:

  A
  [`DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The derivative order, 2 to 4.

## Value

A named numeric vector, keyed by
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, order)`.

## Details

The log-determinant is the sum of the logarithms of the diagonal
entries, so it is a sum of functions of one free value each. Every mixed
component therefore vanishes, and a pure one is the corresponding
derivative of
[`diag_dlog`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
counted once per entry the free value owns.

## See also

[`diagonal_matrix`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`param_logdet`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
