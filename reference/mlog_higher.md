# Third and Fourth Derivatives of a Matrix Exponential

The derivative components of orders three and four of \\M = e^{S}\\ with
respect to the free entries of the symmetric matrix \\S\\.

## Usage

``` r
mlog_higher(s, eta, order)
```

## Arguments

- s:

  A
  [`MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The derivative order, 3 or 4.

## Value

A named list of matrices, keyed by
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, order)`.

## Details

The Frechet derivatives of the exponential contract chains of directions
against divided differences of \\\exp\\ in the eigenvalues, and the sum
runs over every ordering of the directions rather than over the distinct
ones, so a tuple with repeated indices is counted with its multiplicity
rather than corrected for it afterwards. The divided differences come
from the Opitz representation, an exponential of a small bidiagonal
matrix read off its corner, which stays exact where the quotient
recursion cancels catastrophically under near-repeated eigenvalues.

## See also

[`matrix_log`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
