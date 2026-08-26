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
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` symmetric matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension`.

## Details

The Frechet derivatives of the exponential contract chains of directions
against divided differences of \\\exp\\ in the eigenvalues, and the sum
runs over every ordering of the directions, never over the distinct
ones, so a tuple with repeated indices is counted with its multiplicity
instead of being corrected for afterwards. The divided differences come
from the Opitz representation, an exponential of a small bidiagonal
matrix read off its corner, which stays exact where the quotient
recursion cancels catastrophically under near-repeated eigenvalues.

## See also

[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md),
which computes one component,
[`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md)
for the shared tables, and
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
for the cost this order carries.
