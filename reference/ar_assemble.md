# The Matrix and Its Derivatives, From the Packed Arrays

Fills the Toeplitz matrix \\M\_{ij} = \gamma\_{\lvert i - j \rvert}\\
from one column of
[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)'s
packed rows, taking either the value column or one derivative component.
Every derivative of the matrix is Toeplitz too, the Toeplitz structure
being a property of the family, fixed as the point moves, so one
indexing operation serves all five cases.

## Usage

``` r
ar_assemble(s, tay, order = 0L, tuple = NULL)
```

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object, whose `dimension` is read.

- tay:

  The arrays of
  [`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md).

- order:

  The derivative order 1 to 4, or 0 for the value.

- tuple:

  The index tuple of that order, ignored at order 0.

## Value

A symmetric `s@dimension` by `s@dimension` numeric matrix, with no
dimnames.

## See also

[`ar_pack_col()`](https://statmodels7.github.io/parameters7/reference/ar_pack_col.md),
which locates the column, and
[`ar_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_derivative.md),
which loops this over the tuples of an order.
