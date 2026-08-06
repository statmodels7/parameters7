# The Matrix and Its Derivatives, From the Packed Arrays

Fills the Toeplitz matrix of the scaled autocorrelations, taking either
the value column or one derivative component out of the packed rows of
[`ar_taylor`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md).

## Usage

``` r
ar_assemble(s, tay, order = 0L, tuple = NULL)
```

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- tay:

  The arrays of
  [`ar_taylor`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md).

- order:

  The derivative order, or 0 for the value.

- tuple:

  The index tuple of that order, ignored at order 0.

## Value

A symmetric numeric matrix.
