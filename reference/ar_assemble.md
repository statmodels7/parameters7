# The Matrix and Its Derivatives, From the Jets

Fills the Toeplitz matrix of the scaled autocorrelations, taking either
the value or one derivative component from each jet.

## Usage

``` r
ar_assemble(s, j, order = 0L, position = 1L)
```

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- j:

  The jets of
  [`ar_jets`](https://statmodels7.github.io/parameters7/reference/ar_jets.md).

- order:

  The derivative order, or 0 for the value.

- position:

  Which component of that order.

## Value

A symmetric numeric matrix.
