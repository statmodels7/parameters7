# The Blocks' Derivative Components, Keyed by Local Index Tuple

Fetches one block's derivatives of a given order and re-keys them by the
sorted local index tuple, so that a lookup from the composite's
enumeration needs no assumption about how the two orderings correspond.

## Usage

``` r
block_derivs_by_tuple(block, eta, order)
```

## Arguments

- block:

  A
  [`matrix_parameter`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).

- eta:

  The block's own stretch of the free vector.

- order:

  The derivative order, 1 to 4.

## Value

A named list of matrices, keyed as `"1"`, `"1,2"` and so on over the
block's own free indices.
