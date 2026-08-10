# Assemble a Block Diagonal's Derivatives of a Given Order

Places each block's own component in the rows and columns that block
occupies, and returns a zero matrix for any tuple whose indices are not
all owned by one block.

## Usage

``` r
block_diag_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`BlockDiagParam`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md).

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order, 1 to 4.

## Value

A named list of matrices keyed as `param_tuple_names(s, order)`.
