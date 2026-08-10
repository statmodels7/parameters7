# Assemble a Block Diagonal's Log-Determinant Derivatives

The log-determinant is the sum of the blocks', hence separable across
blocks, so a component is the owning block's own and zero for any tuple
spanning two blocks.

## Usage

``` r
block_diag_logdet_derivs(s, eta, order)
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

A named numeric vector keyed as `param_tuple_names(s, order)`.
