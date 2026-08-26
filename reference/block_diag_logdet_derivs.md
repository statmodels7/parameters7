# Assemble a Block Diagonal's Log-Determinant Derivatives

The log-determinant is the sum of the blocks', so it is separable across
blocks: a component is the owning block's own, and 0 for a tuple
spanning two blocks.

## Usage

``` r
block_diag_logdet_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md).

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
as `param_tuple_names(s, order)` and in that order.

## Details

The structure is
[`block_diag_derivs()`](https://statmodels7.github.io/parameters7/reference/block_diag_derivs.md)'s
with numbers in place of matrices, including the per-block cache and the
sorted-tuple key. It is a different separability from the one
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md) or
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
have: there the log-determinant separates over the free values
themselves, here it separates over the blocks and a block's own mixed
components survive.

## See also

[`param_dlogdet.BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.BlockDiagParam.md),
which calls this, and
[`block_diag_derivs()`](https://statmodels7.github.io/parameters7/reference/block_diag_derivs.md),
its twin for the value.
