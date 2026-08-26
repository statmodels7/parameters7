# Assemble a Block Diagonal's Derivatives of a Given Order

Places each block's own component in the rows and columns that block
occupies, and returns a zero matrix for a tuple whose indices are not
all owned by one block.

## Usage

``` r
block_diag_derivs(s, eta, order)
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

A list of `choose(s@n_free + order - 1, order)` symmetric matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.

## Details

Each block's derivatives of the order are fetched at most once and held
in a list for the rest of the call, since the composite's enumeration
visits a block's tuples several times and a block's fourth-order array
is the expensive thing here. The zero matrix is one object shared by
every cross-block component: at order 4 with two blocks of 3 and 2 free
values that is 50 of the 70 components.

## See also

[`block_derivs_by_tuple()`](https://statmodels7.github.io/parameters7/reference/block_derivs_by_tuple.md)
for one block's components, and
[`param_d1.BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.BlockDiagParam.md),
which calls this.
