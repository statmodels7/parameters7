# Derivatives of a Block-Diagonal Parameter

[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for a
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
parameter: each block's own derivatives, placed in the rows and columns
that block occupies, and **exactly zero** for a tuple spanning two
blocks. Nothing is rederived and nothing is differenced, so the
composite is as exact as its blocks are.

## Arguments

- s:

  A
  [`BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

At order 1, a list of `s@n_free` symmetric matrices named by
`s@free_names`; above it, `choose(s@n_free + k - 1, k)` of them keyed as
`param_tuple_names(s, k)` and in that order. Each is `s@dimension` by
`s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.

## Details

The four share
[`block_diag_derivs()`](https://statmodels7.github.io/parameters7/reference/block_diag_derivs.md)
and differ only in the order they pass. The zeros are structural: they
follow from the free values of one block not entering another, and at
order 4 with blocks of 3 and 2 free values they are 50 of the 70
components. See
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
for the argument.

## See also

[`block_diag_derivs()`](https://statmodels7.github.io/parameters7/reference/block_diag_derivs.md),
which assembles them, and
[`param_dlogdet.BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.BlockDiagParam.md)
for the log-determinant's own.
