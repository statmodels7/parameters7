# Log-Determinant Derivatives of a Block-Diagonal Parameter

[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for a
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
parameter: the blocks' own, placed by owner, with every cross-block
component exactly zero because the log-determinant is a sum over the
blocks.

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

A numeric vector: at order 1, `s@n_free` values named by `s@free_names`;
above it, `choose(s@n_free + k - 1, k)` values keyed as
`param_tuple_names(s, k)` and in that order.

## Details

The four share
[`block_diag_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/block_diag_logdet_derivs.md)
and differ only in the order they pass. The first order names its result
by `s@free_names`, one value per free value; the orders above it are
keyed by tuple, as the contract requires.

A block's **own** mixed components are not zero, this separability being
over blocks, never over free values. Compare
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
where it is over free values and every mixed component vanishes.

## See also

[`block_diag_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/block_diag_logdet_derivs.md),
which assembles them, and
[`param_logdet.BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.BlockDiagParam.md)
for the quantity differentiated.
