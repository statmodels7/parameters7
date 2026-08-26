# Free Vector of a Block-Diagonal Parameter

Inverts each diagonal block through its own parameter and concatenates
the results, so the composite is exact wherever its blocks are: measured
on `block_diag(log_cholesky(2), ar1(3))`, the round trip closes to \\7
\times 10^{-17}\\.

## Arguments

- s:

  A
  [`BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md)
  object.

- m:

  A symmetric `s@dimension` by `s@dimension` matrix, block diagonal in
  this parameter's blocks, already checked for shape and symmetry by the
  generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, the blocks' free vectors
concatenated.

## Details

A matrix whose off-diagonal blocks are not zero is rejected with
`'m' is not block diagonal in the blocks of this parameter.`, that being
a matrix this family cannot represent; the tolerance is \\10^{-8}\\
relative to `max(1, max(abs(m)))`. The blocks are inverted first, so a
diagonal block that is outside its own family is reported by that
family's message, which is the more specific of the two.

## See also

[`param_value.BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/param_value.BlockDiagParam.md),
the map this inverts.
