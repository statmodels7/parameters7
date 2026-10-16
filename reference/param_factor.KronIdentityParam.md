# Factor of a Block Replication

The factor has closed form, \\I_m \otimes L\\, which lifts the inner
parameter's Cholesky factor. A block-diagonal matrix of identical blocks
factors blockwise, so one inner factorization serves every block and the
method does not decompose any matrix of side \\md\\.

## Arguments

- s:

  A
  [`KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A `s@dimension` by `s@dimension` lower triangular numeric matrix, block
diagonal with `m` identical blocks, satisfying
`L %*% t(L) == param_value(s, eta)`, with no dimnames.

## See also

[`param_solve.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_solve.KronIdentityParam.md)
for the blockwise solve.
