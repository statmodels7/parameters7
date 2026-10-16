# Log-Determinant of a Block Replication

The log-determinant has closed form, \\\log\|M\|\_{+} = m
\log\|S\|\_{+}\\, which is the inner parameter's
log-(pseudo-)determinant times the number of blocks. A block-diagonal
determinant is the product of the blocks' and the blocks are identical,
so the method does not decompose any matrix of side \\md\\.

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

A single number.

## See also

[`param_dlogdet.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.KronIdentityParam.md)
for its derivatives.
