# The Inner Parameter, the Block Count, and the Lift

The three one-line accessors every
[`KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md)
method uses. `.kron_inner()` returns the per-block parameter,
`.kron_m()` the number of blocks, and `.kron_lift()` places a `d` by `d`
matrix into `m` identical diagonal blocks as `kronecker(diag(m), x)`.

## Arguments

- s:

  A
  [`KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md)
  object.

- m:

  For `.kron_lift()`, the number of blocks.

- x:

  A numeric matrix to replicate.

## Value

`.kron_inner()` a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
`.kron_m()` a single integer, `.kron_lift()` an `m * nrow(x)` by
`m * ncol(x)` numeric matrix with `x` on the diagonal blocks and zeros
elsewhere.

## Details

They exist so that the twelve methods below read as one line each and
the storage of `param_params` appears in one place. `.kron_lift()` is
the whole arithmetic of the composition: the value, every derivative
component and the factor are the inner quantity passed through it.

Note that [`kronecker()`](https://rdrr.io/r/base/kronecker.html) drops
dimnames, so a lifted matrix carries none where the inner one carried
`v1`, `v2`, ...

## See also

[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
which stores what these read.
