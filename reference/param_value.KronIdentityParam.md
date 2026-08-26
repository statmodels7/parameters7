# Value of a Block Replication

\#' Returns \\I_m \otimes S(\eta)\\, the inner parameter's matrix placed
into \\m\\ identical diagonal blocks. One
[`kronecker()`](https://rdrr.io/r/base/kronecker.html) call and no
arithmetic of its own: the value is the inner value lifted.

[`kronecker()`](https://rdrr.io/r/base/kronecker.html) drops the inner
parameter's labels, so the result is relabeled over the composite side,
`v1`, `v2`, ..., `v(md)`.

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

A `s@dimension` by `s@dimension` symmetric numeric matrix, block
diagonal with `m` identical blocks, and labeled `v1`, `v2`, ..., `vp` on
both margins. Positive definite exactly when the inner parameter is.

## See also

[`param_free.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_free.KronIdentityParam.md)
for the inverse and
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
for the construction.
