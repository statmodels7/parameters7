# Log-Determinant Hessian of a Block Replication

Closed form: \\m\\ times the inner parameter's Hessian, with the same
keying. Zero wherever the inner parameter's is, so replicating a
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
gives exact zeros here and replicating an
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
does not.

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

A numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
`param_tuple_names(s)`.

## See also

[`param_dlogdet.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.KronIdentityParam.md)
for the order below.
