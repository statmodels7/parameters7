# Fourth Log-Determinant Derivatives of a Block Replication

Closed form: \\m\\ times the inner parameter's fourth-order vector.
Exact where the inner parameter's is, so the composition never
introduces the accuracy loss a fallback at this order would.

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

A numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
`param_tuple_names(s, 4)`.

## See also

[`param_d3logdet.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.KronIdentityParam.md)
for the order below.
