# Third Log-Determinant Derivatives of a Block Replication

Closed form: \\m\\ times the inner parameter's third-order vector.

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

A numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
`param_tuple_names(s, 3)`.

## See also

[`param_d4logdet.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.KronIdentityParam.md)
for the order above.
