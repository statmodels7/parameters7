# Third Derivatives of a Block Replication

Closed form: the inner parameter's third derivatives, each lifted into
\\m\\ identical blocks.

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

A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
`param_tuple_names(s, 3)` and in that order, each block diagonal with no
dimnames.

## See also

[`param_d4.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.KronIdentityParam.md)
for the order above.
