# Fourth Derivatives of a Block Replication

Closed form: the inner parameter's fourth derivatives, each lifted into
\\m\\ identical blocks. Each component is a dense matrix of side \\md\\,
so the memory grows with \\(md)^2\\ for every component.

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

A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
`param_tuple_names(s, 4)` and in that order, each block diagonal and
labeled `v1`, `v2`, ..., `v(md)` on both margins.

## See also

[`param_d3.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.KronIdentityParam.md)
for the order below.
