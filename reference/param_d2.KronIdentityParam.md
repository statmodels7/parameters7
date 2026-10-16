# Second Derivatives of a Block Replication

The second derivative has closed form, \\\partial\_{kl} M = I_m \otimes
\partial\_{kl} S\\. The blocks share their free values, so a cross-block
component does not arise, and the list consists exactly of the inner
parameter's components, lifted.

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

A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
`param_tuple_names(s)` and in that order, each block diagonal and
labeled `v1`, `v2`, ..., `v(md)` on both margins.

## See also

[`param_d1.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.KronIdentityParam.md)
and
[`param_d3.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.KronIdentityParam.md)
for the neighboring orders.
