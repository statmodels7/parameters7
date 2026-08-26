# Second Derivatives of a Block Replication

Closed form: \\\partial\_{kl} M = I_m \otimes \partial\_{kl} S\\. The
blocks share their free values, so there is no cross-block component to
compute and the list is exactly the inner parameter's, lifted.

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
`param_tuple_names(s)` and in that order, each block diagonal with no
dimnames.

## See also

[`param_d1.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.KronIdentityParam.md)
and
[`param_d3.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.KronIdentityParam.md)
for the neighboring orders.
