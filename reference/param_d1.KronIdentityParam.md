# First Derivatives of a Block Replication

Closed form: \\\partial_k M = I_m \otimes \partial_k S\\, the inner
parameter's first derivatives lifted one at a time. Exact whenever the
inner parameter's are, and numerical whenever they are not, so
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
on the composite reports the inner parameter's own answer.

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

A list of `s@n_free` symmetric matrices named by `s@free_names`, each
`s@dimension` by `s@dimension`, block diagonal and labeled `v1`, `v2`,
..., `vp` on both margins.

## See also

[`param_d2.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.KronIdentityParam.md)
for the order above.
