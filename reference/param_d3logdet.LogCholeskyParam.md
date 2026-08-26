# Third Log-Determinant Derivatives of a Log-Cholesky Parameter

Closed form, and identically zero. \\\log\|M\| = 2\sum_i \eta_i\\ is
linear in the free vector, so every derivative above the first vanishes
at every \\\eta\\. The zeros are exact, so a consumer can drop the term
instead of carrying small numbers through a contraction.

It also means this family cannot exercise the order: a check of a
log-determinant's third derivative against a numerical reference passes
here whatever is missing, both sides being zero.
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
are where this order has content.

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic. Its values do not enter the result.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of `choose(s@n_free + 2, 3)` zeros, keyed as
`param_tuple_names(s, 3)`.

## See also

[`param_d2logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.LogCholeskyParam.md)
for the order below, and
[`param_d4logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.LogCholeskyParam.md)
for the one above.
