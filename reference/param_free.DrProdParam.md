# Free Vector of a Scales-Times-Correlation Parameter

Reads the standard deviations off the diagonal as
\\\sqrt{\Sigma\_{jj}}\\, divides them out, and hands the resulting
correlation matrix to the correlation block, so the round trip is as
exact as the block's.

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- m:

  A symmetric positive definite `s@dimension` by `s@dimension` matrix,
  already checked for shape and symmetry by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`: the linked standard deviations
followed by the correlation block's own free vector, named by
`s@free_names`.

## Details

A matrix with a non-positive diagonal entry is rejected with
`'m' must have a positive diagonal.`; anything else outside the family
is reported by the correlation block's own message, which is the more
specific of the two.

The division always leaves a unit diagonal, so a block that carried a
scale of its own would receive a matrix from which its scale cannot be
recovered, and the round trip would fail without an error.
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
rejects such a block at construction.

## See also

[`param_value.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_value.DrProdParam.md),
the map this inverts.
