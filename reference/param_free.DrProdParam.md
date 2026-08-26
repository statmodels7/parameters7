# Free Vector of a Scales-Times-Correlation Parameter

Reads the standard deviations off the diagonal as
\\\sqrt{\Sigma\_{jj}}\\, divides them out, and hands the resulting
correlation matrix to the correlation block. Exact wherever the block
is: measured at \\p = 3\\ with the default block, the round trip closes
to \\7 \times 10^{-16}\\.

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
followed by the correlation block's own free vector.

## Details

A matrix with a non-positive diagonal entry is rejected with
`'m' must have a positive diagonal.`; anything else outside the family
is reported by the correlation block's own message, which is the more
specific of the two.

This is where a `correlation` block that carries a scale of its own
shows: the division always leaves a unit diagonal, so such a block is
handed a matrix its own scale cannot be recovered from, and the round
trip fails without any error being signaled. See
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
on that requirement.

## See also

[`param_value.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_value.DrProdParam.md),
the map this inverts.
