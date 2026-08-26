# Solve and Factor of a Scales-Times-Correlation Parameter

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for a
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
parameter. Both come from the correlation block with a scaling on either
side: \\\Sigma^{-1} = D^{-1} R^{-1} D^{-1}\\, and \\\Sigma =
(DL)(DL)^\top\\ for \\L\\ the correlation's own factor. Neither forms
\\\Sigma\\ and neither factorizes it.

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- b:

  A numeric matrix with `s@dimension` rows, defaulted to the identity by
  the generic.
  [`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
  takes no `b`.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
returns a numeric matrix with `s@dimension` rows and as many columns as
`b`;
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
a lower triangular `s@dimension` by `s@dimension` matrix.

## Details

The scaling is done by dividing and multiplying rows, which for a
diagonal matrix is what the products amount to. Measured at \\p = 3\\
with the default block, the solve agrees with
[`solve()`](https://rdrr.io/r/base/solve.html) on the assembled matrix
to \\7 \times 10^{-15}\\ and `tcrossprod(param_factor(s, eta))` with the
matrix to \\4 \times 10^{-16}\\. The factor is lower triangular, as
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
requires, because scaling the rows of a lower triangular matrix leaves
it lower triangular.

Neither is ever refused for rank: this family admits no deficient
correlation block, so it is always of full rank.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the two contracts.
