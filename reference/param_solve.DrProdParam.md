# Solve and Factor of a Scales-Times-Correlation Parameter

Defines
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for a
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
parameter. Both come from the correlation block with a scaling on either
side: \\\Sigma^{-1} = D^{-1} R^{-1} D^{-1}\\, and \\\Sigma =
(DL)(DL)^\top\\ for \\L\\ the correlation's own factor. The methods do
not form \\\Sigma\\ and do not factorize it.

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

The methods scale the rows by division and multiplication, which for a
diagonal matrix is equivalent to forming the products. The factor is
lower triangular, as
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
requires, because scaling the rows of a lower triangular matrix leaves
it lower triangular.

The methods do not reject a matrix because of its rank: this family
admits only correlation blocks of full rank, so it always has full rank.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the two generics.
