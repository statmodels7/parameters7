# Factor of a Correlation Parameter

Returns \\L\\ by assembling it from the angles, so no factorization is
taken: the factor is what the parametrization holds. Its rows are unit
vectors, which is where \\R = LL^\top\\ gets its unit diagonal.

It is the natural way to simulate a correlated standard normal vector:
\\Lz\\ with \\z\\ standard normal has correlation matrix \\R\\ and unit
variances.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A `s@dimension` by `s@dimension` lower triangular numeric matrix with a
positive diagonal and unit row norms, satisfying
`L %*% t(L) == param_value(s, eta)`, and carrying no dimnames.

## See also

[`param_value.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_value.CorrelationParam.md)
for the matrix, and
[`param_factor.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_factor.matrix_parameter.md)
for what a family without a closed form pays.
