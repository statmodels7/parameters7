# Log-Determinant of a Scales-Times-Correlation Parameter

The log-determinant is in closed form, \\2\sum_j \log d_j + \log\lvert R
\rvert\\, and the scales contribute twice because they multiply on both
sides. The correlation's term comes from that family by whatever route
that family has, so a closed form there stays a closed form here and the
determinant of the assembled matrix is not taken.

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number.

## See also

[`param_dlogdet.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.DrProdParam.md)
for its derivatives, and
[`param_logdet.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.CorrelationParam.md)
for the term that the default block supplies.
