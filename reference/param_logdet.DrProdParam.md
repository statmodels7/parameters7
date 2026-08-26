# Log-Determinant of a Scales-Times-Correlation Parameter

Closed form, \\2\sum_j \log d_j + \log\lvert R \rvert\\, the scales
contributing twice because they multiply on both sides. The
correlation's term comes from that family by whatever route it has, so a
closed form there stays a closed form here and no determinant of the
assembled matrix is taken.

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
for the term the default block supplies.
