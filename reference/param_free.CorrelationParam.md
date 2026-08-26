# Free Vector of a Correlation Parameter

Returns the angles behind a correlation matrix, read off its Cholesky
factor and carried onto the free scale by the link. Exact:
\\\theta\_{i1} = \arccos L\_{i1}\\, and each subsequent angle divides
out the sines already recovered before taking an arc cosine. A true
inverse of
[`param_value.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_value.CorrelationParam.md),
the round trip closing to \\2 \times 10^{-16}\\.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- m:

  A correlation matrix of side `s@dimension`: symmetric with a unit
  diagonal and positive definite, already checked for shape and symmetry
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

Three rejections. `m` must have a unit diagonal, or it is not a
correlation matrix. It must be positive definite, tested spectrally
through
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md).
And its factor must not reach an angle of exactly 0 or \\\pi\\, where
the link has no finite value: that happens at a correlation of \\\pm
1\\, which is on the boundary of the set without being in it.

The last case is worth knowing before it is met. The parametrization
reaches a correlation of \\-1\\ to every printed digit at a free value
near 10, and the matrix is then singular in double precision, so a
matrix produced by a fit that ran to its boundary cannot be inverted
back.

## See also

[`param_value.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_value.CorrelationParam.md),
the map this inverts, and
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
for the definiteness test.
