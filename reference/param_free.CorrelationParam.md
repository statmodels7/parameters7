# Free Vector of a Correlation Parameter

The angles read off the Cholesky factor, exact: \\\theta\_{i1}\\ is the
arc cosine of \\L\_{i1}\\, and each subsequent angle divides out the
sines already recovered. A matrix that is not a correlation matrix, or
one whose factor reaches an angle of \\0\\ or \\\pi\\, is refused.

## Arguments

- s:

  A
  [`CorrelationParam`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- m:

  A correlation matrix.

- ...:

  Unused.

## Value

A named numeric vector of free values.
