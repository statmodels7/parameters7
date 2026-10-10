# The Scale and the Common Correlation of a Compound Symmetry

Declares two quantities: the marginal variance and the correlation that
every pair shares. The Jacobian is diagonal, through
[`readable_diagonal()`](https://statmodels7.github.io/parameters7/reference/readable_diagonal.md),
each quantity being one link of one free value.

The correlation's interval is built on the inverse hyperbolic tangent.
The family's own link is bounded below by \\-1/(p-1)\\, which is above
\\-1\\ for \\p \> 2\\, so an interval built on this scale can extend
below the bound that the parametrization enforces; it stays inside
\\(-1, 1)\\, the range of a correlation.

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object, whose two links are read.

- eta:

  A numeric vector of two free values.

- ...:

  Ignored.

## Value

A list as described in
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md).
