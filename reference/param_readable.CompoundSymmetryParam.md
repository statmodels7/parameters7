# The Scale and the Common Correlation of a Compound Symmetry

Declares two quantities: the marginal variance and the correlation every
pair shares. The Jacobian is diagonal, through
[`readable_diagonal()`](https://statmodels7.github.io/parameters7/reference/readable_diagonal.md),
each quantity being one link of one free value.

The correlation's interval is built on the inverse hyperbolic tangent,
which is worth knowing here: the family's own link is bounded at
\\-1/(p-1)\\ and never at \\-1\\, so an interval built this way can
reach below the bound the parametrization enforces. It respects \\(-1,
1)\\, which is the interval a reader of a correlation expects.

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
