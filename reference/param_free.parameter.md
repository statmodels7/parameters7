# Refusal to Invert Without a Closed Form

The base class refuses rather than inverting the map numerically: an
optimisation-based inverse would return a plausible \\\eta\\ for a
matrix outside the set the family parametrises.

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- m:

  A symmetric numeric matrix.

- ...:

  Unused.

## Value

Never returns; raises an error.
