# Refusal to Invert Without a Closed Form

The base class refuses rather than inverting the map numerically: an
optimisation-based inverse would return a plausible \\\eta\\ for a
matrix outside the set the family parametrises.

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- m:

  A symmetric numeric matrix.

- ...:

  Unused.

## Value

Never returns; raises an error.
