# Free Vector of an AR(1) Parameter

The scale is the common diagonal entry and the correlation is the first
off-diagonal one divided by it, both exact; the rest of the matrix is
then checked against the pattern those two imply, and a matrix that does
not match is rejected rather than fitted.

## Arguments

- s:

  An
  [`Ar1Param`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- m:

  An AR(1) matrix.

- ...:

  Unused.

## Value

A named numeric vector of two free values.
