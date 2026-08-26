# The Scale and the Correlation of an AR(1)

Declares two quantities: the marginal variance and the correlation at
lag one, the two things an AR(1) covariance is about. Their intervals
are built on the log and the inverse hyperbolic tangent, so a variance
stays positive and a correlation stays inside \\(-1, 1)\\.

The Jacobian is diagonal, through
[`readable_diagonal()`](https://statmodels7.github.io/parameters7/reference/readable_diagonal.md):
each quantity is one link of one free value.

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object, whose two links are read.

- eta:

  A numeric vector of two free values.

- ...:

  Ignored.

## Value

A list as described in
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md).
