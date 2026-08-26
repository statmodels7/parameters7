# Solve of an AR(1) Parameter

Exact, and tridiagonal. The precision of an AR(1) covariance is
\\\\\sigma^2(1-\rho^2)\\^{-1}\\ times the matrix with 1 at the two
corners of the diagonal, \\1+\rho^2\\ elsewhere on it, and \\-\rho\\ on
the two first off-diagonals. Every other entry is exactly zero, an AR(1)
process being Markov, so the inverse is written down and no
factorization is performed: the cost is \\O(p)\\ entries against the
base class's \\O(p^3)\\ Cholesky.

Measured at \\p = 4\\ and \\\rho = 0.6\\: the entries beyond the first
off-diagonal are 0 exactly, and the whole matrix agrees with
[`base::solve()`](https://rdrr.io/r/base/solve.html) on the assembled
covariance to \\1 \times 10^{-16}\\.

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- b:

  A numeric matrix with `s@dimension` rows, already coerced from a
  vector and defaulted to the identity by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric matrix with `s@dimension` rows and as many columns as `b`.

## See also

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
which explains why the tridiagonal precision is not itself AR(1), and
[`param_solve.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_solve.matrix_parameter.md)
for the factorization this avoids.
