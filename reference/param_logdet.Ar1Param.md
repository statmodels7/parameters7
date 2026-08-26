# Log-Determinant of an AR(1) Parameter

Closed form. The determinant of the correlation pattern is
\\(1-\rho^2)^{p-1}\\, so

\$\$\log\|M\| = p\log\sigma^2 + (p-1)\log(1-\rho^2).\$\$

Two logarithms and no factorization, whatever \\p\\ is, against the base
class's \\O(p^3)\\ eigendecomposition. Measured at \\p = 4\\ against the
eigenvalues of the assembled matrix, the two agree to the printed digit.

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number, finite at every free vector, the rhobit link keeping
\\1-\rho^2\\ strictly positive.

## See also

[`param_dlogdet.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.Ar1Param.md)
for its four derivative orders, and
[`ar1_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_terms.md)
for the two terms.
