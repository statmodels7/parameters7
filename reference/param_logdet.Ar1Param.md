# Log-Determinant of an AR(1) Parameter

Closed form. The determinant of the correlation pattern is
\\(1-\rho^2)^{p-1}\\, so

\$\$\log\|M\| = p\log\sigma^2 + (p-1)\log(1-\rho^2).\$\$

Two logarithms and no factorization, whatever \\p\\ is, where the base
method takes an \\O(p^3)\\ eigendecomposition. With \\\rho = \tanh z\\,
the second term is evaluated as \\(p-1)\\2\log 2 - 2\|z\| - 2\log(1 +
e^{-2\|z\|})\\\\, which does not lose accuracy where \\\rho\\ rounds to
\\-1\\ or 1.

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

A single number, finite at every free vector.

## See also

[`param_dlogdet.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.Ar1Param.md)
for its four derivative orders, and
[`ar1_logdet_chain()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_chain.md)
for the derivatives of the correlation's term.
