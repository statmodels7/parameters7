# Solve of a Compound Symmetry Parameter

Exact, by Sherman-Morrison. The inverse of \\\sigma^2\\(1-\rho)I + \rho
J\\\\ is

\$\$\frac{1}{\sigma^2(1-\rho)} \left\[I - \frac{\rho}{1 +
(p-1)\rho}J\right\],\$\$

compound symmetric again, so the inverse is written down and no
factorization is performed: the cost is \\O(p^2)\\ per column of `b`,
against the \\O(p^3)\\ Cholesky factorization of the base method. The
factors \\1 + (p-1)\rho\\ and \\1 - \rho\\ are evaluated from the second
free value, as \\p\omega\\ and \\p(1-\omega)/(p-1)\\ with \\\omega\\ the
logistic function of that value, so the inverse stays accurate where the
correlation rounds to one of its bounds.

An exchangeable covariance has an exchangeable precision, so the family
is closed under the choice of side.

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
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

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for the generic, and
[`param_solve.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_solve.matrix_parameter.md)
for the factorization this avoids.
