# Log-Determinant of a Compound Symmetry Parameter

Closed form, from the two distinct eigenvalues
\\\sigma^2\\1+(p-1)\rho\\\\ and \\\sigma^2(1-\rho)\\:

\$\$\log\|M\| = p\log\sigma^2 + \log\\1+(p-1)\rho\\ +
(p-1)\log(1-\rho).\$\$

Three logarithms and no factorization, whatever \\p\\ is, against the
base class's \\O(p^3)\\ eigendecomposition. Measured at \\p = 4\\
against the eigenvalues of the assembled matrix, the two agree to the
printed digit.

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number, finite at every free vector: both eigenvalues are
strictly positive inside the correlation's link bounds.

## See also

[`param_dlogdet.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CompoundSymmetryParam.md)
for its four derivative orders, and
[`cs_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_terms.md)
for the two terms.
