# Log-Determinant of a Compound Symmetry Parameter

Closed form, from the two distinct eigenvalues
\\\sigma^2\\1+(p-1)\rho\\\\ and \\\sigma^2(1-\rho)\\:

\$\$\log\|M\| = p\log\sigma^2 + \log\\1+(p-1)\rho\\ +
(p-1)\log(1-\rho).\$\$

Three logarithms and no factorization, whatever \\p\\ is, where the base
method takes an \\O(p^3)\\ eigendecomposition.

The two correlation logarithms are evaluated from the second free value,
as \\\log p + \log\omega\\ and \\\log(p/(p-1)) + \log(1-\omega)\\ with
\\\omega\\ the logistic function of that value, through
`plogis(., log.p = TRUE)`, and not from the rounded correlation. The
value is then accurate at every free value, including those where the
correlation rounds to one of its bounds.

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
[`cs_logdet_chain()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_chain.md)
for the derivatives of the correlation's term.
