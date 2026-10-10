# Default Log-Determinant Gradient

The method that every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
of its own. It evaluates the trace identity

\$\$\partial_k \log\|M\| = \mathrm{tr}\\\left(M^{+} \partial_k
M\right),\$\$

with \\M^{+}\\ the Moore-Penrose inverse formed from the directions that
the declared rank keeps, which is the ordinary inverse for a full-rank
family. The trace is computed as `sum(mi * dk)`, the elementwise product
summed, both matrices being symmetric, so no matrix product is formed.

The identity is exact, so the accuracy is entirely the accuracy of
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md):
rounding with an analytic first derivative, and that of a central
difference with a numerical one.

## Arguments

- s:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## See also

[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
for the generic,
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
for the derivative arrays this reads, and
[`spectrum_pinv()`](https://statmodels7.github.io/parameters7/reference/spectrum_pinv.md)
for \\M^{+}\\.
