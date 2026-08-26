# First Derivatives of a Transition Matrix Parameter

Closed form. The free value `alr{i}.{j}` belongs to row \\i\\ alone, so
\\\partial P\\ is zero in every other row and, inside row \\i\\, is the
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
first derivative \\\partial_b \pi_a = \pi_a(\delta\_{ab} - \pi_b)\\: the
covariance structure of a categorical indicator.

Each component therefore sums to zero along its own row, the row sums of
\\P\\ being the constant 1.

## Arguments

- s:

  A
  [`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `s@n_free` \\K \times K\\ matrices named by `s@free_names`,
each supported on one row.

## See also

[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md),
which assembles it,
[`param_d1.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.SimplexParam.md)
for the per-row formula, and
[`param_d2.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.TransitionMatrixParam.md)
for the order above.
