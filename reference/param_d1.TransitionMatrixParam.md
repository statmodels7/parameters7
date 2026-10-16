# First Derivatives of a Transition Matrix Parameter

Computes the first derivatives in closed form. The free value
`alr{i}.{j}` belongs only to row \\i\\, so \\\partial P\\ is zero in
every other row and, inside row \\i\\, equals the
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
first derivative \\\partial_b \pi_a = \pi_a(\delta\_{ab} - \pi_b)\\,
which is the covariance structure of a categorical indicator.

Each component therefore sums to zero along its own row, because the row
sums of \\P\\ are the constant 1.

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

A list of `s@n_free` \\K \times K\\ matrices named by `s@free_names`.
Each matrix is supported on one row.

## See also

[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md),
which assembles it,
[`param_d1.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.SimplexParam.md)
for the per-row formula, and
[`param_d2.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.TransitionMatrixParam.md)
for the order above.
