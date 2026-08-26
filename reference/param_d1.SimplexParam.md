# First Derivatives of a Simplex Parameter

Closed form:

\$\$\partial_b \pi_a = \pi_a(\delta\_{ab} - \pi_b),\$\$

which is the covariance structure of a categorical indicator: the same
array a multinomial score already carries. Every component sums to zero
over the category index, \\\sum_a \pi_a\\ being the constant 1.

## Arguments

- s:

  A
  [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `s@n_free` numeric vectors of length \\K\\, named by
`s@free_names`, each summing to 0.

## See also

[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md),
which builds the array, and
[`param_d2.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.SimplexParam.md)
for the order above.
