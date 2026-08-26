# Second Derivatives of a Simplex Parameter

Closed form, one further application of the product rule to \\\partial_b
\pi_a = \pi_a(\delta\_{ab} - \pi_b)\\. The array is the third cumulant
of a categorical indicator, and every component sums to zero over the
category index.

Unlike a diagonal family, nothing here is separable: the entries of
\\\pi\\ are coupled through the shared normalizing constant, so a
component in two different free values is generally non-zero.

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

A list of `choose(s@n_free + 1, 2)` numeric vectors of length \\K\\,
keyed as `param_tuple_names(s)` and in that order, each summing to 0.

## See also

[`param_d1.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.SimplexParam.md)
and
[`param_d3.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.SimplexParam.md)
for the neighboring orders.
