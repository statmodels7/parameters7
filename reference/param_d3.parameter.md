# Default Third Derivatives

The method every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
of its own. It applies one product stencil per index tuple directly to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
never to a lower-order numerical derivative, with a one-dimensional
factor per distinct component of the order that component's multiplicity
asks for. The step is \\\varepsilon^{1/5}\max(1, \|\eta_k\|)\\, about
\\7.4 \times 10^{-4}\\ near the origin, and the truncation error is of
order \\h^2\\; measured against a closed form on a \\2 \times 2\\
covariance the gap is \\7 \times 10^{-6}\\ on entries of size 12. No
family in this package reaches it.

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 2, 3)` estimates keyed as
`param_tuple_names(s, 3)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## See also

[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md),
which does the work and writes the stencil out, and
[`param_d4.parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4.parameter.md)
for the order above.
