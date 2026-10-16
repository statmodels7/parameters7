# Default Third Derivatives

The method that every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
of its own. It applies one tensor-product stencil per index tuple,
through
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
at its default accuracy 4, to the highest analytic derivative below the
third and to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
where there is none, never to a lower-order numerical derivative. The
families in this package do not reach it.

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
