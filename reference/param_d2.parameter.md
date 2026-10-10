# Default Second Derivatives

The method that every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
of its own. It takes exactly one difference per component, of whichever
quantity the family already supplies: the analytic
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
where there is one, at the order-1 step, and
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
itself where there is not, through a three-point second difference on
the diagonal and a four-point mixed stencil off it, at the order-2 step.
The truncation error is of order \\h^2\\ on every route. The families in
this package do not reach it.

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

A list of `choose(s@n_free + 1, 2)` estimates keyed as
`param_tuple_names(s)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## See also

[`numerical_d2()`](https://statmodels7.github.io/parameters7/reference/numerical_d2.md),
which does the work and writes the three routes out, and
[`param_d1.parameter()`](https://statmodels7.github.io/parameters7/reference/param_d1.parameter.md)
for the order below.
