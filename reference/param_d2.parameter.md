# Default Second Derivatives

This is the method that every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
of its own. It takes exactly one difference per component, through
[`numerical_d2()`](https://statmodels7.github.io/parameters7/reference/numerical_d2.md)
at its default accuracy 4, of the quantity that the family already
supplies: the analytic
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
where there is one, and
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
itself where there is not. The truncation error is of order \\h^4\\ on
every route. The families in this package do not reach it.

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
`param_tuple_names(s)` and in that order, each of which has the shape of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and is symmetrized for a matrix family.

## See also

[`numerical_d2()`](https://statmodels7.github.io/parameters7/reference/numerical_d2.md),
which does the work and describes the two routes, and
[`param_d1.parameter()`](https://statmodels7.github.io/parameters7/reference/param_d1.parameter.md)
for the order below.
