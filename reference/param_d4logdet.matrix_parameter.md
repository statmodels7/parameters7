# Default Fourth Log-Determinant Derivatives

This is the method that every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
of its own. Each component is one tensor-product central difference,
through
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md)
at accuracy 4, of the highest analytic log-determinant derivative below
the fourth: one first difference of
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
where the family registers it, and otherwise, for the tuple \\(k, l, m,
n)\\, one second-order stencil of the \\(k, l)\\ component of
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
in the remaining two components, the second difference where \\m = n\\
and the product of two first differences where they differ.

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

A numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
`param_tuple_names(s, 4)` and in that order.

## Accuracy

With analytic
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
or with a family's own
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
the method applies the difference in a single layer to an exact
quantity, and the result is accurate to about \\10^{-9}\\.

When only
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
is supplied, the base-class
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
would read numerical arrays and is not used: the stencil is a fourth
difference of
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
which is exact, accurate to about \\10^{-6}\\. A family that needs this
order to more digits supplies
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
in closed form, or its own
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md).
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports which components are numerical.

## See also

[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for the generic,
[`param_d3logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
for the order below, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
to find out which components are numerical.
