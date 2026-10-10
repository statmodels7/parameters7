# Default Fourth Log-Determinant Derivatives

The method that every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
of its own. For the tuple \\(k, l, m, n)\\ it takes the \\(k, l)\\
component of
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
and applies one second-order stencil in the remaining two components:
the three-point second difference where \\m = n\\, and the four-point
mixed stencil

\$\$\frac{g(\eta + h_m e_m + h_n e_n) - g(\eta + h_m e_m - h_n e_n) -
g(\eta - h_m e_m + h_n e_n) + g(\eta - h_m e_m - h_n e_n)} {4 h_m
h_n}\$\$

where they differ, both at the order-2 step \\\varepsilon^{1/4}\max(1,
\|\eta_k\|)\\. It costs three or four evaluations of the whole
second-order block per tuple.

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
the differencing is a single layer on an exact quantity, at the order-2
step.

With only
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
supplied, the numerical arrays would feed a numerical second-order block
that is then differenced twice more, and the error can exceed the
quantity itself. The method signals an error instead, through
[`check_analytic_arrays()`](https://statmodels7.github.io/parameters7/reference/check_analytic_arrays.md).
A family that needs a fourth derivative of its log-determinant supplies
at least
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
