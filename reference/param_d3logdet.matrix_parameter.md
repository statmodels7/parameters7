# Default Third Log-Determinant Derivatives

This is the method that every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
of its own. Each component is one tensor-product central difference,
through
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md)
at accuracy 4, of the highest analytic log-determinant derivative below
the third, as
[`param_analytic()`](https://statmodels7.github.io/parameters7/reference/param_analytic.md)
reads it: for the tuple \\(k, l, m)\\, one first difference of the \\(k,
l)\\ component of
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
along \\\eta_m\\, at the step \\\varepsilon^{1/5}\max(1, \|\eta_m\|)\\.

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

A numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
`param_tuple_names(s, 3)` and in that order.

## Accuracy

[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
is an exact identity given the matrix derivative arrays, so it is
differenced where
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
are analytic or where the family writes it out, and the error is of
order \\10^{-11}\\. Where the base-class
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
would read numerical arrays it is not used: the stencil is applied to
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
where
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
is analytic, and otherwise to
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
which is exact, as a third difference. A numerical quantity is therefore
never differenced again, and the third difference of the log-determinant
is accurate to about \\10^{-7}\\.

## See also

[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
for the generic,
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
for the quantity differenced, and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the order above.
