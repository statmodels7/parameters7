# Derivatives of a Parameter by One Tensor Stencil per Component

Computes the derivatives of order `order` of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
or, with `what = "logdet"`, of
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
through
[`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html).
The analytic derivatives of lower order, as
[`param_analytic()`](https://statmodels7.github.io/parameters7/reference/param_analytic.md)
reads them, are declared to it, so each component is one tensor-product
central difference of the highest analytic quantity below it, and of the
value itself where there is none.

## Usage

``` r
param_tensor(s, eta, order, what = c("value", "logdet"), accuracy = 4L)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order, 1 to 4.

- what:

  `"value"` or `"logdet"`.

- accuracy:

  The order of accuracy of the one-dimensional stencils.

## Value

For `"value"`, a list of estimates keyed as
`param_tuple_names(s, order)` and in that order, each of which has the
shape of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result. For `"logdet"`, a named numeric vector in the same order.

## Details

The coordinates are the free values. The value is flattened to a vector,
and each component is returned with the attributes of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and, for a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
symmetrized as `(m + t(m)) / 2`. The step of a component whose
difference is of order \\r\\ is \\\varepsilon^{1/(r + a)}\max(1,
\|\eta_k\|)\\, with \\a\\ the accuracy; the unconstrained scale has no
bounds.

## See also

[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
to
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md),
[`param_d3logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md),
which call this.
