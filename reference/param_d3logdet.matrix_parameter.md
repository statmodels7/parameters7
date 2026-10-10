# Default Third Log-Determinant Derivatives

The method that every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
of its own. For the tuple \\(k, l, m)\\ it takes the \\(k, l)\\
component of
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
and applies one three-point central difference in the remaining
component,

\$\$\partial\_{klm} \log\|M\| \approx
\frac{\partial\_{kl}\log\|M\|\\(\eta + h e_m) -
\partial\_{kl}\log\|M\|\\(\eta - h e_m)}{2h},\$\$

at the order-1 step \\\varepsilon^{1/3}\max(1, \|\eta_m\|)\\. Since the
tuple is sorted, the differenced component is the largest index, and the
pair is looked up by matching sorted indices against
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)'s
order-2 list. It costs two evaluations of the whole second-order block
per tuple.

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
is an exact identity **given** the matrix derivative arrays, so
differencing it is a single numerical layer whenever
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
are analytic.

A family's own
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
is exact, and differencing it is a single layer too. Where the
base-class
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
reads numerical arrays the layers would compound, and the method signals
an error instead of returning a value;
[`check_analytic_arrays()`](https://statmodels7.github.io/parameters7/reference/check_analytic_arrays.md)
is the guard. A family that needs this order writes
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
out, or its own
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md).

## See also

[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
for the generic,
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
for the quantity differenced, and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the order above.
