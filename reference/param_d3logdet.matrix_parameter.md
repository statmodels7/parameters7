# Default Third Log-Determinant Derivatives

The method every
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

## How accurate it is depends on the family, not on this method

[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
is an exact identity **given** the matrix derivative arrays, so
differencing it is a single numerical layer whenever
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
are analytic. Measured on a \\4 \times 4\\ AR(1) covariance with
analytic arrays, the answer agrees with the closed form to \\3 \times
10^{-10}\\ on entries of size 4.5.

Where the arrays are themselves numerical the layers would compound, and
the method refuses instead of answering: the same measurement gives \\8
\times 10^{-3}\\, a relative error near two parts in a thousand, and
nothing downstream could tell that number from the accurate one.
[`check_analytic_arrays()`](https://statmodels7.github.io/parameters7/reference/check_analytic_arrays.md)
is the guard, and a family that needs this order should write
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
out, which recovers the \\10^{-10}\\.

## See also

[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
for the generic,
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
for the quantity differenced, and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the order above.
