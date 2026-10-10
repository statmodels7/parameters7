# Log-Determinant Derivatives of a Sum of Fixed Matrices

[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for a
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
parameter: the cyclic trace expansion in the weights, carried onto the
free scale by a chain rule with a diagonal Jacobian.

## Arguments

- s:

  A
  [`SumStructParam()`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector: at order 1, `s@n_free` values named by `s@free_names`;
above it, `choose(s@n_free + k - 1, k)` values keyed as
`param_tuple_names(s, k)` and in that order.

## Details

The four share
[`sum_struct_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_logdet_derivs.md)
and differ only in the order they pass. Where the family is rank
deficient they differentiate the log pseudo-determinant that
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
returns. Unlike the value's own derivatives, the mixed components are in
general not zero: for the two variance components of the example in
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
the mixed second derivative is \\-0.245\\ against \\0.245\\ for each
pure one. The log-determinant of a sum is not a sum, which is why this
family needs an expansion where
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
and
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
need none.

It is the most expensive quantity of the family.

## See also

[`sum_struct_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_logdet_derivs.md),
which assembles them, and
[`param_logdet.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.SumStructParam.md)
for the quantity differentiated.
