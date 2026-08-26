# Derivatives of a Sum of Fixed Matrices

[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for a
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
parameter. The value being linear in the weights, a component is
**exactly zero** unless every index names the same free value, and is
then that weight's \\m\\-th derivative times its own fixed matrix.
Nothing is differenced and no matrix arithmetic is done.

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

At order 1, a list of `s@n_free` symmetric matrices named by
`s@free_names`; above it, `choose(s@n_free + k - 1, k)` of them keyed as
`param_tuple_names(s, k)` and in that order. Each is `s@dimension` by
`s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.

## Details

The four share
[`sum_struct_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_derivs.md)
and differ only in the order they pass. At \\K = 2\\ the zeros are 1 of
the 3 second-order components, 2 of 4 at third order and 3 of 5 at
fourth. Compare
[`param_dlogdet.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.SumStructParam.md),
where nothing vanishes: the value is linear in the weights and its
log-determinant is not.

## See also

[`sum_struct_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_derivs.md),
which assembles them.
