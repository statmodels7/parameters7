# Value of a Sum of Fixed Matrices

The weighted sum \\\sum_k c_k(\eta_k) P_k\\, accumulated in place over
the components. Positive semidefinite at every free vector, a
non-negative combination of positive semidefinite matrices being one;
positive **definite** only where the components' null spaces meet at the
origin, which is the condition `rank` records.

The value is labeled `v1`, `v2`, ..., `vp` on both margins, the
convention
[`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md)
states and every family in the package follows.

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

A symmetric positive semidefinite `s@dimension` by `s@dimension` numeric
matrix, labeled `v1`, `v2`, ..., `vp` on both margins.

## See also

[`param_free.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_free.SumStructParam.md)
for the inverse, and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for the parametrization.
