# Value of a Sum of Fixed Matrices

Returns the weighted sum \\\sum_k c_k(\eta_k) P_k\\, accumulated in
place over the components. The matrix is positive semidefinite at every
free vector, because a non-negative combination of positive semidefinite
matrices is positive semidefinite, and it is positive definite only
where the components' null spaces meet at the origin, which is the
condition that `rank` records.

The value is labeled `v1`, `v2`, ..., `vp` on both margins, the
convention that
[`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md)
states.

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
