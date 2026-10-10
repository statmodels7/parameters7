# The Scale of a Fixed Matrix

Declares the multiplier \\h(\eta)\\, the one quantity that describes a
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
beyond its fixed matrix, with its interval on the log scale. The fixed
matrix is supplied by the caller and is not reported.

A parameter built with `link = NULL` has no free value and nothing to
declare, so this method returns `NULL` there, as the base method does.

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object, whose `param_params$link` is read.

- eta:

  A numeric vector of at most one free value; `numeric(0)` for a fixed
  parameter.

- ...:

  Ignored.

## Value

A list as described in
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md),
or `NULL` when the matrix is fixed and carries no free value.
