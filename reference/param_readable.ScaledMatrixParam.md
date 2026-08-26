# The Scale of a Fixed Matrix

Declares the multiplier \\h(\eta)\\, the one quantity a
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
is about, with its interval on the log scale. The fixed matrix itself is
the caller's own and needs no reporting.

A parameter built with `link = NULL` has no free value and nothing to
declare, so this returns `NULL` there, which is the same answer the base
class gives and means the same thing.

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
