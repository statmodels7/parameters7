# No Declared Quantities

This method is inherited by every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
that declares no interpretable quantities of its own. It returns `NULL`.
For
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
the composition wrappers and the other families without a method of
their own, the matrix itself is the quantity that a reader reads, and a
consumer that receives `NULL` reports the matrix.

Because the result is `NULL` and not an error, a consumer can call the
generic on every parameter in a model without checking which families
declare something.

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object. The method does not read it.

- eta:

  A numeric vector of free values. The method does not read it.

- ...:

  Ignored.

## Value

`NULL`.

## See also

[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
for the structure of the declaration and the list of families that
declare something.
