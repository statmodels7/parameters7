# No Declared Quantities

The method that every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it declares no interpretable quantities of its own. It
returns `NULL`. For
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
the five composition wrappers and the other families without a method of
their own, the matrix itself is the quantity that a reader reads, and a
consumer that receives `NULL` reports the matrix.

Because the result is `NULL` and not an error, a consumer can call the
generic on every parameter in a model without knowing which families
declare something.

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object. Not read.

- eta:

  A numeric vector of free values. Not read.

- ...:

  Ignored.

## Value

`NULL`.

## See also

[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
for the structure of the declaration and the list of families that
declare something.
