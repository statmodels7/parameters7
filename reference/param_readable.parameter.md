# No Declared Quantities

The method every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it declares no interpretable quantities of its own. It
returns `NULL`, which is a positive statement and no kind of gap: for
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
and the four compositions the matrix itself is what a reader reads, and
a consumer that gets `NULL` reports the matrix as it already would.

Returning `NULL` instead of throwing is what a consumer needs in order
to loop over every parameter in a model and ask each one, without
knowing which kinds declare something.

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
for the contract and the list of families that do declare something.
