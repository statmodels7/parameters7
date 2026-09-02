# Solve and Factor of an Inverse Parameter

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
is the inner family's **value**, no inversion being needed: the inverse
of \\S^{-1}\\ is \\S\\.
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
is the lower triangular \\L\\ with \\LL^\top\\ the value, which is the
package's convention and the transpose of what
[`base::chol()`](https://rdrr.io/r/base/chol.html) returns.

## Arguments

- s:

  An
  [`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- b:

  A numeric matrix or vector with `s@dimension` rows, or `NULL` for the
  whole inverse.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
the inner value, or its product with `b`;
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
a lower triangular matrix `L` with `tcrossprod(L)` the value.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for the family.
