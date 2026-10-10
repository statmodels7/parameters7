# Solve and Factor of a Sum of Fixed Matrices

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for a
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
parameter. Both come from the assembled matrix, by
[`solve()`](https://rdrr.io/r/base/solve.html) and by a transposed
[`chol()`](https://rdrr.io/r/base/chol.html): a sum of fixed matrices
has no structure that a solve could exploit, where
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
has a tridiagonal inverse and
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
assembles its factor directly.

## Arguments

- s:

  A
  [`SumStructParam()`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- b:

  A numeric matrix with `s@dimension` rows, defaulted to the identity by
  the generic.
  [`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
  takes no `b`.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
returns a numeric matrix with `s@dimension` rows and as many columns as
`b`;
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
a lower triangular `s@dimension` by `s@dimension` matrix. Both carry no
dimnames, while the value and the derivative arrays carry `v1`, `v2`,
...; the names of the value are removed because this family computes
both results from its own value.

## Details

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
is the call [`solve()`](https://rdrr.io/r/base/solve.html) on the
assembled matrix; the factor is `t(chol(M))`, lower triangular with \\M
= L L^\top\\ as
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
requires.

Both are rejected by the generic where the family is rank deficient, a
singular matrix having neither an inverse nor a Cholesky factor.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the two generics.
