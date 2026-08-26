# Solve and Factor of a Sum of Fixed Matrices

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for a
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
parameter. Both come from the assembled matrix, by
[`solve()`](https://rdrr.io/r/base/solve.html) and by a transposed
[`chol()`](https://rdrr.io/r/base/chol.html): a sum of fixed matrices
has no structure a solve could exploit, where
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
has a tridiagonal inverse and
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
holds its factor already.

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
a lower triangular `s@dimension` by `s@dimension` matrix. Both are bare,
where the value and the derivative arrays carry `v1`, `v2`, ...: this is
the one family whose two answers are built from its own value rather
than from its structure, so they are unnamed explicitly to match the
twelve families that never label them.

## Details

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
therefore agrees with [`solve()`](https://rdrr.io/r/base/solve.html) on
the assembled matrix exactly, being the same call; the factor is
`t(chol(M))`, lower triangular with \\M = L L^\top\\ as
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
requires, and reproduces the matrix to \\4 \times 10^{-16}\\.

Both are refused by the generic where the family is rank deficient, a
singular matrix having neither an inverse nor a Cholesky factor.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the two contracts.
