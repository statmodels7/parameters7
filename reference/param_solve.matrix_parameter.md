# Default Solve

The method that every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
of its own. It takes the Cholesky factor \\L\\ from
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
and applies it twice, `backsolve(t(l), forwardsolve(l, b))`, which is
\\L^{-\top}L^{-1}B = M^{-1}B\\. The inverse is not formed, and each
triangular system is solved once.

The result is exact up to rounding, as
[`base::solve()`](https://rdrr.io/r/base/solve.html) is, which is why
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
does not list it. The generic has already rejected a rank-deficient
family and one whose value is not a symmetric matrix, and filled `b`
with the identity when the caller left it out, so by the time this runs
`b` is a matrix of the right height.

## Arguments

- s:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object, of full rank.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- b:

  A numeric matrix with `s@dimension` rows, already coerced from a
  vector and defaulted to the identity by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric matrix with `s@dimension` rows and as many columns as `b`.

## Details

Positive definiteness is decided inside
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md),
from the eigenvalues, before any factorization is attempted; see
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
for why the verdict does not come from whether
[`base::chol()`](https://rdrr.io/r/base/chol.html) signals an error.

The cost is one \\O(p^3)\\ factorization plus \\O(p^2)\\ per column of
`b`. The families that override this method are listed under
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md).

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for the generic and the argument handling,
[`param_factor.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_factor.matrix_parameter.md)
for the factor, and
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
for the definiteness test.
