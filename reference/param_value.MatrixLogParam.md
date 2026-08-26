# Value of a Matrix Logarithm Parameter

Returns \\M = \exp(S)\\, through the eigendecomposition of \\S\\: \\Q
\exp(\Lambda) Q^\top\\, the exponential applied to the eigenvalues.
Positive definiteness needs no checking, every \\e^{\lambda}\\ being
positive whatever \\S\\ is. The cost is one \\O(p^3)\\
eigendecomposition.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A symmetric positive definite `s@dimension` by `s@dimension` numeric
matrix with dimnames `v1`, `v2`, ...

## See also

[`param_free.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_free.MatrixLogParam.md)
for the inverse,
[`mlog_s()`](https://statmodels7.github.io/parameters7/reference/mlog_s.md)
for the assembly of \\S\\, and
[`param_solve.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_solve.MatrixLogParam.md),
which is this map at \\-\eta\\.
