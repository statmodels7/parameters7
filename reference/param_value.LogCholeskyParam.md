# Matrix of a Log-Cholesky Parameter

Returns \\M = L L^\top\\, with \\L\\ assembled from the free vector by
[`chol_assemble()`](https://statmodels7.github.io/parameters7/reference/chol_assemble.md)
and the product taken as `tcrossprod(l)`. Positive definiteness needs no
checking: \\L\\ is triangular with a diagonal that is an exponential, so
it has full rank at every finite \\\eta\\, and \\L L^\top\\ is therefore
in the interior of the cone. The cost is one \\O(p^3)\\ product and no
factorization.

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A symmetric positive definite `s@dimension` by `s@dimension` numeric
matrix, with dimnames `v1`, `v2`, ...

## See also

[`param_free.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_free.LogCholeskyParam.md)
for the inverse,
[`param_factor.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_factor.LogCholeskyParam.md)
for \\L\\ itself, and
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for the ordering of `eta`.
