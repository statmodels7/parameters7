# Factor of a Log-Cholesky Parameter

Returns \\L\\ by assembling it from the free vector, which is all this
family's parametrization is: no factorization is taken, because the
factor is what \\\eta\\ holds. It is the only family here whose factor
is free, and it is \\O(p^2)\\ against the base class's \\O(p^3)\\
Cholesky.

The diagonal is `exp(eta[1:p])` and the entries below it are the
remaining free values, so a caller who wants a standard deviation off
the factor can read it from `eta` directly.

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

A `s@dimension` by `s@dimension` lower triangular numeric matrix with a
positive diagonal, satisfying `L %*% t(L) == param_value(s, eta)`
exactly, and carrying no dimnames.

## See also

[`chol_assemble()`](https://statmodels7.github.io/parameters7/reference/chol_assemble.md),
which does the assembly, and
[`param_factor.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_factor.matrix_parameter.md)
for what every other family pays.
