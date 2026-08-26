# Fourth Log-Determinant Derivatives of a Matrix Logarithm Parameter

Closed form, and identically zero, for the reason it is zero at second
and third order: \\\log\|M\| = \mathrm{tr}(S)\\ is linear in the free
vector. The zeros are exact.

This is the order at which a numerical fallback is least usable, so a
family whose log-determinant is not linear gains most from a closed form
here; see
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the measured accuracy of the alternative.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A named numeric vector of zeros.
