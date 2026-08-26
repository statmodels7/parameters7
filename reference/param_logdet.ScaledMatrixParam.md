# Log-Determinant of a Scaled Parameter

Closed form,

\$\$\log\|M\|\_+ = r \log h(\eta) + \log\|P\|\_+,\$\$

with \\r\\ the rank and \\\log\|P\|\_+\\ a constant of the object,
computed once at construction. The scale enters once per non-zero
eigenvalue, which is where the factor \\r\\ comes from. One expression
covers both cases: it is the log-determinant where the family is of full
rank and the log pseudo-determinant otherwise. Under the log link it is
\\r\eta + \log\|P\|\_+\\, so a step of 1 in the free value moves it by
\\r\\.

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic. For a fixed parameter it is `numeric(0)` and the
  result is `logdet_p` alone.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number.

## See also

[`param_dlogdet.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.ScaledMatrixParam.md)
for its gradient, and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for the rank it uses.
