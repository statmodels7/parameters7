# Matrix of a Diagonal Parameter

Returns \\M = \mathrm{diag}(h(\eta_1), \dots, h(\eta_p))\\, the inverse
link applied to each free value and placed on the diagonal. For a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
the one free value is recycled, giving \\h(\eta_1) I\\. Positive
definiteness follows from the link, whose range
[`check_positive_link()`](https://statmodels7.github.io/parameters7/reference/check_positive_link.md)
restricted to the non-negative half line at construction, so nothing is
tested here.

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A diagonal positive definite `s@dimension` by `s@dimension` numeric
matrix with dimnames `v1`, `v2`, ...

## See also

[`diag_entries()`](https://statmodels7.github.io/parameters7/reference/diag_entries.md),
which applies the link, and
[`param_free.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_free.DiagMatrixParam.md)
for the inverse.
