# First Derivatives of a Diagonal Parameter

Closed form, and it is the link's own derivative placed on a diagonal:

\$\$\partial_k M = h'(\eta_k)\\ \textstyle\sum\_{i:\\ o_i = k}
E\_{ii},\$\$

where \\o_i\\ is the free value entry \\i\\ belongs to. For a
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
that is one entry, so \\\partial_k M\\ is \\h'(\eta_k) E\_{kk}\\; for a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
the one free value owns the whole diagonal and \\\partial_1 M =
h'(\eta_1) I\\.

\\h'\\ comes from
[`linkfunctions7::dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.html),
so its accuracy is the link's and nothing is differenced.

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

A list of `s@n_free` diagonal matrices named by `s@free_names`, each
`s@dimension` by `s@dimension` with dimnames `v1`, `v2`, ...

## See also

[`diag_owner()`](https://statmodels7.github.io/parameters7/reference/diag_owner.md)
for the ownership map, and
[`param_d2.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.DiagMatrixParam.md)
for the order above.
