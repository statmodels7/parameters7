# First Derivatives of a Log-Cholesky Parameter

Closed form, by the Leibniz rule on \\M = L L^\top\\. Writing \\L_k\\
for the derivative of the factor in the \\k\\-th free value,

\$\$\partial_k M = L_k L^\top + L L_k^\top.\$\$

The factor's derivative is a single-entry matrix: \\L\_{ii} E\_{ii}\\
for a diagonal free value, because the parametrization holds its
logarithm and \\\partial \exp(\eta_i)/\partial \eta_i = \exp(\eta_i) =
L\_{ii}\\, and \\E\_{ij}\\ for a value below the diagonal, which enters
\\L\\ linearly. Each component therefore costs one row and one column of
\\L\\, never a matrix product.

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

A list of `s@n_free` symmetric matrices named by `s@free_names`, each
`s@dimension` by `s@dimension` with dimnames `v1`, `v2`, ...

## See also

[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
which assembles it,
[`chol_dfactor()`](https://statmodels7.github.io/parameters7/reference/chol_dfactor.md)
for \\\partial^S L\\, and
[`param_d2.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.LogCholeskyParam.md)
for the order above.
