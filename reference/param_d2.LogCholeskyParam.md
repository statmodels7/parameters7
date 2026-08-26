# Second Derivatives of a Log-Cholesky Parameter

Closed form. Differentiating \\\partial_k M = L_k L^\top + L L_k^\top\\
once more,

\$\$\partial\_{kl} M = L\_{kl} L^\top + L_k L_l^\top + L_l L_k^\top + L
L\_{kl}^\top,\$\$

and the factor's second derivative \\L\_{kl}\\ is non-zero only where
\\k = l\\ names a diagonal free value, where it is \\L\_{ii} E\_{ii}\\
again. A below-diagonal value enters \\L\\ linearly, so any second
derivative touching one of those drops the outer terms and leaves the
two cross products.

One consequence is worth knowing when reading a result: \\M\\ is
quadratic in each below-diagonal free value, so the component keyed by
that value twice is a constant in \\\eta\\, and every third and fourth
derivative repeating it is exactly zero.

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

A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
`param_tuple_names(s)` and in that order, each `s@dimension` by
`s@dimension`.

## See also

[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
which assembles it, and
[`param_d1.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.LogCholeskyParam.md)
for the order below.
