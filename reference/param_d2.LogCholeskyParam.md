# Second Derivatives of a Log-Cholesky Parameter

Closed form. Differentiating \\\partial_k M = L_k L^\top + L L_k^\top\\
again gives \\\partial\_{kl} M = L\_{kl} L^\top + L_k L_l^\top + L_l
L_k^\top + L L\_{kl}^\top\\, and the factor's second derivative
\\L\_{kl}\\ is non-zero only when \\k = l\\ is a diagonal value, where
it is \\L\_{ii} E\_{ii}\\ again.

## Arguments

- s:

  A
  [`LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A named list of symmetric matrices.
