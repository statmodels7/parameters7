# First Derivatives of a Log-Cholesky Structure

Closed form. Writing \\L_k\\ for the derivative of the factor in the
\\k\\-th free value, \\\partial_k M = L_k L^\top + L L_k^\top\\. The
factor's derivative is \\L\_{ii} E\_{ii}\\ for a diagonal value, because
the parametrisation is its logarithm, and \\E\_{ij}\\ for a value below
the diagonal.

## Arguments

- s:

  A
  [`LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/LogCholeskyStruct.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A named list of symmetric matrices.
