# Third Derivatives of a Log-Cholesky Parameter

Closed form by the Leibniz rule on \\M = L L^\top\\: each component
distributes its three differentiations over the two factors, and a
factor differentiated more than once survives only in a repeated
diagonal direction, where every derivative of \\e^{\eta_k}\\ is itself.

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
