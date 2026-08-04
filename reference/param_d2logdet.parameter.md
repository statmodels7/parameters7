# Default Log-Determinant Hessian

Fallback: \\\mathrm{tr}(M^{+} \partial\_{kl} M) - \mathrm{tr}(M^{+}
\partial_k M\\ M^{+} \partial_l M)\\, the derivative of the identity
behind
[`param_dlogdet`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md).

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A named numeric vector.
