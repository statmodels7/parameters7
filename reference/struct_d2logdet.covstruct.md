# Default Log-Determinant Hessian

Fallback: \\\mathrm{tr}(M^{+} \partial\_{kl} M) - \mathrm{tr}(M^{+}
\partial_k M\\ M^{+} \partial_l M)\\, the derivative of the identity
behind
[`struct_dlogdet`](https://statmodels7.github.io/covstructs7/reference/struct_dlogdet.md).

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A named numeric vector.
