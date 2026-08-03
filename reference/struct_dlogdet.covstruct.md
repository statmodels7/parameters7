# Default Log-Determinant Gradient

Fallback: \\\mathrm{tr}(M^{+} \partial_k M)\\, with the pseudo-inverse
formed from the directions the declared rank keeps, which is the
ordinary inverse when the family is of full rank.

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
