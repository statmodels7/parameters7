# The Quantities of the Inverted Family

Declares the quantities of the family that an
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md),
[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
or
[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
parameter inverts, with the same Jacobian, and adds `inverse = TRUE`.
The free values of an inverse family are those of the inner family, so
the inner family's quantities are functions of them; they describe the
inverse of the matrix that
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns. For
[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
they are the marginal variance and the lag-one correlation of the AR(1)
process whose precision the matrix is. An inverse of an inverse clears
the mark.

## Arguments

- s:

  An
  [`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md)
  object, whose inner family is read.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- ...:

  Passed to the inner family's method.

## Value

`NULL` when the inner family declares nothing, otherwise a list as
described in
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md),
with `inverse` set.
