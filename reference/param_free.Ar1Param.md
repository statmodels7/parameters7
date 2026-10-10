# Free Vector of an AR(1) Parameter

Returns the two free values behind an AR(1) matrix. The variance is the
common diagonal entry and the correlation is the first off-diagonal
entry divided by it, both read exactly; the remaining entries of the
matrix are then compared with the pattern that those two values imply.

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- m:

  An AR(1) numeric matrix of side `s@dimension`, already checked for
  shape and symmetry by the generic. Its diagonal must be constant, and
  every entry must equal \\\sigma^2\rho^{\|i-j\|}\\ for the two values
  read off the diagonal and the first off-diagonal.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length 2, named by `s@free_names`.

## Details

A matrix that does not match the implied pattern is rejected and is
never projected onto the nearest AR(1) matrix. A Toeplitz matrix whose
lag-2 entry is not \\\sigma^2\rho^2\\ is not in the family; a caller who
wants the nearest AR(1) matrix computes that projection separately and
passes the result to this method.

## See also

[`param_value.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_value.Ar1Param.md),
the map this inverts.
