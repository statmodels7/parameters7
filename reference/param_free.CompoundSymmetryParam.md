# Free Vector of a Compound Symmetry Parameter

Returns the two free values behind a compound-symmetric matrix. The
variance is the common diagonal entry and the correlation is the common
off-diagonal entry divided by it, both read exactly and then carried
onto the free scale by the two links. The round trip closes to \\3
\times 10^{-16}\\.

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object.

- m:

  A compound-symmetric numeric matrix of side `s@dimension`, already
  checked for shape and symmetry by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length 2, named by `s@free_names`.

## Details

A matrix whose diagonal is not constant, or whose off-diagonal entries
are not all equal, is **rejected** with a message naming which, and is
never averaged into the nearest compound-symmetric matrix. A silent
projection would hide the caller's mistake, and a caller who wants the
nearest such matrix can average the entries and invert that.

## See also

[`param_value.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_value.CompoundSymmetryParam.md),
the map this inverts.
