# Free Vector of a Compound Symmetry Parameter

Returns the two free values behind a compound-symmetric matrix. The
variance is the common diagonal entry and the correlation is the common
off-diagonal entry divided by it, both read exactly and then carried
onto the free scale by the two links. The round trip closes up to
rounding.

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
never averaged into the nearest compound-symmetric matrix, because a
silent projection would hide a mistake in the caller's matrix. A caller
who wants the nearest such matrix computes it separately and passes the
result in.

## See also

[`param_value.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_value.CompoundSymmetryParam.md),
the map this inverts.
