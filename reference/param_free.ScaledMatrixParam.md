# Free Vector of a Scaled Parameter

Returns \\g(h)\\ where \\h\\ is the multiple of \\P\\ that `m` is, read
off the ratio at the entry of \\P\\ of largest magnitude, which is where
it is best determined. Exact, and a true inverse of
[`param_value.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.ScaledMatrixParam.md).

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object.

- m:

  A symmetric numeric matrix, a positive multiple of the object's fixed
  matrix, already checked for shape and symmetry by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free` named by `s@free_names`, or
`numeric(0)` for a fixed parameter.

## Details

Three rejections. The ratio must be positive, or `m` is not a positive
multiple. `m` must then agree with `h * p` to \\10^{-8}\\ relative
everywhere, or it is not a multiple of \\P\\ at all: the ratio at one
entry is not enough, since any matrix has *some* ratio there. And for a
fixed parameter, built with `link = NULL`, the multiple must be 1 to
\\10^{-8}\\, the object having no free value to absorb anything else;
the result is then `numeric(0)`.

## See also

[`param_value.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.ScaledMatrixParam.md),
the map this inverts.
