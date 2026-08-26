# Free Vector of a Diagonal Parameter

Returns \\g\\ applied to the diagonal of `m`, the link in the forward
direction, after checking that `m` really is in the family's set. Exact,
and a true inverse of
[`param_value.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.DiagMatrixParam.md),
the link being a bijection.

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- m:

  A diagonal positive definite `s@dimension` by `s@dimension` numeric
  matrix, already checked for shape and symmetry by the generic. It must
  have a constant diagonal when `s` came from
  [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md).

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

Three rejections, each with its own message. `m` must be diagonal,
tested as `max(abs(off)) > 1e-10 * max(1, max(abs(d)))` on the
off-diagonal part, so an asymmetry of rounding size passes and a real
off-diagonal entry does not. Its diagonal entries must all be positive.
And for a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
they must all be **equal**, tested by
`diff(range(d)) > 1e-10 * max(abs(d))`, a diagonal that varies not being
a scalar multiple of the identity; only the first entry is then read.

## See also

[`param_value.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.DiagMatrixParam.md),
the map this inverts.
