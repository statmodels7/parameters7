# Free Vector of a Simplex Parameter

Returns the additive log-ratio, \\\eta_a = \log(\pi_a/\pi_K)\\, exact
and a true inverse of
[`param_value.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_value.SimplexParam.md):
the round trip closes to \\2 \times 10^{-16}\\.

## Arguments

- s:

  A
  [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  object.

- m:

  A probability vector of length \\K\\: strictly positive and summing to
  1.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

Two rejections, both with their own message. A vector with a
non-positive entry is outside the **open** simplex, and \\\log 0\\ is
not finite; a vector that does not sum to 1 is not a probability vector,
and it is **not renormalized**, a silent repair being the kind of thing
that hides a caller's defect for a long time.

The first rejection is the one a fit runs into.
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
saturates a large free value to an exact 0 in the reference category, so
a value produced at the boundary cannot be inverted back.

## See also

[`param_value.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_value.SimplexParam.md),
the map this inverts.
