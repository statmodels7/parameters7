# Free Vector of a Simplex Parameter

Returns the additive log-ratio, \\\eta_a = \log(\pi_a/\pi_K)\\, exact
and a true inverse of
[`param_value.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_value.SimplexParam.md),
up to rounding.

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

Three rejections, each with its own message. An argument that is not a
numeric vector of length \\K\\ is rejected. A vector with a non-positive
or missing entry is outside the **open** simplex, and \\\log 0\\ is not
finite. A vector that does not sum to 1 is not a probability vector, and
it is **not renormalized**, because a silent repair would hide an error
in the caller.

The second rejection is the one that a fit meets: a large positive free
value makes
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
return an exact 0 in the reference category, so a value produced at the
boundary cannot be inverted back.

## See also

[`param_value.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_value.SimplexParam.md),
the map this inverts.
