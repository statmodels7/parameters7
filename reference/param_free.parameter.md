# Rejection to Invert Without a Closed Form

The method every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
of its own. It always signals an error, naming the family, and it is the
one generic whose base method refuses instead of computing. Nothing here
is approximated, which is the point: an inverse obtained by minimizing
\\\lVert V(\eta) - m \rVert\\ would hand back a plausible \\\eta\\ for a
matrix that is nowhere in the family's set, and the caller could not
tell that answer from a correct one. The inverse map is written out
exactly or refused.

All fifteen families in this package write theirs out, so this method is
reached only by a family defined elsewhere.

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, whose `param_name` goes into the message.

- m:

  A value of the family's shape. Never read.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

Never returns; always signals an error.

## See also

[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
for the generic and for what each family's own method rejects, and
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
for the forward map.
