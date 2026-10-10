# Default Inverse Map

The method that every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
of its own. It always signals an error naming the family; it is the one
generic whose base method rejects the call instead of computing. An
inverse obtained by minimizing \\\lVert V(\eta) - m \rVert\\ would
return an \\\eta\\ even for a matrix that lies outside the family's set,
and the caller could not distinguish that result from a correct one, so
the inverse map is either written out exactly or rejected.

Every family in this package writes its inverse out, so this method is
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
