# Third Derivatives of a Simplex Parameter

Closed form, the same product rule applied a third time, which is the
cumulant recursion of a categorical indicator at third order. Exact,
where a family without a closed form would get a product stencil good to
about six digits. Every component sums to zero over the category index.

## Arguments

- s:

  A
  [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 2, 3)` numeric vectors of length \\K\\,
keyed as `param_tuple_names(s, 3)` and in that order, each summing to 0.

## See also

[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md),
which builds the array, and
[`param_d4.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.SimplexParam.md)
for the order above.
