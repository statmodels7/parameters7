# Fourth Derivatives of a Simplex Parameter

Closed form, the cumulant recursion at fourth order, which is where the
contract stops. Exact, and that matters most here: a product stencil at
fourth order keeps about five digits, and this array is what a
fourth-order chain rule through a link reads. Every component sums to
zero over the category index, measured at \\1.7 \times 10^{-17}\\.

The array holds \\K(K-1)^4\\ entries, so it is the largest object the
family builds; at \\K = 4\\ that is 324.

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

A list of `choose(s@n_free + 3, 4)` numeric vectors of length \\K\\,
keyed as `param_tuple_names(s, 4)` and in that order, each summing to 0.

## See also

[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md),
which builds the array,
[`param_d3.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.SimplexParam.md)
for the order below, and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for the alternative.
