# Default Fourth Derivatives

The method every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
of its own. It applies one product stencil per index tuple directly to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
as
[`param_d3.parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3.parameter.md)
does, with the multiplicities summing to four. The step is
\\\varepsilon^{1/6}\max(1, \|\eta_k\|)\\, about \\2.5 \times 10^{-3}\\
near the origin, and rounding amplified by \\h^{-4}\\ leaves about five
digits: measured against a closed form on a \\2 \times 2\\ covariance
the gap is \\1.2 \times 10^{-4}\\ on entries of size 24. Enough to catch
a wrong closed form, not enough to fit with. No family in this package
reaches it.

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 3, 4)` estimates keyed as
`param_tuple_names(s, 4)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## See also

[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md),
which does the work, and
[`param_d3.parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3.parameter.md)
for the order below.
