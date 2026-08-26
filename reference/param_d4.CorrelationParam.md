# Fourth Derivatives of a Correlation Parameter

Closed form, the Leibniz rule with four differentiations distributed
over the two factors of \\R = LL^\top\\. This is the order at which a
numerical route is least usable, keeping about five digits, and the
spherical construction pays nothing for it: every factor derivative is a
product of trigonometric tables that
[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
has already built.

[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)'s
vanishing rules do most of the work. A quadruple spanning three rows of
\\L\\ contributes nothing at all, since a Leibniz term splits the four
indices between two factors and each factor must stay within one row.
The diagonal is exactly zero.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
`param_tuple_names(s, 4)` and in that order, each with a zero diagonal.

## See also

[`corr_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_derivative.md),
which assembles it,
[`param_d3.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.CorrelationParam.md)
for the order below, and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for the alternative.
