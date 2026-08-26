# Derivatives of an AR(1) Parameter

Closed form at all four orders. One page covers
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for this family because the four are the same product at different
orders: \\M = \sigma^2 P(\rho)\\, so a component with \\a\\ scale
indices and \\b\\ correlation indices is the \\a\\-th derivative of the
scale times the \\b\\-th derivative of the pattern.

The pattern's derivatives are the ones that need work here. An entry is
\\\rho^{\|i-j\|}\\, a power where
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)'s
pattern is linear, so each is composed with the rhobit link to the order
asked, through
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md).
Each distinct lag is composed once.

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A named list of `order + 1` symmetric `s@dimension` by `s@dimension`
matrices; see **Details** for the keying of each order.

## Details

The four methods return lists of `order + 1` matrices, keyed by the
tuple names of their own order:

- [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md):
  2 entries, keyed by `free_names`.

- [`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md):
  3 entries, `param_tuple_names(s, 2)`.

- [`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md):
  4 entries, `param_tuple_names(s, 3)`.

- [`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md):
  5 entries, `param_tuple_names(s, 4)`.

Every component with at least one correlation index has a zero diagonal,
the lag-0 entries of the pattern being the constant 1.

## See also

[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md),
which assembles all four,
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md)
for the pattern and its derivatives, and
[`param_d1.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.CompoundSymmetryParam.md)
for the sibling family, whose pattern is linear.
