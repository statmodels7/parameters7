# Derivatives of a Compound Symmetry Parameter

Closed form at all four orders. One page covers
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for this family because the four are the same product taken at different
orders: the value is \\M = \sigma^2 P(\rho)\\ with \\P(\rho) = I +
\rho(J - I)\\ linear in the correlation, so a component with \\a\\ scale
indices and \\b\\ correlation indices is

\$\$\partial^{(a,b)} M = \frac{\mathrm{d}^a
\sigma^2}{\mathrm{d}\eta_1^a} \cdot \frac{\mathrm{d}^b
P}{\mathrm{d}\eta_2^b},\$\$

one derivative of each chain. Only the counts \\a\\ and \\b\\ matter, so
a tuple is fully described by them and no order is a special case.

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
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

The lists are this short because there are two free values: at \\d = 2\\
the count \\\binom{d+k-1}{k}\\ is \\k+1\\.

Every component with at least one correlation index has a zero diagonal,
the diagonal of \\P\\ being the constant 1; and any component whose
correlation indices number two or more carries only the link's own
higher derivative, \\P\\ being linear in \\\rho\\.

## See also

[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md),
which assembles all four,
[`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md)
for the pattern and its derivatives, and
[`param_dlogdet.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CompoundSymmetryParam.md)
for the log-determinant's orders.
