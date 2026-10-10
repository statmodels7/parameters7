# Derivative Components of a Transition Matrix Parameter

Assembles a whole derivative order from the row-wise simplex tensors. A
component whose free values span two rows is the zero matrix; one whose
free values all belong to row \\i\\ is that row's
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
component embedded in row \\i\\ of a matrix of zeros. The four
derivative methods of the family are one call each to this function.

## Usage

``` r
tm_derivative(s, eta, order)
```

## Arguments

- s:

  A
  [`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` \\K \times K\\ matrices
keyed as `param_tuple_names(s, order)` and in that order, most of them
exactly zero. Each non-zero one is supported on a single row.

## Details

The rows are parametrized independently, so
[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
is evaluated once per row, and each slice of a row's tensor is embedded
in that row of a \\K \times K\\ matrix of zeros. The cross-row
components are not computed; each is returned as a matrix of zeros.

## See also

[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md),
which does the per-row work, and
[`tm_positions()`](https://statmodels7.github.io/parameters7/reference/tm_positions.md)
for the row map.
