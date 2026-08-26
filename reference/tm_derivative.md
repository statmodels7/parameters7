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
is evaluated once per row and
[`simplex_components()`](https://statmodels7.github.io/parameters7/reference/simplex_components.md)
slices it with a `wrap` that does the embedding. Nothing of size \\K^2
(K(K-1))^{\text{order}}\\ is built: the cross-row components are never
computed, only skipped.

## See also

[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
and
[`simplex_components()`](https://statmodels7.github.io/parameters7/reference/simplex_components.md),
which do the per-row work, and
[`tm_positions()`](https://statmodels7.github.io/parameters7/reference/tm_positions.md)
for the row map.
