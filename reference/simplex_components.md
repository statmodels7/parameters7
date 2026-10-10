# Extract Named Components From Softmax Tensors

Slices the tensors of
[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
into the named list that a derivative generic returns, one entry per
distinct index tuple, keyed as `param_tuple_names(s, order)`.

## Usage

``` r
simplex_components(s, tens, order, wrap = identity)
```

## Arguments

- s:

  The
  [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  to which the tuples belong.

- tens:

  The tensor list from
  [`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md).

- order:

  The order to extract: 1, 2, 3 or 4. `tens` must carry at least this
  order.

- wrap:

  A function applied to each raw slice before it is stored, `identity`
  by default. Every caller in the package uses the default.

## Value

A named list of `choose(s@n_free + order - 1, order)` entries, keyed as
`param_tuple_names(s, order)` and in that order, each entry `wrap()`'s
result.

## Details

The `wrap` argument is a function applied to each slice before it is
stored, `identity` by default. A simplex component is a vector of length
\\K\\ and is returned as it stands.
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
does not call this function: its derivative method slices the tensors
itself and places each slice in one row of a \\K \times K\\ matrix of
zeros.

## See also

[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
for the arrays sliced, and
[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md),
which slices the same arrays row by row for a transition matrix.
