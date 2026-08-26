# Extract Named Components From Softmax Tensors

Slices the tensors of
[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
into the named list a derivative generic returns, one entry per distinct
index tuple, keyed as `param_tuple_names(s, order)`.

## Usage

``` r
simplex_components(s, tens, order, wrap = identity)
```

## Arguments

- s:

  The parameter the tuples belong to, a
  [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  or, through `wrap`, one row's worth of a
  [`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md).

- tens:

  The tensor list from
  [`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md).

- order:

  The order to extract: 1, 2, 3 or 4. `tens` must carry at least this
  order.

- wrap:

  A function applied to each raw slice before it is stored. `identity`
  by default, which is the simplex's own case;
  [`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
  passes a function embedding the slice in a row.

## Value

A named list of `choose(s@n_free + order - 1, order)` entries, keyed as
`param_tuple_names(s, order)` and in that order, each entry `wrap()`'s
result.

## Details

The `wrap` argument is how
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
reuses this. A simplex component is a vector of length \\K\\ and is
returned as it stands; a transition matrix's is that vector placed into
one row of a \\K \times K\\ matrix of zeros, and the caller passes the
function that does the placing.

## See also

[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
for the arrays sliced, and
[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md)
for the `wrap` that embeds a row.
