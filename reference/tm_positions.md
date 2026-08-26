# Row and Chart Coordinate of Each Free Value

Returns, for each of the \\K(K-1)\\ free values, the row it belongs to
and its position inside that row's simplex chart, in the order the free
vector uses. Every method of the family reads it to split the free
vector into rows.

## Usage

``` r
tm_positions(s)
```

## Arguments

- s:

  A
  [`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md)
  object, whose `param_params$n_state` is read.

## Value

A list with two integer vectors of length `s@n_free`: `row`, with values
in `1:K`, and `coord`, with values in `1:(K-1)`.

## Details

The free vector runs row by row, so at \\K = 3\\ the rows are
`1 1 2 2 3 3` and the coordinates `1 2 1 2 1 2`, which labels the free
values `alr1.1 alr1.2 alr2.1 alr2.2 alr3.1 alr3.2`.

## See also

[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
which turns this into `free_names`, and
[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md),
which uses it to place a row's component.
