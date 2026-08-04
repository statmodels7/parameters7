# Extract Named Components From Softmax Tensors

Slices the tensors of
[`simplex_tensors`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
into the named lists the derivative generics return.

## Usage

``` r
simplex_components(s, tens, order, wrap = identity)
```

## Arguments

- s:

  The parameter the tuples belong to.

- tens:

  The tensor list.

- order:

  The order to extract.

- wrap:

  A function applied to each raw slice, for the transition matrix to
  embed a row; the identity here.

## Value

A named list keyed as `param_tuple_names(s, order)`.
