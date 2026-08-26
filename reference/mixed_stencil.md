# One Product Stencil for a Mixed Partial Derivative

Differentiates `f` once in each component an index tuple names, using a
central factor per distinct component of the order that component's
multiplicity asks for. The whole tensor product is summed in one pass,
so the result is a single stencil, never a composition of lower-order
numerical derivatives.

## Usage

``` r
mixed_stencil(f, eta, tuple)
```

## Arguments

- f:

  A function of the free vector, returning anything that supports `+`
  and scalar multiplication.

- eta:

  The point, a numeric vector.

- tuple:

  An integer vector of component indices, possibly with repetitions, as
  [`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
  returns. Its length is the derivative order.

## Value

The stencil's value, shaped like `f(eta)`.

## Details

Each factor's nodes and weights come from
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
and
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
at accuracy 2: two points at order one, three at order two, five at
orders three and four. They are read from there instead of being written
out here, because a table of stencil coefficients kept in a second place
is a table that can come to disagree with the first.

The step for each factor is
[`fd_step()`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
at that component's value and at the **total** order of the tuple, so
all factors of one component share a step and the division at the end is
by \\\prod_j h_j^{m_j}\\. Nodes whose weight is zero are skipped, so a
repeated-index factor costs one evaluation fewer than its node count
suggests.

## See also

[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
for the weights,
[`fd_along()`](https://statmodels7.github.io/parameters7/reference/fd_along.md)
for the single-component case, and
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md),
the two callers.
