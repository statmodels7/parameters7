# One Product Stencil for a Mixed Partial Derivative

Differentiates `f` once in each component the index tuple names, a
central factor per distinct component of the order its multiplicity
asks. The result is a single product stencil, not a composition of
lower-order numerical derivatives.

## Usage

``` r
mixed_stencil(f, eta, tuple)
```

## Arguments

- f:

  A function of the free vector.

- eta:

  The point.

- tuple:

  An integer vector of component indices, possibly repeated.

## Value

The stencil's value, shaped like `f(eta)`.

## Details

Each factor's nodes and weights are numericals7's, read from
[`fd_offsets`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
and
[`fd_weights`](https://statmodels7.github.io/numericals7/reference/fd_weights.html).
They were transcribed here once – two points at order one, three at two,
the five-point forms at three and four – and a table of stencil
coefficients written out in a second place is a table that can disagree
with the first.

## See also

[`fd_weights`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
