# One Product Stencil for a Mixed Partial Derivative

Differentiates `f` once in each component the index tuple names, a
central factor per distinct component of the order its multiplicity
asks: two points at order one, three at two, the five-point forms at
three and four. The result is a single product stencil, not a
composition of lower-order numerical derivatives.

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
