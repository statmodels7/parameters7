# One Stencil Along One Free Value

Applies one central stencil to a function of the free vector, along a
single component of it, and returns the derivative in whatever shape `f`
produces. The nodes and weights come from
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
and
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
at accuracy 2, so a difference taken here uses the same stencils as the
rest of the toolkit.

## Usage

``` r
fd_along(f, eta, k, order = 1L, h = NULL)
```

## Arguments

- f:

  A function of the free vector, returning anything that supports `+`
  and scalar multiplication.

- eta:

  The free vector, numeric.

- k:

  The component to differentiate along: a position in `1:length(eta)`.

- order:

  The derivative order: 1, 2, 3 or 4.

- h:

  The step. `NULL`, the default, takes
  [`fd_step()`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
  at `eta[k]` and this order. Pass a value when several stencils have to
  share one step, as a mixed derivative does.

## Value

Whatever `f` returns, differentiated `order` times in component `k`.

## Details

The sum \\h^{-k}\sum_j w_j f(\eta + s_j h e_k)\\ is accumulated term by
term, skipping the nodes whose weight is zero, so `f` is called once per
non-zero weight: three times for a first or second derivative, five for
a third or fourth. The accumulator starts at `NULL` and takes the shape
of the first term, which is why `f` may return a matrix, a vector or a
scalar without this function knowing which.

[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
cannot serve here. Its `f` maps a vector of points to the values at
those points, and what is differentiated here maps a whole free vector
to a matrix, so the nodes have to be placed along one coordinate by
hand.

## See also

[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
for the weights,
[`fd_step()`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
for the step, and
[`mixed_stencil()`](https://statmodels7.github.io/parameters7/reference/mixed_stencil.md),
which composes factors across several components.
