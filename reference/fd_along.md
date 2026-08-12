# One Stencil Along One Free Value

Applies a numericals7 central stencil to a function of the free vector,
along one of its components.

## Usage

``` r
fd_along(f, eta, k, order = 1L, h = NULL)
```

## Arguments

- f:

  A function of the free vector.

- eta:

  The free vector.

- k:

  The component to differentiate along.

- order:

  The derivative order.

- h:

  The step;
  [`fd_step`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
  when missing.

## Value

Whatever `f` returns, differentiated.

## Details

The nodes and the weights come from
[`fd_offsets`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
and
[`fd_weights`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
rather than being written out here, so a difference in this package is
the same object the rest of the toolkit differentiates with.
[`fd_derivative`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
cannot be used directly: its `f` maps a vector of points to the values
at those points, and what is differentiated here maps a whole free
vector to a matrix.

## See also

[`fd_weights`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
