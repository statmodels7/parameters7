# Derivative Components of a Transition Matrix Parameter

Assembles a derivative order from the row-wise simplex tensors: a
component whose free values span two rows is zero, and one inside a row
embeds that row's simplex component.

## Usage

``` r
tm_derivative(s, eta, order)
```

## Arguments

- s:

  A
  [`TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The derivative order, 1 to 4.

## Value

A named list of \\K \times K\\ matrices.
