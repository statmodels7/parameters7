# Assemble a Sum of Fixed Matrices' Log-Determinant Derivatives

Carries the trace expansion in the weights onto the free scale. The map
is diagonal, so the chain rule groups the tuple by index and takes one
set partition per group, each block contributing a derivative of the
inverse link and one differentiation in that weight.

## Usage

``` r
sum_struct_logdet_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`SumStructParam`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order, 1 to 4.

## Value

A named numeric vector keyed as `param_tuple_names(s, order)`.
