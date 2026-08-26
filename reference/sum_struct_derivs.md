# Assemble a Sum of Fixed Matrices' Derivatives of a Given Order

The value is linear in the weights, so a component is zero unless every
index of the tuple names the same free value, and it is then that
weight's \\m\\-th derivative times its own fixed matrix.

## Usage

``` r
sum_struct_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`SumStructParam()`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` symmetric matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.

## Details

No arithmetic is done on the matrices at all: a surviving component is
one scalar times a component the object has held since construction, and
the rest share a single zero matrix. At \\K = 2\\ that is 1 of the 3
second-order components, 2 of 4 at third order and 3 of 5 at fourth.

## See also

[`sum_struct_weight_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_weight_derivs.md)
for the scalars, and
[`param_d1.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.SumStructParam.md),
which calls this.
