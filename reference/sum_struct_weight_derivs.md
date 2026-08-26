# Derivatives of the Weights of a Sum of Fixed Matrices

Returns each weight and its first four derivatives in the free value
that carries it, for every component at once, as a matrix with one row
per order.

## Usage

``` r
sum_struct_weight_derivs(s, eta)
```

## Arguments

- s:

  A
  [`SumStructParam()`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

## Value

A 5 by \\K\\ numeric matrix, row \\m+1\\ holding the \\m\\-th derivative
of the inverse link at each free value, so row 1 is the weights
themselves.

## Details

Each weight depends on one free value only, so the table is complete:
there are no cross-derivatives between weights to record, and that is
the whole reason a derivative of the value naming two weights is zero.
It is computed once per call and read for every component of the order.

## See also

[`sum_struct_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_derivs.md)
and
[`sum_struct_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_logdet_derivs.md),
the two callers.
