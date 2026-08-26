# Assemble a Sum of Fixed Matrices' Log-Determinant Derivatives

Carries
[`sum_struct_trace_term()`](https://statmodels7.github.io/parameters7/reference/sum_struct_trace_term.md)'s
expansion in the weights onto the free scale. The map from free values
to weights is diagonal, so the chain rule groups the tuple by index and
takes one set partition per group, each block contributing a derivative
of the inverse link and one differentiation in that weight.

## Usage

``` r
sum_struct_logdet_derivs(s, eta, order)
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

A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
as `param_tuple_names(s, order)` and in that order.

## Details

This is Faa di Bruno with a diagonal inner map, written out here instead
of going through
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md),
because the outer function is a derivative in several weights at once,
never a univariate composition. The product over the groups' partitions
is the grid the loop walks; the partitions come from
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html),
the one enumeration the toolkit keeps.

It is the dearest quantity this family computes: at side 6 and \\K = 5\\
the fourth order costs 0.013 s against 0.0007 s for the value's own
fourth derivatives, the expansion evaluating a chain of matrix products
per ordering. Against one stencil on the analytic order below, the four
orders agree to \\7 \times 10^{-11}\\, \\1 \times 10^{-11}\\, \\3 \times
10^{-12}\\ and \\3 \times 10^{-13}\\.

## See also

[`sum_struct_trace_term()`](https://statmodels7.github.io/parameters7/reference/sum_struct_trace_term.md)
for the inner expansion, and
[`param_dlogdet.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.SumStructParam.md),
which calls this.
