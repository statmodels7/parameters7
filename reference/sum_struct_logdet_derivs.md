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
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md),
because the outer function is a derivative in several weights at once,
never a univariate composition. The product over the groups' partitions
is the grid the loop walks; the partitions come from
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html),
the one enumeration the toolkit keeps.

Where the family is rank deficient, the derivatives are those of the log
pseudo-determinant: the inverse in the expansion is replaced by \\(M +
ZZ^\top)^{-1} - ZZ^\top\\, with \\Z\\ the orthonormal basis of the
declared null space, which does not move with the free vector.

It is the most expensive quantity of this family, the expansion
evaluating a chain of matrix products for each ordering.

## See also

[`sum_struct_trace_term()`](https://statmodels7.github.io/parameters7/reference/sum_struct_trace_term.md)
for the inner expansion, and
[`param_dlogdet.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.SumStructParam.md),
which calls this.
