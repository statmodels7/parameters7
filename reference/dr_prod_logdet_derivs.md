# Assemble a Scales-Times-Correlation Log-Determinant Derivative

The log-determinant is \\2\sum_j \log d_j + \log\lvert R\rvert\\, a sum
of one term per scale and one term for the correlation. A component is
therefore the scale's own where every index names **one** scale
coordinate, the correlation's own where every index is a correlation
coordinate, and 0 otherwise.

## Usage

``` r
dr_prod_logdet_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
as `param_tuple_names(s, order)` and in that order.

## Details

Two separabilities at once, and it is worth keeping them apart: the
scales separate from each other, so a component naming two different
scales is zero, and the correlation separates from the scales, so a
mixed component is zero. What survives is \\2\\\partial^m \log d_k\\ on
the diagonal of the scale block, through
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md),
and whatever the correlation family answers.

## See also

[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for the scale terms, and
[`param_dlogdet.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.DrProdParam.md),
which calls this.
